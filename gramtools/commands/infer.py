## @file
# From quasimap output (allele coverages), infer a personalised reference.
# Personalised reference can be produced as fasta and/or vcf.
# Reads can subsequently be mapped onto the inferred reference for variant calling: see `discover.py`.

import json
import os
import logging

import vcf
import collections

from .. import genotyper
from .. import paths
from .. import version
from .. import prg_local_parser
# For reference of vcf record attributes: cf https://pyvcf.readthedocs.io/en/latest/API.html#vcf-model-record

log = logging.getLogger('gramtools')


def parse_args(common_parser, subparsers):
    parser = subparsers.add_parser('infer',
                                   parents=[common_parser])
    parser.add_argument('--run-dir', '--run-directory',
                        help='Common directory for gramtools running commands. Will contain infer outputs.'
                        'The outputs are: the personalised reference as a (haploid) fasta,'
                        'and the personalised reference as a vcf. '
                             'Note: the sample name used in the vcf output is the --run-dir directory name.',
                        dest='run_dir',
                        type=str,
                        required=True)

    parser.add_argument('--vcf-mode',
                        help='If "single", the output vcf will not contain alleles with zero-coverage in successfully genotyped sites.'
                             'If "population", all alleles will be represented for all sites.'
                             'Default is "single".',
                        type=str,
                        required=False,
                        choices=['single','population'],
                        default='single')

    parser.add_argument('--mean-depth',
                        help='The average read depth.'
                             'By default, the estimation from `quasimap` is used.',
                        type=int,
                        required=False)


    parser.add_argument('--error-rate',
                        help='The estimated per-base error probability in the reads.'
                             'By default, the estimation from \`quasimap\` is used.',
                        type=float,
                        required=False)


    parser.add_argument('--haploid',
                        help='',
                        action='store_true',
                        required=False)


def run(args):
    log.info('Start process: infer')
    _paths = paths.generate_infer_paths(args)

    ## Read in the genotyping parameters from the read stats file.
    if not os.path.exists(_paths['read_stats']) and (args.error_rate is None or args.mean_depth is None):
        log.error("One of --error-rate or --mean-depth genotyping parameter needs to be read from file {} "
                  "but it does not exist."
                  "This file should have been generated at `quasimap` stage.".format(_paths['read_stats']))

    with open(_paths["read_stats"]) as f:
        read_stats = json.load(f)

    mean_depth, error_rate = args.mean_depth, args.error_rate
    if args.mean_depth is None:
        mean_depth = int(read_stats["Read_depth"]["Mean"])

    if args.error_rate is None:
        error_rate = float(read_stats["Quality"]["Error_rate_mean"])


    # Load coverage stats
    all_per_base_coverage = _load_per_base_coverage(_paths['allele_base_coverage'])
    allele_groups, all_groups_site_counts = \
        _load_grouped_allele_coverage(_paths['grouped_allele_counts_coverage'])

    # Set up the genotyping (iterator).
    genotypers = _max_likelihood_alleles(mean_depth,
                                             error_rate,
                                             all_per_base_coverage,
                                             all_groups_site_counts,
                                             allele_groups)

    # Make the vcf
    log.info("Inferring personalised vcf. Output at {}".format(_paths["inferred_vcf"]))
    allele_indexes = _dump_vcf(genotypers, _paths, args)


    # Make the fasta
    log.info("Generating personalised fasta. Output at {}".format(_paths["inferred_fasta"]))
    allele_indexes = iter(allele_indexes)

    # Get the reference ID. This is needed in the personalised fasta so that vcf produced against it in `discover`
    # Has CHROM column corresponding to the original reference fasta, for joining (eg vcf clustering) purposes.
    # Note: perl_generated_fa writes the [first word] of the original fasta to perl_generated_fa header.
    with open(_paths["perl_generated_fa"]) as prg_fasta:
        ref_ID = prg_fasta.readline().strip().replace(">","")

    Prg_Parser = prg_local_parser.Prg_Local_Parser(_paths["prg"], _paths["inferred_fasta"],
                                                   "{} personalised reference from gramtools `infer`".format(ref_ID),
                                                   allele_indexes)
    ref_size = Prg_Parser.parse()

    with open(_paths["inferred_ref_size"], "w") as f:
        f.write(str(ref_size))

    log.info('End process: infer')


## From the vcf describing the original prg, write a new vcf containing the ALT allele inferred at each variant site.
def _dump_vcf(genotypers, _paths, args):
    if _paths['perl_generated_vcf'] is '' or _paths['perl_generated_vcf'] is None:
        log.error('VCF file must be used as part of gramtools build command. Exiting.')
    
    file_handle = open(_paths['perl_generated_vcf'], 'r')
    vcf_reader = vcf.Reader(file_handle)

    # Genotyping information we will output
    format_keys = ['GT', 'DP', 'COV', 'GT_CONF']  # Note of caution: PyVcf (and the VCF format spec) requires 'GT' to be the first field.
    FORMAT = ":".join(format_keys)

    # These namedtuples are defined here for pyVCF operations
    CallData = collections.namedtuple('CallData',format_keys)
    Format = collections.namedtuple('Format','id num type desc')

    ## Update Metadata headers ##
    headers = collections.OrderedDict({
        'fileformat':'VCF4.2', # Some metadata lines, including fileformat, are expected to be 'singular', ie not a list
        'source': ['gramtools, version {}, last commit {}'.format(version.version.version_number,version.version.last_git_commit_hash[:6])],
    })

    if vcf_reader.metadata.get('fileformat',None) is None:
        vcf_reader.metadata['fileformat'] = headers['fileformat']
    vcf_reader.metadata['source'] = headers['source']

    ## Update Format headers ##
    formats = collections.OrderedDict({
        'GT' : Format(id='GT', num='1', type='String', desc='Genotype'),
        'DP' : Format(id='DP', num='1', type='Integer', desc='Total read depth on variant site, from gramtools'),
        'COV' : Format(id='COV', num='R', type='Integer', desc='Number of reads on ref and alt allele(s)'),
        'GT_CONF' : Format(id='GT_CONF', num='1', type='Float', desc='Genotype confidence. Difference in log likelihood of most likely and next most likely genotype'),
    })

    vcf_reader.formats = formats

    ## Name the sample ##
    sample_name = os.path.basename(os.path.normpath(_paths['run_dir']))
    vcf_reader.samples = [sample_name]

    # Build an updater object. It will update the vcf records differently based on infer mode specified at command-line.
    VcfUpdater = _VcfUpdater.factory(args.vcf_mode, CallData, FORMAT, sample_name)

    with open(_paths["inferred_vcf"], 'w') as file_handle:
        vcf_writer = vcf.Writer(file_handle, vcf_reader) # Takes the reader as template

        for vcf_record, gtyper in zip(vcf_reader, genotypers):

            VcfUpdater.update_vcf_record(vcf_record, gtyper)

            vcf_writer.write_record(vcf_record)

    return VcfUpdater.allele_indexes



class _VcfUpdater(object):

    def __init__(self,CallData, FORMAT, sample_name):
        self.CallData = CallData
        self.FORMAT = FORMAT
        self.sample_name = sample_name
        self.allele_indexes = []

    @staticmethod
    def factory(type, CallData, FORMAT, sample_name):
        if type=="single":
            return _SingleUpdater(CallData, FORMAT, sample_name)
        elif type=="population":
            return _PopulationUpdater(CallData, FORMAT, sample_name)


class _SingleUpdater(_VcfUpdater):
    def __init__(self,CallData, FORMAT, sample_name):
        _VcfUpdater.__init__(self, CallData,FORMAT,sample_name)

    ## Computes and records depth, genotype, coverage, genotype confidence.
    # Also removes non-zero coverage alleles if site has been genotyped.
    def update_vcf_record(self, vcf_record, gtyper):
        most_likely_single_allele = _extract_allele_index(gtyper)

        #Case: we had no coverage anywhere in the variant site. We will keep all REF and ALTs in the vcf output.
        if '.' in gtyper.genotype:
            depth, gt_conf = '0', '0.0'
            cov_string = ','.join(['0' for x in range(1 + len(vcf_record.ALT))])
            genotype = './.'

        #Case: we have a genotype. Remove all non-zero coverage alleles in output vcf.
        else:
            # Mark which alleles have non-zero coverage (we are counting reads mapping **only** to a given allele here)
            non_zero_covs = [x for x in range(1 + len(vcf_record.ALT)) if x in gtyper.singleton_alleles_cov]

            #If all coverage is on REF, add a 'dummy' ALT which will get 0 coverage.
            if non_zero_covs == [0]:
                vcf_record.ALT = [vcf_record.ALT[0]]
                non_zero_covs.append(1)

            # Make new ALT from all non-zero coverages, excluding REF (which is at index 0)
            else:
                vcf_record.ALT = [vcf_record.ALT[x - 1] for x in non_zero_covs if x > 0]


            genotype_indexes = sorted(list(gtyper.genotype)) # The indexes of called genotype

            # Add called indexes in case zero coverage allele got called (TODO: is this possible?)
            relevant_covs = sorted(list(set(non_zero_covs + genotype_indexes)))

            if 0 not in gtyper.singleton_alleles_cov: # If no coverage on REF, provide REF index (0)- REF will get 0 cov.
                relevant_covs = [0] + relevant_covs

            cov_values = [gtyper.singleton_alleles_cov.get(x,0) for x in relevant_covs]
            cov_string = ','.join([str(x) for x in cov_values])

            depth = str(sum(gtyper.allele_combination_cov.values()))
            gt_conf = str(gtyper.genotype_confidence)

            # We need the GT indexes to refer to the alleles we have selected (ref + all non-zero coverage alts).
            # So we rebase the indexes of the called alleles
            rebased_most_likely_single_allele = None # Will stay None if most likely haploid call does not appear in most likely diploid call.
            for i in range(len(genotype_indexes)):
                gt = genotype_indexes[i]
                # Convert the allele index to a position relative to the alleles that will get written in new vcf record.
                genotype_indexes[i] = relevant_covs.index(gt)
                # Also rebase the most likely single allele
                if most_likely_single_allele == gt:
                    rebased_most_likely_single_allele = genotype_indexes[i]

            # Case: haploid call
            if len(genotype_indexes) == 1:
                genotype_index = genotype_indexes[0]
                genotype = str(genotype_index) + '/' + str(genotype_index)
            # Case: diploid call. (Can generalise to multiploid call)
            else:
                # As a matter of convention, we will output the most likely 'haploid call' allele first.
                # We will not do so,however, if the most likely 'haploid call' does not appear in the most likely diploid call.
                if rebased_most_likely_single_allele is not None:
                    other_genotype_indexes = set(genotype_indexes)
                    other_genotype_indexes.discard(rebased_most_likely_single_allele)
                    genotype_indexes = [rebased_most_likely_single_allele] + sorted(list(other_genotype_indexes))

                genotype = '/'.join([str(x) for x in genotype_indexes])

        sample = vcf.model._Call(site = vcf_record, sample = self.sample_name, data=self.CallData(GT=genotype, DP=depth, COV=cov_string, GT_CONF=gt_conf))
        vcf_record.samples = [sample]  # Removes pre-existing genotype calls if there were any.
        vcf_record.FORMAT = self.FORMAT

        self.allele_indexes.append(most_likely_single_allele)


class _PopulationUpdater(_VcfUpdater):
    def __init__(self,CallData, FORMAT, sample_name):
        _VcfUpdater.__init__(self, CallData,FORMAT,sample_name)

    ## Computes and records depth, genotype, coverage, genotype confidence.
    # Keeps all alleles at all sites even if they have no coverage.
    def update_vcf_record(self, vcf_record, gtyper):
        most_likely_single_allele = _extract_allele_index(gtyper)

        if '.' in gtyper.genotype: #Case: we had no coverage anywhere in the variant site.
            depth, gt_conf = '0', '0.0'
            cov_string = ','.join(['0' for _ in range(1 + len(vcf_record.ALT))])
            genotype = './.'

        else: #Case: we have a genotype. Remove all non-zero coverage alleles in output vcf.
            cov_values = [gtyper.singleton_alleles_cov.get(x,0) for x in range(1 + len(vcf_record.ALT))]
            cov_string = ','.join([str(x) for x in cov_values])

            depth = str(sum(gtyper.allele_combination_cov.values()))
            gt_conf = str(gtyper.genotype_confidence)

            genotype_indexes = sorted(list(gtyper.genotype))

            if len(genotype_indexes) == 1:
                genotype_index = genotype_indexes[0]
                genotype = str(genotype_index) + '/' + str(genotype_index)
            else:
                # Make sure we output the **most likely single** allele first
                other_genotype_indexes = gtyper.genotype.copy()
                other_genotype_indexes.discard(most_likely_single_allele)

                ordered_genotyped_indexes = [most_likely_single_allele] + sorted(list(other_genotype_indexes))
                genotype = '/'.join([str(x) for x in ordered_genotyped_indexes])

        sample = vcf.model._Call(site = vcf_record, sample = self.sample_name, data=self.CallData(GT=genotype, DP=depth, COV=cov_string, GT_CONF=gt_conf))
        vcf_record.samples = [sample]  # Removes pre-existing genotype calls if there were any.
        vcf_record.FORMAT = self.FORMAT

        self.allele_indexes.append(most_likely_single_allele)


## Generate max. likelihood call for each allele.
def _max_likelihood_alleles(mean_depth,
                            error_rate,
                            all_per_base_coverage,
                            all_groups_site_counts,
                            allele_groups):
    sites_coverage = iter(zip(all_per_base_coverage, all_groups_site_counts))
    for per_base_coverage, groups_site_counts in sites_coverage:
        gtyper = genotyper.Genotyper(mean_depth,
                                     error_rate,
                                     groups_site_counts,
                                     per_base_coverage,
                                     allele_groups)
        gtyper.run()

        yield gtyper



## Returns the most likely single allele id
def _extract_allele_index(gtyper):
    if '.' in gtyper.genotype:  # If no genotyping information, return REF allele index (is 0)
        return 0

    allele_index = None

    # "alleles" is allele groups, we ignore groups which consist of more than just one allele id
    for alleles, likelihood in gtyper.likelihoods:
        if len(alleles) != 1:
           continue
        allele_index = list(alleles)[0]
        break
    if allele_index is None:
        raise ValueError('Genotyper likelihoods does not have any valid haploid data')

    return allele_index


def _load_grouped_allele_coverage(fpath):
    with open(fpath, 'r') as fhandle:
        data = json.load(fhandle)
    groups_coverage = data['grouped_allele_counts']

    allele_groups = groups_coverage['allele_groups']
    for key, value in allele_groups.items():
        allele_groups[key] = set(value)

    groups_site_counts = groups_coverage['site_counts']
    return allele_groups, groups_site_counts


def _load_per_base_coverage(fpath):
    with open(fpath, 'r') as fhandle:
        data = json.load(fhandle)
    data = data['allele_base_counts']
    return data


