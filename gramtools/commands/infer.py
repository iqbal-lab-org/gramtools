## @file
# From quasimap output (allele coverages), infer a personalised reference.
# Personalised reference can be produced as fasta and/or vcf.
# Reads can subsequently be mapped onto the inferred reference for variant calling: see `discover.py`.

import json
import logging
import sys

import vcf
import collections

from .. import genotyper
from .. import paths
from .. import version
from .. import prg_local_parser

log = logging.getLogger('gramtools')


def parse_args(common_parser, subparsers):
    parser = subparsers.add_parser('infer',
                                   parents=[common_parser])
    parser.add_argument('--gram-dir',
                        help='Directory where the output of `build` is stored',
                        dest='gram_directory',
                        type=str,
                        required=True)

    parser.add_argument('--quasimap-dir',
                        help='Directory where the output of `quasimap` is stored',
                        dest='quasimap_directory',
                        type=str,
                        required=True)

    parser.add_argument('--mean-depth',
                        help='The average read depth of the reads mapped to the prg using `quasimap`',
                        type=int,
                        required=True)

    parser.add_argument('--error-rate',
                        help='The estimated per-base error probability in the reads mapped to the prg using `quasimap`',
                        type=float,
                        required=True)

    parser.add_argument('--haploid',
                        help='',
                        action='store_true',
                        required=False)

    parser.add_argument('--output-fasta',
                        help='Output the inferred personalised reference as a fasta',
                        type=str,
                        required=False)

    parser.add_argument('--output-vcf',
                        help='Output the inferred personalised reference as a vcf.'
                             'If specified, You will need to specify whether to output vcf in "single" or "population" mode',
                        type=str,
                        required=False)

    parser.add_argument('--vcf_mode',
                        help='If "single", the output vcf will not contain alleles with zero-coverage in successfully genotyped sites.'
                             'If "population", all alleles will be represented for all sites.',
                        type=str,
                        required='--output-vcf' in sys.argv,
                        choices=['single','population'],)


def run(args):
    log.info('Start process: infer')
    _paths = paths.generate_infer_paths(args)

    if args.output_fasta is None and args.output_vcf is None:
        log.error("Either --output-fasta or --output-vcf must be specified")
        exit(1)

    all_per_base_coverage = _load_per_base_coverage(_paths['allele_base_coverage'])
    allele_groups, all_groups_site_counts = \
        _load_grouped_allele_coverage(_paths['grouped_allele_counts_coverage'])

    genotypers = _max_likelihood_alleles(args.mean_depth,
                                             args.error_rate,
                                             all_per_base_coverage,
                                             all_groups_site_counts,
                                             allele_groups)

    if args.output_fasta is not None:
        prg_parser = prg_local_parser.parse(_paths['prg'])
        writer = prg_local_parser.FastaWriter(args.output_fasta)
        _dump_fasta(prg_parser, genotypers, writer)
        writer.close()

    if args.output_vcf is not None:
        _dump_vcf(genotypers, _paths, args)

    log.info('End process: infer')


## From the vcf describing the original prg, write a new vcf containing the ALT allele inferred at each variant site.
def _dump_vcf(genotypers, _paths, args):
    if _paths['perl_generated_vcf'] is '' or _paths['perl_generated_vcf'] is None:
        log.error('VCF file must be used as part of gramtools build command. Exiting.')
    
    file_handle = open(_paths['perl_generated_vcf'], 'r')
    vcf_reader = vcf.Reader(file_handle)

    # Get some metadata
    sample_name = vcf_reader.samples[0]

    # Genotyping information we will output
    format_keys = ['GT', 'DP', 'COV', 'GT_CONF']  # Note of caution: PyVcf (and the VCF format spec) requires 'GT' to be the first field.
    FORMAT = ":".join(format_keys)

    # These namedtuples are defined here for pyVCF operations
    CallData = collections.namedtuple('CallData',format_keys)
    Format = collections.namedtuple('Format','id num type desc')

    ## Make vcf header specifications
    headers = collections.OrderedDict({
        'fileformat':'VCF4.2', # Some metadata lines, including fileformat, are expected to be 'singular', ie not a list
        'source': ['gramtools, version {}, last commit {}'.format(version.version.version_number,version.version.last_git_commit_hash)],
    })

    formats = collections.OrderedDict({
        'GT' : Format(id='GT', num='1', type='String', desc='Genotype'),
        'DP' : Format(id='DP', num='1', type='Integer', desc='Total read depth from gramtools'),
        'COV' : Format(id='COV', num='R', type='Integer', desc='Number of reads on ref and alt allele(s)'),
        'GT_CONF' : Format(id='GT_CONF', num='1', type='Float', desc='Genotype confidence. Difference in log likelihood of most likely and next most likely genotype'),
    })

    vcf_reader.metadata = headers
    vcf_reader.formats = formats



    # What we want to capture from the original vcf describing the prg
    record_attributes = [
        'CHROM',
        'POS',
        'ID',
        'REF',
        'ALT',
        'QUAL',
        'FILTER',
        'INFO',
        'FORMAT',
        'sample_indexes',
        'samples',
    ]


    # Build an updater object. It will update the vcf records differently based on infer mode specified at command-line.
    VcfUpdater = _VcfUpdater.factory(args.vcf_mode, CallData, FORMAT, sample_name)

    with open(args.output_vcf, 'w') as file_handle:
        vcf_writer = vcf.Writer(file_handle, vcf_reader) # Takes the reader as template

        for vcf_record, gtyper in zip(vcf_reader, genotypers):
            # TODO: might want to use copy.deepcopy instead of bothering with record_attributes. Or, would pyvcf expose a copy constructor?
            attributes = {name: getattr(vcf_record, name, None) for name in record_attributes}
            new_record = vcf.model._Record(**attributes)

            VcfUpdater.update_vcf_record(new_record, gtyper)

            vcf_writer.write_record(new_record)


class _VcfUpdater(object):

    @staticmethod
    def factory(type, CallData, FORMAT, sample_name):
        if type=="single":
            return _SingleUpdater(CallData, FORMAT, sample_name)
        elif type=="population":
            return _PopulationUpdater(CallData, FORMAT, sample_name)


class _SingleUpdater(_VcfUpdater):
    def __init__(self,CallData, FORMAT, sample_name):
        self.CallData = CallData
        self.FORMAT = FORMAT
        self.sample_name = sample_name

    ## Computes and records depth, genotype, coverage, genotype confidence.
    # Also removes non-zero coverage alleles if site has been genotyped.
    def update_vcf_record(self, new_record, gtyper):
        if '.' in gtyper.genotype: #Case: we had no coverage anywhere in the variant site. We will keep all REF and ALTs in the vcf output.
            depth, gt_conf = '0', '0.0'
            cov_string = ','.join(['0' for x in range(1 + len(new_record.ALT))])
            genotype = './.'

        else: #Case: we have a genotype. Remove all non-zero coverage alleles in output vcf.
            non_zero_covs = [x for x in range(1 + len(new_record.ALT)) if x in gtyper.singleton_alleles_cov] # Mark which alleles have non-zero coverage (we are counting reads mapping **only** to a given allele here)

            if non_zero_covs == [0]:  #If all coverage is on REF, add a 'dummy' ALT which will get 0 coverage.
                new_record.ALT = [new_record.ALT[0]]
                non_zero_covs.append(1)
            else:
                new_record.ALT = [new_record.ALT[x - 1] for x in non_zero_covs if x > 0] # Exclude the REF from consideration here. It is marked as '0' allele.


            genotype_indexes = sorted(list(gtyper.genotype))

            #TODO: is it actually possible that a zero coverage allele gets called?
            relevant_covs = sorted(list(set(non_zero_covs + genotype_indexes)))

            if 0 not in gtyper.singleton_alleles_cov: # Add what will be 0 coverage on the REF allele.
                relevant_covs = [0] + relevant_covs

            cov_values = [gtyper.singleton_alleles_cov.get(x,0) for x in relevant_covs]
            cov_string = ','.join([str(x) for x in cov_values])

            depth = str(sum(gtyper.allele_combination_cov.values()))
            gt_conf = str(gtyper.genotype_confidence)

            # Rebase the indexes so that they refer to the called alleles relative to the selected (non-zero coverage) alleles.
            genotype_pos = 0; index_length = len(genotype_indexes)
            for rebased_index, allele_index in enumerate(relevant_covs):
                if allele_index == genotype_indexes[genotype_pos]:
                    genotype_indexes[genotype_pos] = rebased_index
                    genotype_pos+=1
                    if genotype_pos == index_length:
                        break

            assert genotype_pos == index_length

            if index_length == 1:
                genotype_index = genotype_indexes[0]
                genotype = str(genotype_index) + '/' + str(genotype_index)
            else:
                genotype = '/'.join([str(x) for x in genotype_indexes])

        sample = vcf.model._Call(site = new_record, sample = self.sample_name, data=self.CallData(GT=genotype, DP=depth, COV=cov_string,GT_CONF=gt_conf))
        new_record.samples = [sample]  # Removes pre-existing genotype calls if there were any.
        new_record.FORMAT = self.FORMAT


class _PopulationUpdater(_VcfUpdater):
    def __init__(self,CallData, FORMAT, sample_name):
        self.CallData = CallData
        self.FORMAT = FORMAT
        self.sample_name = sample_name

    ## Computes and records depth, genotype, coverage, genotype confidence.
    # Keeps all alleles at all sites even if they have no coverage.
    def update_vcf_record(self, new_record, gtyper):
        if '.' in gtyper.genotype: #Case: we had no coverage anywhere in the variant site.
            depth, gt_conf = '0', '0.0'
            cov_string = ','.join(['0' for _ in range(1 + len(new_record.ALT))])
            genotype = './.'

        else: #Case: we have a genotype. Remove all non-zero coverage alleles in output vcf.
            cov_values = [gtyper.singleton_alleles_cov.get(x,0) for x in range(1 + len(new_record.ALT))]
            cov_string = ','.join([str(x) for x in cov_values])

            depth = str(sum(gtyper.allele_combination_cov.values()))
            gt_conf = str(gtyper.genotype_confidence)

            genotype_indexes = sorted(list(gtyper.genotype))

            if len(genotype_indexes) == 1:
                genotype_index = genotype_indexes[0]
                genotype = str(genotype_index) + '/' + str(genotype_index)
            else:
                genotype = '/'.join([str(x) for x in genotype_indexes])

        sample = vcf.model._Call(site = new_record, sample = self.sample_name, data=self.CallData(GT=genotype, DP=depth, COV=cov_string,GT_CONF=gt_conf))
        new_record.samples = [sample]  # Removes pre-existing genotype calls if there were any.
        new_record.FORMAT = self.FORMAT


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



def _dump_fasta(prg_parser, genotypers, writer):
    gtyper = next(genotypers)
    allele_index =_extract_allele_index(gtyper)

    for cursor in prg_parser:
        if cursor.just_left_site:
            try:
                allele_index = _extract_allele_index(next(genotypers))
            except StopIteration:
                allele_index = 0

        if cursor.on_marker:
            continue

        if not cursor.within_allele:
            writer.append(cursor.char)
            continue

        if cursor.allele_id == allele_index:
            writer.append(cursor.char)


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


