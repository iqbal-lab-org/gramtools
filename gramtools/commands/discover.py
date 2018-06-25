import os
import glob
import shutil
import logging

import vcf
import cortex
from Bio import SeqIO, Seq

from .. import paths
from .. import prg

log = logging.getLogger('gramtools')


def parse_args(common_parser, subparsers):
    parser = subparsers.add_parser('discover',
                                   parents=[common_parser])
    parser.add_argument('--gram-directory',
                        help='',
                        type=str,
                        required=True)
    parser.add_argument('--inferred-vcf',
                        help='',
                        type=str,
                        required=True)
    parser.add_argument('--reads',
                        help='',
                        type=str,
                        required=True)
    parser.add_argument('--output-vcf',
                        help='',
                        type=str,
                        required=False)


def run(args):
    log.info('Start process: discover')
    _paths = paths.generate_discover_paths(args)
    _dump_inferred_fasta(args.inferred_vcf, _paths)

    try:
        cortex.calls(_paths['inferred_reference'],
                     args.reads,
                     _paths['cortex_directory'],
                     "sample_name")
    except RuntimeError:
        pass

    cortex_vcf_path = _discover_cortex_vcf_file_path(_paths)
    shutil.copyfile(cortex_vcf_path, _paths['cortex_vcf_path'])
    shutil.rmtree(_paths['cortex_directory'])

    log.info('End process: discover')


def _discover_cortex_vcf_file_path(_paths):
    path_pattern = os.path.join(_paths['cortex_directory'],
                                'cortex_output/vcfs/*_wk_*FINAL*raw.vcf')
    found = list(glob.glob(path_pattern, recursive=True))
    assert len(found) == 1, "Multiple possible output cortex VCF files found"
    return found[0]


def _dump_inferred_fasta(inferred_vcf_file_path, _paths):
    with open(_paths['prg'], 'r') as file_handle:
        prg_seq = file_handle.read()
    reference = _get_reference(prg_seq)

    file_handle = open(inferred_vcf_file_path, 'r')
    vcf_reader = vcf.Reader(file_handle)

    inferred_reference = _get_inferred_reference(reference, vcf_reader)
    inferred_reference = ''.join(inferred_reference)
    record = [
        SeqIO.SeqRecord(
            Seq.Seq(inferred_reference, Seq.IUPAC.unambiguous_dna),
            "",
            description="inferred personal reference generated using gramtools")
    ]

    log.debug("Writing inferred fasta reference:\n%s", _paths['inferred_reference'])
    SeqIO.write(record, _paths['inferred_reference'], "fasta")


def _get_inferred_reference(reference, vcf_reader):
    inferred_reference = []
    record = next(vcf_reader, None)

    start_ref_idx = -1
    end_ref_idx = -1

    for idx, x in enumerate(reference):
        within_replaced_ref = start_ref_idx < idx <= end_ref_idx
        if within_replaced_ref:
            continue

        vcf_finished = record is None
        if vcf_finished:
            inferred_reference.append(x)
            continue

        start_ref_idx = record.POS - 1
        end_ref_idx = start_ref_idx + len(record.REF) - 1

        at_allele_insert_idx = idx == start_ref_idx
        if at_allele_insert_idx:
            inferred_reference += list(str(record.ALT[0]))
            record = next(vcf_reader, None)
            continue

        inferred_reference.append(x)
    return inferred_reference


def _get_reference(prg_seq):
    """Get reference by parsing a PRG sequence."""
    log.debug('Parsing PRG for reference')
    regions = prg.parse(prg_seq)

    reference = []
    for region in regions:
        region = region.alleles[0]
        reference += list(region)
    return reference
