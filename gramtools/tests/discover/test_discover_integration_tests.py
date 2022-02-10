from tempfile import mkdtemp
from unittest import mock, TestCase
import shutil
from pathlib import Path

from pysam import VariantFile

from gramtools.commands.discover import discover
import gramtools.tests as gram_testing


class DiscoverRunner:
    def __init__(self, base_dir: Path):
        self.disco_dir = Path(mkdtemp())

        mock_disco_paths = mock.patch(
            "gramtools.commands.discover.discover.DiscoverPaths",
            spec=True,
            reads_files=list(base_dir.glob("*.fq.gz")),
            pers_ref=base_dir / "pers_ref.fa",
            geno_vcf=base_dir / "geno.vcf.gz",
            discov_vcf=self.disco_dir / "cortex.vcf",
            final_vcf=self.disco_dir / "rebased.vcf",
        )
        mock_enforce_ploidy = mock.patch(
            "gramtools.commands.discover.discover.enforce_genotyping_was_haploid"
        )

        mock_disco_paths.start()
        mock_enforce_ploidy.start()

        args = gram_testing.gramtools_main.root_parser.parse_args(
            f"discover -i {base_dir} -o {self.disco_dir} --force".split()
        )
        args.mem_height = 3
        discover.run(args)
        mock_disco_paths.stop()
        mock_enforce_ploidy.stop()

        self.disco_paths = mock_disco_paths.start()

    def __del__(self):
        shutil.rmtree(self.disco_dir)


class TestDiscoveredVariants(TestCase):
    def test_discover_one_variant(self):
        """
        Generated reads header:
            ##ART_Illumina	read_length	50
            @CM	art_illumina -ss HS25 -i pers_ref_with_var.fa -l 50 -f 50 -o r -rs 1589365086
            @SQ	chr1	375
            ##Header End
        The reads are then split in two files to test we can use several
        """
        runner = DiscoverRunner(gram_testing.data_dir / "IT4")
        self.assertEqual(2, len(runner.disco_paths.reads_files))

        cortex_records = list(VariantFile(runner.disco_paths.discov_vcf).fetch())
        self.assertEqual(1, len(cortex_records))
        rec = cortex_records[0]
        self.assertEqual("chr1", rec.chrom)
        self.assertEqual(72, rec.pos)
        self.assertEqual("G", rec.ref)
        self.assertEqual(1, len(rec.alts))
        self.assertEqual("GCCAAACC", rec.alts[0])

        rebased_records = list(VariantFile(runner.disco_paths.final_vcf).fetch())
        self.assertEqual(1, len(rebased_records))
        rebased_rec = rebased_records[0]
        self.assertEqual("chr1", rebased_rec.chrom)
        self.assertEqual(74, rebased_rec.pos)
        self.assertEqual("T", rebased_rec.ref)
        self.assertEqual(1, len(rebased_rec.alts))
        self.assertEqual("GCCAAACC", rebased_rec.alts[0])
