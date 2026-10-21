"""Tests for the cluster-side build of one population's member-only dataset.

Run on a toy VCF small enough to read by hand. Skipped if bcftools is absent.
"""

import json
import os
import shutil
import subprocess
import tempfile
import unittest

HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=1000000>\n"
    '##INFO=<ID=AC,Number=A,Type=Integer,Description="stale">\n'
    '##INFO=<ID=AN,Number=1,Type=Integer,Description="stale">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    "\ts1\ts2\ts3\ts4\ts5\ts6\n"
)

# members: s1 s2 s3     controls: s4 s5 s6
# INFO is deliberately WRONG (AC=99): the build must recompute it, not trust it.
SITES = [
    # carried by a member only
    ("chr1", 100, "A", "G", ["0|1", "0|0", "0|0", "0|0", "0|0", "0|0"]),
    # carried by a control only: must NOT reach the member dataset
    ("chr1", 200, "C", "T", ["0|0", "0|0", "0|0", "1|0", "0|0", "0|0"]),
    # carried by both; s1 is homozygous, s2 and s5 are heterozygous
    ("chr1", 300, "A", "G", ["1|1", "0|1", "0|0", "0|0", "0|1", "0|0"]),
    # multi-allelic: dropped, so an AC=0 ALT allele can never slip through
    ("chr1", 400, "G", "A,C", ["1|2", "0|0", "0|0", "0|0", "0|0", "0|0"]),
    # indel: dropped
    ("chr1", 500, "AC", "A", ["0|1", "0|0", "0|0", "0|0", "0|0", "0|0"]),
    # carried by nobody at all
    ("chr1", 600, "T", "C", ["0|0", "0|0", "0|0", "0|0", "0|0", "0|0"]),
    # one member and one control
    ("chr1", 700, "G", "A", ["0|0", "0|0", "0|1", "0|0", "0|0", "1|0"]),
]


def _vcf(sites=SITES):
    rows = [
        f"{c}\t{p}\t.\t{r}\t{a}\t50\tPASS\tAC=99;AN=99\tGT\t" + "\t".join(g)
        for c, p, r, a, g in sites
    ]
    return HEADER + "\n".join(rows) + "\n"


@unittest.skipUnless(shutil.which("bcftools"), "bcftools not installed")
class TestBuild(unittest.TestCase):
    script = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "build_population_dataset.sh"
    )

    def build(
        self,
        vcf_text=None,
        members=("s1", "s2", "s3"),
        controls=("s4", "s5", "s6"),
        **env,
    ):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        src = os.path.join(self.tmp, "src.vcf")
        with open(src, "w") as fh:
            fh.write(vcf_text if vcf_text is not None else _vcf())
        subprocess.run(["bgzip", "-f", src], check=True)
        subprocess.run(["bcftools", "index", "-t", src + ".gz"], check=True)
        for name, group in (("members", members), ("controls", controls)):
            with open(os.path.join(self.tmp, f"{name}.txt"), "w") as fh:
                fh.write("\n".join(group) + "\n")
        self.out = os.path.join(self.tmp, "out")
        self.result = subprocess.run(
            [
                self.script,
                "YRI",
                src + ".gz",
                os.path.join(self.tmp, "members.txt"),
                os.path.join(self.tmp, "controls.txt"),
                self.out,
            ],
            capture_output=True,
            text=True,
            env={**os.environ, **env},
        )
        return self.result

    def sites(self, name):
        path = os.path.join(self.out, name)
        with open(path) as fh:
            return [tuple(line.split()) for line in fh.read().splitlines()]

    def member_records(self):
        out = subprocess.run(
            [
                "bcftools",
                "query",
                "-f",
                "%POS\t%REF\t%ALT\t%AC\t%AN\n",
                os.path.join(self.out, "members.vcf.gz"),
            ],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        return [tuple(line.split("\t")) for line in out.splitlines()]

    def test_the_member_vcf_holds_only_sites_a_member_carries(self):
        result = self.build()
        self.assertEqual(result.returncode, 0, result.stderr)
        # 200 (control only), 400 (multi-allelic), 500 (indel) and 600
        # (nobody) are all gone.
        self.assertEqual([r[0] for r in self.member_records()], ["100", "300", "700"])

    def test_allele_counts_are_recomputed_for_the_members_not_inherited(self):
        self.build()
        # Every INFO value in the source says AC=99, AN=99.
        self.assertEqual(
            self.member_records(),
            [
                ("100", "A", "G", "1", "6"),
                ("300", "A", "G", "3", "6"),
                ("700", "G", "A", "1", "6"),
            ],
        )

    def test_sites_no_member_carries_are_listed_as_known_absent(self):
        self.build()
        self.assertEqual(
            self.sites("absent_sites.tsv"),
            [("chr1", "200", "C", "T"), ("chr1", "600", "T", "C")],
        )

    def test_sites_a_member_carries_are_listed_as_known_present(self):
        self.build()
        self.assertEqual(
            self.sites("present_sites.tsv"),
            [
                ("chr1", "100", "A", "G"),
                ("chr1", "300", "A", "G"),
                ("chr1", "700", "G", "A"),
            ],
        )

    def test_verification_lists_are_capped_and_reproducible(self):
        self.build(N_VERIFY="2", SEED="7")
        first = (self.sites("present_sites.tsv"), self.sites("absent_sites.tsv"))
        self.build(N_VERIFY="2", SEED="7")
        second = (self.sites("present_sites.tsv"), self.sites("absent_sites.tsv"))
        self.assertEqual(len(first[0]), 2)
        self.assertEqual(first, second)

    def test_every_sample_gets_a_het_table_of_its_own(self):
        self.build()
        for sample in ("s1", "s2", "s3", "s4", "s5", "s6"):
            self.assertTrue(
                os.path.exists(os.path.join(self.out, "het", f"{sample}.tsv")), sample
            )

    def test_het_tables_hold_only_heterozygous_biallelic_snps(self):
        self.build()
        het = {
            s: self.sites(f"het/{s}.tsv") for s in ("s1", "s2", "s3", "s4", "s5", "s6")
        }
        self.assertEqual(het["s1"], [("chr1", "100", "A", "G")])
        # s1 is homozygous at 300; s2 is heterozygous there.
        self.assertEqual(het["s2"], [("chr1", "300", "A", "G")])
        # a control's het site that no member carries is still listed: it is
        # exactly the query that should answer NO against the member beacon.
        self.assertEqual(het["s4"], [("chr1", "200", "C", "T")])
        self.assertEqual(het["s6"], [("chr1", "700", "G", "A")])

    def test_a_sample_with_no_het_sites_still_gets_an_empty_file(self):
        vcf = _vcf(
            [("chr1", 100, "A", "G", ["1|1", "0|1", "0|0", "0|0", "0|0", "0|0"])]
        )
        self.build(vcf)
        self.assertEqual(self.sites("het/s1.tsv"), [])

    def test_copies_the_two_lists_so_the_output_is_self_contained(self):
        self.build()
        for name in ("members.txt", "controls.txt"):
            with open(os.path.join(self.out, name)) as a, \
                    open(os.path.join(self.tmp, name)) as b:
                self.assertEqual(a.read(), b.read())

    def test_writes_a_dataset_file_naming_the_split(self):
        self.build(SPLIT_TAG="S1")
        with open(os.path.join(self.out, "dataset.yaml")) as fh:
            text = fh.read()
        self.assertIn("id: YRI_M3_S1", text)
        self.assertIn("assembly: GRCh38", text)
        self.assertIn("3 members", text)

    def test_writes_a_summary_a_later_gate_can_read(self):
        self.build()
        with open(os.path.join(self.out, "build_summary.json")) as fh:
            summary = json.load(fh)
        self.assertEqual(summary["member_sites"], 3)
        self.assertEqual(summary["n_members"], 3)
        self.assertEqual(summary["n_controls"], 3)
        self.assertEqual(summary["dataset_id"], "YRI_M3_S1")

    def test_skips_the_beacon_pipeline_when_no_repo_is_given(self):
        result = self.build()
        self.assertIn("skipping the beacon pipeline", result.stderr)

    def test_refuses_a_sample_that_is_both_member_and_control(self):
        result = self.build(controls=("s3", "s5", "s6"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("both", result.stderr)

    def test_refuses_missing_genotypes_rather_than_mixing_allele_numbers(self):
        # AN would then differ by site, and a dataset built that way has no
        # single N. A reference panel should be complete: say so, do not guess.
        vcf = _vcf(
            [
                ("chr1", 100, "A", "G", ["0|.", "0|1", "0|0", "0|0", "0|0", "0|0"]),
                ("chr1", 300, "A", "G", ["1|1", "0|1", "0|0", "0|0", "0|1", "0|0"]),
            ]
        )
        result = self.build(vcf)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("allele number", result.stderr)

    def test_refuses_when_no_site_varies_in_the_members(self):
        vcf = _vcf(
            [("chr1", 200, "C", "T", ["0|0", "0|0", "0|0", "1|0", "0|0", "0|0"])]
        )
        result = self.build(vcf)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no site varies", result.stderr)

    def test_refuses_an_existing_output_directory(self):
        self.build()
        again = subprocess.run(
            [
                self.script,
                "YRI",
                os.path.join(self.tmp, "src.vcf.gz"),
                os.path.join(self.tmp, "members.txt"),
                os.path.join(self.tmp, "controls.txt"),
                self.out,
            ],
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(again.returncode, 0)
        self.assertIn("already exists", again.stderr)


if __name__ == "__main__":
    unittest.main()
