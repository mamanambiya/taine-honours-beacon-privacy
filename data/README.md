# Data

**Nothing in this directory is committed.** `.gitignore` excludes `data/**`
except this file. Record provenance here; keep the bytes out of git.

## Cohorts needed

| Arm | Source | Build | Notes |
| --- | --- | --- | --- |
| African (live beacon) | H3A V6 African Reference Panel, `V6HC-S_AFR` | GRCh38 | 1,895 samples, 42,331,932 variants — the dataset behind `beacon.afrigen-d.org` |
| African (experimental) | 1000 Genomes AFR | GRCh38 | For the size-matched replication |
| European (experimental) | 1000 Genomes EUR / CEU | GRCh38 | The arm Shringarpure & Bustamante measured |

## V6HC-S_AFR on ILIFU

```text
/cbio/dbs/refpanels/h3a_reference_panels/version_6/v6hc_s_african/
    V6HC-S_AFR/vcfs/V6HC-S_AFR_chr{1..22}_all.vcf.gz
```

Autosomes only, GRCh38. This is the same data the production beacon was loaded
from, so an allele-frequency spectrum computed here is the spectrum that governs
the beacon's `D_N`.

**Never SSH to the ILIFU login node.** Work on a compute node — see the
`hpc-clusters` guidance for how to get an interactive session.

## 1000 Genomes build

The beacon answers `501` for any build other than GRCh38, so both arms must be
on GRCh38. The published Phase 1 data that Shringarpure & Bustamante used is
GRCh37, and Phase 3 as originally released is GRCh37 too.

Use the GRCh38-aligned high-coverage (30x) call set, or lift over and document
the liftover — but do not mix builds between arms, because a coordinate
mismatch produces exactly the signal this study is trying to measure.

## A population cut from a larger panel

A VCF cut to a subset of samples keeps every site of the parent panel, including
the ones no sample in the subset carries. Load that into a beacon and it answers
yes at those sites: the dataset then says "this site is in the panel", not "a
member carries this allele", and attack power measures nothing.

Cut the members' VCF with `bcftools view -S members.txt -c1` (restricted to
biallelic SNPs, with `AN`, `AC` and `AF` recomputed for the members), and check
the result: no record with `AC = 0`, and one allele number, `2N`, at every site.
`analysis/lrt/build_population_dataset.sh` does this and refuses to continue
if it cannot show it.

## Recording what you used

For every dataset actually analysed, note here: the source URL or path, the
release/version string, the date retrieved, the sample count after any
filtering, and the checksum if one is published. A power comparison between
cohorts is only interpretable if both arms' provenance is pinned.

## Analysed datasets

### V6HC-S_AFR allele-count histogram (2026-10-23)

| | |
| --- | --- |
| Source | `/cbio/dbs/refpanels/h3a_reference_panels/version_6/v6hc_s_african/V6HC-S_AFR/vcfs/V6HC-S_AFR_chr{1..22}_all.vcf.gz`, GRCh38, autosomes |
| Release | `version_6` (the path's own label; no other release string was read) |
| Samples | 1,895 (`N`; allele number `2N = 3790`) |
| Counted | alternate alleles of biallelic and multi-allelic **SNPs** (`VTYPES=snps`), one entry per ALT allele, counts recomputed from genotypes with `bcftools +fill-tags`; no `-f PASS` filter |
| Result | 42,331,932 ALT alleles in total, 83,379 table rows across the 22 chromosomes (counts repeat per chromosome; `sfs_fit.py` sums them) |
| Tool | `analysis/simulation/sfs_histogram.sh` at `d5920f9`, `bcftools` from `beacon-tools.sif`, one Slurm array task per chromosome (array job 14452681, merge job 14452682) on ILIFU compute nodes |
| Table | `/cbio/users/mamana/sfs/v6hc_s_afr.tsv`, copied to the author's Mac; 667,077 bytes, `sha256 df1febf6173d4c771694636fee30754d4b25a70da6fda6bb1b6a8dac11bcc992` |
| Fit record | `results/tables/sfs_fit_v6hc_s_afr.json` (gitignored; regenerate with the command in `analysis/simulation/docs/fitting-the-sfs.md`) |

**Checked 2026-10-23 (merge job log).** All 22 chromosomes report
`skipped 0 at sites with AN != 3790` and `0 with AC = 0`, so no site was dropped
for missing genotypes and no record without a carrier was counted. The 22 per-chromosome
`kept` counts (chr1 3,346,708 through chr21 560,238) sum to 42,331,932, the same
total the fit read from the merged table. A slice test on
`chr22:16,000,000-17,000,000` kept 9,383 ALT alleles.

**Still not checked.** `bcftools` warned that the `.tbi` index is older than the
VCF it indexes; region reads worked and every chromosome's counts are consistent,
but the warning was not traced to a cause. The 42,331,932 total equals the variant
count in the table at the top of this file; whether the two were counted the same
way (SNP ALT alleles against all records) was not checked.

**What the fit says.** The single-Beta fit sits on its search bound
(`a' = exp(-7)`, `at_boundary = true`), so its parameters and its query budget are
not to be quoted. `D_N` taken straight from the histogram
(`analysis/simulation/empirical_dn.py`, PR #9), target heterozygous, `delta = 1e-6`,
`alpha = 0.05`, 95% power:

| Beacon size | `D_N` | Queries |
| --- | --- | --- |
| 20 | 0.0824 | 31 |
| 65 | 0.0338 | 78 |
| 500 | 0.0067 | 404 |
| 1,000 | 0.0032 | 851 |

It assumes Hardy-Weinberg pairing, independent sites and unrelated panel samples.
