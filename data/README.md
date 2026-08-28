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

## Recording what you used

For every dataset actually analysed, note here: the source URL or path, the
release/version string, the date retrieved, the sample count after any
filtering, and the checksum if one is published. A power comparison between
cohorts is only interpretable if both arms' provenance is pinned.
