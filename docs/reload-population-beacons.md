# Reloading the 1000G population beacons

Why the first load was unusable, what replaces it, and how to run it. The plan
this follows is
[`superpowers/plans/2026-10-21-population-beacon-rebuild.md`](superpowers/plans/2026-10-21-population-beacon-rebuild.md).

**Status.** The code and its tests are done. Nothing here has been run on
ILIFU or on `afrigend-all-dev`; every command below is untested against the
cluster and the VM. Phase 0 is the gate that comes first.

## What went wrong

The ten population beacons on `afrigend-all-dev` (five African, five European,
chromosomes 20 to 22) each list 6,240,864 records. A beacon answers "does anyone
in this dataset carry this allele?", so a record must exist only where a member
carries the alternate allele. These datasets list every site of the parent panel
instead. Measured on 2026-10-21: 20 of 20 sampled YRI sites are listed in all ten
datasets; about 74% of YRI records and 83% of CEU records have an alternate
allele count of zero; the API answered yes at 3 of 3 YRI sites with allele
frequency 0. The oracle says "this site is in the panel", so attack power would
be near zero for every population and an African versus European comparison
would measure nothing.

## What the tutorial and upload missed

Read from `afrigen-beacon-v2` at `c5283d2`:

- The transformer (`vcf_to_beacon.py`) reads no genotypes. It writes one record
  per ALT allele of every site that passes QUAL and DP.
- The API's `exists` means "some record matches in the dataset", with no allele
  count guard.
- The ILIFU guide's only slicing recipe is a region cut. Nothing says a sample or
  population subset must drop the sites that do not vary in those samples.
- The right recipe, `bcftools view -S samples.txt -c1`, is in
  `examples/README.md` only. The example loads correctly; the guide that uses it
  never teaches why.
- The guide's and tutorial's checks are positive-only: a known variant answers
  yes, an empty region answers no. Nothing asks for an allele nobody carries to
  answer no. Ten datasets with the same variant count, visible in `/api/datasets`,
  is not something the guide tells you to compare.
- Our own `worklist-for-taine.md` (task 2) has the same gap.

Proposals for the beacon repository's owners (not made here): a guide section
for datasets cut from a subset of samples, with `-S`, `-c1` and a negative
control; a transformer warning when many records have `AC = 0`; a cross-dataset
count check; a tutorial step where an in-panel site your samples do not carry
must return false.

## The design

Per population: 40 members, 40 held-out controls, chosen by a seeded split. The
member dataset holds only biallelic SNPs that vary in those 40 people, with allele
counts recomputed for them. Controls never enter the beacon, so each population
has same-population non-members for the attack. Dataset ids are `<POP>_M40_S1`.

## The pieces

| | |
| --- | --- |
| `analysis/splits/make_split.py` | Seeded, order-independent member/control split; refuses to overwrite a split |
| `analysis/lrt/build_population_dataset.sh` | One population: member VCF, verification lists, het tables, dataset file |
| `analysis/lrt/beacon_check.py` | Verifies a loaded dataset against that output: positive control, negative control, dataset scoping, members and controls |

Scope limits, stated rather than hidden: biallelic SNPs only; chromosomes 20 to
22 only; one split (`S1`); `N = 40`, which is not comparable with the 1,895 of
V6HC-S_AFR (the closed form and the fitted spectrum cover that).

## Run order

**Phase 0 (ILIFU, read-only), done 2026-10-21.** Gate G0 passed: the YRI source
is the NYGC 30x high-coverage call set (`CCDG_13607_B01_GRM_WGS_2019-02-19`, on
the GRCh38 full analysis set reference), not a lift-over. Its header records the
cut as `view -v snps --force-samples -S extract_YRI.txt` with no `-c1`, which
confirms the cause. Two properties of that source shaped the build script: it
carries sites that fail VQSR (FILTER is not PASS), and about 3% of sites have
missing genotypes (allele number below 2N; 96.7% of the first 400,000 loaded YRI
records had AN = 216 and the rest ran down to 42). The command that read the
header, for the record (if a source ever turns out to be the IGSR GRCh38
re-mapping rather than the NYGC 30x set the project chose, stop and report):

```bash
$BCFTOOLS view -h /cbio/users/arctai001/vcf_samples/snps_YRI.vcf.gz | grep -E '^##(source|reference|fileDate|bcftools_viewCommand)'
```

**Split (your machine).** One sample list per population (`bcftools query -l`),
then:

```bash
python3 analysis/splits/make_split.py --population YRI --samples data/YRI_samples.txt --members 40 --controls 40 --seed 1 --out-dir analysis/splits
```

**Build (ILIFU compute node).** Copy the lists and the script across, then one
job per population:

```bash
sbatch --account=b56-cbio-001-ag --job-name=popbuild --array=0-9 --cpus-per-task=2 --mem=8G --time=1-0:0:0 --output=/cbio/users/mamana/pop/build-%A_%a.out --wrap 'P=$(echo CEU ESN FIN GBR GWD IBS LWK MSL TSI YRI | cut -d" " -f$((SLURM_ARRAY_TASK_ID+1))); bash /cbio/users/mamana/pop/build_population_dataset.sh $P /cbio/users/arctai001/vcf_samples/snps_$P.vcf.gz /cbio/users/mamana/pop/${P}_members.txt /cbio/users/mamana/pop/${P}_controls.txt /cbio/users/mamana/pop/build/$P'
```

The build keeps only sites whose FILTER is PASS (or missing) and that have a called
genotype for every member and control (`PASS_ONLY=1`, `REQUIRE_COMPLETE=1`, the
defaults), so the allele number is 2N at every site. Switch either off with `0`;
switching off `REQUIRE_COMPLETE` brings back the refusal, because a dataset whose
allele number varies by site has no single N. Everything after the first step reads
one cleaned file (`clean.vcf.gz`), so the member VCF, the verification lists and
the het tables describe the same sites.

`BCFTOOLS` must be exported first (see `analysis/simulation/docs/fitting-the-sfs.md`).
Set `BEACON_REPO` and `BEACON_CONTAINER` to also run the beacon transform; leave
them unset to build the member VCF and the lists only.

**Load and verify (the dev VM).** Load each transformed dataset as the ILIFU
guide's Route A describes, register it with `--samples 40`, flush the cache, then
for each population:

```bash
python3 analysis/lrt/beacon_check.py --base-url http://127.0.0.1:8086 --build-dir <build>/<POP> --out results/tables/beacon_check_<POP>.json
```

Exit code 0 means every control passed. **No power number is quoted until it is
0 for all ten populations.**
