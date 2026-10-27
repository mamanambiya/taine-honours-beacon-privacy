# Rebuilding the twelve population beacons: the steps

> **What this describes (2026-10-27).** The rebuild that is loaded on
> `afrigend-all-dev`: 12 datasets `<POP>_M65_S1`, each 65 members with 20
> held-out controls: the ten 1000 Genomes populations plus pooled `AFR` and
> `EUR`, cut from the NYGC 30x GRCh38 panel (chr 20 to 22). Result: 9,859,803
> variants (the sum of each `member_sites`) and `beacon_check.py` PASS on all
> 12. The pooled sets are drawn only from samples that are not members of a
> population dataset, because the beacon importer keys individuals by sample
> id (beacon issue 130). Every query runs through `localhost` on the VM; no
> public URL is needed.

Every command to run, in order, with what a correct result looks like. The
design and the reasons are in
[`reload-population-beacons.md`](reload-population-beacons.md); this file is
only the steps.

**Status.** These are the commands that were run on 2026-10-27, with the
mistakes of the first draft corrected. Each step says what to expect. If a
result is different, stop and say which step it was; do not carry on.

## Who runs what

| Part | Where | Who |
| --- | --- | --- |
| A. Get the scripts onto ILIFU | your own computer | Taine |
| B. Build the twelve datasets | ILIFU | Taine |
| C. Bring the results to your computer | your own computer | Taine |
| D. Load them into the beacon | `afrigend-all-dev` | whoever has VM access |
| E. Verify | a computer with VM access | whoever has VM access |
| F. The V6HC-S_AFR counts | ILIFU | Taine |

Parts A to C and F need only an ILIFU account. Parts D and E need the SSH alias
`afrigend-all-dev`, which Taine may not have; if so, send the folder from
Part C to the person who does. Work on `/cbio/users/$USER` can fail with `Disk
quota exceeded` when the account is over its quota; use `/scratch3/users/$USER`
for scratch files then.

Conventions: run each command on its own, one at a time. `$IU` is your ILIFU
username; on your own computer set it once per terminal window. On ILIFU,
`$USER` already is your username. Never run the build on the login node.

## Part A. Get the scripts onto ILIFU

Your own computer. Clone the repository (you are a collaborator, so GitHub
will accept your own login):

```bash
git clone https://github.com/mamanambiya/taine-honours-beacon-privacy.git
```

Your own computer:

```bash
cd taine-honours-beacon-privacy
```

Your own computer. Expect `5f6eab9` or a later commit:

```bash
git log -1 --format=%h
```

Your own computer. Set your ILIFU username (replace the name with yours, once
per window):

```bash
IU=your-ilifu-username
```

Your own computer. Make the folder on ILIFU:

```bash
ssh $IU@transfer.ilifu.ac.za "mkdir -p /cbio/users/$IU/pop"
```

Your own computer. Copy the three scripts and the sample lists:

```bash
scp analysis/splits/make_split.py analysis/lrt/build_population_dataset.sh analysis/lrt/beacon_check.py $IU@transfer.ilifu.ac.za:/cbio/users/$IU/pop/
```

```bash
scp analysis/splits/samples/*.txt $IU@transfer.ilifu.ac.za:/cbio/users/$IU/pop/
```

Your own computer. Check what arrived:

```bash
ssh $IU@transfer.ilifu.ac.za "ls -l /cbio/users/$IU/pop/"
```

Expect the three scripts of 4263, 10238 and 11350 bytes, and 14 `.txt` sample
lists. If you cannot clone, open
each file on GitHub, click **Raw**, and paste it into a new file with the same
name on ILIFU; then run `wc -l` on it (226 lines for
`build_population_dataset.sh`) and `bash -n` on it (no output means it parses).

## Part B. Build the twelve datasets (ILIFU)

Log in and get a compute node. Never run any of this on the login node:

```bash
ssh $IU@slurm.ilifu.ac.za
```

ILIFU login node:

```bash
srun --pty --nodes=1 --ntasks-per-node=1 --cpus-per-task=2 --mem=16G --time=1-0:0:0 bash
```

ILIFU compute node. Set the names, once per shell:

```bash
POP=/cbio/users/$USER/pop
```

```bash
SRC=/cbio/users/$USER/vcf_samples
```

```bash
SIF=/cbio/users/$USER/beacon-tools.sif
```

**B1. Cut the one source VCF.** All twelve datasets are built from one file:
the NYGC 30x GRCh38 phased panel, chr 20 to 22, restricted to the 1,004 AFR and
EUR samples. Check the panel is there (expect three files):

```bash
ls -l /cbio/dbs/1kg/b38/nygc_30x/1kGP_high_coverage_Illumina.chr2[012].filtered.SNV_INDEL_SV_phased_panel.vcf.gz
```

Set the name for bcftools (B4 sets it again; it is needed now):

```bash
export BCFTOOLS="singularity exec --bind /cbio,/scratch3 $SIF bcftools"
```

Make the list of all 1,004 samples. Expect `1004`:

```bash
cat $POP/afr_samples.txt $POP/eur_samples.txt | sort -u > $POP/union_samples.txt; wc -l < $POP/union_samples.txt
```

Cut each chromosome, then join them (about 13 minutes in all). Expect three
files, then a joined file of about 830 MB:

```bash
for C in 20 21 22; do $BCFTOOLS view --threads 8 -S $POP/union_samples.txt /cbio/dbs/1kg/b38/nygc_30x/1kGP_high_coverage_Illumina.chr$C.filtered.SNV_INDEL_SV_phased_panel.vcf.gz -Oz -o $POP/part_chr$C.vcf.gz; $BCFTOOLS index -t $POP/part_chr$C.vcf.gz; done
```

```bash
$BCFTOOLS concat $POP/part_chr20.vcf.gz $POP/part_chr21.vcf.gz $POP/part_chr22.vcf.gz -Oz -o $POP/source_afr_eur_chr20-22.vcf.gz
```

```bash
$BCFTOOLS index -t $POP/source_afr_eur_chr20-22.vcf.gz
```

Check. Expect `3713694` records and `1004` samples:

```bash
$BCFTOOLS view -H $POP/source_afr_eur_chr20-22.vcf.gz | wc -l
```

```bash
$BCFTOOLS query -l $POP/source_afr_eur_chr20-22.vcf.gz | wc -l
```

The index warning "older than the data file" on the panel is harmless here.
Delete the three `part_chr*` files afterwards to save space.

**B2. Update the tools image.** It must contain the carrier-check warning
(beacon PR #133):

```bash
cd /cbio/users/$USER
```

```bash
singularity pull --force beacon-tools.sif docker://ghcr.io/afrigen-d/variant-checker-beacon/beacon-tools:latest
```

Check. Expect one line that mentions `--genotype-check-every`:

```bash
singularity exec $SIF python /opt/beacon-tools/vcf_transform/vcf_to_beacon.py --help | grep genotype-check
```

**B3. Get the beacon repository** (it holds the pipeline script the build
calls). If `/cbio/users/$USER/variant-checker-beacon` already exists, run the
second command instead of the first:

```bash
git clone https://github.com/AfriGen-D/variant-checker-beacon.git /cbio/users/$USER/variant-checker-beacon
```

```bash
git -C /cbio/users/$USER/variant-checker-beacon pull
```

Check. Expect `1` or more:

```bash
grep -c skip-import /cbio/users/$USER/variant-checker-beacon/scripts/run_ilifu_pipeline.sh
```

**B4. Say how to run bcftools and where the beacon tools are:**

```bash
export BCFTOOLS="singularity exec --bind /cbio,/scratch3 $SIF bcftools"
```

```bash
export BEACON_REPO=/cbio/users/$USER/variant-checker-beacon
```

```bash
export BEACON_CONTAINER=$SIF
```

**B5. Check the sample lists** (one per population, copied in Part A):

```bash
wc -l $POP/CEU_samples.txt $POP/ESN_samples.txt $POP/FIN_samples.txt $POP/GBR_samples.txt $POP/GWD_samples.txt $POP/IBS_samples.txt $POP/LWK_samples.txt $POP/MSL_samples.txt $POP/TSI_samples.txt $POP/YRI_samples.txt
```

Expect 99, 98, 99, 91, 111, 107, 99, 85, 107, 108 in that order. MSL, with 85,
is the limit: 65 members plus 20 controls is exactly 85.

**B6. Choose 65 members and 20 held-out controls per dataset.** The same seed
always gives the same split. Check `python3` works first (expect a version
number):

```bash
python3 --version
```

The ten populations first:

```bash
for P in CEU ESN FIN GBR GWD IBS LWK MSL TSI YRI; do python3 $POP/make_split.py --population $P --samples $POP/${P}_samples.txt --members 65 --controls 20 --seed 1 --out-dir $POP; done
```

Now the pooled `AFR` and `EUR`. They must not contain anyone who is a member of
a population dataset, so they draw from the `*_disjoint_samples.txt` lists
(176 AFR and 178 EUR samples; these are the pooled lists minus the 325 members
of the five populations on that side):

```bash
for P in AFR EUR; do python3 $POP/make_split.py --population $P --samples $POP/${P}_disjoint_samples.txt --members 65 --controls 20 --seed 1 --out-dir $POP; done
```

Check. Expect `36`:

```bash
ls $POP/*_members.txt $POP/*_controls.txt $POP/*_split.json | wc -l
```

Check there is no overlap. Expect `0` and `0`:

```bash
cat $POP/ESN_members.txt $POP/GWD_members.txt $POP/LWK_members.txt $POP/MSL_members.txt $POP/YRI_members.txt | sort > $POP/afr_popmembers.txt; sort $POP/AFR_members.txt | comm -12 - $POP/afr_popmembers.txt | wc -l
```

```bash
cat $POP/CEU_members.txt $POP/FIN_members.txt $POP/GBR_members.txt $POP/IBS_members.txt $POP/TSI_members.txt | sort > $POP/eur_popmembers.txt; sort $POP/EUR_members.txt | comm -12 - $POP/eur_popmembers.txt | wc -l
```

**B7. Submit the twelve builds**, one job per dataset. All use the one source
file:

```bash
sbatch --account=b56-cbio-001-ag --job-name=popbuild --array=0-11 --cpus-per-task=2 --mem=8G --time=1-0:0:0 --output=/cbio/users/$USER/pop/build-%A_%a.out --wrap 'P=$(echo CEU ESN FIN GBR GWD IBS LWK MSL TSI YRI AFR EUR | cut -d" " -f$((SLURM_ARRAY_TASK_ID+1))); bash /cbio/users/$USER/pop/build_population_dataset.sh $P /cbio/users/$USER/pop/source_afr_eur_chr20-22.vcf.gz /cbio/users/$USER/pop/${P}_members.txt /cbio/users/$USER/pop/${P}_controls.txt /cbio/users/$USER/pop/build/$P'
```

Expect `Submitted batch job <number>`. If it says the account is invalid, list
yours and use it in place of `b56-cbio-001-ag`:

```bash
sacctmgr -nP show assoc user=$USER format=account
```

Watch the jobs (twelve tasks, running or finished):

```bash
squeue -u $USER
```

(The 2026-10-27 run used two interactive compute-node sessions with two builds
at a time instead of Slurm; both give the same files. Each build took 46 to 68
minutes and the twelve took about 3 hours 12 minutes in all.)

**B8. Check the results** when the queue is empty. Each log's last line should
read `<POP>: N sites vary in 65 members -> <POP>_M65_S1`:

```bash
tail -n 2 $POP/build-*.out
```

Expect these values of `member_sites` (all twelve share one cleaned universe of
3,232,835 sites and allele number 130):

| Dataset | member_sites | Dataset | member_sites |
| --- | --- | --- | --- |
| CEU | 610,835 | LWK | 1,038,569 |
| ESN | 1,002,769 | MSL | 1,032,590 |
| FIN | 565,814 | TSI | 609,724 |
| GBR | 584,895 | YRI | 1,017,091 |
| GWD | 1,027,561 | AFR | 1,077,584 |
| IBS | 664,903 | EUR | 627,468 |

The African datasets are above the European ones, and the pooled AFR above any
single African population:

```bash
grep -H member_sites $POP/build/*/build_summary.json
```

If a build stopped with `allele number is not 2N`, `no site varies` or any
other `build_population_dataset:` message, do not fix it by hand: tell Mamana that
line.

**B9. Find the transformed files** (no compression is needed; the load in
Part D reads them with `cat`):

```bash
ls -l $POP/build/YRI/beacon/vcf_output
```

Expect `variants_batch.jsonl`, `individuals.json`, `variant_genotypes.json` and
`transformation_summary.json` under `$POP/build/YRI/beacon/vcf_output`. If the
names or the folder differ, stop and tell Mamana.

**ROLLBACK:** not part of the steps above. It removes one population's build so
it can be run again (the script refuses to write into an existing folder). Use
the population in place of `YRI`:

```bash
rm -rf $POP/build/YRI
```

## Part C. Bring the results to your computer

Your own computer. Do not copy the two big intermediate VCFs. (The 12 folders
are about 8.8 GB; the heterozygous tables that `beacon_check.py` needs are 170
to 240 MB each. If the machine that runs Part E is the VM, stream them there
instead: `ssh ... "tar cf - ..." | ssh afrigend-all-dev "tar xf -"`.)

```bash
rsync -av --exclude 'clean.vcf.gz*' --exclude 'members.vcf.gz*' $IU@transfer.ilifu.ac.za:/cbio/users/$IU/pop/build/ ~/pop-build/
```

Your own computer. Expect twelve folders:

```bash
ls ~/pop-build
```

Your own computer. Keep the split files in the repository, because they are
what makes the analysis reproducible:

```bash
scp $IU@transfer.ilifu.ac.za:/cbio/users/$IU/pop/*_members.txt $IU@transfer.ilifu.ac.za:/cbio/users/$IU/pop/*_controls.txt $IU@transfer.ilifu.ac.za:/cbio/users/$IU/pop/*_split.json analysis/splits/
```

## Part D. Load them into the beacon (afrigend-all-dev)

Needs the SSH alias `afrigend-all-dev`. The new ids are `<POP>_M65_S1`. On
2026-10-27 the earlier 20-member datasets (`<POP>_M20_S1`) were deleted first,
after a backup (`mongodump` of their variants and dataset records to
`~/backup-m20-2026-10-27/`); do that, or keep them, before loading.

Your own computer. Copy the folders to the VM:

```bash
rsync -av ~/pop-build/ afrigend-all-dev:~/pop_build/
```

The VM had 12 GB free. The twelve load files are about 3.8 GB in all, so stream
them straight from ILIFU into MongoDB in D5 rather than keeping copies.

Your own computer. Open a shell on the VM and stay in it for the rest of Part D:

```bash
ssh afrigend-all-dev
```

On the VM. Check there is room (expect several GB free; MongoDB needs it):

```bash
df -h / | tail -1
```

**D1. Start the beacon** if it is stopped (the data volume is kept; never add
`-v` to `down`):

```bash
cd ~/variant-checker-beacon
```

```bash
docker compose -f compose/docker-compose.dev.yml up -d --build mongodb redis beacon-api
```

The API container runs as the user `beacon` and writes its log to `./logs`. If
`beacon-api-boolean` stays `unhealthy` and `docker logs beacon-api-boolean`
shows `PermissionError: ... /app/logs/beacon.log`, that folder is not writable
by it. The folder and the existing `beacon.log` inside it both need to be
writable by it. Fix them, then restart the API:

```bash
chmod -R a+rwX logs
```

```bash
docker compose -f compose/docker-compose.dev.yml restart beacon-api
```

Expect `beacon-dev-mongodb`, `beacon-dev-redis` and `beacon-api-boolean`, each
`(healthy)` after about a minute:

```bash
docker ps --format '{{.Names}}  {{.Status}}' | grep beacon
```

Expect `200`:

```bash
curl -s -o /dev/null -w "%{http_code}\n" http://localhost:8086/api/health
```

**D2. Raise the query limit on this instance only.** The API is bound to the
VM's loopback, so the raised limit is not reachable from outside. Back up the
file first:

```bash
cp .env.boolean .env.boolean.bak-2026-10-21
```

See whether the setting is already there:

```bash
grep -n RATELIMIT .env.boolean
```

If it printed nothing, add it. If it printed a line, edit that line to read the
same value instead:

```bash
echo 'RATELIMIT_QUERY_ENDPOINT=1000000/hour' >> .env.boolean
```

There is a second limit, `RATELIMIT_DEFAULT` (100/hour), for every other
endpoint (`/api/datasets` and the rest). Without raising it too, a plain
`curl http://localhost:8086/api/datasets` after a verification run answers
`{"detail": "Request was throttled."}`. Edit that line the same way:

```bash
sed -i 's|^RATELIMIT_DEFAULT=.*|RATELIMIT_DEFAULT=1000000/hour|' .env.boolean
```

This setting only lifts the limit if the API was built from beacon repo pull
request 134 (or later). Before that fix the query endpoint kept a hardcoded
`50/hour` and answered HTTP 429 after 50 queries whatever this variable said
("Expected available in 3593 seconds"). Check that the file you are building
from has it, and rebuild the image:

```bash
grep -n 'def get_rate' beacon_api/views_boolean.py
```

Expect one line. Then recreate only the API container:

```bash
docker compose -f compose/docker-compose.dev.yml up -d --build --force-recreate --no-deps beacon-api
```

Expect `1000000/hour`:

```bash
docker exec beacon-api-boolean printenv RATELIMIT_QUERY_ENDPOINT
```

**D3. Names used below**, once per shell:

```bash
M=beacon-dev-mongodb
```

```bash
R=beacon-dev-redis
```

```bash
A=beacon-api-boolean
```

**D4. Check nothing is loaded yet under the new ids.** Loading a dataset twice
doubles it and is hard to undo. Expect twelve lines, all `0`:

```bash
for P in CEU ESN FIN GBR GWD IBS LWK MSL TSI YRI AFR EUR; do echo -n "$P "; docker exec $M mongosh beacon_db --quiet --eval "print(db.variants.countDocuments({dataset_ids: ['${P}_M65_S1']}))"; done
```

**D5. Load each dataset once.** After each, the count in MongoDB must equal that
dataset's `member_sites` from B8. Run it in `tmux` or detached (`setsid nohup`)
so a lost connection does not stop it. The 2026-10-27 run streamed each file
from ILIFU through your computer, with `ServerAliveInterval` so a stalled link
fails instead of hanging, and a 15-minute limit per dataset. Your own computer,
one dataset (replace `CEU`, and `compute-NNN` with your node alias):

```bash
ssh -o ServerAliveInterval=20 -o ServerAliveCountMax=3 ILIFU-compute-NNN_slurm.ilifu.ac.za "cat /cbio/users/$IU/pop/build/CEU/beacon/vcf_output/variants_batch.jsonl" | ssh -o ServerAliveInterval=20 -o ServerAliveCountMax=3 afrigend-all-dev 'docker exec -i beacon-dev-mongodb mongoimport --db beacon_db --collection variants --type json --file /dev/stdin'
```

Expect `N document(s) imported successfully. 0 document(s) failed to import.`
Check the count on the VM (expect the `member_sites` of that dataset, here
610835):

```bash
docker exec $M mongosh beacon_db --quiet --eval "print(db.variants.countDocuments({dataset_ids: ['CEU_M65_S1']}))"
```

Repeat for `ESN`, `FIN`, `GBR`, `GWD`, `IBS`, `LWK`, `MSL`, `TSI`, `YRI`, `AFR`
and `EUR`. **Never run the same dataset twice.** If an import stops partway
(the count is non-zero but below `member_sites`, and does not grow), delete that
dataset's partial documents before trying again, and only those:

```bash
docker exec $M mongosh beacon_db --quiet --eval "printjson(db.variants.deleteMany({dataset_ids: ['GWD_M65_S1']}))"
```

When all twelve are in, the total must be the sum of the twelve `member_sites`.
Expect `9859803`:

```bash
docker exec $M mongosh beacon_db --quiet --eval "print(db.variants.countDocuments({}))"
```

**D6. Register each dataset with its sample count.** Define the names once:

```bash
DATA=$HOME/pop_build; TOOLS=/opt/beacon-tools; MONGO_URI="mongodb://localhost:27017/"; IMAGE=ghcr.io/afrigen-d/variant-checker-beacon/beacon-tools:latest
```

```bash
beacon-docker() { docker run --rm --network host -u "$(id -u):$(id -g)" -e HOME=/tmp -v "$DATA:$DATA" -w "$DATA" --entrypoint= $IMAGE "$@"; }
```

```bash
for P in CEU ESN FIN GBR GWD IBS LWK MSL TSI YRI AFR EUR; do beacon-docker python $TOOLS/data_import/register_dataset.py $DATA/$P/dataset.yaml --db beacon_db --mongo-uri "$MONGO_URI" --samples 65; done
```

Each `dataset.yaml` (and `build_summary.json`) must be on the VM under
`$DATA/<POP>/`; copy those two small files with `tar` if they are not. Expect,
for each, `Registered dataset '<POP>_M65_S1' (...): N variants, 65 samples`, with
N equal to that dataset's `member_sites`.

**D7. Rebuild the indexes** (safe to repeat) and clear the cache:

```bash
docker exec $A python manage.py create_indexes --settings=beacon_project.settings_boolean
```

Expect a last line starting `Done.`.

The API serves `/api/datasets` from a page cache in Redis, so a list read after
a change can show the old datasets. Clear the cache. `FLUSHDB` also removes the
rate-limit counters, which is harmless here; deleting only the cache keys is the
narrower choice, and what the 2026-10-27 run did. Expect `OK`:

```bash
docker exec $R redis-cli FLUSHDB
```

**D8. Check the list.** Expect the twelve new ids, each with `sampleCount` 65
and `variantCount` equal to its `member_sites` (and none of the old
`<POP>_M20_S1` ids if you deleted them). A dataset record's key is `_id`, not
`id`, if you query it in MongoDB:

```bash
curl -s http://localhost:8086/api/datasets | python3 -m json.tool
```

**ROLLBACK:** not part of the steps above. To remove one new dataset (use its
id; this touches only that dataset) and then clear the cache:

```bash
docker exec $M mongosh beacon_db --quiet --eval "printjson(db.variants.deleteMany({dataset_ids: ['YRI_M65_S1']}))"
```

```bash
docker exec $M mongosh beacon_db --quiet --eval "printjson(db.datasets.deleteOne({_id: 'YRI_M65_S1'}))"
```

```bash
docker exec $R redis-cli FLUSHDB
```

To put the query limit back, restore the file and recreate the API:

```bash
cp .env.boolean.bak-2026-10-21 .env.boolean
```

```bash
docker compose -f compose/docker-compose.dev.yml up -d --force-recreate --no-deps beacon-api
```

## Part E. Verify against the data

Run the check on the VM, against `localhost`. No public URL and no tunnel are
needed. The check reads each dataset's ground truth (`present_sites.tsv`,
`absent_sites.tsv`, `members.txt`, `controls.txt` and the `het/` tables, 170 to
240 MB per dataset, about 2.4 GB for twelve). Put them under
`~/pop65_build/<POP>/` on the VM, for example by streaming from ILIFU:

```bash
ssh ILIFU-compute-NNN_slurm.ilifu.ac.za "cd /cbio/users/$IU/pop/build/CEU && tar cf - present_sites.tsv absent_sites.tsv members.txt controls.txt het" | ssh afrigend-all-dev "mkdir -p ~/pop65_build/CEU && cd ~/pop65_build/CEU && tar xf -"
```

Repeat for the other eleven. Copy `beacon_check.py` to `~/pop65_build/`. On
the VM, run all twelve detached (about 4 to 5 minutes each):

```bash
cd ~/pop65_build && mkdir -p checks && for P in CEU ESN FIN GBR GWD IBS LWK MSL TSI YRI AFR EUR; do python3 beacon_check.py --base-url http://127.0.0.1:8086 --build-dir ~/pop65_build/$P --dataset ${P}_M65_S1 --out checks/beacon_check_$P.json > checks/beacon_check_$P.log 2>&1; echo "$P exit $? $(tail -1 checks/beacon_check_$P.log)"; done
```

(Wrap it in a script run with `setsid nohup bash script.sh &` so it survives a
closed terminal.) Expect `exit 0` twelve times and, in every log, `PASS (4451
queries)`. The check asks a known site first (it must answer yes), then 100 sites
a member carries (all yes), 100 sites no member carries (all no), a
dataset-scoping control, and the heterozygous sites of every member and control.
**No attack power number is quoted until all twelve exit 0.** The results from
2026-10-27 are in `results/tables/beacon_check_<POP>.json` (git-ignored; the
outcome is also recorded in `data/README.md`).

## Part F. The V6HC-S_AFR allele counts (ILIFU)

This feeds the fitted spectrum. Your own computer. Copy the script:

```bash
scp analysis/simulation/sfs_histogram.sh $IU@transfer.ilifu.ac.za:/cbio/users/$IU/pop/
```

ILIFU compute node (as in Part B):

```bash
mkdir -p /cbio/users/$USER/sfs
```

```bash
export BCFTOOLS="singularity exec --bind /cbio,/scratch3 /cbio/users/$USER/beacon-tools.sif bcftools"
```

Try one small slice first. The last line it prints says how many ALT alleles it
kept and how many sites it skipped; `skipped ... AN != 3790` should be 0 or
tiny:

```bash
VIEW_ARGS="-r chr22:16000000-17000000" bash /cbio/users/$USER/pop/sfs_histogram.sh /cbio/dbs/refpanels/h3a_reference_panels/version_6/v6hc_s_african/V6HC-S_AFR/vcfs/V6HC-S_AFR_chr22_all.vcf.gz 1895 /cbio/users/$USER/sfs/slice.tsv
```

Then one job per chromosome:

```bash
sbatch --account=b56-cbio-001-ag --job-name=sfs --array=1-22 --cpus-per-task=2 --mem=8G --time=1-0:0:0 --output=/cbio/users/$USER/sfs/sfs-%A_%a.out --wrap 'bash /cbio/users/$USER/pop/sfs_histogram.sh /cbio/dbs/refpanels/h3a_reference_panels/version_6/v6hc_s_african/V6HC-S_AFR/vcfs/V6HC-S_AFR_chr${SLURM_ARRAY_TASK_ID}_all.vcf.gz 1895 /cbio/users/$USER/sfs/v6hc_s_afr_chr${SLURM_ARRAY_TASK_ID}.tsv'
```

When the queue is empty, join the 22 files. Expect `22` first:

```bash
ls /cbio/users/$USER/sfs/v6hc_s_afr_chr*.tsv | wc -l
```

```bash
cat /cbio/users/$USER/sfs/v6hc_s_afr_chr*.tsv > /cbio/users/$USER/sfs/v6hc_s_afr.tsv
```

Your own computer. Bring the small table back and put it in the repository:

```bash
scp $IU@transfer.ilifu.ac.za:/cbio/users/$IU/sfs/v6hc_s_afr.tsv data/v6hc_s_afr.tsv
```

Then the fit, which needs numpy only:

```bash
python3 analysis/simulation/sfs_fit.py data/v6hc_s_afr.tsv --n-individuals 1895 --label "V6HC-S_AFR" --delta 1e-6 --out results/tables/sfs_fit_v6hc_s_afr.json
```

Expect `converged` true and `at_boundary` false in the output and in the JSON.
If either is not, the fit is a bound, not a measurement: tell Mamana and
paste the output.

**ROLLBACK:** not part of the steps above. It removes only the counts:

```bash
rm -f /cbio/users/$USER/sfs/v6hc_s_afr_chr*.tsv /cbio/users/$USER/sfs/v6hc_s_afr.tsv /cbio/users/$USER/sfs/slice.tsv
```
