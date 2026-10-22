# Rebuilding the ten population beacons: the steps

Every command to run, in order, with what a correct result looks like. The
design and the reasons are in
[`reload-population-beacons.md`](reload-population-beacons.md); this file is
only the steps.

**Status.** None of these commands has been run on ILIFU or on
`afrigend-all-dev`. Each step says what to expect. If a result is different,
stop and say which step it was; do not carry on.

## Who runs what

| Part | Where | Who |
| --- | --- | --- |
| A. Get the scripts onto ILIFU | your own computer | Taine |
| B. Build the ten datasets | ILIFU | Taine |
| C. Bring the results to your computer | your own computer | Taine |
| D. Load them into the beacon | `afrigend-all-dev` | whoever has VM access |
| E. Verify | a computer with VM access | whoever has VM access |
| F. The V6HC-S_AFR counts | ILIFU | Taine |

Parts A to C and F need only an ILIFU account. Parts D and E need the SSH alias
`afrigend-all-dev`, which Taine may not have; if so, send the folder from
Part C to the person who does.

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

Your own computer. Copy the three scripts:

```bash
scp analysis/splits/make_split.py analysis/lrt/build_population_dataset.sh analysis/lrt/beacon_check.py $IU@transfer.ilifu.ac.za:/cbio/users/$IU/pop/
```

Your own computer. Check what arrived:

```bash
ssh $IU@transfer.ilifu.ac.za "ls -l /cbio/users/$IU/pop/"
```

Expect three files of 4263, 10238 and 11350 bytes. If you cannot clone, open
each file on GitHub, click **Raw**, and paste it into a new file with the same
name on ILIFU; then run `wc -l` on it (226 lines for
`build_population_dataset.sh`) and `bash -n` on it (no output means it parses).

## Part B. Build the ten datasets (ILIFU)

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

**B1. Check the ten source files exist.** Expect ten lines, one per
population:

```bash
ls -l $SRC/snps_CEU.vcf.gz $SRC/snps_ESN.vcf.gz $SRC/snps_FIN.vcf.gz $SRC/snps_GBR.vcf.gz $SRC/snps_GWD.vcf.gz $SRC/snps_IBS.vcf.gz $SRC/snps_LWK.vcf.gz $SRC/snps_MSL.vcf.gz $SRC/snps_TSI.vcf.gz $SRC/snps_YRI.vcf.gz
```

If they are in another account's folder, set `SRC` to that folder, provided you
can read it.

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

**B5. Make the sample lists** (one per population):

```bash
for P in CEU ESN FIN GBR GWD IBS LWK MSL TSI YRI; do $BCFTOOLS query -l $SRC/snps_$P.vcf.gz > $POP/${P}_samples.txt; done
```

Check. Expect ten lines, each at least 40:

```bash
wc -l $POP/*_samples.txt
```

**B6. Choose 20 members and 20 held-out controls per population.** The same
seed always gives the same split. Check `python3` works first (expect a
version number):

```bash
python3 --version
```

```bash
for P in CEU ESN FIN GBR GWD IBS LWK MSL TSI YRI; do python3 $POP/make_split.py --population $P --samples $POP/${P}_samples.txt --members 20 --controls 20 --seed 1 --out-dir $POP; done
```

Check. Expect `30`:

```bash
ls $POP/*_members.txt $POP/*_controls.txt $POP/*_split.json | wc -l
```

**B7. Submit the ten builds**, one job per population:

```bash
sbatch --account=b56-cbio-001-ag --job-name=popbuild --array=0-9 --cpus-per-task=2 --mem=8G --time=1-0:0:0 --output=/cbio/users/$USER/pop/build-%A_%a.out --wrap 'P=$(echo CEU ESN FIN GBR GWD IBS LWK MSL TSI YRI | cut -d" " -f$((SLURM_ARRAY_TASK_ID+1))); bash /cbio/users/$USER/pop/build_population_dataset.sh $P /cbio/users/$USER/vcf_samples/snps_$P.vcf.gz /cbio/users/$USER/pop/${P}_members.txt /cbio/users/$USER/pop/${P}_controls.txt /cbio/users/$USER/pop/build/$P'
```

Expect `Submitted batch job <number>`. If it says the account is invalid, list
yours and use it in place of `b56-cbio-001-ag`:

```bash
sacctmgr -nP show assoc user=$USER format=account
```

Watch the jobs (ten tasks, running or finished):

```bash
squeue -u $USER
```

**B8. Check the results** when the queue is empty. Each log's last line should
read `<POP>: N sites vary in 20 members -> <POP>_M20_S1`:

```bash
tail -n 2 $POP/build-*.out
```

Expect ten different values of `member_sites`, each far below 6,240,864, with
the African populations above the European ones:

```bash
grep -H member_sites $POP/build/*/build_summary.json
```

If a build stopped with `allele number is not 2N`, `no site varies` or any
other `build_population_dataset:` message, do not fix it by hand: tell Mamana that
line.

**B9. Find the transformed files and compress them** (the load file is big):

```bash
ls -l $POP/build/YRI/beacon/vcf_output
```

Expect `variants_batch.jsonl` and `individuals.json` under
`$POP/build/YRI/beacon/vcf_output`. If the names or the folder differ, stop
and tell Mamana.
Otherwise:

```bash
for P in CEU ESN FIN GBR GWD IBS LWK MSL TSI YRI; do gzip -v $POP/build/$P/beacon/vcf_output/variants_batch.jsonl; done
```

**ROLLBACK:** not part of the steps above. It removes one population's build so
it can be run again (the script refuses to write into an existing folder). Use
the population in place of `YRI`:

```bash
rm -rf $POP/build/YRI
```

## Part C. Bring the results to your computer

Your own computer. Do not copy the two big intermediate VCFs:

```bash
rsync -av --exclude 'clean.vcf.gz*' --exclude 'members.vcf.gz*' $IU@transfer.ilifu.ac.za:/cbio/users/$IU/pop/build/ ~/pop-build/
```

Your own computer. Expect ten folders:

```bash
ls ~/pop-build
```

Your own computer. Keep the split files in the repository, because they are
what makes the analysis reproducible:

```bash
scp $IU@transfer.ilifu.ac.za:/cbio/users/$IU/pop/*_members.txt $IU@transfer.ilifu.ac.za:/cbio/users/$IU/pop/*_controls.txt $IU@transfer.ilifu.ac.za:/cbio/users/$IU/pop/*_split.json analysis/splits/
```

## Part D. Load them into the beacon (afrigend-all-dev)

Needs the SSH alias `afrigend-all-dev`. The existing ten datasets stay where
they are: the new ones have different ids (`<POP>_M20_S1`), and deleting the
old ones is a separate decision.

Your own computer. Copy the folders to the VM:

```bash
rsync -av ~/pop-build/ afrigend-all-dev:~/pop_build/
```

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
doubles it and is hard to undo. Expect ten lines, all `0`:

```bash
for P in CEU ESN FIN GBR GWD IBS LWK MSL TSI YRI; do echo -n "$P "; docker exec $M mongosh beacon_db --quiet --eval "print(db.variants.countDocuments({dataset_ids: ['${P}_M20_S1']}))"; done
```

**D5. Load each population once.** Use `tmux` so a lost connection does not
stop it. Run one population at a time; after each, the number imported must
equal that population's `member_sites` from Part B (see
`~/pop_build/<POP>/build_summary.json`):

```bash
zcat ~/pop_build/YRI/beacon/vcf_output/variants_batch.jsonl.gz | docker exec -i $M mongoimport --db beacon_db --collection variants --type json --file /dev/stdin
```

If that population's file was not compressed (the name has no `.gz`), read it
with `cat` instead of `zcat`:

```bash
cat ~/pop_build/YRI/beacon/vcf_output/variants_batch.jsonl | docker exec -i $M mongoimport --db beacon_db --collection variants --type json --file /dev/stdin
```

Expect `N document(s) imported successfully. 0 document(s) failed to import.`
Repeat with `CEU`, `ESN`, `FIN`, `GBR`, `GWD`, `IBS`, `LWK`, `MSL` and `TSI`
in place of `YRI`. **Never run the same population twice.**

**D6. Register each dataset with its sample count.** Define the names once:

```bash
DATA=$HOME/pop_build; TOOLS=/opt/beacon-tools; MONGO_URI="mongodb://localhost:27017/"; IMAGE=ghcr.io/afrigen-d/variant-checker-beacon/beacon-tools:latest
```

```bash
beacon-docker() { docker run --rm --network host -u "$(id -u):$(id -g)" -e HOME=/tmp -v "$DATA:$DATA" -w "$DATA" --entrypoint= $IMAGE "$@"; }
```

```bash
for P in CEU ESN FIN GBR GWD IBS LWK MSL TSI YRI; do beacon-docker python $TOOLS/data_import/register_dataset.py $DATA/$P/dataset.yaml --db beacon_db --mongo-uri "$MONGO_URI" --samples 20; done
```

Expect, for each, `Registered dataset '<POP>_M20_S1' (...): N variants, 20
samples`, with N equal to that population's `member_sites`.

**D7. Rebuild the indexes** (safe to repeat) and clear the cache:

```bash
docker exec $A python manage.py create_indexes --settings=beacon_project.settings_boolean
```

Expect a last line starting `Done.`:

```bash
docker exec $R redis-cli FLUSHDB
```

Expect `OK`.

**D8. Check the list.** Expect the ten new ids, each with `sampleCount` 20 and
`variantCount` equal to its `member_sites`; the ten old datasets are still
listed:

```bash
curl -s http://localhost:8086/api/datasets | python3 -m json.tool
```

**ROLLBACK:** not part of the steps above. To remove one new dataset (use its
id; this touches only that dataset) and then clear the cache:

```bash
docker exec $M mongosh beacon_db --quiet --eval "printjson(db.variants.deleteMany({dataset_ids: ['YRI_M20_S1']}))"
```

```bash
docker exec $M mongosh beacon_db --quiet --eval "printjson(db.datasets.deleteOne({_id: 'YRI_M20_S1'}))"
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

On a computer with VM access, from the repository folder. Open the tunnel in a
second window and leave it running:

```bash
ssh -N -L 8086:127.0.0.1:8086 afrigend-all-dev
```

First window:

```bash
mkdir -p results/tables
```

```bash
for P in CEU ESN FIN GBR GWD IBS LWK MSL TSI YRI; do python3 analysis/lrt/beacon_check.py --base-url http://127.0.0.1:8086 --build-dir ~/pop-build/$P --per-individual 10 --out results/tables/beacon_check_$P.json; echo "$P exit $?"; done
```

Expect `exit 0` ten times. The check asks a known site first (it must answer
yes), then 100 sites a member carries (all yes), 100 sites no member carries
(all no), a dataset-scoping control, and every member's and control's
heterozygous sites. **No attack power number is quoted until all ten exit 0.**

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
