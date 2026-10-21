#!/usr/bin/env bash
# Build one population's member-only beacon dataset, plus the held-out material
# the attack and its verification need. Runs on an ILIFU compute node.
#
# Usage:
#   build_population_dataset.sh POP SOURCE.vcf.gz MEMBERS.txt CONTROLS.txt OUTDIR
#
# Why this exists. A beacon answers "does anyone in this dataset carry this
# allele?". That is only true if every record is a site some member carries. A
# VCF cut to a subset of samples keeps every site of the parent panel, so a
# dataset built from it says "yes" at sites no member carries and the attack has
# nothing to measure. The cut must drop the sites that do not vary in the
# members; this script does that, and refuses to continue if it cannot show it.
#
# Writes, under OUTDIR (which must not exist):
#   members.txt, controls.txt  copies of the two lists, so OUTDIR is self-contained
#   clean.vcf.gz(.tbi)     the source cut to the members and controls, biallelic SNPs,
#                          filtered as above: ONE site universe every later file is
#                          taken from (an intermediate; safe to delete once verified)
#   members.vcf.gz(.tbi)   biallelic SNPs that vary in the members; INFO AC/AN/AF
#                          RECOMPUTED for the members, never inherited
#   present_sites.tsv      up to N_VERIFY sites a member carries   (beacon must say yes)
#   absent_sites.tsv       up to N_VERIFY sites no member carries  (beacon must say no)
#   het/<sample>.tsv       heterozygous biallelic SNPs of every member and control:
#                          the positions an adversary holding that genome would query
#   dataset.yaml           the dataset description, naming the split
#   build_summary.json     counts a later gate can check
#   beacon/                the beacon transform output, when BEACON_REPO is set
#
# Scope, stated rather than hidden: biallelic SNPs only. Multi-allelic sites are
# dropped because the beacon writes one record per ALT allele, and one ALT can
# have no carrier at a site another ALT varies at. Indels are dropped because
# the attack and its spectrum are about SNPs.
#
# Environment:
#   BCFTOOLS          how to run bcftools (default: bcftools). On ILIFU:
#                       BCFTOOLS="singularity exec --bind /cbio,/scratch3 \
#                                 /cbio/users/$USER/beacon-tools.sif bcftools"
#   SEED              seed for the verification samples (default 1). The sample is
#                     reproducible on one machine; awk differs between machines.
#   N_VERIFY          sites per verification list (default 100)
#   SPLIT_TAG         label in the dataset id (default S1): POP_M<members>_S1
#   DATASET_ID        overrides the generated id
#   PASS_ONLY         1 (default) keeps only sites whose FILTER is PASS or missing;
#                     0 keeps every site, including ones that failed a filter.
#   REQUIRE_COMPLETE  1 (default) keeps only sites with a called genotype for EVERY
#                     member and control, so the allele number is 2N at every site;
#                     0 keeps sites with missing genotypes and then refuses to
#                     continue, because the dataset would have no single N.
#                     Both defaults reflect the real 1000G (NYGC 30x) source: it has
#                     sites that fail VQSR, and about 3% of sites with missing calls.
#   BEACON_REPO       a clone of variant-checker-beacon; when set, its
#                     scripts/run_ilifu_pipeline.sh is run with --skip-import
#   BEACON_CONTAINER  the tools .sif or image, passed to that script
set -euo pipefail

die() { echo "build_population_dataset: $*" >&2; exit 3; }

if [ "$#" -ne 5 ]; then
    sed -n '2,/^set -euo/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//' >&2
    exit 2
fi

pop=$1
src=$2
members=$3
controls=$4
out=$5

B=${BCFTOOLS:-bcftools}
seed=${SEED:-1}
n_verify=${N_VERIFY:-100}
tag=${SPLIT_TAG:-S1}
pass_only=${PASS_ONLY:-1}
require_complete=${REQUIRE_COMPLETE:-1}

for flag in PASS_ONLY:"$pass_only" REQUIRE_COMPLETE:"$require_complete"; do
    case ${flag#*:} in
        0|1) ;;
        *) echo "${flag%%:*} must be 0 or 1, got '${flag#*:}'" >&2; exit 2 ;;
    esac
done

[ -e "$out" ] && die "output directory $out already exists; refusing to overwrite a dataset that may have been verified"
[ -s "$members" ] || die "members list $members is missing or empty"
[ -s "$controls" ] || die "controls list $controls is missing or empty"

overlap=$(sort "$members" | comm -12 - <(sort "$controls") | tr '\n' ' ')
[ -z "$overlap" ] || die "sample(s) in both members and controls: $overlap"

n_members=$(grep -c . "$members")
n_controls=$(grep -c . "$controls")
dataset_id=${DATASET_ID:-${pop}_M${n_members}_${tag}}

mkdir -p "$out/het"
all="$out/.all_samples.txt"
cat "$members" "$controls" > "$all"
cp "$members" "$out/members.txt"
cp "$controls" "$out/controls.txt"

# Reservoir sample of n lines from stdin, then ordered by position.
sample_sites() {
    awk -v n="$n_verify" -v seed="$seed" '
        BEGIN { srand(seed) }
        { if (NR <= n) r[NR] = $0
          else { j = int(rand() * NR) + 1; if (j <= n) r[j] = $0 } }
        END { for (i = 1; i <= n && i <= NR; i++) print r[i] }' \
        | sort -k1,1 -k2,2n
}

# 0. One cleaned site universe. Everything below reads this file, never the source,
#    so the member VCF, the verification lists and the het tables describe exactly the
#    same sites. Filters are applied to the members AND controls together: a site with
#    a missing call in a control is dropped too, so no output mixes universes.
pass_args=()
if [ "$pass_only" = 1 ]; then
    pass_args=(-f 'PASS,.')
fi
clean="$out/clean.vcf.gz"
if [ "$require_complete" = 1 ]; then
    # Two steps on purpose: F_MISSING is evaluated on the samples of the file it
    # reads, and `view -S ... -i` would test it BEFORE the sample cut.
    $B view -S "$all" -m2 -M2 -v snps ${pass_args[@]+"${pass_args[@]}"} "$src" -Ou \
        | $B view -i 'F_MISSING=0' -Oz -o "$clean"
else
    $B view -S "$all" -m2 -M2 -v snps ${pass_args[@]+"${pass_args[@]}"} "$src" -Oz -o "$clean"
fi
clean_sites=$($B view -H "$clean" | wc -l | tr -d ' ')
[ "$clean_sites" -gt 0 ] || die "no site of ${pop} survives the filters (PASS_ONLY=${pass_only}, REQUIRE_COMPLETE=${require_complete})"
$B index -t "$clean"

# 1. The member VCF. Counts are recomputed for the members BEFORE the filter, so
#    -c1 acts on the members' own allele counts, not on whatever the parent
#    panel's INFO says.
$B view -S "$members" "$clean" -Ou \
    | $B +fill-tags -Ou -- -t AN,AC,AF \
    | $B view -c1 -Oz -o "$out/members.vcf.gz"

member_sites=$($B view -H "$out/members.vcf.gz" | wc -l | tr -d ' ')
[ "$member_sites" -gt 0 ] || die "no site varies in the ${n_members} members of ${pop}"
$B index -t "$out/members.vcf.gz"

# 2. Gate: one allele number, equal to 2N, at every site; no AC = 0 survives.
an_found=$($B query -f '%AN\n' "$out/members.vcf.gz" | sort -u | tr '\n' ' ')
[ "$an_found" = "$((2 * n_members)) " ] \
    || die "allele number is not 2N = $((2 * n_members)) at every site (found: ${an_found}); the panel has missing genotypes"
zero=$($B query -f '%AC\n' "$out/members.vcf.gz" | awk '$1 == 0' | wc -l | tr -d ' ')
[ "$zero" -eq 0 ] || die "$zero records with AC = 0 survived the filter"

# 3. Verification lists.
$B query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$out/members.vcf.gz" \
    | sample_sites > "$out/present_sites.tsv"

$B view -S "$members" "$clean" -Ou \
    | $B +fill-tags -Ou -- -t AC \
    | $B view -e 'AC>0' -Ou \
    | $B query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
    | sample_sites > "$out/absent_sites.tsv"
[ -s "$out/absent_sites.tsv" ] \
    || echo "warning: no site is absent from every member; the negative control has nothing to test" >&2

# 4. Het tables, one pass over the source for members and controls together. The
#    sample order of the genotype columns is the order of that same stream.
names=$($B view -h "$clean" \
    | awk -F'\t' '/^#CHROM/ { s = $10; for (i = 11; i <= NF; i++) s = s "\t" $i; print s }')
$B view "$clean" \
    | $B query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
    | awk -F'\t' -v OFS='\t' -v dir="$out/het" -v names="$names" '
        BEGIN {
            n = split(names, nm, "\t")
            for (i = 1; i <= n; i++) { f = dir "/" nm[i] ".tsv"; printf "" > f; close(f) }
        }
        {
            for (i = 1; i <= n; i++) {
                g = $(i + 4)
                if (g == "0|1" || g == "1|0" || g == "0/1" || g == "1/0")
                    print $1, $2, $3, $4 >> (dir "/" nm[i] ".tsv")
            }
        }'
rm -f "$all"

# 5. The dataset description.
cat > "$out/dataset.yaml" <<EOF
dataset:
  id: ${dataset_id}
  name: 1000 Genomes ${pop}, ${n_members} members (split ${tag}), GRCh38
  description: >-
    ${n_members} members of ${pop}, chosen by a seeded split with ${n_controls}
    held-out controls. Only biallelic SNPs that vary in these ${n_members}
    people are listed, so a record means a member carries the alternate allele.
    Source: $(basename "$src").
  assembly: GRCh38
  contact:
    name: AfriGen-D Beacon team
    email: beacon@afrigen-d.org
EOF

# 6. The beacon transform, when a clone of the beacon repository is given.
if [ -n "${BEACON_REPO:-}" ]; then
    "$BEACON_REPO/scripts/run_ilifu_pipeline.sh" "$out/members.vcf.gz" \
        --dataset-file "$out/dataset.yaml" \
        ${BEACON_CONTAINER:+--container "$BEACON_CONTAINER"} \
        --skip-import --output "$out/beacon"
else
    echo "skipping the beacon pipeline: BEACON_REPO is not set" >&2
fi

version=$($B --version 2>&1 | sed -n 1p)
cat > "$out/build_summary.json" <<EOF
{
  "population": "${pop}",
  "dataset_id": "${dataset_id}",
  "source_vcf": "${src}",
  "n_members": ${n_members},
  "n_controls": ${n_controls},
  "pass_only": $([ "$pass_only" = 1 ] && echo true || echo false),
  "require_complete_genotypes": $([ "$require_complete" = 1 ] && echo true || echo false),
  "clean_sites": ${clean_sites},
  "member_sites": ${member_sites},
  "allele_number": $((2 * n_members)),
  "verification_sample_size": ${n_verify},
  "seed": ${seed},
  "bcftools": "${version}"
}
EOF
echo "${pop}: ${member_sites} sites vary in ${n_members} members -> ${dataset_id}" >&2
