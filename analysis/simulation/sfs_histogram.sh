#!/usr/bin/env bash
# Allele-count histogram of one VCF: the only thing that has to leave the cluster
# for the SFS fit (sfs_fit.py).
#
# Usage:
#   sfs_histogram.sh VCF.gz N_INDIVIDUALS [OUT.tsv]
#
# Writes two tab-separated columns, `alternate-allele count` and `variants`, one
# row per count that occurs, to OUT.tsv (or stdout). Tables for several VCFs can
# simply be concatenated; sfs_fit.read_histogram sums repeated counts.
#
# Counts are per ALT ALLELE, not per site: a multi-allelic site contributes one
# entry for each of its ALT alleles, which is how the beacon loader stores it.
# Counts are recomputed from the genotypes with `bcftools +fill-tags`, never
# read from INFO, because a panel's INFO AC may describe a larger parent cohort.
#
# Sites whose allele number is not 2N (missing genotypes) are skipped and
# counted on stderr. A site with a count above 2N would otherwise be a wrong-N
# error that the fit refuses; skipping and reporting is the honest alternative.
#
# Environment:
#   BCFTOOLS   how to run bcftools (default: bcftools). On ILIFU the tool lives in
#              the container, so give the whole command, for example
#                BCFTOOLS="singularity exec --bind /cbio,/scratch3 \
#                          /cbio/users/$USER/beacon-tools.sif bcftools"
#              (a shell function such as beacon-sing is not visible to a script).
#   VTYPES     variant types to keep: snps (default), indels, snps,indels, any.
#              The 2015 attack and the SFS it assumes are about SNPs; widen it
#              only if the beacon's query sites include indels.
#   VIEW_ARGS  extra `bcftools view` arguments, for example "-f PASS" or
#              "-r chr22:16000000-17000000" for a test slice.
set -euo pipefail

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    sed -n '2,/^set -euo/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//' >&2
    exit 2
fi

vcf=$1
n_individuals=$2
out=${3:-/dev/stdout}
bcftools_cmd=${BCFTOOLS:-bcftools}
vtypes=${VTYPES:-snps}
view_args=${VIEW_ARGS:-}

case $n_individuals in
    ''|*[!0-9]*) echo "N_INDIVIDUALS must be a positive integer, got '$n_individuals'" >&2
                 exit 2 ;;
esac

type_args=()
if [ "$vtypes" != "any" ]; then
    type_args=(-v "$vtypes")
fi

# shellcheck disable=SC2086  # BCFTOOLS and VIEW_ARGS are deliberately word-split
$bcftools_cmd view ${type_args[@]+"${type_args[@]}"} $view_args "$vcf" \
    | $bcftools_cmd +fill-tags -- -t AN,AC \
    | $bcftools_cmd query -f '%AN\t%AC\n' \
    | awk -F'\t' -v two_n=$((2 * n_individuals)) '
        BEGIN { OFS = "\t" }
        {
            if ($1 + 0 != two_n) { wrong_an++; next }
            m = split($2, ac, ",")
            for (i = 1; i <= m; i++) {
                if (ac[i] == "." || ac[i] + 0 == 0) { zero++; continue }
                if (ac[i] + 0 > two_n) { too_high++; continue }
                hist[ac[i] + 0]++
                kept++
            }
        }
        END {
            for (k in hist) print k, hist[k] | "sort -n"
            close("sort -n")
            printf "kept %d ALT alleles; skipped %d at sites with AN != %d, " \
                   "%d with AC = 0\n", kept, wrong_an, two_n, zero > "/dev/stderr"
            if (kept == 0) {
                print "sfs_histogram: no ALT allele was counted, so there is " \
                      "nothing to fit (wrong N, an empty region, or bcftools " \
                      "produced no records; check the messages above)" > "/dev/stderr"
                exit 1
            }
        }' > "$out"
