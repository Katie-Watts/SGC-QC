#!/usr/bin/env bash
# Summarise PLINK 1.9 --clump-verbose output, with MHC/HLA collapsing.
# Two tables:
#   1. <out_prefix>_nclumps.csv -> file, n_clumps
#      = non-MHC clumps + 1 if any clump's index is in the MHC (collapsed once).
#   2. <out_prefix>_clumps.csv  -> file, index_variant, n_gwsig_in_clump
#      = all non-MHC clumps, PLUS a single MHC representative: the MHC clump with
#        the MOST GW-significant variants (ties broken by strongest index p).
#
# MHC/HLA region (hg38): chr6:25,000,000-34,000,000, judged by index position.
# n_gwsig_in_clump = variants in the clump (index + members) with p < SIG.
#
# Usage:   ./verbose_summary.sh <verbose_dir> <out_prefix> [sig_threshold]
# Example: ./verbose_summary.sh clumped_1mb clumped_1mb/verbose_summary 5e-8

set -uo pipefail

DIR="${1:-.}"
OUTPRE="${2:-verbose_summary}"
SIG="${3:-5e-8}"

MHC_CHR=6; MHC_START=25000000; MHC_END=34000000

NCLUMPS="${OUTPRE}_nclumps.csv"
CLUMPS="${OUTPRE}_clumps.csv"
echo "file,n_clumps,n_clumps_ge2_gwsig,n_clumps_ge5_gwsig" > "$NCLUMPS"
echo "file,index_variant,index_pvalue,n_gwsig_in_clump" > "$CLUMPS"

shopt -s nullglob
files=("$DIR"/*.verbose.clumped)
if [[ ${#files[@]} -eq 0 ]]; then
  echo "No *.verbose.clumped files in $DIR" >&2; exit 1
fi

for f in "${files[@]}"; do
  name="$(basename "$f" .verbose.clumped)"
  awk -v NAME="$name" -v SIG="$SIG" \
      -v mc="$MHC_CHR" -v ms="$MHC_START" -v me="$MHC_END" '
    function is_mhc(snp,   a) {
      # snp is chr:pos:ref:alt ; split to get chr and pos
      split(snp, a, ":")
      return (a[1]==mc && a[2]+0>=ms && a[2]+0<=me)
    }
    function flush(   inmhc) {
      if (!have_index) return
      inmhc = is_mhc(idx_snp)
      if (inmhc) {
        mhc_present=1
        # keep the MHC clump with the MOST SIGNIFICANT index SNP (lowest index
        # p); tie-break on the richer clump (more GW-sig variants). The reported
        # n_gwsig is that chosen clump own count.
        if (!mhc_seen || idx_p < mhc_best_p || \
            (idx_p == mhc_best_p && gwsig > mhc_best_gwsig)) {
          mhc_seen=1
          mhc_best_gwsig=gwsig; mhc_best_p=idx_p
          mhc_best_snp=idx_snp; mhc_best_praw=idx_praw
        }
      } else {
        print "CLUMP\t" NAME "," idx_snp "," idx_praw "," gwsig
        nclumps_nonmhc++
        if (gwsig >= 2) ge2_nonmhc++
        if (gwsig >= 5) ge5_nonmhc++
      }
    }
    # clump summary header -> next data line is the summary (index) line
    /^ *CHR +F +SNP +BP +P +TOTAL/ { expect_summary=1; next }
    expect_summary && $1 ~ /^[0-9XY]+$/ {
      flush(); have_index=1; gwsig=0
      idx_snp=$3; idx_praw=$5; idx_p=$5+0
      if (idx_p>0 && idx_p<SIG+0) gwsig++    # index p from summary line
      expect_summary=0; next
    }
    /^ *\(INDEX\)/ { next }                   # redundant with summary line
    have_index && /^ +[0-9XY]+:[0-9]+:/ {
      p=$NF
      if (p ~ /^[0-9.eE+-]+$/ && p+0>0 && p+0<SIG+0) gwsig++
      next
    }
    END {
      flush()
      if (mhc_present) {
        print "CLUMP\t" NAME "," mhc_best_snp "," mhc_best_praw "," mhc_best_gwsig
      }
      total = nclumps_nonmhc + (mhc_present ? 1 : 0)
      ge2 = ge2_nonmhc + ((mhc_present && mhc_best_gwsig >= 2) ? 1 : 0)
      ge5 = ge5_nonmhc + ((mhc_present && mhc_best_gwsig >= 5) ? 1 : 0)
      print "COUNT\t" NAME "," total "," ge2 "," ge5
    }
  ' "$f" | while IFS=$'\t' read -r tag payload; do
      case "$tag" in
        CLUMP) echo "$payload" >> "$CLUMPS" ;;
        COUNT) echo "$payload" >> "$NCLUMPS" ;;
      esac
    done
done

echo "Wrote:"
echo "  $NCLUMPS  (file, n_clumps; MHC collapsed to 1)"
echo "  $CLUMPS   (file, index_variant, n_gwsig_in_clump; one MHC representative)"
