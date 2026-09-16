#!/usr/bin/env bash
# Verbose-ONLY clumping pass: for each GWAS, produce the --clump-verbose output
# (<name>.verbose.clumped) if it's missing. Does NOT run or touch the standard
# clump. Independent skip logic keyed on the verbose output existing, so it fills
# in only the verbose files a half-finished run hasn't made yet.
#
# Uses .verbose.done markers (separate from the standard run's .done) so it can
# also skip no-signal phenotypes (which legitimately produce no .verbose.clumped)
# without re-running them each restart.
#
# Usage:
#   ./clump_verbose_only.sh <gwas_dir> <panel_prefix> <out_dir> [n_parallel] [plink_bin]
# Example:
#   ./clump_verbose_only.sh . panel_subset clumped 2 ./plink

set -uo pipefail

GWAS_DIR="${1:-.}"
PANEL="${2:-panel_subset}"
OUT_DIR="${3:-clumped}"
NPAR="${4:-2}"
PLINK="${5:-./plink}"

# --- clumping parameters (match your standard run) ---
P1=5e-8; P2=1e-5; R2=0.01; KB=500
C_CHR="chromosome"; C_POS="position"; C_REF="ref"; C_ALT="alt"; C_P="pvalue"
# -----------------------------------------------------

mkdir -p "$OUT_DIR" "$OUT_DIR/tmp"

if [[ "$PLINK" == */* ]]; then
  PLINK="$(cd "$(dirname "$PLINK")" && pwd)/$(basename "$PLINK")"
fi
[[ -x "$PLINK" ]] || { echo "ERROR: plink not executable at: $PLINK" >&2; exit 1; }
PANEL="$(cd "$(dirname "$PANEL")" && pwd)/$(basename "$PANEL")"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"

process_one() {
  local f="$1"
  local name; name="$(basename "$f" .tsv.gz)"
  local outv="$OUT_DIR/${name}.verbose"      # verbose output prefix
  local tmp="$OUT_DIR/tmp/${name}.verbose.idp.tsv"

  # skip if the verbose run already completed, judged from its log:
  #   success   -> log contains "Results written to"
  #   no signal -> log contains "No significant --clump results"
  # An interrupted run has a log with neither (or no log) and is (re)run.
  # (.verbose.done markers aren't used here -- the log is the source of truth.)
  if [[ -f "${outv}.log" ]] && \
     grep -qE "Results written to|No significant --clump results" "${outv}.log" 2>/dev/null; then
    echo "[skip] $name (verbose already done)"
    return 0
  fi

  # generate SNP<TAB>P (id = chr:pos:ref:alt)
  gzip -dc "$f" | awk -F'\t' -v OFS='\t' \
    -v cc="$C_CHR" -v cp="$C_POS" -v cr="$C_REF" -v ca="$C_ALT" -v cpv="$C_P" '
    NR==1 {
      for (i=1;i<=NF;i++) h[$i]=i
      ic=h[cc]; ip=h[cp]; ir=h[cr]; ia=h[ca]; iv=h[cpv]
      if (!ic||!ip||!ir||!ia||!iv) { print "MISSING_COLUMN" > "/dev/stderr"; exit 2 }
      print "SNP","P"; next
    }
    { chr=$ic; sub(/^chr/,"",chr)
      print chr":"$ip":"toupper($ir)":"toupper($ia), $iv }' > "$tmp"
  if [[ $? -ne 0 || ! -s "$tmp" ]]; then
    echo "[FAIL] $name (id/p extraction)"; rm -f "$tmp"; return 1
  fi

  # verbose clump only
  "$PLINK" --bfile "$PANEL" --clump "$tmp" \
           --clump-p1 "$P1" --clump-p2 "$P2" --clump-r2 "$R2" --clump-kb "$KB" \
           --clump-snp-field SNP --clump-field P \
           --clump-verbose \
           --out "$outv" >"${outv}.log" 2>&1

  # report outcome from the verbose log (skip on next run is log-based)
  if grep -qE "Results written to|No significant --clump results" "${outv}.log" 2>/dev/null; then
    if [[ -f "${outv}.clumped" ]]; then echo "[done] $name (verbose)"
    else echo "[done] $name (verbose, no clumps)"; fi
  else
    echo "[incomplete] $name: verbose log has no completion line — will retry next run"
    rm -f "${outv}.clumped"
  fi
  rm -f "$tmp"
}
export -f process_one
export OUT_DIR PANEL PLINK P1 P2 R2 KB C_CHR C_POS C_REF C_ALT C_P

find "$GWAS_DIR" -maxdepth 1 -name '*.tsv.gz' ! -name '._*' -print0 \
  | xargs -0 -P "$NPAR" -I{} bash -c 'process_one "$@"' _ {}

echo "All done (verbose). Results in: $OUT_DIR"
