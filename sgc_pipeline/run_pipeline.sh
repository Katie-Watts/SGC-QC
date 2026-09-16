#!/usr/bin/env bash
###############################################################################
#                                                                             #
#   SGC skin-GWAS post-meta-analysis pipeline  —  RUNBOOK                      #
#                                                                             #
#   One place that says WHAT to run, in WHAT ORDER, with WHICH parameters.    #
#   Not a fire-and-forget orchestrator: each STAGE is a function you can run   #
#   on its own, so you can redo one piece with different settings without      #
#   rerunning everything.                                                     #
#                                                                             #
#   USAGE                                                                     #
#     ./run_pipeline.sh <stage> [<stage> ...]     run named stage(s), in order#
#     ./run_pipeline.sh all                        run the whole thing        #
#     ./run_pipeline.sh list                       list stages and exit       #
#                                                                             #
#   Edit the PARAMETERS block below, then call the stage you want, e.g.       #
#     ./run_pipeline.sh clump clean summary                              #
#     ./run_pipeline.sh overview meta_results                                 #
#     ./run_pipeline.sh workbook                                              #
#                                                                             #
#   Every stage just prints and runs the underlying script(s); copy a single  #
#   command out of here and run it by hand any time you want to tweak a flag.  #
#                                                                             #
###############################################################################
set -uo pipefail

# ===========================================================================
#  PARAMETERS  — edit these, they flow into every stage
# ===========================================================================

# --- locations -------------------------------------------------------------
GWAS_DIR="."                       # folder of <PHENO>_<STRATUM>.tsv.gz sumstats
SCRIPTS="."                        # folder holding all the scripts
OUT="results"                      # everything this pipeline writes goes here
PANEL="panel_subset"               # PLINK bfile prefix for the LD reference panel
PLINK="./plink"                    # PLINK 1.9 binary
THREADS=4                          # parallel workers (scale to your machine)

# --- the known-loci catalogue (novelty source) ----------------------------
# Same columns as the "Previous Known loci- variants" sheet:
#   Phenotype  GWAS  Study  Lead variant  Chr  BP  P-value  Clump start  Clump end
NOVELTY="novelty_info.txt"

# (mapping.xlsx is no longer needed — it was only for cohort-level files, which
#  aren't part of this meta-analysis pipeline. Meta files have fixed columns.)

# --- clumping parameters ---------------------------------------------------
#   (These live inside clump_verbose_only.sh as the single source of truth;
#    repeated here only so the runbook documents them. Edit them THERE.)
#   P1=5e-8  P2=1e-5  R2=0.01  KB=500   (500 kb radius = 1 Mb window)

# --- thresholds ------------------------------------------------------------
SIG="5e-8"                         # genome-wide significance
MIN_GWSIG=5                        # "robust locus" = clump with >= this many GW-sig
                                   #   (MHC is always collapsed to ONE clump/locus)

# --- file extension --------------------------------------------------------
SUFFIX=".tsv.gz"                   # sumstat extension

# ===========================================================================
#  derived paths (usually leave alone)
# ===========================================================================
CLUMPDIR="$OUT/clumped_1mb"        # PLINK writes raw *.verbose.clumped here
RESULTS="$CLUMPDIR"                 # cleaned files stay HERE (main content)
UNCLEAN="$CLUMPDIR/unclean"         # raw (pre-clean) copies are kept here
VSUM="$CLUMPDIR/verbose_summary"   # verbose_summary.sh output prefix
mkdir -p "$OUT" "$CLUMPDIR"

say()  { printf '\n\033[1;36m== %s ==\033[0m\n' "$*"; }
run()  { printf '  $ %s\n' "$*"; eval "$@"; }

# ===========================================================================
#  STAGE 0  (OPTIONAL)  reference-panel prep
#  Only needed the FIRST time, or when the panel/GWAS build changes. If your
#  panel_subset.{bed,bim,fam} already exists and matches the GWAS build, SKIP.
# ===========================================================================
stage_panel() {
  say "STAGE 0  reference-panel prep (optional)"
  echo "  Only run if the panel isn't built yet. Uncomment the pieces you need."
  # Rewrite the panel .bim ids to chr:pos:ref:alt in the GWAS orientation so a
  # GWAS-built id matches the panel (uses one representative GWAS):
  # run "Rscript $SCRIPTS/harmonize_bim.R all_hg38_ids.bim ATOPIC_DERM_ALL$SUFFIX all_hg38_matched"
  #
  # Build a keep-list to shrink the panel to variants present in a GWAS:
  # run "Rscript $SCRIPTS/make_keep_list.R ATOPIC_DERM_ALL$SUFFIX keep_ids.txt"
  # run "$PLINK --bfile all_hg38_matched --extract keep_ids.txt --make-bed --out $PANEL"
  echo "  (left commented — edit stage_panel in the runbook if you need it)"
}

# ===========================================================================
#  STAGE 1  clump  (PLINK --clump-verbose on every GWAS)
#  Verbose-only pass. Resumable (skips files already clumped). MHC handled
#  downstream (collapsed to one locus everywhere).
# ===========================================================================
stage_clump() {
  say "STAGE 1  verbose clumping"
  run "$SCRIPTS/clump_verbose_only.sh '$GWAS_DIR' '$PANEL' '$CLUMPDIR' $THREADS '$PLINK'"
}

# ===========================================================================
#  STAGE 2  clean  (keep a raw backup, then clean the files in place)
#  PLINK wrote *.verbose.clumped straight into clumped_1mb/. We copy the raw
#  originals into clumped_1mb/unclean/ (so nothing is lost), then run
#  clean_verbose.sh IN PLACE on clumped_1mb/ -- so the cleaned files remain the
#  main content of clumped_1mb/ and every downstream stage reads them there.
#  (Replaces the old move-then-clean; there is no separate 'move' stage now.)
# ===========================================================================
stage_clean() {
  say "STAGE 2  back up raw copies -> unclean/, then clean in place"
  run "mkdir -p '$UNCLEAN'"
  # copy (not move) the raw files so the cleaned versions stay in clumped_1mb/
  run "cp -n '$CLUMPDIR'/*.verbose.clumped '$UNCLEAN'/ 2>/dev/null || true"
  run "$SCRIPTS/clean_verbose.sh '$CLUMPDIR'"
}

# ===========================================================================
#  STAGE 4  summary  (per-GWAS clump counts; MHC collapsed to 1)
#  Produces  <VSUM>_nclumps.csv  and  <VSUM>_clumps.csv
#  _clumps.csv is the seed for the Clump_index_variants tab.
# ===========================================================================
stage_summary() {
  say "STAGE 4  verbose_summary (clump counts + index list)"
  run "$SCRIPTS/verbose_summary.sh '$RESULTS' '$VSUM' '$SIG'"
}

# ===========================================================================
#  STAGE 5  min_p  (lowest p per GWAS -> GWAS Minimum p tabs)
# ===========================================================================
stage_minp() {
  say "STAGE 5  minimum p-value per GWAS"
  run "Rscript $SCRIPTS/min_pvalue.R '$GWAS_DIR' '$OUT/min_pvalues.csv'"
}

# ===========================================================================
#  WORKBOOK TABS
# ===========================================================================

# --- Overview: phenotype x stratum presence grid ---------------------------
stage_overview() {
  say "TAB  Overview"
  run "python3 $SCRIPTS/build_overview.py --gwas-dir '$GWAS_DIR' --suffix '$SUFFIX' --out '$OUT/overview.tsv'"
}

# --- Meta-Analysis Results: lambda / N / n_sig from sumstats + clump counts -
stage_meta_results() {
  say "TAB  Meta-Analysis Results"
  run "python3 $SCRIPTS/build_meta_results.py --gwas-dir '$GWAS_DIR' --suffix '$SUFFIX' \
        --nclumps '${VSUM}_nclumps.csv' --out '$OUT/meta_results.tsv' --threads $THREADS"
}

# --- Clump_index_variants (.A): index list + novelty + QC metrics ----------
stage_clump_index() {
  say "TAB  Clump_index_variants (.A)"
  # 1) seed from verbose_summary _clumps.csv -> GWAS, index_variant, p, n_gwsig
  run "python3 $SCRIPTS/make_clump_index.py --clumps '${VSUM}_clumps.csv' \
        --out '$OUT/clump_index.tsv'"
  # 2) novelty vs the catalogue (ancestry-aware; GWAS col = PHENO_STRATUM)
  run "python3 $SCRIPTS/annotate_novelty.py --variants '$OUT/clump_index.tsv' \
        --known '$NOVELTY' --pheno-col GWAS --variant-col index_variant \
        --stratum-from-gwas --related --out '$OUT/clump_index_novelty.tsv'"
  # 3) expected-vs-observed QC (needs the verbose files)
  run "python3 $SCRIPTS/clump_expected_gwsig.py --verbose-dir '$RESULTS' \
        --out '$OUT/clump_expected_vs_observed.tsv' --sig '$SIG'"
  echo "  rsid / AF / n_cohorts / i2 / dir-concordance / cases-controls are"
  echo "  added later: add_meta_metrics.py + add_n_cohort.py (need meta files),"
  echo "  then rsid & AF pasted manually."
}

# --- Cross-disease variants -------------------------------------------------
#  Two steps: (1) build the loci + phenotype pairs from the verbose files,
#  (2) add the novelty rollup + direction-difference to make the full sheet.
#  --related-map is an optional TSV (child<TAB>parent) for "novel but related",
#  e.g.  PSOR_VULGARIS <TAB> PSOR .  Leave RELATED_MAP empty to skip it.
RELATED_MAP=""                     # e.g. "$SCRIPTS/related_map.tsv"
stage_cross() {
  say "TAB  Cross-disease variants"
  run "python3 $SCRIPTS/cross_disease_verbose.py '$RESULTS' '$OUT/cross_disease.tsv' '$SIG'"
  local rm_arg=""
  [[ -n "$RELATED_MAP" ]] && rm_arg="--related-map '$RELATED_MAP'"
  run "python3 $SCRIPTS/build_cross_disease.py --cross '$OUT/cross_disease.tsv' \
        --known '$NOVELTY' $rm_arg \
        --gwas-dir '$GWAS_DIR' --suffix '$SUFFIX' \
        --out '$OUT/cross_disease_full.tsv' --threads $THREADS"
}

# --- suggestive-hits track  (feeds sex & ancestry specific tabs) -----------
#  Extract p<1e-5 across all GWAS -> combined file -> per-stratum analyses.
stage_suggestive() {
  say "TAB-prep  suggestive-hits extraction + strata split"
  run "Rscript $SCRIPTS/extract_suggestive.R '$GWAS_DIR' '$OUT/suggestive_hits.csv' $THREADS"
  run "Rscript $SCRIPTS/analyse_strata.R '$OUT/suggestive_hits.csv' '$OUT/strata_results'"
  # drop MHC then distance-prune each stratum set:
  run "Rscript $SCRIPTS/remove_mhc.R '$OUT/strata_results_sex_specific.csv' '$OUT/strata_results_ancestry_specific.csv'"
  run "Rscript $SCRIPTS/prune_250kb.R '$OUT/strata_results_sex_specific_noMHC.csv' '$OUT/sex_pruned.csv'"
  run "Rscript $SCRIPTS/prune_250kb.R '$OUT/strata_results_ancestry_specific_noMHC.csv' '$OUT/ancestry_pruned.csv'"
}

# --- Sex specific loci (.A) -------------------------------------------------
stage_sex() {
  say "TAB  Sex specific loci (.A)"
  run "Rscript $SCRIPTS/sex_lookup.R '$OUT/sex_pruned.csv' '$GWAS_DIR' '$OUT/sex_effect_comparison.csv' $THREADS"
  run "Rscript $SCRIPTS/sex_diff_test.R '$OUT/sex_effect_comparison.csv' '$OUT/sex_diff_test.csv'"
  # novelty for the sex loci (variant col + Phenotype col; no stratum gate):
  run "python3 $SCRIPTS/annotate_novelty.py --variants '$OUT/sex_diff_test.csv' \
        --known '$NOVELTY' --pheno-col Phenotype --chr-col chromosome --pos-col position \
        --related --out '$OUT/sex_specific_novelty.tsv'"
}

# --- Ancestry specific loci (.A) -------------------------------------------
stage_ancestry() {
  say "TAB  Ancestry specific loci (.A)"
  run "Rscript $SCRIPTS/ancestry_lookup.R '$OUT/ancestry_pruned.csv' '$GWAS_DIR' '$OUT/ancestry_effect_comparison.csv' $THREADS"
  run "Rscript $SCRIPTS/ancestry_diff_test.R '$OUT/ancestry_effect_comparison.csv' '$OUT/ancestry_diff_test.csv'"
  run "python3 $SCRIPTS/annotate_novelty.py --variants '$OUT/ancestry_diff_test.csv' \
        --known '$NOVELTY' --pheno-col Phenotype --chr-col chromosome --pos-col position \
        --stratum-col Ancestry --related --out '$OUT/ancestry_specific_novelty.tsv'"
}

# --- Known-loci recovery (Known loci in our GWAS tab) ----------------------
stage_recovery() {
  say "TAB  Known loci recovery"
  run "python3 $SCRIPTS/known_loci_recovery.py --xlsx '$NOVELTY_XLSX' \
        --gwas-dir '$GWAS_DIR' --out-prefix '$OUT/recovery' --mode exact --threads $THREADS"
  echo "  (needs the catalogue as an xlsx sheet; point NOVELTY_XLSX at it)"
}

# ===========================================================================
#  ASSEMBLE THE WORKBOOK  — collate every tab TSV into one .xlsx
# ===========================================================================
stage_workbook() {
  say "ASSEMBLE  final workbook"
  run "python3 $SCRIPTS/assemble_workbook.py --out-dir '$OUT' \
        --novelty '$NOVELTY' --xlsx '$OUT/SGC_Meta_Analysis_Results.xlsx'"
}

# ===========================================================================
#  PLOTS
# ===========================================================================
stage_plots() {
  say "PLOTS"
  WB="$OUT/SGC_Meta_Analysis_Results.xlsx"
  # unique loci per phenotype
  run "python3 $SCRIPTS/count_unique_clumps.py --verbose-dir '$RESULTS' \
        --out '$OUT/clumps_per_phenotype.tsv' --min-gwsig $MIN_GWSIG --collapse-mhc"
  run "python3 $SCRIPTS/plot_unique_clumps.py '$OUT/clumps_per_phenotype.tsv' --out '$OUT/unique_loci.png'"
  # cross-disease distribution (>=5 GW-sig)
  run "python3 $SCRIPTS/plot_pheno_distribution.py '$WB' --min-gwsig $MIN_GWSIG --out '$OUT/pheno_distribution.png'"
  # trait x locus heatmap (>=5), clustered
  run "python3 $SCRIPTS/heatmap_trait_locus.py --verbose-dir '$RESULTS' \
        --out '$OUT/trait_locus_heatmap.png' --min-gwsig $MIN_GWSIG --cap 50"
  # highlights: loci spanning most phenotypes
  run "python3 $SCRIPTS/highlights_top_loci.py --verbose-dir '$RESULTS' \
        --out '$OUT/highlights_top_loci.png' --min-gwsig $MIN_GWSIG --top 30"
  # QQ + Manhattan per GWAS (heavy; comment out if not needed)
  # run "Rscript $SCRIPTS/gwas_plots.R '$GWAS_DIR' '$OUT/plots' $THREADS"
}

# ===========================================================================
#  stage table + dispatcher
# ===========================================================================
ALL_STAGES=(panel clump clean summary minp overview meta_results \
            clump_index cross suggestive sex ancestry recovery workbook plots)

usage() {
  cat <<EOF
Stages (run in this order for a full build):
  panel         (optional) reference-panel prep
  clump         PLINK --clump-verbose on every GWAS
  clean         back up raw copies to clumped_1mb/unclean/, then clean in place
  summary       verbose_summary.sh -> clump counts + index list
  minp          min p-value per GWAS
  overview      Overview tab
  meta_results  Meta-Analysis Results tab
  clump_index   Clump_index_variants (.A): index + novelty + QC
  cross         Cross-disease variants (.A)
  suggestive    extract p<1e-5, split sex/ancestry, MHC-drop, prune
  sex           Sex specific loci (.A)
  ancestry      Ancestry specific loci (.A)
  recovery      Known-loci recovery tab
  workbook      assemble all tabs into the .xlsx
  plots         summary plots

  all           every stage above, in order
  list          print this list

Examples:
  ./run_pipeline.sh clump clean summary
  ./run_pipeline.sh overview meta_results workbook
  ./run_pipeline.sh all
EOF
}

[[ $# -eq 0 ]] && { usage; exit 1; }
if [[ "$1" == "list" || "$1" == "-h" || "$1" == "--help" ]]; then usage; exit 0; fi
if [[ "$1" == "all" ]]; then set -- "${ALL_STAGES[@]}"; fi

for stage in "$@"; do
  fn="stage_$stage"
  if declare -F "$fn" >/dev/null; then
    "$fn"
  else
    echo "Unknown stage: $stage" >&2; usage; exit 2
  fi
done
say "done"
