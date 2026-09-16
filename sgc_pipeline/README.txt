SGC meta-analysis pipeline
===================================
Start with run_pipeline.sh: edit the PARAMETERS block at the top, then
  ./run_pipeline.sh list         # see all stages
  ./run_pipeline.sh all          # run everything in order
  ./run_pipeline.sh clean summary # or run individual stages

Expected layout (defaults; set OUT="." for a flat layout):
  working_dir/  scripts + <PHENO>_<STRATUM>.tsv.gz + novelty_info.txt +
                panel_subset.{bed,bim,fam} + plink
  clumped_1mb/           cleaned .verbose.clumped (main content)
  clumped_1mb/unclean/   raw pre-clean copies

SCRIPTS BY STAGE
----------------
runbook            run_pipeline.sh

panel (optional)   harmonize_bim.R  make_keep_list.R
clump              clump_verbose_only.sh
clean              clean_verbose.sh
summary            verbose_summary.sh
minp               min_pvalue.R

overview           build_overview.py
meta_results       build_meta_results.py
clump_index        make_clump_index.py  annotate_novelty.py  clump_expected_gwsig.py
                   (later: add_meta_metrics.py  add_n_cohort.py; rsid/AF manual)
cross              cross_disease_verbose.py  build_cross_disease.py
suggestive         extract_suggestive.R  analyse_strata.R  remove_mhc.R  prune_250kb.R
sex                sex_lookup.R  sex_diff_test.R  annotate_novelty.py
ancestry           ancestry_lookup.R  ancestry_diff_test.R  annotate_novelty.py
recovery           known_loci_recovery.py
workbook           assemble_workbook.py
plots              count_unique_clumps.py  plot_unique_clumps.py
                   plot_pheno_distribution.py  heatmap_trait_locus.py
                   highlights_top_loci.py  (gwas_plots.R optional)

ON DEMAND (per-variant, not in a stage)
  extract_variant.py          pull a variant from any GWAS file(s)

INPUTS NEEDED
  <PHENO>_<STRATUM>.tsv.gz   the GWAS sumstats (cols: chromosome position ref alt beta se pvalue)
  novelty_info.txt           known-loci catalogue (Previous Known loci- variants format)
  panel_subset.{bed,bim,fam} LD reference panel   +   plink (1.9)
