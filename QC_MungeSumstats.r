# Usage: Rscript gwas_pipeline.R manifest.tsv outdir
library(MungeSumstats)
library(data.table)

args       <- commandArgs(trailingOnly = TRUE)
manifest_f <- args[1]      # one GWAS path per line
outdir     <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh38")
BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")

paths <- data.table::fread(manifest_f, header = FALSE)[[1]]

results <- lapply(seq_along(paths), function(i) {
  gwas_path <- paths[i]
  out_f     <- file.path(outdir,
                         paste0(tools::file_path_sans_ext(basename(gwas_path),
                                                          compression = TRUE),
                                "_munged.tsv.gz"))
  message("[", i, "/", length(paths), "] ", gwas_path)

 tryCatch({
    MungeSumstats::format_sumstats(
      path                   = gwas_path,
      save_path              = out_f,
      INFO_filter            = 0.3,
      FRQ_filter             = 0.005,
      N_dropNA               = FALSE,
      snp_ids_are_rs_ids     = FALSE,
      log_folder_ind         = TRUE,
      log_folder             = file.path(outdir, "logs"),
      log_mungesumstats_msgs = TRUE,
      force_new              = TRUE
    )
  }, error = function(e) {
    message("FAILED: ", gwas_path, " -- ", conditionMessage(e))
    NULL
  })
})

names(results) <- paths
saveRDS(results, file.path(outdir, "munge_results.rds"))