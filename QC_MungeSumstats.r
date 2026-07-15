# Usage: Rscript gwas_pipeline.R path_list.tsv outdir
library(MungeSumstats)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
path_list_f <- args[1]      # one GWAS path per line
outdir <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38")
BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh37")
BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")

paths <- data.table::fread(path_list_f, header = FALSE)[[1]]

results <- lapply(seq_along(paths), function(i) {
  gwas_path <- paths[i]
  out_f <- file.path(outdir,paste0(tools::file_path_sans_ext(basename(gwas_path),
                     compression = TRUE), "_munged.tsv.gz"))
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
   
#add CHR:POS(GRCh38):EFFECT:OTHER as first column
  munged_f <- if (is.list(res)) res$sumstats else res
  dt <- data.table::fread(munged_f)
  dt[, SNP_ID := paste(CHR, BP, A2, A1, sep = ":")] # A2 = effect, A1 = other
  data.table::setcolorder(dt, c("SNP_ID", setdiff(names(dt), "SNP_ID")))
  data.table::fwrite(dt, munged_f, sep = "\t", compress = "gzip")
   
    res
  }, error = function(e) {
    message("FAILED: ", gwas_path, " -- ", conditionMessage(e))
    NULL
  })
})

names(results) <- paths
saveRDS(results, file.path(outdir, "munge_results.rds"))
