#!/usr/bin/env Rscript
# Extract the lowest p-value from each .tsv.gz GWAS file in a folder.
# Usage: Rscript min_pvalue.R <in_dir> <out_csv>
# Example: Rscript min_pvalue.R . min_pvalues.csv

suppressPackageStartupMessages(library(data.table))

args    <- commandArgs(trailingOnly = TRUE)
in_dir  <- if (length(args) >= 1) args[1] else "."
out_csv <- if (length(args) >= 2) args[2] else "min_pvalues.csv"

col_p <- "pvalue"   # p-value column name

files <- list.files(in_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)

res <- rbindlist(lapply(files, function(f) {
  name <- sub("\\.tsv\\.gz$", "", basename(f))
  # read only the p-value column for speed
  d <- tryCatch(fread(f, select = col_p, showProgress = FALSE),
                error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(data.table(file = name, min_p = NA_real_))
  p <- suppressWarnings(as.numeric(d[[col_p]]))
  p <- p[is.finite(p) & p > 0 & p <= 1]
  data.table(file = name, min_p = if (length(p)) min(p) else NA_real_)
}), fill = TRUE)

setorder(res, min_p)
fwrite(res, out_csv)
cat("Wrote", out_csv, "with", nrow(res), "rows\n")
