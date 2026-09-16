#!/usr/bin/env Rscript
# Remove variants in the MHC/HLA region (chr6, hg38) from the two strata
# results files. Writes new *_noMHC.csv files and reports how many were removed.
#
# Default region: chr6:25,000,000-34,000,000 (25-34 Mb) -- the conservative
# window widely used for MHC exclusion in GWAS (hg19/hg38). Change with args.
#
# Usage:   Rscript remove_mhc.R <file1.csv> [file2.csv ...] [--start=N --end=N]
# Example: Rscript remove_mhc.R strata_results_sex_specific.csv strata_results_ancestry_specific.csv
#
# Alternative regions (hg38):
#   extended MHC (gene-based): --start=25726063 --end=33400644
#   classical MHC (NCBI):      --start=28510120 --end=33480577

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)

# defaults
mhc_chr   <- 6
mhc_start <- 25000000
mhc_end   <- 34000000

# parse optional --start= / --end= / --chr= flags; the rest are file paths
files <- character(0)
for (a in args) {
  if (grepl("^--start=", a))      mhc_start <- as.numeric(sub("^--start=", "", a))
  else if (grepl("^--end=", a))   mhc_end   <- as.numeric(sub("^--end=", "", a))
  else if (grepl("^--chr=", a))   mhc_chr   <- sub("^--chr=", "", a)
  else                            files     <- c(files, a)
}

if (length(files) == 0) {
  files <- c("strata_results_sex_specific.csv",
             "strata_results_ancestry_specific.csv")
}

cat(sprintf("MHC region (hg38): chr%s:%s-%s\n\n",
            mhc_chr, format(mhc_start, scientific = FALSE),
            format(mhc_end, scientific = FALSE)))

for (f in files) {
  if (!file.exists(f)) { cat("SKIP (not found):", f, "\n"); next }

  d <- fread(f, showProgress = FALSE)
  n0 <- nrow(d)
  if (n0 == 0) {
    out <- sub("\\.csv$", "_noMHC.csv", f)
    fwrite(d, out)
    cat(sprintf("%s: empty, wrote %s\n", basename(f), basename(out)))
    next
  }

  # normalize chromosome for comparison (strip any "chr" prefix)
  chr_col <- as.character(d$chromosome)
  chr_col <- sub("^chr", "", chr_col)

  in_mhc <- chr_col == as.character(mhc_chr) &
            d$position >= mhc_start &
            d$position <= mhc_end

  d_keep <- d[!in_mhc]
  out <- sub("\\.csv$", "_noMHC.csv", f)
  fwrite(d_keep, out)

  cat(sprintf("%s: %d rows, removed %d in MHC, kept %d -> %s\n",
              basename(f), n0, sum(in_mhc), nrow(d_keep), basename(out)))
}
