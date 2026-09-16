#!/usr/bin/env Rscript
# Greedy distance pruning: within each source_file and chromosome, keep only the
# most significant variant in any 250 kb window. Iteratively pick the lowest-p
# variant, remove everything within +/- 250 kb of it, repeat.
#
# Usage:   Rscript prune_250kb.R <in_csv> <out_csv> [window_kb]
# Example: Rscript prune_250kb.R strata_results_sex_specific_noMHC.csv sex_pruned.csv
#          Rscript prune_250kb.R hits.csv hits_pruned.csv 500   # 500 kb window

suppressPackageStartupMessages(library(data.table))

args    <- commandArgs(trailingOnly = TRUE)
in_csv  <- if (length(args) >= 1) args[1] else stop("need input csv")
out_csv <- if (length(args) >= 2) args[2] else sub("\\.csv$", "_pruned.csv", in_csv)
win     <- (if (length(args) >= 3) as.numeric(args[3]) else 250) * 1000  # kb -> bp

d <- fread(in_csv, showProgress = FALSE)
if (nrow(d) == 0) { fwrite(d, out_csv); cat("empty input; wrote", out_csv, "\n"); quit() }

# greedy pruning within one group (already single source_file + chromosome)
prune_group <- function(dt) {
  dt   <- dt[order(pvalue)]          # strongest first
  pos  <- dt$position
  kept <- logical(nrow(dt))
  taken <- rep(FALSE, nrow(dt))      # removed because near an already-kept lead
  for (i in seq_len(nrow(dt))) {
    if (taken[i]) next               # already within 250kb of a stronger lead
    kept[i] <- TRUE                  # this is a lead SNP
    # mark all not-yet-taken variants within the window of this lead
    near <- abs(pos - pos[i]) <= win & !taken
    taken[near] <- TRUE
    taken[i]    <- TRUE
  }
  dt[kept]
}

# apply per source_file + chromosome
res <- d[, prune_group(.SD), by = .(source_file, chromosome)]

# tidy ordering
setorder(res, source_file, chromosome, pvalue)
fwrite(res, out_csv)

cat(sprintf("Input: %d variants -> kept %d lead variants across %d source_file(s)\n",
            nrow(d), nrow(res), uniqueN(res$source_file)))
cat("Wrote:", out_csv, "\n")
