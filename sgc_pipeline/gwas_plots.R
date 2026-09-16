#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Batch QQ + Manhattan plots for GWAS summary stats (.tsv.gz) using qqman.
# Runs in parallel via parallel::mclapply (forking; works on macOS/Linux).
#
# Usage:
#   Rscript gwas_plots.R <in_dir> <out_dir> [n_cores]
#
# Examples:
#   Rscript gwas_plots.R /path/to/folder /path/to/plots
#   Rscript gwas_plots.R /path/to/folder /path/to/plots 6
#
# Requires: data.table, qqman  (install.packages(c("data.table","qqman")))
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(qqman)
  library(parallel)
})

# --- command-line args: in_dir, out_dir, [n_cores] -------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript gwas_plots.R <in_dir> <out_dir> [n_cores]")
}
in_dir  <- args[1]
out_dir <- args[2]
n_cores <- if (length(args) >= 3) as.integer(args[3]) else max(1L, detectCores() - 1L)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- set these to YOUR files' column names ---------------------------------
col_chr <- "chromosome"
col_bp  <- "position"
col_p   <- "pvalue"
col_snp <- NA                 # no SNP id column; set to a name if you have one
# ---------------------------------------------------------------------------

# Keep each fork single-threaded so data.table's own threading doesn't
# oversubscribe against mclapply. Comment out to let fread multithread.
setDTthreads(1)

gws   <- 5e-8                  # genome-wide significance threshold
files <- list.files(in_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)
message("Found ", length(files), " files; using ", n_cores, " cores.")

process_one <- function(f) {
  name <- sub("\\.tsv\\.gz$", "", basename(f))

  # skip files whose plots already exist (resume an interrupted run)
  qq_png  <- file.path(out_dir, paste0(name, "_qq.png"))
  man_png <- file.path(out_dir, paste0(name, "_manhattan.png"))
  if (file.exists(qq_png) && file.exists(man_png)) {
    message(sprintf("[%s] skip (exists): %s", format(Sys.time(), "%H:%M:%S"), name))
    return(data.table(file = name, status = "skipped_exists",
                      n_snps = NA, lambda_gc = NA, min_p = NA, n_gws = NA))
  }

  d <- tryCatch(fread(f, showProgress = FALSE), error = function(e) NULL)
  if (is.null(d)) {
    return(data.table(file = name, status = "read_error",
                      n_snps = NA, lambda_gc = NA, min_p = NA, n_gws = NA))
  }

  # guard against missing expected columns
  if (!all(c(col_chr, col_bp, col_p) %in% names(d))) {
    return(data.table(file = name, status = "missing_columns",
                      n_snps = NA, lambda_gc = NA, min_p = NA, n_gws = NA))
  }

  setnames(d, col_chr, "CHR")
  setnames(d, col_bp,  "BP")
  setnames(d, col_p,   "P")
  if (!is.na(col_snp) && col_snp %in% names(d)) setnames(d, col_snp, "SNP")

  # normalize chromosome: strip "chr", map X/Y/MT to numbers qqman accepts
  d[, CHR := sub("^chr", "", as.character(CHR))]
  d[, CHR := fifelse(CHR == "X", "23",
              fifelse(CHR == "Y", "24",
              fifelse(CHR %in% c("MT","M"), "25", CHR)))]
  d[, CHR := suppressWarnings(as.integer(CHR))]
  d[, BP  := suppressWarnings(as.integer(BP))]
  d[, P   := suppressWarnings(as.numeric(P))]

  d <- d[is.finite(CHR) & is.finite(BP) & is.finite(P) & P > 0 & P <= 1]
  if (nrow(d) == 0) {
    return(data.table(file = name, status = "no_valid_rows",
                      n_snps = 0, lambda_gc = NA, min_p = NA, n_gws = NA))
  }

  # qqman's manhattan() errors when no SNP column exists (broken fallback in
  # some versions -> tapply length mismatch). Synthesize one from chr:pos.
  if (!("SNP" %in% names(d))) d[, SNP := paste(CHR, BP, sep = ":")]

  # genomic inflation factor (lambda GC)
  chisq  <- qchisq(1 - d$P, df = 1)
  lambda <- median(chisq) / qchisq(0.5, df = 1)

  n_snps <- nrow(d); min_p <- min(d$P); n_gws <- sum(d$P < gws)

  # y-axis top: extend past the strongest point (and the sig line) so the
  # tallest peaks aren't clipped by manhattan()'s default ylim.
  ymax <- max(-log10(min_p), -log10(gws)) * 1.05
  ymax <- ceiling(ymax)

  # QQ plot (lambda in title)
  png(qq_png, width = 1000, height = 1000, res = 150)
  qq(d$P, main = sprintf("%s  (lambda = %.3f)", name, lambda))
  dev.off()

  # Manhattan plot
  png(man_png, width = 1600, height = 800, res = 150)
  if ("SNP" %in% names(d)) {
    manhattan(d, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
              main = name, genomewideline = -log10(gws), ylim = c(0, ymax))
  } else {
    manhattan(d, chr = "CHR", bp = "BP", p = "P",
              main = name, genomewideline = -log10(gws), ylim = c(0, ymax))
  }
  dev.off()

  message(sprintf("[%s] done: %s", format(Sys.time(), "%H:%M:%S"), name))
  data.table(file = name, status = "ok",
             n_snps = n_snps, lambda_gc = round(lambda, 4),
             min_p = min_p, n_gws = n_gws)
}

results <- mclapply(files, process_one, mc.cores = n_cores)

summary_dt <- rbindlist(results, fill = TRUE)
setorder(summary_dt, -lambda_gc)
fwrite(summary_dt, file.path(out_dir, "gwas_summary.csv"))

message("Done. Plots and gwas_summary.csv in: ", out_dir)
