#!/usr/bin/env Rscript
# For each variant in the ancestry-pruned file, pull beta/se/pvalue from its OWN
# non-EUR stratum and from the EUR file of the same phenotype, so the effect-size
# difference vs Europeans is visible.
#
# Memory-safe: processes ONE phenotype at a time (reads only that phenotype's
# needed non-EUR strata + EUR, does its lookups, discards them), so peak memory
# holds at most n_cores phenotypes' files -- independent of total phenotype count.
#
# Input columns: source_file, chromosome, position, ref, alt (prune_250kb.R out)
#                source_file = <PHENOTYPE>_<ANCESTRY>, e.g. ATOPIC_DERM_AFR.
# Source files:  <src_dir>/<PHENOTYPE>_<ANCESTRY>.tsv.gz and _EUR.tsv.gz
#
# Usage:   Rscript ancestry_lookup.R <variants_csv> <src_dir> <out_csv> [n_cores]
# Example: Rscript ancestry_lookup.R ancestry_pruned.csv . ancestry_effect_comparison.csv 3

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

args    <- commandArgs(trailingOnly = TRUE)
var_csv <- if (length(args) >= 1) args[1] else stop("need variants csv")
src_dir <- if (length(args) >= 2) args[2] else "."
out_csv <- if (length(args) >= 3) args[3] else "ancestry_effect_comparison.csv"
n_cores <- if (length(args) >= 4) as.integer(args[4]) else 3L
setDTthreads(1)

# per-phenotype partial results (resume support): a phenotype with a .done
# marker is skipped on rerun; partials are combined into out_csv at the end.
part_dir <- paste0(out_csv, ".parts")
dir.create(part_dir, showWarnings = FALSE, recursive = TRUE)

v <- fread(var_csv, showProgress = FALSE)
if (nrow(v) == 0) stop("empty input: ", var_csv)

v[, phenotype := sub("_[^_]+$", "", source_file)]
v[, ancestry  := sub("^.*_",   "", source_file)]
v[, variant   := paste(chromosome, position, ref, alt, sep = ":")]

targets <- unique(v[, .(phenotype, ancestry, chromosome, position, ref, alt, variant)])
phenos  <- unique(targets$phenotype)

read_stratum <- function(pheno, strat) {
  f <- file.path(src_dir, paste0(pheno, "_", strat, ".tsv.gz"))
  if (!file.exists(f)) return(NULL)
  d <- tryCatch(fread(f, showProgress = FALSE), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  needed <- c("chromosome","position","ref","alt","beta","se","pvalue")
  if (!all(needed %in% names(d))) return(NULL)
  d <- d[, ..needed]
  d[, key := paste(chromosome, position, ref, alt, sep = ":")]
  setkey(d, key); d
}

pick <- function(d, key) {
  if (is.null(d)) return(list(beta = NA_real_, se = NA_real_, pvalue = NA_real_))
  hit <- d[key]
  if (nrow(hit) == 0 || is.na(hit$beta[1]))
    return(list(beta = NA_real_, se = NA_real_, pvalue = NA_real_))
  list(beta = hit$beta[1], se = hit$se[1], pvalue = hit$pvalue[1])
}

process_pheno <- function(p) {
  part  <- file.path(part_dir, paste0(gsub("[^A-Za-z0-9_.-]", "_", p), ".csv"))
  donef <- paste0(part, ".done")
  if (file.exists(donef)) {   # resume: already completed
    message(sprintf("[%s] skip (done): %s", format(Sys.time(), "%H:%M:%S"), p))
    return(invisible(NULL))
  }
  tp   <- targets[phenotype == p]
  ancs <- unique(tp$ancestry)                 # only the non-EUR strata we need

  # read this phenotype's needed files: each non-EUR ancestry + EUR
  anc_cache <- setNames(lapply(ancs, function(a) {
    d <- read_stratum(p, a)
    if (is.null(d)) message("missing/unreadable: ",
                            file.path(src_dir, paste0(p, "_", a, ".tsv.gz")))
    d
  }), ancs)
  eur <- read_stratum(p, "EUR")
  if (is.null(eur)) message("missing/unreadable: ", file.path(src_dir, paste0(p, "_EUR.tsv.gz")))

  out <- rbindlist(lapply(seq_len(nrow(tp)), function(i) {
    k <- tp$variant[i]
    a <- pick(anc_cache[[tp$ancestry[i]]], k)
    e <- pick(eur, k)
    data.table(
      phenotype = p, ancestry = tp$ancestry[i], variant = k,
      chromosome = tp$chromosome[i], position = tp$position[i],
      ref = tp$ref[i], alt = tp$alt[i],
      anc_beta = a$beta, anc_se = a$se, anc_p = a$pvalue,
      eur_beta = e$beta, eur_se = e$se, eur_p = e$pvalue,
      beta_diff = a$beta - e$beta
    )
  }))
  fwrite(out, part)          # write partial first
  file.create(donef)          # mark done only after the write succeeds
  message(sprintf("[%s] done: %s (%d variants)",
                  format(Sys.time(), "%H:%M:%S"), p, nrow(tp)))
  invisible(NULL)             # anc_cache/eur dropped when this returns
}

invisible(mclapply(phenos, process_pheno, mc.cores = n_cores))

# combine all per-phenotype partials into the final output
parts <- list.files(part_dir, pattern = "\\.csv$", full.names = TRUE)
res <- rbindlist(lapply(parts, function(pf) {
  x <- tryCatch(fread(pf, showProgress = FALSE), error = function(e) NULL)
  if (is.null(x) || nrow(x) == 0) NULL else x
}), fill = TRUE)

if (is.null(res) || nrow(res) == 0) {
  fwrite(data.table(), out_csv)
  cat("No results.\n")
} else {
  setorder(res, phenotype, ancestry, chromosome, position)
  fwrite(res, out_csv)
  cat(sprintf("Looked up %d variants across %d phenotype(s) -> %s\n",
              nrow(res), uniqueN(res$phenotype), out_csv))
}
cat("Partial results kept in:", part_dir, "(delete to force full re-run)\n")
