#!/usr/bin/env Rscript
# For each variant in a pruned results file, pull beta/se/pvalue for BOTH sexes
# (FEMALE and MALE) of that variant's phenotype from the original source .tsv.gz
# files, so the actual effect-size difference is visible.
#
# Memory-safe: processes ONE phenotype at a time (reads its two sex files, does
# its lookups, discards them), so peak memory holds at most n_cores phenotypes'
# files -- independent of how many phenotypes the pruned set spans.
#
# Input columns: source_file, chromosome, position, ref, alt (prune_250kb.R out)
# Source files:  <src_dir>/<PHENOTYPE>_FEMALE.tsv.gz , _MALE.tsv.gz
#                with columns chromosome position ref alt beta se pvalue ...
#
# Usage:   Rscript sex_lookup.R <variants_csv> <src_dir> <out_csv> [n_cores]
# Example: Rscript sex_lookup.R sex_pruned.csv . sex_effect_comparison.csv 3

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

args    <- commandArgs(trailingOnly = TRUE)
var_csv <- if (length(args) >= 1) args[1] else stop("need variants csv")
src_dir <- if (length(args) >= 2) args[2] else "."
out_csv <- if (length(args) >= 3) args[3] else "sex_effect_comparison.csv"
n_cores <- if (length(args) >= 4) as.integer(args[4]) else 3L
setDTthreads(1)

# per-phenotype partial results (resume support): a phenotype with a .done
# marker is skipped on rerun; partials are combined into out_csv at the end.
part_dir <- paste0(out_csv, ".parts")
dir.create(part_dir, showWarnings = FALSE, recursive = TRUE)

v <- fread(var_csv, showProgress = FALSE)
if (nrow(v) == 0) stop("empty input: ", var_csv)

v[, phenotype := sub("_[^_]+$", "", source_file)]
v[, variant   := paste(chromosome, position, ref, alt, sep = ":")]

targets <- unique(v[, .(phenotype, chromosome, position, ref, alt, variant)])
phenos  <- unique(targets$phenotype)
sexes   <- c("FEMALE", "MALE")

# read one stratum file keyed on chr:pos:ref:alt (NULL if missing/malformed)
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

# process a SINGLE phenotype end-to-end, then its files go out of scope
process_pheno <- function(p) {
  part  <- file.path(part_dir, paste0(gsub("[^A-Za-z0-9_.-]", "_", p), ".csv"))
  donef <- paste0(part, ".done")
  if (file.exists(donef)) {   # resume: already completed
    message(sprintf("[%s] skip (done): %s", format(Sys.time(), "%H:%M:%S"), p))
    return(invisible(NULL))
  }
  fem <- read_stratum(p, "FEMALE")
  mal <- read_stratum(p, "MALE")
  if (is.null(fem)) message("missing/unreadable: ", file.path(src_dir, paste0(p, "_FEMALE.tsv.gz")))
  if (is.null(mal)) message("missing/unreadable: ", file.path(src_dir, paste0(p, "_MALE.tsv.gz")))

  tp <- targets[phenotype == p]
  out <- rbindlist(lapply(seq_len(nrow(tp)), function(i) {
    k <- tp$variant[i]
    f <- pick(fem, k); m <- pick(mal, k)
    data.table(
      phenotype = p, variant = k,
      chromosome = tp$chromosome[i], position = tp$position[i],
      ref = tp$ref[i], alt = tp$alt[i],
      female_beta = f$beta, female_se = f$se, female_p = f$pvalue,
      male_beta   = m$beta, male_se   = m$se, male_p   = m$pvalue,
      beta_diff   = f$beta - m$beta
    )
  }))
  fwrite(out, part)          # write partial first
  file.create(donef)          # mark done only after the write succeeds
  message(sprintf("[%s] done: %s (%d variants)",
                  format(Sys.time(), "%H:%M:%S"), p, nrow(tp)))
  invisible(NULL)             # fem/mal dropped when this function returns
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
  setorder(res, phenotype, chromosome, position)
  fwrite(res, out_csv)
  cat(sprintf("Looked up %d variants across %d phenotype(s) -> %s\n",
              nrow(res), uniqueN(res$phenotype), out_csv))
}
cat("Partial results kept in:", part_dir, "(delete to force full re-run)\n")
