#!/usr/bin/env Rscript
# Extract all variants with pvalue < threshold from every .tsv.gz file in a
# folder into ONE combined CSV, tagging each row with its source filename.
# Resumable: each input's hits are written to a per-file partial as soon as it
# finishes, so an interrupted run resumes and skips already-processed files.
#
# Usage:   Rscript extract_suggestive.R <in_dir> <out_csv> [n_cores]
# Example: Rscript extract_suggestive.R . suggestive_hits.csv 3

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

args    <- commandArgs(trailingOnly = TRUE)
in_dir  <- if (length(args) >= 1) args[1] else "."
out_csv <- if (length(args) >= 2) args[2] else "suggestive_hits.csv"
n_cores <- if (length(args) >= 3) as.integer(args[3]) else 3L

col_p     <- "pvalue"     # p-value column name
threshold <- 1e-5         # suggestive significance cutoff

setDTthreads(1)

# per-file partial results live here; used for resume and combined at the end
part_dir <- paste0(out_csv, ".parts")
dir.create(part_dir, showWarnings = FALSE, recursive = TRUE)

files <- list.files(in_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)
message("Found ", length(files), " files; using ", n_cores, " cores.")

process_one <- function(f) {
  name  <- sub("\\.tsv\\.gz$", "", basename(f))
  part  <- file.path(part_dir, paste0(name, ".csv"))
  donef <- file.path(part_dir, paste0(name, ".done"))

  # resume: skip if this file was already completed
  if (file.exists(donef)) {
    message(sprintf("[%s] skip (done): %s", format(Sys.time(), "%H:%M:%S"), name))
    return(invisible(NULL))
  }

  d <- tryCatch(fread(f, showProgress = FALSE), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0 || !(col_p %in% names(d))) {
    fwrite(data.table(), part)                 # empty partial
    file.create(donef)                          # mark done
    message(sprintf("[%s] done (no hits): %s", format(Sys.time(), "%H:%M:%S"), name))
    return(invisible(NULL))
  }

  d[, (col_p) := suppressWarnings(as.numeric(get(col_p)))]
  hits <- d[get(col_p) > 0 & get(col_p) < threshold]

  if (nrow(hits) > 0) {
    hits[, source_file := name]
    setcolorder(hits, c("source_file", setdiff(names(hits), "source_file")))
    fwrite(hits, part)
  } else {
    fwrite(data.table(), part)
  }
  file.create(donef)                            # mark done only after write succeeds
  message(sprintf("[%s] done (%d hits): %s",
                  format(Sys.time(), "%H:%M:%S"), nrow(hits), name))
  invisible(NULL)
}

invisible(mclapply(files, process_one, mc.cores = n_cores))

# --- combine all non-empty partials into the final output ------------------
parts <- list.files(part_dir, pattern = "\\.csv$", full.names = TRUE)
combined <- rbindlist(lapply(parts, function(p) {
  x <- tryCatch(fread(p, showProgress = FALSE), error = function(e) NULL)
  if (is.null(x) || nrow(x) == 0) NULL else x
}), fill = TRUE)

if (is.null(combined) || nrow(combined) == 0) {
  message("No variants below ", threshold, " in any file.")
  fwrite(data.table(source_file = character()), out_csv)
} else {
  setorder(combined, pvalue)
  fwrite(combined, out_csv)
  message("Wrote ", out_csv, " with ", nrow(combined),
          " variants from ", uniqueN(combined$source_file), " files.")
}

message("Partial results kept in: ", part_dir,
        "  (delete this folder to force a full re-run)")
