#!/usr/bin/env Rscript
# Formal test for an effect-size difference between a variant's non-EUR ancestry
# and EUR, using the two-sample z-test on betas from independent GWAS:
#
#     z = (beta_anc - beta_eur) / sqrt(se_anc^2 + se_eur^2)
#     p = 2 * pnorm(-abs(z))            (two-sided)
#
# Valid because the ancestry and EUR GWAS share no individuals, so the
# covariance term is zero. Raw p-values only (no multiple-testing correction).
#
# Input: the ancestry comparison file from ancestry_lookup.R, with columns
#        anc_beta, anc_se, eur_beta, eur_se (+ identifiers incl. ancestry).
#
# Usage:   Rscript ancestry_diff_test.R <comparison_csv> <out_csv>
# Example: Rscript ancestry_diff_test.R ancestry_effect_comparison.csv ancestry_diff_test.csv

suppressPackageStartupMessages(library(data.table))

args    <- commandArgs(trailingOnly = TRUE)
in_csv  <- if (length(args) >= 1) args[1] else "ancestry_effect_comparison.csv"
out_csv <- if (length(args) >= 2) args[2] else "ancestry_diff_test.csv"

d <- fread(in_csv, showProgress = FALSE)
if (nrow(d) == 0) stop("empty input: ", in_csv)

need <- c("anc_beta","anc_se","eur_beta","eur_se")
if (!all(need %in% names(d)))
  stop("input missing required columns: ", paste(setdiff(need, names(d)), collapse=", "))

# two-sample z-test on the beta difference (non-EUR minus EUR)
d[, beta_diff := anc_beta - eur_beta]
d[, se_diff   := sqrt(anc_se^2 + eur_se^2)]
d[, z_diff    := beta_diff / se_diff]
d[, p_diff    := 2 * pnorm(-abs(z_diff))]

# rows with any missing beta/se can't be tested -> NA (kept, flagged)
d[, testable := is.finite(anc_beta) & is.finite(anc_se) &
                is.finite(eur_beta) & is.finite(eur_se) & se_diff > 0]
d[testable == FALSE, c("z_diff","p_diff") := NA_real_]

setorder(d, p_diff, na.last = TRUE)
fwrite(d, out_csv)

n_test <- sum(d$testable)
n_sig  <- sum(d$p_diff < 0.05, na.rm = TRUE)
cat(sprintf("Tested %d of %d variants; %d with raw p_diff < 0.05\n",
            n_test, nrow(d), n_sig))
cat("Wrote:", out_csv, "\n")
