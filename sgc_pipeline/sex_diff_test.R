#!/usr/bin/env Rscript
# Formal test for a sex difference in effect size at each variant, using the
# two-sample z-test on betas from independent female/male GWAS:
#
#     z = (beta_female - beta_male) / sqrt(se_female^2 + se_male^2)
#     p = 2 * pnorm(-abs(z))            (two-sided)
#
# Valid because the two GWAS share no individuals (sex-stratified), so the
# covariance term is zero. Raw p-values only (no multiple-testing correction).
#
# Input: the sex comparison file from sex_lookup.R, with columns
#        female_beta, female_se, male_beta, male_se (+ identifiers).
#
# Usage:   Rscript sex_diff_test.R <comparison_csv> <out_csv>
# Example: Rscript sex_diff_test.R sex_effect_comparison.csv sex_diff_test.csv

suppressPackageStartupMessages(library(data.table))

args    <- commandArgs(trailingOnly = TRUE)
in_csv  <- if (length(args) >= 1) args[1] else "sex_effect_comparison.csv"
out_csv <- if (length(args) >= 2) args[2] else "sex_diff_test.csv"

d <- fread(in_csv, showProgress = FALSE)
if (nrow(d) == 0) stop("empty input: ", in_csv)

need <- c("female_beta","female_se","male_beta","male_se")
if (!all(need %in% names(d)))
  stop("input missing required columns: ", paste(setdiff(need, names(d)), collapse=", "))

# two-sample z-test on the beta difference
d[, beta_diff := female_beta - male_beta]
d[, se_diff   := sqrt(female_se^2 + male_se^2)]
d[, z_diff    := beta_diff / se_diff]
d[, p_diff    := 2 * pnorm(-abs(z_diff))]

# rows with any missing beta/se can't be tested -> NA (kept, flagged)
d[, testable := is.finite(female_beta) & is.finite(female_se) &
                is.finite(male_beta)   & is.finite(male_se)   & se_diff > 0]
d[testable == FALSE, c("z_diff","p_diff") := NA_real_]

setorder(d, p_diff, na.last = TRUE)
fwrite(d, out_csv)

n_test <- sum(d$testable)
n_sig  <- sum(d$p_diff < 0.05, na.rm = TRUE)
cat(sprintf("Tested %d of %d variants; %d with raw p_diff < 0.05\n",
            n_test, nrow(d), n_sig))
cat("Wrote:", out_csv, "\n")
