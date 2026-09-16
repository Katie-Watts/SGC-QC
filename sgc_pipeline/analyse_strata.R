#!/usr/bin/env Rscript
# Analyse the combined suggestive-hits file two ways.
#
#  (1) SEX-SPECIFIC: variants genome-wide significant (p < 5e-8) for a
#      phenotype in FEMALE but absent from the file for MALE (same phenotype),
#      and vice versa.
#  (2) ANCESTRY-SPECIFIC: variants genome-wide significant (p < 5e-8) in any of
#      AFR / AMR / EAS / SAS for a phenotype, but absent from the file for EUR.
#
# "Absent from the file" = not present at p < 1e-5 in that stratum (the input
# only contains suggestive rows), so this captures signal that is genome-wide
# significant in one stratum and not even suggestive in the comparator.
#
# Matching across strata is by chromosome + position + ref + alt.
#
# Usage:   Rscript analyse_strata.R <hits_csv> <out_prefix>
# Example: Rscript analyse_strata.R suggestive_hits.csv strata_results

suppressPackageStartupMessages(library(data.table))

args       <- commandArgs(trailingOnly = TRUE)
hits_csv   <- if (length(args) >= 1) args[1] else "suggestive_hits.csv"
out_prefix <- if (length(args) >= 2) args[2] else "strata_results"

gws <- 5e-8   # genome-wide significance threshold

d <- fread(hits_csv, showProgress = FALSE)
if (nrow(d) == 0) stop("Input file is empty: ", hits_csv)

# --- split source_file into phenotype (before last _) and stratum (after) ---
d[, phenotype := sub("_[^_]+$", "", source_file)]
d[, stratum   := sub("^.*_",   "", source_file)]

# variant key: chr:pos:ref:alt
d[, variant := paste(chromosome, position, ref, alt, sep = ":")]

# --- report strata found; anything not in these sets is ignored (e.g. ALL) --
sex_levels    <- c("FEMALE", "MALE")
anc_levels    <- c("AFR", "AMR", "EAS", "SAS", "EUR")
used_levels   <- c(sex_levels, anc_levels)
found_strata  <- sort(unique(d$stratum))
ignored       <- setdiff(found_strata, used_levels)
cat("Strata found:", paste(found_strata, collapse = ", "), "\n")
if (length(ignored) > 0)
  cat("Ignoring strata not in the sex/ancestry sets:",
      paste(ignored, collapse = ", "), "\n")
cat("\n")

# =====================================================================
# (1) SEX-SPECIFIC
# =====================================================================
ds <- d[stratum %in% sex_levels]

# genome-wide significant rows, per stratum
sig_f <- ds[stratum == "FEMALE" & pvalue < gws]
sig_m <- ds[stratum == "MALE"   & pvalue < gws]

# ALL rows present in the file for each stratum (any suggestive p), for the
# "absent from the file" test -- keyed by phenotype + variant
present_f <- unique(ds[stratum == "FEMALE", .(phenotype, variant)])
present_m <- unique(ds[stratum == "MALE",   .(phenotype, variant)])

# female-significant but not present at all in male file (same phenotype)
female_only <- sig_f[!present_m, on = .(phenotype, variant)]
female_only[, group := "FEMALE_only"]

# male-significant but not present at all in female file
male_only <- sig_m[!present_f, on = .(phenotype, variant)]
male_only[, group := "MALE_only"]

sex_out <- rbindlist(list(female_only, male_only), fill = TRUE)
setorder(sex_out, phenotype, pvalue)
fwrite(sex_out, paste0(out_prefix, "_sex_specific.csv"))

# =====================================================================
# (2) ANCESTRY-SPECIFIC: sig in any of AFR/AMR/EAS/SAS, absent in EUR
# =====================================================================
non_eur <- c("AFR", "AMR", "EAS", "SAS")
da <- d[stratum %in% anc_levels]

# genome-wide significant rows in any non-EUR ancestry
sig_noneur <- da[stratum %in% non_eur & pvalue < gws]

# everything present in the EUR file (any suggestive p) per phenotype+variant
present_eur <- unique(da[stratum == "EUR", .(phenotype, variant)])

# non-EUR-significant but not present at all in the EUR file
ancestry_only <- sig_noneur[!present_eur, on = .(phenotype, variant)]
setorder(ancestry_only, phenotype, pvalue)
fwrite(ancestry_only, paste0(out_prefix, "_ancestry_specific.csv"))

# =====================================================================
# summary to stdout
# =====================================================================
cat("=== Sex-specific ===\n")
cat("  FEMALE-significant, absent in MALE file:", nrow(female_only), "variants\n")
cat("  MALE-significant, absent in FEMALE file:", nrow(male_only),   "variants\n")
cat("  written:", paste0(out_prefix, "_sex_specific.csv"), "\n\n")

cat("=== Ancestry-specific (non-EUR sig, absent in EUR) ===\n")
cat("  total:", nrow(ancestry_only), "variants\n")
if (nrow(ancestry_only) > 0) {
  print(ancestry_only[, .N, by = stratum])
}
cat("  written:", paste0(out_prefix, "_ancestry_specific.csv"), "\n")
