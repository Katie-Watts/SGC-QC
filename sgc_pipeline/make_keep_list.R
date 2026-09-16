#!/usr/bin/env Rscript
# Build a variant-id keep-list from a GWAS, as chr:pos:ref:alt (the same
# convention harmonize_bim.R rewrote the panel ids to). Use the output with
# PLINK --extract to shrink the reference panel before clumping.
#
# Usage:   Rscript make_keep_list.R <gwas(.tsv/.tsv.gz)> <out.txt>
# Example: Rscript make_keep_list.R ATOPIC_DERM_ALL.tsv keep_ids.txt
#
# GWAS columns needed: chromosome position ref alt.

suppressPackageStartupMessages(library(data.table))
setDTthreads(0)

args    <- commandArgs(trailingOnly = TRUE)
gwas_f  <- if (length(args) >= 1) args[1] else stop("need GWAS file")
out_f   <- if (length(args) >= 2) args[2] else "keep_ids.txt"

g <- fread(gwas_f, showProgress = FALSE,
           select = c("chromosome","position","ref","alt"))
g[, chromosome := sub("^chr","", as.character(chromosome))]
g[, id := paste(chromosome, position, toupper(ref), toupper(alt), sep = ":")]

# unique ids (a GWAS shouldn't have dups, but guard anyway)
ids <- unique(g$id)
writeLines(ids, out_f)

cat("GWAS variants:", nrow(g), "\n")
cat("Unique ids written:", length(ids), "->", out_f, "\n")
