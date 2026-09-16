#!/usr/bin/env Rscript
# Rewrite a reference-panel .bim so each variant's ID becomes chr:pos:ref:alt in
# the GWAS's allele orientation, for variants present in BOTH panel and GWAS
# (match by chr+pos and allele pair; allow ref/alt swap and, for unambiguous
# SNPs, strand flip). Lets a GWAS-generated chr:pos:ref:alt id match the panel.
#
# STREAMED / MEMORY-BOUNDED: the GWAS is reduced to a compact key->gid map once;
# the .bim is then processed in CHUNKS (read -> key -> match -> append output),
# so the full panel table never coexists with its output. Peak memory is ~one
# chunk + the GWAS map, independent of panel size -- safe on a loaded machine.
#
# CLUMPING only: only column 2 (ID) is rewritten; .bed and A1/A2 are untouched.
#
# Usage:   Rscript harmonize_bim.R <panel.bim> <gwas(.tsv/.tsv.gz)> <out_prefix> [chunk_rows]
# Example: Rscript harmonize_bim.R all_hg38_ids.bim gwas.tsv all_hg38_matched
#          Rscript harmonize_bim.R all_hg38_ids.bim gwas.tsv all_hg38_matched 2000000
#
# GWAS cols: chromosome position ref alt (+others). bim: chr id cM pos A1 A2.

suppressPackageStartupMessages(library(data.table))
setDTthreads(0)

args     <- commandArgs(trailingOnly = TRUE)
bim_f    <- if (length(args) >= 1) args[1] else stop("need panel .bim")
gwas_f   <- if (length(args) >= 2) args[2] else stop("need GWAS file")
out_pre  <- if (length(args) >= 3) args[3] else "harmonized"
chunk_n  <- if (length(args) >= 4) as.integer(args[4]) else 2000000L

enc   <- c(A=1L, T=2L, C=3L, G=4L)
compi <- c(2L, 1L, 4L, 3L)
code  <- function(v) { i <- enc[v]; i[is.na(i)] <- 0L; i }
comp  <- function(i) fifelse(i >= 1L & i <= 4L, compi[pmax(i,1L)], 0L)
pair_int <- function(a, b) {
  lo <- pmin(a, b); hi <- pmax(a, b)
  fifelse(lo == 0L | hi == 0L, NA_integer_, lo * 10L + hi)
}

# ---- 1. GWAS -> compact key->gid map, then free the big table ----
g <- fread(gwas_f, showProgress = FALSE,
           select = c("chromosome","position","ref","alt"))
g[, chr := sub("^chr","", as.character(chromosome))]
g[, gid := paste(chr, position, toupper(ref), toupper(alt), sep = ":")]
g[, rc := code(toupper(ref))][, ac := code(toupper(alt))]
g[, pk := pair_int(rc, ac)]
g[, k  := fifelse(!is.na(pk),
                  paste(chr, position, pk, sep = ":"),
                  paste(chr, position,
                        pmin(toupper(ref),toupper(alt)),
                        pmax(toupper(ref),toupper(alt)), sep=":"))]
gk <- g[, .(gid = gid[1]), by = k]
rm(g); gc()
setkey(gk, k)

# ---- 2. stream the bim in chunks ----
bim_out  <- paste0(out_pre, ".bim")
match_out<- paste0(out_pre, ".matched.txt")
if (file.exists(bim_out))   file.remove(bim_out)
if (file.exists(match_out)) file.remove(match_out)

status_tot <- list()   # named counts accumulated across chunks
n_total <- 0L; n_matched <- 0L

process_chunk <- function(b, first_chunk) {
  setnames(b, c("chr","id","cm","pos","A1","A2"))
  b[, chr := sub("^chr","", as.character(chr))]
  b[, A1 := toupper(A1)][, A2 := toupper(A2)]
  b[, c1 := code(A1)][, c2 := code(A2)]

  # direct key + lookup
  b[, pk := pair_int(c1, c2)]
  b[, kd := fifelse(!is.na(pk),
                    paste(chr, pos, pk, sep = ":"),
                    paste(chr, pos, pmin(A1,A2), pmax(A1,A2), sep=":"))]
  b[, gid_direct := gk[.(kd), gid]]

  # flipped key + lookup (ACGT SNPs only)
  b[, cf1 := comp(c1)][, cf2 := comp(c2)]
  b[, pkf := pair_int(cf1, cf2)]
  b[, kf := fifelse(!is.na(pkf), paste(chr, pos, pkf, sep = ":"), NA_character_)]
  b[, gid_flip := NA_character_]
  ix <- which(!is.na(b$kf))
  if (length(ix)) b[ix, gid_flip := gk[.(kf[ix]), gid]]

  b[, palindromic := c1 >= 1L & c1 <= 4L & c1 == comp(c2)]

  b[, status := "no_pos_or_allele_match"]
  b[, new_id := NA_character_]
  b[!is.na(gid_direct) & !palindromic,
    `:=`(new_id = gid_direct, status = "match_same_strand")]
  b[is.na(new_id) & !is.na(gid_flip) & !palindromic,
    `:=`(new_id = gid_flip, status = "match_flipped_strand")]
  b[is.na(new_id) & palindromic & (!is.na(gid_direct) | !is.na(gid_flip)),
    status := "ambiguous_palindromic"]

  b[, id_out := fifelse(!is.na(new_id), new_id, id)]

  # append this chunk's bim rows and matched ids
  fwrite(b[, .(chr, id_out, cm, pos, A1, A2)], bim_out,
         sep = "\t", col.names = FALSE, quote = FALSE, append = !first_chunk)
  mm <- b[!is.na(new_id), new_id]
  if (length(mm)) {
    con_m <- file(match_out, open = if (first_chunk) "w" else "a")
    writeLines(mm, con_m)
    close(con_m)
  }

  list(n = nrow(b), nmatch = length(mm),
       status = table(b$status))
}

# read the bim in chunks via an open connection
con <- file(bim_f, "r")
first <- TRUE
repeat {
  lines <- readLines(con, n = chunk_n)
  if (length(lines) == 0) break
  b <- fread(text = lines, header = FALSE, sep = "\t", showProgress = FALSE)
  r <- process_chunk(b, first)
  first <- FALSE
  n_total   <- n_total + r$n
  n_matched <- n_matched + r$nmatch
  for (nm in names(r$status)) {
    prev <- if (nm %in% names(status_tot)) status_tot[[nm]] else 0L
    status_tot[[nm]] <- prev + r$status[[nm]]
  }
  rm(b, lines); gc(FALSE)
  cat(sprintf("[%s] processed %d variants (%d matched so far)\n",
              format(Sys.time(), "%H:%M:%S"), n_total, n_matched))
}
close(con)

# ---- report ----
rep <- data.table(status = names(status_tot),
                  N = as.integer(unlist(status_tot)))
setorder(rep, -N)
fwrite(rep, paste0(out_pre, ".report.txt"), sep = "\t")

cat("\nPanel variants:", n_total, "\n")
cat("Matched to GWAS (id rewritten):", n_matched, "\n")
print(rep)
cat("\nWrote:", bim_out, "/", match_out, "/", paste0(out_pre, ".report.txt"), "\n")
