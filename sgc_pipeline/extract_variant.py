#!/usr/bin/env python3
"""
Extract a variant (by chromosome + position) from a GWAS file, with the header.

Finds the chromosome and position columns by NAME (case-insensitive), so it
works across cohorts with different headers. Handles gzip and plain text.
Position-only match: alleles are ignored, so it returns every row at that
chr:pos (useful for spotting multi-allelics).

Usage:
    python3 extract_variant.py FILE CHR POS
    python3 extract_variant.py BBJ1_ATOPIC_DERM_EAS_ALL.txt.gz 10 39140359

    # several files at once (header printed once per file, with a filename tag):
    python3 extract_variant.py 10 39140359 --files cohort_files/*_ALL*.txt.gz

    # restrict to one phenotype's cohort files by name:
    python3 extract_variant.py 10 39140359 --pheno MELANOMA_MALIGNT \\
        --files cohort_files/*_ALL*.txt.gz

Prints the header row then matching data row(s), tab-separated.
"""

import argparse
import gzip
import io
import os
import re
import sys

CHR_ALIASES = {"CHR", "CHROMOSOME", "CHROM", "#CHROM", "CHR_ID"}
POS_ALIASES = {"POS", "POSITION", "GENPOS", "BP", "BP_HG38", "BASE_PAIR_LOCATION"}


def norm_chr(c):
    c = re.sub(r"^chr", "", str(c).strip(), flags=re.I).upper()
    return {"23": "X", "24": "Y", "25": "X", "26": "MT", "M": "MT"}.get(c, c)


def open_text(path):
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8",
                                errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def sniff_delim(line):
    if "\t" in line:
        return "\t"
    if "," in line:
        return ","
    return None


def find_cols(header):
    ic = ip = None
    for i, h in enumerate(header):
        u = h.strip().upper()
        if ic is None and u in CHR_ALIASES:
            ic = i
        if ip is None and u in POS_ALIASES:
            ip = i
    return ic, ip


def extract(path, chrom, pos, tag=False):
    want_chr, want_pos = norm_chr(chrom), str(int(pos))
    with open_text(path) as fh:
        first = fh.readline()
        delim = sniff_delim(first)
        header = first.rstrip("\n").split(delim) if delim else first.split()
        ic, ip = find_cols(header)
        if ic is None or ip is None:
            sys.stderr.write(f"WARN [{path}] could not find chr/pos columns "
                             f"in header: {header[:8]}...\n")
            return 0
        pre = f"{path}\t" if tag else ""
        # header (print once)
        sys.stdout.write(pre + (delim if delim else " ").join(header) + "\n")
        n = 0
        for line in fh:
            core = line.rstrip("\n")
            f = core.split(delim) if delim else core.split()
            if len(f) <= max(ic, ip):
                continue
            if norm_chr(f[ic]) == want_chr and f[ip].strip() == want_pos:
                sys.stdout.write(pre + core + "\n")
                n += 1
        if n == 0:
            sys.stderr.write(f"  [{path}] no row at {want_chr}:{want_pos}\n")
        return n


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("args", nargs="+",
                    help="FILE CHR POS  (or  CHR POS --files F1 F2 ...)")
    ap.add_argument("--files", nargs="+", help="one or more GWAS files")
    ap.add_argument("--pheno", default=None,
                    help="only search files whose name contains this string "
                         "(e.g. ATOPIC_DERM); applied to --files")
    a = ap.parse_args()

    if a.files:
        # form: CHR POS --files ...
        if len(a.args) < 2:
            ap.error("need CHR and POS before --files")
        chrom, pos = a.args[0], a.args[1]
        files = a.files
        if a.pheno:
            kept = [f for f in files if a.pheno in os.path.basename(f)]
            sys.stderr.write(f"--pheno {a.pheno!r}: {len(kept)}/{len(files)} "
                             f"files match\n")
            files = kept
            if not files:
                sys.exit(f"No files contain {a.pheno!r} in their name.")
        multi = len(files) > 1
        total = 0
        for f in files:
            total += extract(f, chrom, pos, tag=multi)
        sys.stderr.write(f"\n{total} matching row(s) across {len(files)} file(s)\n")
    else:
        # form: FILE CHR POS
        if a.pheno:
            sys.stderr.write("WARN: --pheno only applies with --files; ignored\n")
        if len(a.args) != 3:
            ap.error("usage: FILE CHR POS   (or  CHR POS --files ...)")
        path, chrom, pos = a.args
        extract(path, chrom, pos)


if __name__ == "__main__":
    main()
