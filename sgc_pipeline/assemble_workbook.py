#!/usr/bin/env python3
"""
Assemble the final workbook from the pipeline's per-tab TSV outputs.

Collects whatever tab files exist in --out-dir and lays each on its sheet with
the workbook's column order (including the manually-filled rsid / AF columns,
left blank for you to paste). Also embeds the novelty catalogue so the file is
self-contained. Missing inputs are skipped with a note, so you can assemble a
partial workbook and re-run as more stages complete.

Only the CURRENT sheets are produced (the ".A" versions and the ones with no
duplicate) -- no "_old" or plain-duplicate tabs.

Inputs looked for in --out-dir (all optional):
    overview.tsv                       -> Overview
    meta_results.tsv                   -> Meta-Analysis Results
    min_pvalues.csv                    -> GWAS Minimum p_ALL
    clump_index_novelty.tsv            -> Clump_index_variants   (+ blank rsid/AF)
    cross_disease.tsv                  -> Cross-disease variants (+ blank rsid/AF)
    sex_specific_novelty.tsv           -> Sex specific loci      (+ blank rsid/AF)
    ancestry_specific_novelty.tsv      -> Ancestry specific loci (+ blank rsid/AF)
    recovery_recovery_exact.tsv        -> Known loci in our GWAS
--novelty (the catalogue TSV)          -> Previous Known loci- variants

Usage:
    python3 assemble_workbook.py --out-dir results \\
        --novelty novelty_info.txt --xlsx results/SGC_Meta_Analysis_Results.xlsx
"""

import argparse
import csv
import os
import sys


def read_tsv(path):
    if not path or not os.path.exists(path):
        return None
    with open(path, encoding="utf-8-sig", errors="replace") as fh:
        first = fh.readline()
        delim = "\t" if "\t" in first else ("," if "," in first else "\t")
        fh.seek(0)
        rows = list(csv.reader(fh, delimiter=delim))
    return rows if rows else None


def insert_blank_cols(rows, after_col_name, new_names):
    """Insert empty columns (new_names) immediately AFTER the column called
    after_col_name. Returns modified rows. If the column isn't found, prepend."""
    hdr = rows[0]
    try:
        i = hdr.index(after_col_name) + 1
    except ValueError:
        i = 1
    for r_i, r in enumerate(rows):
        fill = new_names if r_i == 0 else [""] * len(new_names)
        rows[r_i] = r[:i] + fill + r[i:]
    return rows


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default="results")
    ap.add_argument("--novelty", default=None,
                    help="catalogue TSV to embed as 'Previous Known loci- variants'")
    ap.add_argument("--xlsx", required=True)
    args = ap.parse_args()

    from openpyxl import Workbook
    from openpyxl.styles import Font, PatternFill, Alignment

    d = args.out_dir
    P = lambda f: os.path.join(d, f)

    # (sheet name, source file, [(after_col, [blank cols to insert]) ...])
    plan = [
        ("Overview", P("overview.tsv"), []),
        ("Meta-Analysis Results", P("meta_results.tsv"), []),
        ("GWAS Minimum p_ALL", P("min_pvalues.csv"), []),
        ("Clump_index_variants", P("clump_index_novelty.tsv"),
         [("index_variant", ["rsid", "AF"])]),
        ("Cross-disease variants", P("cross_disease_full.tsv"),
         [("Locus_variant", ["rsid", "AF"])]),
        ("Ancestry specific loci", P("ancestry_specific_novelty.tsv"),
         [("variant", ["rsid", "AF"])]),
        ("Sex specific loci", P("sex_specific_novelty.tsv"),
         [("Variant", ["rsid", "AF"])]),
        ("Known loci in our GWAS", P("recovery_recovery_exact.tsv"), []),
    ]
    if args.novelty:
        plan.append(("Previous Known loci- variants", args.novelty, []))

    wb = Workbook()
    wb.remove(wb.active)
    hfont = Font(bold=True, color="FFFFFF")
    hfill = PatternFill("solid", fgColor="4472C4")
    halign = Alignment(horizontal="left", vertical="center")

    made, skipped = [], []
    for sheet, path, inserts in plan:
        rows = read_tsv(path)
        if rows is None:
            skipped.append(sheet)
            continue
        for after, cols in inserts:
            rows = insert_blank_cols(rows, after, cols)
        ws = wb.create_sheet(title=sheet[:31])
        for r in rows:
            ws.append(r)
        for c in ws[1]:
            c.font = hfont; c.fill = hfill; c.alignment = halign
        ws.freeze_panes = "A2"
        if ws.max_row > 1 and ws.max_column >= 1:
            ws.auto_filter.ref = ws.dimensions
        made.append((sheet, ws.max_row - 1))

    if not made:
        sys.exit("Nothing to assemble — no tab files found in " + d)

    wb.save(args.xlsx)
    sys.stderr.write(f"Wrote {args.xlsx}\n  sheets:\n")
    for s, n in made:
        sys.stderr.write(f"    {s}  ({n} rows)\n")
    if skipped:
        sys.stderr.write("  skipped (no input yet): " + ", ".join(skipped) + "\n")


if __name__ == "__main__":
    main()
