#!/usr/bin/env python3
"""
Seed the Clump_index_variants table from verbose_summary.sh's _clumps.csv.

verbose_summary.sh writes  <pre>_clumps.csv  with columns
    file, index_variant, index_pvalue, n_gwsig_in_clump
(one row per clump; MHC already collapsed to a single representative). This
script just renames/reshapes it into the sheet's leading columns:
    GWAS, index_variant, index_variant_pvalue, n_gwsig_variants_in_clump

Everything else on the sheet (rsid, AF, n_cohorts_variant, i2,
beta_direction_concordance, Novelty/Matched/Details, cohorts/cases/controls,
the clump_* QC columns) is added by later steps or by hand -- this only lays
down the four columns that come straight from clumping.

Usage:
    python3 make_clump_index.py --clumps verbose_summary_clumps.csv \\
        --out clump_index.tsv
"""

import argparse
import csv
import sys


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--clumps", required=True,
                    help="verbose_summary _clumps.csv")
    ap.add_argument("--out", default="clump_index.tsv")
    args = ap.parse_args()

    n = 0
    with open(args.clumps, newline="", encoding="utf-8-sig") as fh, \
         open(args.out, "w") as out:
        r = csv.DictReader(fh)
        out.write("GWAS\tindex_variant\tindex_variant_pvalue\t"
                  "n_gwsig_variants_in_clump\n")
        for row in r:
            gwas = row.get("file") or row.get("File")
            var = row.get("index_variant")
            p = row.get("index_pvalue", "")
            ng = row.get("n_gwsig_in_clump", "")
            if not gwas or not var:
                continue
            out.write(f"{gwas.strip()}\t{var.strip()}\t{p}\t{ng}\n")
            n += 1
    sys.stderr.write(f"Wrote {args.out}: {n} clump index variants\n")


if __name__ == "__main__":
    main()
