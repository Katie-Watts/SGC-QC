#!/usr/bin/env python3
"""
Build the OVERVIEW table: one row per phenotype, a Yes/blank cell per stratum
showing which <PHENOTYPE>_<STRATUM>.tsv.gz GWAS files exist, and a count.

Scans --gwas-dir for *.tsv.gz, splits each name into PHENOTYPE + STRATUM (the
last underscore token, matched against the known strata), and tallies presence.

Layout matches the workbook's Overview sheet:
    row 1:  Phenotype | Meta-analysis run ?
    row 2:  (blank)   | ALL (combined/all individuals) | Female | Male | AFR |
            AMR | EAS | EUR | SAS | (count)
    row 3+: <pheno>   | Yes/blank per stratum ...                     | N

Output is a TSV (phenotype, then one column per stratum, then n_strata); the
workbook assembler places it on the Overview sheet with the two-row header.

Usage:
    python3 build_overview.py --gwas-dir . --out overview.tsv
"""

import argparse
import os
import re
import sys
from collections import defaultdict

# display order + labels exactly as the workbook Overview sheet uses them
STRATA = ["ALL", "FEMALE", "MALE", "AFR", "AMR", "EAS", "EUR", "SAS"]
LABELS = {"ALL": "ALL (combined/all individuals)", "FEMALE": "Female",
          "MALE": "Male", "AFR": "AFR", "AMR": "AMR", "EAS": "EAS",
          "EUR": "EUR", "SAS": "SAS"}
KNOWN = set(STRATA)


def split_pheno_stratum(stem):
    m = re.match(r"^(.*)_([A-Za-z]+)$", stem)
    if m and m.group(2).upper() in KNOWN:
        return m.group(1), m.group(2).upper()
    return stem, ""


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gwas-dir", default=".")
    ap.add_argument("--suffix", default=".tsv.gz")
    ap.add_argument("--out", default="overview.tsv")
    args = ap.parse_args()

    present = defaultdict(set)          # phenotype -> {strata}
    unmatched = []
    for fn in sorted(os.listdir(args.gwas_dir)):
        if not fn.endswith(args.suffix) or fn.startswith("._"):
            continue
        stem = fn[:-len(args.suffix)]
        pheno, stratum = split_pheno_stratum(stem)
        if not stratum:
            unmatched.append(fn)
            continue
        present[pheno].add(stratum)

    with open(args.out, "w") as out:
        # header rows mirror the sheet (labels row + count col)
        out.write("Phenotype\t" + "\t".join(LABELS[s] for s in STRATA)
                  + "\tn_strata\n")
        for pheno in sorted(present):
            cells = ["Yes" if s in present[pheno] else "" for s in STRATA]
            n = sum(1 for c in cells if c)
            out.write(pheno + "\t" + "\t".join(cells) + f"\t{n}\n")

    sys.stderr.write(
        f"Overview: {len(present)} phenotypes, "
        f"{sum(len(v) for v in present.values())} GWAS files -> {args.out}\n")
    if unmatched:
        sys.stderr.write(
            f"NOTE: {len(unmatched)} file(s) had no recognised stratum suffix "
            f"and were skipped, e.g. {unmatched[:5]}\n")


if __name__ == "__main__":
    main()
