#!/usr/bin/env python3
"""
Plot n_loci per phenotype from count_unique_clumps.py output.

Reads the TSV (phenotype, n_loci, n_clumps_raw, strata, per_stratum_raw) and
draws a sorted horizontal bar chart of n_loci, one bar per phenotype, value
labelled at the end of each bar. Figure height scales with the number of
phenotypes so labels never overlap.

Usage:
    python3 plot_unique_clumps.py clumps_per_phenotype.tsv
    python3 plot_unique_clumps.py clumps_per_phenotype.tsv --out loci.png --top 40
"""

import argparse
import csv
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("tsv", help="count_unique_clumps.py output")
    ap.add_argument("--out", default=None,
                    help="output image (default: <tsv stem>_n_loci.png)")
    ap.add_argument("--top", type=int, default=None,
                    help="show only the N phenotypes with most loci")
    ap.add_argument("--title", default="Unique loci per phenotype")
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()

    rows = []
    with open(args.tsv, encoding="utf-8") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            try:
                rows.append((r["phenotype"], int(r["n_loci"])))
            except (KeyError, ValueError):
                continue
    if not rows:
        sys.exit(f"ERROR: no usable rows in {args.tsv} "
                 f"(need columns phenotype, n_loci)")

    rows.sort(key=lambda x: x[1], reverse=True)   # descending -> tallest on left
    if args.top:
        rows = rows[:args.top]
    names = [r[0] for r in rows]
    vals = [r[1] for r in rows]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    w = max(6.0, 0.28 * len(rows) + 1.5)
    fig, ax = plt.subplots(figsize=(w, 6))
    bars = ax.bar(range(len(rows)), vals, color="#3b6ea5", edgecolor="none")
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels(names, fontsize=11, rotation=90)
    ax.set_ylabel("Number of unique loci", fontsize=13)
    ax.set_title(args.title, fontweight="bold", fontsize=14)
    ax.tick_params(axis="y", labelsize=11)
    ax.margins(x=0.005)

    ymax = max(vals) if vals else 1
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2,
                b.get_height() + ymax * 0.01,
                str(v), ha="center", va="bottom", fontsize=9)
    ax.set_ylim(0, ymax * 1.08)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()

    out = args.out or (args.tsv.rsplit(".", 1)[0] + "_n_loci.png")
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight")
    sys.stderr.write(f"Wrote {out}  ({len(rows)} phenotypes, "
                     f"n_loci {min(vals)}-{max(vals)})\n")


if __name__ == "__main__":
    main()
