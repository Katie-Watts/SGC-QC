#!/usr/bin/env python3
"""
Plot the distribution of how many phenotypes each cross-disease locus spans,
from the "Cross-disease variants" sheet.

It reads the per-locus phenotype count -- the column named n_associated_phenotypes
(or n_phenotypes, whichever is present) -- and if neither exists, counts the
non-empty phenotype_k cells per row instead. It draws a bar for each value
(2 phenotypes, 3, 4, ...), labels counts on top, and leaves a caption panel
below the axis where you can drop annotation bullets (pass --notes or edit the
text file it writes).

Usage:
    python3 plot_pheno_distribution.py SGC_Meta_Analysis_Results.xlsx
    python3 plot_pheno_distribution.py SGC_Meta_Analysis_Results.xlsx \
        --out pheno_dist.png \
        --notes "Most loci are shared by exactly 2 phenotypes" \
                "The long tail is dominated by the MHC and pigmentation loci"

If --notes is omitted the caption area is left blank for you to annotate.
"""

import argparse
import sys
from collections import Counter

SHEET = "Cross-disease variants"
COUNT_CANDIDATES = ("n_associated_phenotypes", "n_phenotypes",
                    "n_associated_phenotype", "num_phenotypes")
GWSIG_CANDIDATES = ("most_gwsig_in_clump", "n_gwsig_in_clump",
                    "max_gwsig_in_clump")


def load_counts(xlsx, sheet):
    from openpyxl import load_workbook
def load_counts(xlsx, sheet, min_gwsig=1):
    from openpyxl import load_workbook
    ws = load_workbook(xlsx, read_only=True, data_only=True)[sheet]
    rows = list(ws.iter_rows(values_only=True))
    hdr = [str(h) if h is not None else "" for h in rows[0]]

    # optional GW-sig filter column
    gcol = next((hdr.index(c) for c in GWSIG_CANDIDATES if c in hdr), None)
    if min_gwsig > 1 and gcol is None:
        sys.exit(f"ERROR: --min-gwsig needs a column {GWSIG_CANDIDATES} in "
                 f"'{sheet}'; not found (have {hdr[:5]}...)")

    def passes(r):
        if min_gwsig <= 1:
            return True
        return (gcol < len(r) and isinstance(r[gcol], (int, float))
                and r[gcol] >= min_gwsig)

    # prefer an explicit count column
    col = next((hdr.index(c) for c in COUNT_CANDIDATES if c in hdr), None)
    counts, n_total, n_kept = [], 0, 0
    if col is not None:
        for r in rows[1:]:
            if col < len(r) and isinstance(r[col], (int, float)):
                n_total += 1
                if passes(r):
                    counts.append(int(r[col])); n_kept += 1
        src = hdr[col]
    else:
        # fall back to counting non-empty phenotype_k cells
        pidx = [i for i, h in enumerate(hdr) if h.startswith("phenotype_")]
        if not pidx:
            sys.exit(f"ERROR: no count column {COUNT_CANDIDATES} and no "
                     f"phenotype_* columns in '{sheet}'")
        for r in rows[1:]:
            if not r or not r[0]:
                continue
            n = sum(1 for i in pidx if i < len(r) and r[i]
                    and str(r[i]).strip())
            if n:
                n_total += 1
                if passes(r):
                    counts.append(n); n_kept += 1
        src = f"count of phenotype_* cells ({len(pidx)} slots)"
    if min_gwsig > 1:
        src += f"  [filtered to >={min_gwsig} GW-sig in clump: {n_kept}/{n_total} loci]"
    return counts, src


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("xlsx")
    ap.add_argument("--sheet", default=SHEET)
    ap.add_argument("--out", default=None)
    ap.add_argument("--title",
                    default="Distribution of cross-disease loci by number of phenotypes")
    ap.add_argument("--notes", nargs="*", default=None,
                    help="annotation bullets to print below the plot")
    ap.add_argument("--min-gwsig", type=int, default=1,
                    help="keep only loci whose richest clump has >= this many "
                         "GW-sig variants (uses most_gwsig_in_clump column)")
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()

    counts, src = load_counts(args.xlsx, args.sheet, args.min_gwsig)
    if not counts:
        sys.exit("ERROR: no per-locus phenotype counts found")
    dist = Counter(counts)
    lo, hi = min(dist), max(dist)
    xs = list(range(lo, hi + 1))
    ys = [dist.get(x, 0) for x in xs]
    total = sum(ys)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    notes = args.notes
    # reserve a caption panel at the bottom if there are notes (or always leave a gap)
    n_note_lines = len(notes) if notes else 3
    fig_h = 5.0 + 0.28 * n_note_lines
    fig, (ax, cap) = plt.subplots(
        2, 1, figsize=(9, fig_h),
        gridspec_kw={"height_ratios": [5, max(1.0, 0.28 * n_note_lines)]})

    bars = ax.bar(xs, ys, color="#3b6ea5", edgecolor="none", width=0.8)
    ax.set_xticks(xs)
    ax.set_xlabel("Number of phenotypes sharing the locus")
    ax.set_ylabel("Number of loci")
    title = args.title
    if args.min_gwsig > 1:
        title += f"\n(loci with \u2265{args.min_gwsig} GW-sig variants in a clump)"
    ax.set_title(title, fontweight="bold")
    ymax = max(ys)
    for b, y in zip(bars, ys):
        if y:
            ax.text(b.get_x() + b.get_width() / 2, y + ymax * 0.012,
                    str(y), ha="center", va="bottom", fontsize=8)
    ax.set_ylim(0, ymax * 1.10)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # caption panel
    cap.axis("off")
    if notes:
        txt = "\n".join(f"\u2022 {n}" for n in notes)
        cap.text(0.01, 0.95, txt, va="top", ha="left", fontsize=9,
                 wrap=True, family="sans-serif")

    fig.tight_layout()
    out = args.out or (args.xlsx.rsplit(".", 1)[0] + "_pheno_distribution.png")
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight")

    # also dump the numbers so you can write annotations from them
    sys.stderr.write(f"Source column: {src}\n")
    sys.stderr.write(f"Wrote {out}\nDistribution ({total} loci):\n")
    for x in xs:
        sys.stderr.write(f"  {x} phenotypes: {dist.get(x,0)}\n")


if __name__ == "__main__":
    main()
