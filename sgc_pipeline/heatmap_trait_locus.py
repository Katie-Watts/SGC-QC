#!/usr/bin/env python3
"""
Heatmap of trait x locus, coloured by association strength (-log10 p, capped).

Pipeline:
  1. parse every *.verbose.clumped, keep p<SIG variants
  2. build loci as union-find components of clumps sharing a GW-sig variant
     (same construction as cross_disease_verbose.py); MHC collapsed to one locus
  3. KEEP ONLY loci where at least one clump has >= MIN_GWSIG (default 5) GW-sig
     variants -- the "qualifying" loci
  4. for each (trait, locus) take the trait's best (smallest) p anywhere in that
     locus, across all its strata; colour = -log10(p), capped at --cap
  5. draw the heatmap: traits (rows) x loci (columns)

All association p-values come straight from the verbose files -- no summary-stat
scan needed, since the clumped output already carries them.

Traits = phenotypes (strata pooled). A cell is blank (shown grey) where the
trait has no GW-sig variant in that locus.

Usage:
  python3 heatmap_trait_locus.py \\
      --verbose-dir clumped_1mb/results \\
      --out trait_locus_heatmap.png \\
      --min-gwsig 5 --cap 50

  # order/limit:
  --min-traits N   only show loci associated with >= N traits (default 1)
  --cap X          -log10(p) colour ceiling (default 50)
  --label-loci     write locus variant ids as x labels (busy; off by default)
"""

import argparse
import math
import os
import re
import sys
from collections import defaultdict

SUMMARY_HDR = re.compile(r"^\s*CHR\s+F\s+SNP\s+BP\s+P\s+TOTAL")
INDEX_LINE = re.compile(r"^\s*\(INDEX\)")
MEMBER_LINE = re.compile(r"^\s+[0-9XYMT]+:[0-9]+:")
CHR_TOKEN = re.compile(r"^[0-9XYMT]+$")
KNOWN_STRATA = {"ALL", "EUR", "AFR", "AMR", "EAS", "SAS", "MALE", "FEMALE"}
MHC_CHR, MHC_START, MHC_END = "6", 25_000_000, 34_000_000


def norm_chr(c):
    c = re.sub(r"^chr", "", str(c).strip(), flags=re.I).upper()
    return {"23": "X", "24": "Y", "25": "X", "26": "MT", "M": "MT"}.get(c, c)


def to_float(t):
    try:
        return float(t)
    except (TypeError, ValueError):
        return None


def split_pheno_stratum(name):
    m = re.match(r"^(.*)_([A-Za-z]+)$", name)
    if m and m.group(2).upper() in KNOWN_STRATA:
        return m.group(1), m.group(2).upper()
    return name, ""


def in_mhc(v):
    p = str(v).split(":")
    if len(p) < 2:
        return False
    try:
        pos = int(p[1])
    except ValueError:
        return False
    return norm_chr(p[0]) == MHC_CHR and MHC_START <= pos <= MHC_END


def parse_verbose(path, sig):
    """clumps: list of dict(index, index_p, sig={variant: p})."""
    clumps, cur, expect = [], None, False
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if SUMMARY_HDR.match(line):
                expect = True
                continue
            f = line.split()
            if expect and f and CHR_TOKEN.match(f[0]):
                if cur:
                    clumps.append(cur)
                cur, expect = None, False
                if len(f) < 5:
                    continue
                p = to_float(f[4])
                cur = {"index": f[2], "index_p": p if p is not None else 1.0,
                       "sig": {}}
                if p is not None and p < sig:
                    cur["sig"][f[2]] = p
                continue
            if INDEX_LINE.match(line):
                continue
            if cur is not None and MEMBER_LINE.match(line):
                if not f:
                    continue
                p = to_float(f[-1])
                if p is None or p >= sig:
                    continue
                prev = cur["sig"].get(f[0])
                if prev is None or p < prev:
                    cur["sig"][f[0]] = p
    if cur:
        clumps.append(cur)
    return clumps


def build(verbose_dir, sig, min_gwsig):
    files = sorted(f for f in os.listdir(verbose_dir)
                   if f.endswith(".verbose.clumped"))
    all_clumps = []
    for fn in files:
        pheno, stratum = split_pheno_stratum(fn[:-len(".verbose.clumped")])
        for c in parse_verbose(os.path.join(verbose_dir, fn), sig):
            if c["sig"]:
                all_clumps.append({"pheno": pheno, "stratum": stratum,
                                   "index": c["index"], "index_p": c["index_p"],
                                   "sig": c["sig"], "n_gwsig": len(c["sig"])})
    sys.stderr.write(f"Parsed {len(files)} files -> {len(all_clumps)} clumps\n")

    parent = list(range(len(all_clumps)))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    v2c = defaultdict(list)
    for i, c in enumerate(all_clumps):
        for v in c["sig"]:
            v2c[v].append(i)
    for idxs in v2c.values():
        for j in idxs[1:]:
            union(idxs[0], j)

    comp = defaultdict(list)
    for i in range(len(all_clumps)):
        comp[find(i)].append(i)

    loci = []
    for members in comp.values():
        rep_count, rep_p = defaultdict(int), {}
        trait_p = {}            # pheno -> best p in this component
        max_ngw = 0
        for i in members:
            c = all_clumps[i]
            rep_count[c["index"]] += 1
            if c["index"] not in rep_p or c["index_p"] < rep_p[c["index"]]:
                rep_p[c["index"]] = c["index_p"]
            max_ngw = max(max_ngw, c["n_gwsig"])
            for v, p in c["sig"].items():
                if c["pheno"] not in trait_p or p < trait_p[c["pheno"]]:
                    trait_p[c["pheno"]] = p
        rep = sorted(rep_count, key=lambda v: (-rep_count[v], rep_p[v]))[0]
        loci.append({"rep": rep, "trait_p": trait_p, "max_ngw": max_ngw,
                     "mhc": in_mhc(rep)})

    # collapse MHC
    mhc = [l for l in loci if l["mhc"]]
    rest = [l for l in loci if not l["mhc"]]
    if mhc:
        trait_p, max_ngw = {}, 0
        rc, rp = defaultdict(int), {}
        for l in mhc:
            max_ngw = max(max_ngw, l["max_ngw"])
            for ph, p in l["trait_p"].items():
                if ph not in trait_p or p < trait_p[ph]:
                    trait_p[ph] = p
            rc[l["rep"]] += 1
            rp[l["rep"]] = min(rp.get(l["rep"], 1.0), min(l["trait_p"].values()))
        rep = sorted(rc, key=lambda v: (-rc[v], rp[v]))[0]
        rest.append({"rep": rep, "trait_p": trait_p, "max_ngw": max_ngw,
                     "mhc": True})
        sys.stderr.write(f"MHC collapsed {len(mhc)} loci into 1\n")

    # keep only qualifying loci: some clump had >= min_gwsig GW-sig variants
    qual = [l for l in rest if l["max_ngw"] >= min_gwsig]
    sys.stderr.write(f"{len(rest)} loci total; {len(qual)} qualify "
                     f"(>= {min_gwsig} GW-sig variants in a clump)\n")
    return qual


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--verbose-dir", required=True)
    ap.add_argument("--out", default="trait_locus_heatmap.png")
    ap.add_argument("--sig", type=float, default=5e-8)
    ap.add_argument("--min-gwsig", type=int, default=5)
    ap.add_argument("--min-traits", type=int, default=1)
    ap.add_argument("--cap", type=float, default=50.0,
                    help="-log10(p) colour ceiling (default 50)")
    ap.add_argument("--no-cluster", action="store_true",
                    help="keep genomic order instead of clustering rows/cols")
    ap.add_argument("--no-xlabels", action="store_true",
                    help="hide locus x-axis labels (busy with many loci)")
    ap.add_argument("--no-grid", action="store_true",
                    help="omit the gridlines between cells")
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()

    loci = build(args.verbose_dir, args.sig, args.min_gwsig)
    if args.min_traits > 1:
        loci = [l for l in loci if len(l["trait_p"]) >= args.min_traits]
        sys.stderr.write(f"{len(loci)} loci with >= {args.min_traits} traits\n")
    if not loci:
        sys.exit("No qualifying loci to plot.")

    traits = sorted({t for l in loci for t in l["trait_p"]})

    # genomic order is the default column order (used when not clustering)
    def sortkey(l):
        p = str(l["rep"]).split(":")
        c = norm_chr(p[0])
        cn = {"X": 23, "Y": 24, "MT": 25}.get(c, 999)
        try: cn = int(c)
        except ValueError: pass
        try: pos = int(p[1])
        except (ValueError, IndexError): pos = 0
        return (cn, pos)
    loci.sort(key=sortkey)

    ti = {t: i for i, t in enumerate(traits)}
    import numpy as np
    M = np.full((len(traits), len(loci)), np.nan)
    for j, l in enumerate(loci):
        for t, p in l["trait_p"].items():
            val = -math.log10(p) if p > 0 else args.cap
            M[ti[t], j] = min(val, args.cap)

    # ---- cluster rows (traits) and columns (loci) ----
    ordering = "genomic position"
    if not args.no_cluster and len(traits) > 2 and len(loci) > 2:
        try:
            from scipy.cluster.hierarchy import linkage, leaves_list
            from scipy.spatial.distance import pdist
            # fill NaN with 0 (no association) for the distance computation only
            F = np.nan_to_num(M, nan=0.0)
            # rows
            row_order = leaves_list(linkage(pdist(F, metric="euclidean"),
                                            method="average"))
            # cols
            col_order = leaves_list(linkage(pdist(F.T, metric="euclidean"),
                                            method="average"))
            M = M[np.ix_(row_order, col_order)]
            traits = [traits[i] for i in row_order]
            loci = [loci[i] for i in col_order]
            ordering = "hierarchical clustering"
        except Exception as e:
            sys.stderr.write(f"WARN clustering failed ({e}); keeping genomic order\n")
    sys.stderr.write(f"Axis order: {ordering}\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad(color="#e8e8e8")           # non-significant / missing = light grey

    # x labels on by default now; wider figure to fit them
    fig_w = max(8, 0.24 * len(loci) + 3)
    fig_h = max(5, 0.30 * len(traits) + 2.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(M, aspect="auto", cmap=cmap, vmin=0, vmax=args.cap,
                   interpolation="nearest")

    ax.set_yticks(range(len(traits)))
    ax.set_yticklabels(traits, fontsize=7)
    ax.set_xticks(range(len(loci)))
    ax.set_xticklabels([l["rep"] for l in loci], fontsize=6, rotation=90)
    ax.set_xlabel(f"{len(loci)} loci (\u2265{args.min_gwsig} GW-sig variants in a "
                  f"clump), ordered by {ordering}", fontsize=10)
    if args.no_xlabels:
        ax.set_xticks([])
    ax.set_title(f"Trait \u00d7 locus association  (-log10 p, capped at {args.cap:g}; "
                 f"grey = not significant)",
                 fontweight="bold", fontsize=12)

    # thin gridlines between cells (minor ticks at half-offsets)
    if not args.no_grid:
        import numpy as np
        ax.set_xticks(np.arange(-0.5, len(loci), 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(traits), 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=0.5)
        ax.tick_params(which="minor", length=0)

    cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
    cb.set_label(f"-log10(p), capped at {args.cap:g}")
    fig.tight_layout()
    fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
    sys.stderr.write(f"Wrote {args.out}  ({len(traits)} traits x {len(loci)} loci)\n")


if __name__ == "__main__":
    main()
