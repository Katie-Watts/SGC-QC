#!/usr/bin/env python3
"""
Highlights plot: the loci that span the MOST phenotypes, restricted to loci
that contain a clump with >= MIN_GWSIG (default 5) GW-sig variants.

Builds loci fresh from the *.verbose.clumped files (union-find over clumps
sharing a GW-sig variant, MHC collapsed -- same construction as
cross_disease_verbose.py), counts how many distinct phenotypes each locus is
associated with, keeps only loci whose richest clump has >= MIN_GWSIG GW-sig
variants, and draws a horizontal bar chart of the top N by phenotype count.

Each bar is one locus, labelled by its representative variant; bar length is the
number of phenotypes; the richest clump's GW-sig count is annotated at the end.

Usage:
  python3 highlights_top_loci.py \\
      --verbose-dir clumped_1mb/results \\
      --out highlights_top_loci.png \\
      --min-gwsig 5 --top 30

  --min-traits N   also require >= N phenotypes (default 2; cross-disease)
  --list           also print the ranked loci and their phenotypes to stderr
"""

import argparse
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


def build(verbose_dir, sig):
    files = sorted(f for f in os.listdir(verbose_dir)
                   if f.endswith(".verbose.clumped"))
    all_clumps = []
    for fn in files:
        pheno, stratum = split_pheno_stratum(fn[:-len(".verbose.clumped")])
        for c in parse_verbose(os.path.join(verbose_dir, fn), sig):
            if c["sig"]:
                all_clumps.append({"pheno": pheno, "index": c["index"],
                                   "index_p": c["index_p"], "sig": c["sig"],
                                   "n_gwsig": len(c["sig"])})
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
        rc, rp = defaultdict(int), {}
        traits = {}
        max_ngw = 0
        for i in members:
            c = all_clumps[i]
            rc[c["index"]] += 1
            if c["index"] not in rp or c["index_p"] < rp[c["index"]]:
                rp[c["index"]] = c["index_p"]
            max_ngw = max(max_ngw, c["n_gwsig"])
            for v, p in c["sig"].items():
                if c["pheno"] not in traits or p < traits[c["pheno"]]:
                    traits[c["pheno"]] = p
        rep = sorted(rc, key=lambda v: (-rc[v], rp[v]))[0]
        loci.append({"rep": rep, "traits": traits, "max_ngw": max_ngw,
                     "mhc": in_mhc(rep)})

    mhc = [l for l in loci if l["mhc"]]
    rest = [l for l in loci if not l["mhc"]]
    if mhc:
        traits, max_ngw = {}, 0
        rc, rp = defaultdict(int), {}
        for l in mhc:
            max_ngw = max(max_ngw, l["max_ngw"])
            for ph, p in l["traits"].items():
                if ph not in traits or p < traits[ph]:
                    traits[ph] = p
            rc[l["rep"]] += 1
            rp[l["rep"]] = min(rp.get(l["rep"], 1.0), min(l["traits"].values()))
        rep = sorted(rc, key=lambda v: (-rc[v], rp[v]))[0]
        rest.append({"rep": rep, "traits": traits, "max_ngw": max_ngw,
                     "mhc": True})
        sys.stderr.write(f"MHC collapsed {len(mhc)} loci into 1\n")
    return rest


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--verbose-dir", required=True)
    ap.add_argument("--out", default="highlights_top_loci.png")
    ap.add_argument("--sig", type=float, default=5e-8)
    ap.add_argument("--min-gwsig", type=int, default=5)
    ap.add_argument("--min-traits", type=int, default=2)
    ap.add_argument("--top", type=int, default=30)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()

    loci = build(args.verbose_dir, args.sig)
    total = len(loci)
    # filter: richest clump has >= min_gwsig GW-sig variants
    loci = [l for l in loci if l["max_ngw"] >= args.min_gwsig]
    after_gw = len(loci)
    loci = [l for l in loci if len(l["traits"]) >= args.min_traits]
    sys.stderr.write(
        f"{total} loci -> {after_gw} with >= {args.min_gwsig} GW-sig variants "
        f"-> {len(loci)} also spanning >= {args.min_traits} phenotypes\n")
    if not loci:
        sys.exit("No loci pass the filter.")

    loci.sort(key=lambda l: (len(l["traits"]), l["max_ngw"]), reverse=True)
    shown = loci[:args.top]

    if args.list:
        sys.stderr.write("\nRanked loci:\n")
        for l in shown:
            phs = sorted(l["traits"], key=lambda t: l["traits"][t])
            sys.stderr.write(f"  {l['rep']:<20} {len(l['traits'])} phenotypes, "
                             f"{l['max_ngw']} GW-sig: {', '.join(phs[:8])}"
                             f"{' ...' if len(phs) > 8 else ''}\n")

    shown.reverse()   # largest at top after barh
    labels = [l["rep"] for l in shown]
    counts = [len(l["traits"]) for l in shown]
    ngw = [l["max_ngw"] for l in shown]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    h = max(4, 0.32 * len(shown) + 1.2)
    fig, ax = plt.subplots(figsize=(10, h))
    bars = ax.barh(range(len(shown)), counts, color="#c1553b", edgecolor="none")
    ax.set_yticks(range(len(shown)))
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Number of phenotypes", fontsize=12)
    ax.tick_params(axis="x", labelsize=10)
    ax.set_title(
        f"Loci spanning the most phenotypes "
        f"(\u2265{args.min_gwsig} GW-sig variants in a clump)",
        fontweight="bold", fontsize=13)

    xmax = max(counts)
    for b, c, g in zip(bars, counts, ngw):
        ax.text(b.get_width() + xmax * 0.01, b.get_y() + b.get_height() / 2,
                f"{c}  ({g} GW-sig)", va="center", ha="left", fontsize=8)
    ax.set_xlim(0, xmax * 1.18)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
    sys.stderr.write(f"Wrote {args.out}  (showing top {len(shown)} of "
                     f"{len(loci)} passing loci)\n")


if __name__ == "__main__":
    main()
