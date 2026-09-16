#!/usr/bin/env python3
"""
Build a cross-disease locus table from verbose clumping results.

A LOCUS is defined by the clumps themselves: two clumps (from different
phenotype_stratum files) belong to the same locus if they share ANY
genome-wide-significant variant (p < SIG; index or member). Clumps are linked
transitively (union-find), so a locus is a connected component of clumps.

A cross-disease locus is a component spanning >= 2 distinct PHENOTYPES
(phenotype = file name with the stratum suffix stripped).

Output (tab-separated), matching the prior cross-disease format:
    Variant <TAB> GWAS
where Variant is the representative (the variant that is a clump INDEX in the
most phenotypes at that locus; ties -> lowest p), and GWAS is
    PHENO1 (STRAT; STRAT; ...); PHENO2 (STRAT; ...); ...
listing, per phenotype, the strata whose clumps fall in that locus.

Usage:
    python3 cross_disease_verbose.py <verbose_dir> <out.tsv> [sig] [strata_csv]//
      strata_csv (optional) restricts which suffixes count, e.g. ALL,EUR,EAS,...
"""
import sys, os, glob, re
from collections import defaultdict

verbose_dir = sys.argv[1] if len(sys.argv) > 1 else "."
out_path    = sys.argv[2] if len(sys.argv) > 2 else "cross_disease.tsv"
SIG         = float(sys.argv[3]) if len(sys.argv) > 3 else 5e-8

KNOWN_STRATA = {"ALL","EUR","AFR","AMR","EAS","SAS","MALE","FEMALE"}

def split_pheno_stratum(name):
    # name like ATOPIC_DERM_AFR -> (ATOPIC_DERM, AFR); fall back to whole name
    m = re.match(r"^(.*)_([A-Za-z]+)$", name)
    if m and m.group(2).upper() in KNOWN_STRATA:
        return m.group(1), m.group(2).upper()
    return name, ""

# ---- parse one verbose file into a list of clumps ----
# each clump: {"index": snp, "index_p": float, "sig": {variant: p} for p<SIG}
def parse_verbose(path):
    clumps = []
    cur = None
    expect_summary = False
    with open(path) as fh:
        for line in fh:
            if re.match(r"^\s*CHR\s+F\s+SNP\s+BP\s+P\s+TOTAL", line):
                expect_summary = True
                continue
            f = line.split()
            if expect_summary and f and re.match(r"^[0-9XY]+$", f[0]):
                # summary line: CHR F SNP BP P TOTAL ...
                if cur: clumps.append(cur)
                snp = f[2]; p = float(f[4])
                cur = {"index": snp, "index_p": p, "sig": {}}
                if p < SIG: cur["sig"][snp] = p
                expect_summary = False
                continue
            if re.match(r"^\s*\(INDEX\)", line):
                continue  # redundant with summary line
            if cur and re.match(r"^\s+[0-9XY]+:[0-9]+:", line):
                snp = f[0]
                try: p = float(f[-1])
                except ValueError: p = 1.0
                if p < SIG:
                    # keep the best (lowest) p if a variant recurs in a clump
                    if snp not in cur["sig"] or p < cur["sig"][snp]:
                        cur["sig"][snp] = p
                continue
    if cur: clumps.append(cur)
    return clumps

# ---- gather all clumps across all files ----
all_clumps = []   # list of dicts with pheno, stratum, index, index_p, sig
files = sorted(glob.glob(os.path.join(verbose_dir, "*.verbose.clumped")))
if not files:
    sys.stderr.write(f"No *.verbose.clumped files in {verbose_dir}\n"); sys.exit(1)

for path in files:
    base = os.path.basename(path)[:-len(".verbose.clumped")]
    pheno, stratum = split_pheno_stratum(base)
    for c in parse_verbose(path):
        if not c["sig"]:
            continue  # no GW-sig variant -> can't link cross-disease
        all_clumps.append({
            "pheno": pheno, "stratum": stratum,
            "index": c["index"], "index_p": c["index_p"], "sig": c["sig"],
        })

# ---- union-find over clumps, linked by shared GW-sig variants ----
parent = list(range(len(all_clumps)))
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb: parent[ra] = rb

# map variant -> clump indices that contain it (as GW-sig), then union them
var_to_clumps = defaultdict(list)
for i, c in enumerate(all_clumps):
    for v in c["sig"]:
        var_to_clumps[v].append(i)
for v, idxs in var_to_clumps.items():
    for j in idxs[1:]:
        union(idxs[0], j)

# ---- assemble components ----
comp = defaultdict(list)
for i in range(len(all_clumps)):
    comp[find(i)].append(i)

rows = []
for root, members in comp.items():
    phenos = {}      # pheno -> set(strata)
    pheno_bestp = {} # pheno -> best (lowest) p among its GW-sig clump variants here
    for i in members:
        c = all_clumps[i]
        phenos.setdefault(c["pheno"], set())
        if c["stratum"]: phenos[c["pheno"]].add(c["stratum"])
        # best p for this phenotype = min p over all its sig variants in the locus
        cmin = min(c["sig"].values())
        if c["pheno"] not in pheno_bestp or cmin < pheno_bestp[c["pheno"]]:
            pheno_bestp[c["pheno"]] = cmin
    if len(phenos) < 2:
        continue  # not cross-disease

    # representative variant = the variant that is a clump INDEX most often
    # across the clumps in this locus (most common index); tie-break lowest p.
    idx_count = defaultdict(int)
    idx_best_p = {}
    for i in members:
        c = all_clumps[i]
        idx_count[c["index"]] += 1
        if c["index"] not in idx_best_p or c["index_p"] < idx_best_p[c["index"]]:
            idx_best_p[c["index"]] = c["index_p"]
    rep = sorted(idx_count.keys(), key=lambda v: (-idx_count[v], idx_best_p[v]))[0]

    # n_gwsig for the row = number of GW-sig variants in the clump containing the
    # representative variant. If the rep appears in several clumps (across
    # phenotypes/strata), use the clump where the rep's own p is lowest.
    rep_clump_gwsig = 0
    rep_best_p_in_clump = None
    for i in members:
        c = all_clumps[i]
        if rep in c["sig"]:
            rp = c["sig"][rep]
            if rep_best_p_in_clump is None or rp < rep_best_p_in_clump:
                rep_best_p_in_clump = rp
                rep_clump_gwsig = len(c["sig"])

    # row carries: rep variant, structured {pheno: (strata_set, best_p)}
    pdata = {p: (phenos[p], pheno_bestp[p]) for p in phenos}
    rows.append({"rep": rep, "pdata": pdata, "n_gwsig": rep_clump_gwsig,
                 "npheno": len(phenos), "best_p": min(pheno_bestp.values())})

# ---- collapse MHC: merge all cross-disease loci whose representative variant
# falls in chr6:25-34Mb (hg38) into a SINGLE row. The merged row unions the
# phenotypes/strata (best p per phenotype = min across merged loci); the
# representative is the most common representative across the merged MHC loci. ----
MHC_CHR, MHC_START, MHC_END = "6", 25_000_000, 34_000_000
def in_mhc(variant):
    parts = variant.split(":")
    if len(parts) < 2:
        return False
    try:
        pos = int(parts[1])
    except ValueError:
        return False
    return parts[0] == MHC_CHR and MHC_START <= pos <= MHC_END

mhc_rows = [r for r in rows if in_mhc(r["rep"])]
non_mhc  = [r for r in rows if not in_mhc(r["rep"])]

if mhc_rows:
    merged = {}   # pheno -> [strata_set, best_p]
    for r in mhc_rows:
        for p, (strata, bp) in r["pdata"].items():
            if p not in merged:
                merged[p] = [set(strata), bp]
            else:
                merged[p][0] |= strata
                merged[p][1] = min(merged[p][1], bp)
    rep_count = defaultdict(int); rep_bestp = {}
    for r in mhc_rows:
        rep_count[r["rep"]] += 1
        if r["rep"] not in rep_bestp or r["best_p"] < rep_bestp[r["rep"]]:
            rep_bestp[r["rep"]] = r["best_p"]
    mhc_rep = sorted(rep_count, key=lambda v: (-rep_count[v], rep_bestp[v]))[0]
    # n_gwsig for the collapsed row = the n_gwsig of the MHC row whose rep is the
    # chosen mhc_rep, taking the largest (richest clump) if several share it.
    mhc_ngwsig = max((r["n_gwsig"] for r in mhc_rows if r["rep"] == mhc_rep),
                     default=0)
    pdata = {p: (s, bp) for p, (s, bp) in merged.items()}
    collapsed = {"rep": mhc_rep, "pdata": pdata, "n_gwsig": mhc_ngwsig,
                 "npheno": len(merged), "best_p": min(v[1] for v in merged.values())}
    rows = non_mhc + [collapsed]
else:
    rows = non_mhc

# sort output: most phenotypes first, then strongest p
rows.sort(key=lambda r: (-r["npheno"], r["best_p"]))

# ---- write wide format: Variant, n_phenotypes, phenotype_1, pvalue_1, ... ----
STRAT_ORD = {"ALL":0,"EUR":1,"AFR":2,"AMR":3,"EAS":4,"SAS":5,"MALE":6,"FEMALE":7}
def order_strata(s):
    return sorted(s, key=lambda x: STRAT_ORD.get(x, 99))

def fmt_p(p):
    return f"{p:.3g}"

# determine the max number of phenotypes to size the header
max_np = max((r["npheno"] for r in rows), default=0)
header = ["Variant", "n_phenotypes", "n_gwsig_in_clump"]
for k in range(1, max_np + 1):
    header += [f"phenotype_{k}", f"pvalue_{k}"]

with open(out_path, "w") as out:
    out.write("\t".join(header) + "\n")
    for r in rows:
        # phenotypes ordered by best p (strongest first)
        items = sorted(r["pdata"].items(), key=lambda kv: kv[1][1])
        cells = [r["rep"], str(r["npheno"]), str(r["n_gwsig"])]
        for pheno, (strata, bp) in items:
            label = f"{pheno} ({'; '.join(order_strata(strata))})" if strata else pheno
            cells += [label, fmt_p(bp)]
        # pad to header width
        cells += [""] * (len(header) - len(cells))
        out.write("\t".join(cells) + "\n")

sys.stderr.write(f"Wrote {out_path}: {len(rows)} cross-disease loci "
                 f"(>=2 phenotypes) from {len(files)} files\n")
