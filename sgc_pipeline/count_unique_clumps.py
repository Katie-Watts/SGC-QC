#!/usr/bin/env python3
"""
Count unique GWAS clumps per phenotype, summed across strata, WITHOUT double
counting the same locus when more than one stratum finds it.

Rule
----
Within a phenotype, two clumps (from any of its strata) are the SAME locus if
they share at least one variant. Clumps are linked transitively (single-linkage
union-find), so a locus is a connected component of clumps. The phenotype's
count is its number of components: a locus found in three strata counts once, a
locus found in one stratum counts once, and they sum.

    <PHENOTYPE>_<STRATUM>.verbose.clumped   ->   phenotype, count

What "share a variant" means is controlled by --link:
    gwsig  (default) link on shared variants with p < SIG. Matches the linking
           rule in cross_disease_verbose.py. A clump with no GW-sig variant
           still counts as its own locus but cannot merge with anything.
    any    link on ANY shared clumped variant (index or member), regardless of
           p. Merges more aggressively.

Variant ids are chr:pos:ref:alt. MHC (chr6:25-34Mb, hg38) is NOT collapsed here
unless you pass --collapse-mhc (then all MHC clumps in a phenotype count as one
locus, matching verbose_summary.sh's single-representative convention).

By default a clump must have >= 2 GW-sig (p<SIG) variants to count at all;
single-hit clumps are dropped before linking (change with --min-gwsig).

Output (TSV): <out>
    phenotype   n_loci   n_clumps_raw   strata   per_stratum_raw
  n_loci          unique loci after cross-strata merge   <- the number you want
  n_clumps_raw    kept clumps before merging (sum over strata), i.e. after the
                  --min-gwsig filter
  strata          which strata contributed
  per_stratum_raw kept clump count in each, e.g. ALL:12; EUR:9; ...

Usage
-----
  python3 count_unique_clumps.py --verbose-dir clumped_1mb/results \\
      --out clumps_per_phenotype.tsv

  # link on any shared variant, and collapse the MHC to one locus
  python3 count_unique_clumps.py --verbose-dir clumped_1mb/results \\
      --out clumps_per_phenotype.tsv --link any --collapse-mhc
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
STRAT_ORD = {s: i for i, s in enumerate(
    ["ALL", "EUR", "AFR", "AMR", "EAS", "SAS", "MALE", "FEMALE"])}
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


def parse_verbose(path, sig, link_any, min_gwsig):
    """Return a list of clumps; each clump is the SET of variant ids used for
    linking. With link_any=False only p<sig variants are kept for linking, but
    a clump with none still appears (as an empty-linkable set) so it is counted.
    Also returns whether the clump touches the MHC (any variant in region).

    n_gwsig (count of p<sig variants) is tracked per clump independently of the
    link set, so a clump can be dropped when it has fewer than min_gwsig GW-sig
    variants even under --link any (where the link set also holds sub-threshold
    members)."""
    clumps = []          # list of (set_of_link_variants, touches_mhc)
    cur_link, cur_mhc, cur_ngw, have = set(), False, 0, False
    expect = False

    def close():
        nonlocal cur_link, cur_mhc, cur_ngw, have
        if have and cur_ngw >= min_gwsig:
            clumps.append((cur_link, cur_mhc))
        cur_link, cur_mhc, cur_ngw, have = set(), False, 0, False

    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if SUMMARY_HDR.match(line):
                expect = True
                continue
            f = line.split()
            if expect and f and CHR_TOKEN.match(f[0]):
                close()
                expect = False
                if len(f) < 5:
                    continue
                have = True
                snp = f[2]
                p = to_float(f[4])
                is_gw = p is not None and p < sig
                if is_gw:
                    cur_ngw += 1
                if link_any or is_gw:
                    cur_link.add(snp)
                if in_mhc(snp):
                    cur_mhc = True
                continue
            if INDEX_LINE.match(line):
                continue
            if have and MEMBER_LINE.match(line):
                if not f:
                    continue
                snp = f[0]
                p = to_float(f[-1])
                if p is None:
                    continue          # e.g. 'not found in dataset'
                is_gw = p < sig
                if is_gw:
                    cur_ngw += 1
                if link_any or is_gw:
                    cur_link.add(snp)
                if in_mhc(snp):
                    cur_mhc = True
                continue
    close()
    return clumps


def count_phenotype(clumps, collapse_mhc):
    """clumps: list of (link_set, touches_mhc). Return n_loci."""
    n = len(clumps)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    var_to = defaultdict(list)
    for i, (vs, _) in enumerate(clumps):
        for v in vs:
            var_to[v].append(i)
    for idxs in var_to.values():
        for j in idxs[1:]:
            union(idxs[0], j)

    if collapse_mhc:
        mhc = [i for i, (_, m) in enumerate(clumps) if m]
        for j in mhc[1:]:
            union(mhc[0], j)

    return len({find(i) for i in range(n)})


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--verbose-dir", required=True)
    ap.add_argument("--out", default="clumps_per_phenotype.tsv")
    ap.add_argument("--sig", type=float, default=5e-8)
    ap.add_argument("--link", choices=["gwsig", "any"], default="gwsig",
                    help="link clumps on shared GW-sig variants (default) or "
                         "on any shared variant")
    ap.add_argument("--min-gwsig", type=int, default=2,
                    help="drop any clump with fewer than this many GW-sig "
                         "(p<SIG) variants before counting/linking (default 2, "
                         "i.e. exclude single-hit clumps; set 1 to keep them)")
    ap.add_argument("--collapse-mhc", action="store_true",
                    help="count all MHC clumps in a phenotype as one locus")
    ap.add_argument("--exclude", nargs="*", default=[],
                    help="verbose file stems to skip (rarely needed)")
    args = ap.parse_args()

    files = sorted(f for f in os.listdir(args.verbose_dir)
                   if f.endswith(".verbose.clumped"))
    if not files:
        sys.exit(f"ERROR: no *.verbose.clumped in {args.verbose_dir}")

    # phenotype -> list of clumps (pooled across strata); and raw per-stratum
    pheno_clumps = defaultdict(list)
    pheno_raw = defaultdict(lambda: defaultdict(int))
    link_any = args.link == "any"
    excluded = set(args.exclude)

    for fn in files:
        base = fn[:-len(".verbose.clumped")]
        if base in excluded:
            continue
        pheno, stratum = split_pheno_stratum(base)
        clumps = parse_verbose(os.path.join(args.verbose_dir, fn),
                               args.sig, link_any, args.min_gwsig)
        pheno_clumps[pheno].extend(clumps)
        pheno_raw[pheno][stratum] += len(clumps)

    rows = []
    for pheno in sorted(pheno_clumps):
        clumps = pheno_clumps[pheno]
        n_loci = count_phenotype(clumps, args.collapse_mhc)
        raw = pheno_raw[pheno]
        strata = sorted(raw, key=lambda s: STRAT_ORD.get(s, 99))
        rows.append([pheno, n_loci, sum(raw.values()), strata,
                     "; ".join(f"{s}:{raw[s]}" for s in strata)])

    with open(args.out, "w") as out:
        out.write("phenotype\tn_loci\tn_clumps_raw\tstrata\tper_stratum_raw\n")
        for pheno, n_loci, raw_tot, strata, per in rows:
            out.write(f"{pheno}\t{n_loci}\t{raw_tot}\t{';'.join(strata)}\t{per}\n")

    tot_loci = sum(r[1] for r in rows)
    tot_raw = sum(r[2] for r in rows)
    sys.stderr.write(
        f"Wrote {args.out}: {len(rows)} phenotypes\n"
        f"Link mode: {args.link}; min GW-sig per clump: {args.min_gwsig}"
        + (", MHC collapsed" if args.collapse_mhc else "") + "\n"
        f"Total unique loci: {tot_loci}  (from {tot_raw} kept clumps; "
        f"{tot_raw - tot_loci} merged as cross-strata duplicates)\n")


if __name__ == "__main__":
    main()
