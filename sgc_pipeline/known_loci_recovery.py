#!/usr/bin/env python3
"""
Recovery of previously-known GW-significant loci by the new analysis.

Question answered:
  "Using this approach we recover X% of known GW-sig loci at p<5e-8,
   and Y% at p<1e-5."

Known loci come from the "Previous Known loci- variants" sheet (one row per
known lead variant: Phenotype, Chr, BP, and a Clump start/end window). For each
known variant the script looks up the NEW summary statistics for that phenotype
and takes the p-value at that locus, then reports what fraction of known loci
are recovered at each threshold.

"Recovered" is judged by EXACT position: the new p AT THE KNOWN VARIANT'S
position, across all strata (best/smallest p wins). A known variant whose exact
position is NOT PRESENT in any of the phenotype's stratum files is EXCLUDED from
the percentage entirely -- it is neither recovered nor a miss, just untestable.
So the rate is: of known variants present in our data, what fraction reach the
threshold. (A window best-p is still reported per locus for reference but does
not drive the percentage.)

Stratum: recovery is judged across ALL strata of the matching phenotype -- a
locus counts as recovered if ANY stratum (ALL, EUR, AFR, ... MALE, FEMALE)
reaches the threshold. The reported p for a locus is the BEST (smallest) p over
every stratum file that exists for that phenotype, and stratum_used names which
stratum gave it. The known list's ancestry labels are ignored; only the
phenotype has to match.

Deduplication: a phenotype can list the same locus several times (multiple
source studies). Loci are collapsed to unique (phenotype, chr, rounded BP)
before rates are computed, so X% is over distinct known loci.

Outputs:
  <prefix>_recovery_per_locus.tsv   every known locus, its new p, recovered flags
  <prefix>_recovery_summary.tsv     per-phenotype and overall X% / Y%
  prints the headline X% / Y% to stderr

Usage:
  python3 known_loci_recovery.py \\
      --xlsx SGC_Meta_Analysis_Results.xlsx \\
      --gwas-dir . \\
      --out-prefix recovery \\
      --threads 8
"""

import argparse
import gzip
import io
import os
import re
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

SHEET = "Previous Known loci- variants"
KNOWN_STRATA = ["ALL", "EUR", "AFR", "AMR", "EAS", "SAS", "MALE", "FEMALE"]
# known-loci ancestry label (lowercased) -> our stratum code, for exact-ancestry
# matching in the per-stratum breakdown. Labels not listed (Any Ancestry,
# Multiancestry, Greater Middle Eastern, Other) map to no specific stratum:
# they count only toward the any-ancestry strata below.
ANC_MAP = {
    "european": "EUR",
    "east asian": "EAS",
    "south asian": "SAS",
    "african american or afro-caribbean": "AFR",
    "african unspecified": "AFR",
    "hispanic or latin american": "AMR",
}
# strata judged against known loci of ANY ancestry (no exact-ancestry filter):
# the combined ALL analysis plus the sex splits (which are not ancestry-specific)
ANY_ANCESTRY_STRATA = {"ALL", "MALE", "FEMALE"}
# strata judged against ONLY their exact-ancestry known loci
EXACT_ANCESTRY_STRATA = {"EUR", "EAS", "SAS", "AFR", "AMR"}


def norm_chr(c):
    c = re.sub(r"^chr", "", str(c).strip(), flags=re.I).upper()
    return {"23": "X", "24": "Y", "25": "X", "26": "MT", "M": "MT"}.get(c, c)


def thr_label(t):
    """Compact label for a threshold, for column names: 5e-08 -> '5e-8'."""
    s = f"{t:g}"
    return s.replace("e-0", "e-").replace("e+0", "e+")


def to_float(t):
    try:
        return float(t)
    except (TypeError, ValueError):
        return None


def scan_sumstat(job):
    """Read one gzipped sumstat file once; return {(chr,pos): p} for exact hits
    and, per known locus in this file, the min p within its window."""
    path, stem, exact_keys, windows, cols = job
    c_chr, c_pos, c_p = cols
    exact = {}
    win_best = {i: None for i in range(len(windows))}
    # index windows by chr for quick membership
    by_chr = defaultdict(list)
    for i, (ch, lo, hi) in enumerate(windows):
        by_chr[ch].append((lo, hi, i))
    try:
        with io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8",
                              errors="replace") as fh:
            header = [h.strip().lstrip("\ufeff")
                      for h in fh.readline().rstrip("\r\n").split("\t")]
            idx = {h: i for i, h in enumerate(header)}
            for n in (c_chr, c_pos, c_p):
                if n not in idx:
                    return (stem, None, None, f"column {n!r} missing")
            ic, ip, ipv = idx[c_chr], idx[c_pos], idx[c_p]
            need = max(ic, ip, ipv)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) <= need:
                    continue
                ch = norm_chr(f[ic])
                try:
                    pos = int(f[ip])
                except ValueError:
                    continue
                p = to_float(f[ipv])
                if p is None:
                    continue
                if (ch, pos) in exact_keys:
                    prev = exact.get((ch, pos))
                    if prev is None or p < prev:
                        exact[(ch, pos)] = p
                for lo, hi, wi in by_chr.get(ch, ()):
                    if lo <= pos <= hi:
                        b = win_best[wi]
                        if b is None or p < b:
                            win_best[wi] = p
    except OSError as e:
        return (stem, None, None, f"could not read: {e}")
    return (stem, exact, win_best, None)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xlsx", required=True)
    ap.add_argument("--sheet", default=SHEET)
    ap.add_argument("--gwas-dir", default=".")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--suffix", default=".tsv.gz")
    ap.add_argument("--sig", type=float, default=5e-8,
                    help="the headline threshold reported to stderr (default 5e-8)")
    ap.add_argument("--thresholds", type=float, nargs="+",
                    default=[5e-8, 1e-5, 1e-4, 1e-3, 0.05],
                    help="p-value thresholds to count recovery at (default: "
                         "5e-8 1e-5 1e-4 0.001 0.05)")
    ap.add_argument("--mode", choices=["exact", "window", "both"],
                    default="exact",
                    help="recovery by exact known position, by best p in the "
                         "Clump start-end window, or both (default exact)")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--chr-col", default="chromosome")
    ap.add_argument("--pos-col", default="position")
    ap.add_argument("--pval-col", default="pvalue")
    args = ap.parse_args()

    from openpyxl import load_workbook
    ws = load_workbook(args.xlsx, read_only=True, data_only=True)[args.sheet]
    rows = list(ws.iter_rows(values_only=True))
    hdr = [str(h) if h is not None else "" for h in rows[0]]
    H = {h: i for i, h in enumerate(hdr)}
    for need in ("Phenotype", "Chr", "BP"):
        if need not in H:
            sys.exit(f"ERROR: column {need!r} not in sheet; have {hdr}")
    ci_ph, ci_chr, ci_bp = H["Phenotype"], H["Chr"], H["BP"]
    ci_lo = H.get("Clump start"); ci_hi = H.get("Clump end")
    ci_anc = H.get("GWAS")            # ancestry label column

    # collect distinct known loci. A locus can be reported for several
    # ancestries; we keep the SET of ancestry-strata it maps to.
    known = {}     # (pheno, chr, bp) -> {chr,bp,lo,hi,pheno, anc_strata:set}
    for r in rows[1:]:
        if not r or r[ci_ph] is None or r[ci_chr] is None or r[ci_bp] is None:
            continue
        try:
            ch = norm_chr(r[ci_chr]); bp = int(r[ci_bp])
        except (ValueError, TypeError):
            continue
        pheno = str(r[ci_ph]).strip()
        key = (pheno, ch, bp)
        lo = int(r[ci_lo]) if ci_lo is not None and isinstance(r[ci_lo], (int, float)) else bp
        hi = int(r[ci_hi]) if ci_hi is not None and isinstance(r[ci_hi], (int, float)) else bp
        anc = str(r[ci_anc]).strip().lower() if ci_anc is not None and r[ci_anc] else ""
        st = ANC_MAP.get(anc)          # None if it maps to no specific stratum
        if key not in known:
            known[key] = {"pheno": pheno, "chr": ch, "bp": bp,
                          "lo": min(lo, hi), "hi": max(lo, hi),
                          "anc_strata": set()}
        if st:
            known[key]["anc_strata"].add(st)
    sys.stderr.write(f"{len(known)} distinct known loci across "
                     f"{len(set(k[0] for k in known))} phenotypes\n")

    # For each phenotype, scan EVERY stratum file that exists. A locus is
    # recovered if any stratum reaches the threshold; we keep the best p.
    pheno_keys = defaultdict(list)      # pheno -> list of known-keys
    for key, d in known.items():
        pheno_keys[d["pheno"]].append(key)

    per_file = []                       # (file, stem, stratum, [keys])
    phenos_with_no_file = set()
    for pheno, keys in pheno_keys.items():
        found_any = False
        for st in KNOWN_STRATA:
            f = os.path.join(args.gwas_dir, f"{pheno}_{st}{args.suffix}")
            if os.path.exists(f):
                per_file.append((f, f"{pheno}_{st}", st, keys))
                found_any = True
        if not found_any:
            phenos_with_no_file.add(pheno)
    no_file = [k for k in known if k[0] in phenos_with_no_file]
    sys.stderr.write(f"{len(per_file)} sumstat files to scan across "
                     f"{len(pheno_keys) - len(phenos_with_no_file)} phenotypes; "
                     f"{len(no_file)} loci have no GWAS file at all\n")

    cols = (args.chr_col, args.pos_col, args.pval_col)
    jobs, job_keys = [], []
    for f, stem, st, keys in per_file:
        exact_keys = set((known[k]["chr"], known[k]["bp"]) for k in keys)
        windows = [(known[k]["chr"], known[k]["lo"], known[k]["hi"]) for k in keys]
        jobs.append((f, stem, exact_keys, windows, cols))
        job_keys.append((keys, st))

    # results[key]      = (best_p_exact, stratum_of_best, best_p_window, stratum_win)
    # per_stratum[key][stratum] = (p_exact, p_window) in that specific stratum
    results = {}
    per_stratum = defaultdict(dict)
    done = 0
    def handle(job_i, res):
        nonlocal done
        stem, exact, win_best, err = res
        done += 1
        if err:
            sys.stderr.write(f"WARN [{stem}] {err}\n")
            return
        keys, st = job_keys[job_i]
        for wi, k in enumerate(keys):
            d = known[k]
            pe = exact.get((d["chr"], d["bp"]))
            pw = win_best.get(wi)
            per_stratum[k][st] = (pe, pw)
            cur = results.get(k, (None, None, None, None))
            be, bse, bw, bsw = cur
            if pe is not None and (be is None or pe < be):
                be, bse = pe, st
            if pw is not None and (bw is None or pw < bw):
                bw, bsw = pw, st
            results[k] = (be, bse, bw, bsw)
        if done % 50 == 0 or done == len(jobs):
            sys.stderr.write(f"  scanned {done}/{len(jobs)} files\n")

    if args.threads > 1:
        with ProcessPoolExecutor(max_workers=args.threads) as ex:
            for i, res in enumerate(ex.map(scan_sumstat, jobs)):
                handle(i, res)
    else:
        for i, job in enumerate(jobs):
            handle(i, scan_sumstat(job))

    # per-locus table
    thr = sorted(set(args.thresholds), reverse=True)   # loosest -> strictest cols
    rec_cols = [f"recovered_{thr_label(t)}" for t in thr]
    fL = f"{args.out_prefix}_recovery_per_locus.tsv"
    with open(fL, "w") as out:
        out.write("\t".join([
            "phenotype", "chr", "bp", "window_lo", "window_hi",
            "present", "best_p_exact", "stratum_exact",
            "best_p_window", "stratum_window",
            *rec_cols, "status"]) + "\n")
        for key, d in known.items():
            if key in results:
                pe, se, pw, sw = results[key]
                # 'present' = exact position seen in >=1 stratum
                if pe is not None:
                    present, status = "1", "tested"
                else:
                    present, status = "0", "absent_from_data"
            else:
                pe = se = pw = sw = None
                present = "0"
                status = "no_gwas_file" if key in no_file else "absent_from_data"
            def rec(p, t):
                # blank when not present -> excluded from %, not a 0
                return "" if p is None else ("1" if p < t else "0")
            out.write("\t".join([
                d["pheno"], d["chr"], str(d["bp"]), str(d["lo"]), str(d["hi"]),
                present,
                "" if pe is None else f"{pe:.3g}", se or "",
                "" if pw is None else f"{pw:.3g}", sw or "",
                *[rec(pe, t) for t in thr],
                status]) + "\n")

    # summary: per phenotype and overall
    def summarise(keys, use_window):
        # denominator = loci testable under the chosen mode:
        #   exact  -> exact position present in >=1 stratum (p_exact not None)
        #   window -> any variant present in the locus window (p_window not None)
        # loci not testable under the mode are excluded from the rate.
        # returns (n_known, n_testable, [recovered_count per threshold in `thr`])
        n_known = len(keys)
        testable = 0
        rec_n = [0] * len(thr)
        for k in keys:
            pe, _se, pw, _sw = results.get(k, (None, None, None, None))
            p = pw if use_window else pe
            if p is None:
                continue                 # not testable in this mode -> excluded
            testable += 1
            for i, t in enumerate(thr):
                if p < t:
                    rec_n[i] += 1
        return n_known, testable, rec_n

    by_ph = defaultdict(list)
    for k in known:
        by_ph[k[0]].append(k)

    modes = {"exact": False, "window": True}
    if args.mode != "both":
        modes = {args.mode: modes[args.mode]}

    # header columns: recovered_<t> and pct_<t> per threshold
    tcols = []
    for t in thr:
        tcols += [f"recovered_{thr_label(t)}", f"pct_{thr_label(t)}"]

    headline = {}
    for mode_name, use_window in modes.items():
        fS = f"{args.out_prefix}_recovery_{mode_name}.tsv"
        with open(fS, "w") as out:
            out.write("\t".join(["phenotype", "n_known", "n_testable", *tcols])
                      + "\n")
            def line(label, keys):
                n_known, testable, rec_n = summarise(keys, use_window)
                pct = lambda x: f"{100*x/testable:.1f}" if testable else "NA"
                cells = []
                for c in rec_n:
                    cells += [str(c), pct(c)]
                out.write("\t".join([label, str(n_known), str(testable), *cells])
                          + "\n")
                return n_known, testable, rec_n
            for ph in sorted(by_ph):
                line(ph, by_ph[ph])
            headline[mode_name] = line("__ALL__", list(known.keys()))
            headline[mode_name] = (*headline[mode_name], fS)

    # ---- second breakdown: per (phenotype, stratum), ancestry-aware ----
    # A known locus is eligible for a stratum's recovery if:
    #   stratum in ANY_ANCESTRY_STRATA (ALL, MALE, FEMALE)  -> any ancestry
    #   stratum in EXACT_ANCESTRY_STRATA (EUR/EAS/...)       -> the locus maps
    #                                                          to that stratum
    # Recovery p is the p IN THAT STRATUM'S file (not best across strata).
    def eligible(k, st):
        if st in ANY_ANCESTRY_STRATA:
            return True
        return st in known[k]["anc_strata"]

    def summarise_stratum(keys, st, use_window):
        n_elig = testable = 0
        rec_n = [0] * len(thr)
        for k in keys:
            if not eligible(k, st):
                continue
            n_elig += 1
            pe, pw = per_stratum.get(k, {}).get(st, (None, None))
            p = pw if use_window else pe
            if p is None:
                continue            # eligible but not present in this stratum
            testable += 1
            for i, t in enumerate(thr):
                if p < t:
                    rec_n[i] += 1
        return n_elig, testable, rec_n

    # which strata actually have data (a file was scanned)
    scanned_strata = sorted({st for _, _, st, _ in per_file},
                            key=lambda s: KNOWN_STRATA.index(s)
                            if s in KNOWN_STRATA else 99)

    for mode_name, use_window in modes.items():
        fB = f"{args.out_prefix}_recovery_{mode_name}_by_stratum.tsv"
        with open(fB, "w") as out:
            out.write("\t".join([
                "phenotype", "stratum", "ancestry_rule",
                "n_eligible_known", "n_testable", *tcols]) + "\n")
            for ph in sorted(by_ph):
                keys = by_ph[ph]
                for st in scanned_strata:
                    ne, nt, rec_n = summarise_stratum(keys, st, use_window)
                    if ne == 0:
                        continue     # no eligible known loci for this stratum
                    rule = "any ancestry" if st in ANY_ANCESTRY_STRATA \
                           else f"exact ({st})"
                    pct = lambda x: f"{100*x/nt:.1f}" if nt else "NA"
                    cells = []
                    for c in rec_n:
                        cells += [str(c), pct(c)]
                    out.write("\t".join([ph, st, rule, str(ne), str(nt), *cells])
                              + "\n")
        sys.stderr.write(f"Per (phenotype, stratum) breakdown -> {fB}\n")

    sys.stderr.write("\n=== Recovery of known loci ===\n")
    for mode_name in modes:
        NK, NT, rec_n, fS = headline[mode_name]
        def pc(x): return f"{100*x/NT:.1f}%" if NT else "NA"
        label = "exact position" if mode_name == "exact" else \
                "window (Clump start-end)"
        sys.stderr.write(
            f"[{label}]  {NK} known; {NT} testable; {NK - NT} excluded\n")
        for t, c in zip(thr, rec_n):
            sys.stderr.write(
                f"    p<{t:g}\trecovered {pc(c)}  ({c}/{NT})\n")
        sys.stderr.write(f"    -> {fS}\n")
    sys.stderr.write(f"Per-locus detail: {fL}\n")


if __name__ == "__main__":
    main()
