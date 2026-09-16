#!/usr/bin/env python3
"""
Build the META-ANALYSIS RESULTS table, one row per <PHENOTYPE>_<STRATUM> GWAS.

Columns (matching the workbook sheet):
    Phenotype, GWAS(=stratum), lambda (meta), N meta variants,
    n variants p<=5e-8, loci_new, n_clumps >=2 GW sig variants,
    n_clumps >=5 GW sig variants, Cohorts used, Total cases, Total controls,
    Notes, Total sample size

COMPUTED FROM SUMSTATS (per file, one streaming pass):
    lambda (meta)      genomic-inflation lambda_GC = median(chisq) / 0.4549,
                       chisq derived from the p-value column
    N meta variants    number of variants with a usable p
    n variants p<=5e-8 genome-wide-significant variant count

FROM THE CLUMPING SUMMARY (verbose_summary.sh <pre>_nclumps.csv), joined by name:
    loci_new                    = n_clumps (MHC already collapsed to 1)
    n_clumps >=2 GW sig variants
    n_clumps >=5 GW sig variants

LEFT BLANK (external, added manually as before):
    Cohorts used, Total cases, Total controls, Notes, Total sample size

Usage:
    python3 build_meta_results.py \\
        --gwas-dir . \\
        --nclumps clumped_1mb/verbose_summary_nclumps.csv \\
        --out meta_results.tsv --threads 8
"""

import argparse
import csv
import gzip
import io
import math
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor

KNOWN = {"ALL", "EUR", "AFR", "AMR", "EAS", "SAS", "MALE", "FEMALE"}
SIG = 5e-8
# median of a 1-df chi-square (the expected median under the null)
CHI2_MEDIAN_1DF = 0.4549364231195724


def split_pheno_stratum(stem):
    m = re.match(r"^(.*)_([A-Za-z]+)$", stem)
    if m and m.group(2).upper() in KNOWN:
        return m.group(1), m.group(2).upper()
    return stem, ""


def open_text(path):
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8",
                                errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def sniff_delim(line):
    if "\t" in line:
        return "\t"
    if "," in line:
        return ","
    return None


def p_to_chi2(p):
    """p-value -> 1-df chi-square via the normal quantile (z^2)."""
    if p is None or p <= 0 or p >= 1:
        return None
    if p < 1e-300:
        p = 1e-300
    try:
        from scipy.stats import norm
        z = norm.isf(p / 2.0)
        return z * z
    except Exception:
        # Acklam inverse-normal fallback
        a=[-3.969683028665376e+01,2.209460984245205e+02,-2.759285104469687e+02,
           1.383577518672690e+02,-3.066479806614716e+01,2.506628277459239e+00]
        b=[-5.447609879822406e+01,1.615858368580409e+02,-1.556989798598866e+02,
           6.680131188771972e+01,-1.328068155288572e+01]
        c=[-7.784894002430293e-03,-3.223964580411365e-01,-2.400758277161838e+00,
           -2.549732539343734e+00,4.374664141464968e+00,2.938163982698783e+00]
        d=[7.784695709041462e-03,3.224671290700398e-01,2.445134137142996e+00,
           3.754408661907416e+00]
        pp = p / 2.0; pl = 0.02425
        if pp < pl:
            q=math.sqrt(-2*math.log(pp))
            z=(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
        elif pp<=1-pl:
            q=pp-0.5; r=q*q
            z=(((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q/(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)
        else:
            q=math.sqrt(-2*math.log(1-pp))
            z=-(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
        return z*z


def scan_file(job):
    """One streaming pass: return (stem, n_variants, n_sig, lambda, err).
    lambda from the median chisq of a sample of p-values (memory-bounded)."""
    path, pcol, stem = job
    try:
        import statistics
        chi2s = []
        n_var = 0
        n_sig = 0
        with open_text(path) as fh:
            first = fh.readline()
            delim = sniff_delim(first)
            header = (first.rstrip("\n").split(delim) if delim
                      else first.split())
            idx = {h.strip().lstrip("\ufeff"): i for i, h in enumerate(header)}
            if pcol not in idx:
                # try case-insensitive
                low = {h.lower(): i for h, i in idx.items()}
                if pcol.lower() in low:
                    ip = low[pcol.lower()]
                else:
                    return (stem, 0, 0, None, f"no {pcol!r} column; header={header[:8]}")
            else:
                ip = idx[pcol]
            for line in fh:
                f = line.rstrip("\n").split(delim) if delim else line.split()
                if len(f) <= ip:
                    continue
                try:
                    p = float(f[ip])
                except ValueError:
                    continue
                if not (0 < p <= 1):
                    continue
                n_var += 1
                if p <= SIG:
                    n_sig += 1
                c = p_to_chi2(p)
                if c is not None:
                    chi2s.append(c)
        lam = None
        if chi2s:
            med = statistics.median(chi2s)
            lam = med / CHI2_MEDIAN_1DF
        return (stem, n_var, n_sig, lam, None)
    except OSError as e:
        return (stem, 0, 0, None, f"read error: {e}")


def load_nclumps(path):
    """verbose_summary_nclumps.csv -> {stem: (n_clumps, ge2, ge5)}."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path, newline="") as fh:
        r = csv.DictReader(fh)
        for row in r:
            stem = row.get("file") or row.get("File")
            if not stem:
                continue
            def gi(*names):
                for n in names:
                    if n in row and row[n] not in (None, ""):
                        try: return int(float(row[n]))
                        except ValueError: return None
                return None
            out[stem.strip()] = (
                gi("n_clumps"),
                gi("n_clumps_ge2_gwsig", "n_clumps_gt1_sig"),
                gi("n_clumps_ge5_gwsig", "n_clumps_gt5_sig"))
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gwas-dir", default=".")
    ap.add_argument("--suffix", default=".tsv.gz")
    ap.add_argument("--pval-col", default="pvalue")
    ap.add_argument("--nclumps", default=None,
                    help="verbose_summary_*_nclumps.csv for clump counts")
    ap.add_argument("--out", default="meta_results.tsv")
    ap.add_argument("--threads", type=int, default=0)
    args = ap.parse_args()

    files = sorted(f for f in os.listdir(args.gwas_dir)
                   if f.endswith(args.suffix) and not f.startswith("._"))
    if not files:
        sys.exit(f"No {args.suffix} files in {args.gwas_dir}")
    jobs = [(os.path.join(args.gwas_dir, f), args.pval_col, f[:-len(args.suffix)])
            for f in files]

    nclumps = load_nclumps(args.nclumps)
    if args.nclumps and not nclumps:
        sys.stderr.write(f"WARN: no rows read from {args.nclumps}\n")

    nt = args.threads if args.threads and args.threads > 0 else (os.cpu_count() or 1)
    nt = min(nt, len(jobs))
    sys.stderr.write(f"Scanning {len(jobs)} GWAS on {nt} process(es)\n")

    results = {}
    done = 0
    def take(res):
        nonlocal done
        stem, nv, ns, lam, err = res
        done += 1
        if err:
            sys.stderr.write(f"WARN [{stem}] {err}\n")
        results[stem] = (nv, ns, lam)
        if done % 25 == 0 or done == len(jobs):
            sys.stderr.write(f"  scanned {done}/{len(jobs)}\n")

    if nt == 1:
        for j in jobs: take(scan_file(j))
    else:
        with ProcessPoolExecutor(max_workers=nt) as ex:
            for res in ex.map(scan_file, jobs):
                take(res)

    hdr = ["Phenotype", "GWAS", "lambda (meta)", "N meta variants",
           "n variants p<=5e-8", "loci_new", "n_clumps >=2 GW sig variants",
           "n_clumps >=5 GW sig variants", "Cohorts used", "Total cases",
           "Total controls", "Notes", "Total sample size"]
    with open(args.out, "w") as out:
        out.write("\t".join(hdr) + "\n")
        for f in files:
            stem = f[:-len(args.suffix)]
            pheno, stratum = split_pheno_stratum(stem)
            nv, ns, lam = results.get(stem, (0, 0, None))
            nc = nclumps.get(stem, (None, None, None))
            out.write("\t".join([
                pheno, stratum,
                "" if lam is None else f"{lam:.3f}",
                str(nv), str(ns),
                "" if nc[0] is None else str(nc[0]),
                "" if nc[1] is None else str(nc[1]),
                "" if nc[2] is None else str(nc[2]),
                "", "", "", "", ""]) + "\n")     # cohorts/cases/controls/notes/N blank

    sys.stderr.write(f"Wrote {args.out} ({len(files)} GWAS)\n"
                     "Blank columns (external): Cohorts used, Total cases, "
                     "Total controls, Notes, Total sample size\n")


if __name__ == "__main__":
    main()
