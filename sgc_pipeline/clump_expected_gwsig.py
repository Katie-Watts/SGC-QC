#!/usr/bin/env python3
"""
Expected vs observed genome-wide-significant variants per clump.

Rationale: under a single causal signal, a clump member's association Z-score is,
in expectation, the lead variant's Z scaled by its LD (correlation r) with the
lead:   E[Z_j] = r_ij * Z_lead.  A member is therefore *expected* to be
genome-wide significant (|Z| > z_sig, where z_sig corresponds to --sig) when

        |r_ij| > z_sig / |Z_lead|      i.e.    r2_ij > (z_sig / Z_lead)^2

PLINK's --clump-verbose output already gives, per clump, the index variant's p
(=> Z_lead) and each member's RSQ (r^2 with the index). So we can predict how
many members SHOULD be GW-sig and compare with how many actually are -- a QC
check on each clump.

Reading the output:
  n_expected_gwsig  index + members with RSQ > (z_sig/Z_index)^2
  n_observed_gwsig  index + members with p < --sig  (what clumping found)
  ratio             observed / expected
  flag              'fewer_than_expected' when observed << expected
                    'more_than_expected'  when observed >> expected
                    (tolerances via --low-ratio / --high-ratio)

Caveats printed once at the end: RSQ comes from the LD reference panel used for
clumping, not necessarily in-sample LD; a panel/ancestry mismatch will widen the
expected-vs-observed gap independently of data quality. The prediction is on |Z|
(RSQ is unsigned), which is what a significance-count check needs.

Usage:
  python3 clump_expected_gwsig.py \\
      --verbose-dir clumped_1mb/results \\
      --out clump_expected_vs_observed.tsv \\
      --sig 5e-8
"""

import argparse
import math
import os
import re
import sys

SUMMARY_HDR = re.compile(r"^\s*CHR\s+F\s+SNP\s+BP\s+P\s+TOTAL")
INDEX_LINE = re.compile(r"^\s*\(INDEX\)")
MEMBER_LINE = re.compile(r"^\s+[0-9XYMT]+:[0-9]+:")
CHR_TOKEN = re.compile(r"^[0-9XYMT]+$")
KNOWN_STRATA = {"ALL", "EUR", "AFR", "AMR", "EAS", "SAS", "MALE", "FEMALE"}
MHC_CHR, MHC_START, MHC_END = "6", 25_000_000, 34_000_000


def to_float(t):
    try:
        v = float(t)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def norm_chr(c):
    c = re.sub(r"^chr", "", str(c).strip(), flags=re.I).upper()
    return {"23": "X", "24": "Y", "25": "X", "26": "MT", "M": "MT"}.get(c, c)


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


def p_to_z(p):
    """Two-sided p -> |Z|. Uses the survival function inverse; clamps tiny p."""
    if p is None or p <= 0:
        return None
    if p < 1e-300:
        p = 1e-300
    try:
        from scipy.stats import norm
        return abs(norm.isf(p / 2.0))
    except Exception:
        # fallback: rational approximation of the normal quantile
        # (Acklam); good to ~1e-9 relative error
        a=[-3.969683028665376e+01,2.209460984245205e+02,-2.759285104469687e+02,
           1.383577518672690e+02,-3.066479806614716e+01,2.506628277459239e+00]
        b=[-5.447609879822406e+01,1.615858368580409e+02,-1.556989798598866e+02,
           6.680131188771972e+01,-1.328068155288572e+01]
        c=[-7.784894002430293e-03,-3.223964580411365e-01,-2.400758277161838e+00,
           -2.549732539343734e+00,4.374664141464968e+00,2.938163982698783e+00]
        d=[7.784695709041462e-03,3.224671290700398e-01,2.445134137142996e+00,
           3.754408661907416e+00]
        pp = p / 2.0
        pl = 0.02425
        if pp < pl:
            q = math.sqrt(-2*math.log(pp))
            z = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
                ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
        elif pp <= 1-pl:
            q = pp-0.5; r=q*q
            z = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
                (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)
        else:
            q = math.sqrt(-2*math.log(1-pp))
            z = -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
                 ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
        return abs(z)


def parse_verbose(path, sig):
    """Yield clump dicts: index, index_p, index_z, members=[(snp,p,rsq)]."""
    clumps, cur = [], None
    expect = False
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if SUMMARY_HDR.match(line):
                expect = True
                continue
            f = line.split()
            if expect and f and CHR_TOKEN.match(f[0]):
                if cur:
                    clumps.append(cur)
                expect = False
                if len(f) < 5:
                    cur = None
                    continue
                p = to_float(f[4])
                cur = {"index": f[2], "index_p": p,
                       "index_z": p_to_z(p), "members": []}
                continue
            if INDEX_LINE.match(line):
                continue
            if cur is not None and MEMBER_LINE.match(line):
                if not f:
                    continue
                snp = f[0]
                p = to_float(f[-1])          # last col = member p
                # RSQ is the 3rd column in the member block: KB RSQ ALLELES F P
                # member line fields: SNP KB RSQ ALLELES F P  (SNP is f[0])
                rsq = to_float(f[2]) if len(f) >= 3 else None
                cur["members"].append((snp, p, rsq))
    if cur:
        clumps.append(cur)
    return clumps


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--verbose-dir", required=True)
    ap.add_argument("--out", default="clump_expected_vs_observed.tsv")
    ap.add_argument("--sig", type=float, default=5e-8)
    ap.add_argument("--low-ratio", type=float, default=0.5,
                    help="flag 'fewer_than_expected' when observed/expected < this "
                         "(default 0.5)")
    ap.add_argument("--high-ratio", type=float, default=2.0,
                    help="flag 'more_than_expected' when observed/expected > this "
                         "(default 2.0)")
    ap.add_argument("--min-expected", type=int, default=2,
                    help="only evaluate clumps expected to have >= this many "
                         "GW-sig variants (default 2; avoids noise on singletons)")
    args = ap.parse_args()

    z_sig = p_to_z(args.sig)
    sys.stderr.write(f"Significance p<{args.sig:g}  =>  |Z| > {z_sig:.3f}\n")

    files = sorted(f for f in os.listdir(args.verbose_dir)
                   if f.endswith(".verbose.clumped")
                   and not f.startswith("._"))
    if not files:
        sys.exit(f"No .verbose.clumped files in {args.verbose_dir}")

    rows = []
    n_missing_rsq = 0
    n_clumps = 0
    for fn in files:
        pheno, stratum = split_pheno_stratum(fn[:-len(".verbose.clumped")])
        for c in parse_verbose(os.path.join(args.verbose_dir, fn), args.sig):
            zc = c["index_z"]
            if zc is None or zc <= 0:
                continue
            n_clumps += 1
            # r2 threshold above which a member is expected GW-sig
            r2_thr = (z_sig / zc) ** 2
            exp = 1                     # the index itself is GW-sig by construction
            obs = 1
            members_with_rsq = 0
            for snp, p, rsq in c["members"]:
                if rsq is not None:
                    members_with_rsq += 1
                    if rsq > r2_thr:
                        exp += 1
                if p is not None and p < args.sig:
                    obs += 1
            if members_with_rsq == 0 and c["members"]:
                n_missing_rsq += 1
            ratio = obs / exp if exp else float("nan")
            flag = ""
            # fewer-than-expected only meaningful when we expected a few
            if exp >= args.min_expected and ratio < args.low_ratio:
                flag = "fewer_than_expected"
            # more-than-expected: interesting even when expected is small
            # (e.g. many GW-sig members despite near-zero LD with the lead)
            elif ratio > args.high_ratio and obs >= args.min_expected:
                flag = "more_than_expected"
            rows.append((pheno, stratum, c["index"],
                         c["index_p"], zc, r2_thr,
                         exp, obs, ratio, len(c["members"]),
                         "MHC" if in_mhc(c["index"]) else "", flag))

    with open(args.out, "w") as out:
        out.write("\t".join([
            "phenotype", "stratum", "index_variant",
            "index_p", "index_z", "r2_threshold",
            "n_expected_gwsig", "n_observed_gwsig", "obs_over_exp",
            "n_members", "region", "flag"]) + "\n")
        for r in rows:
            out.write("\t".join([
                r[0], r[1], r[2],
                "" if r[3] is None else f"{r[3]:.3g}",
                f"{r[4]:.2f}", f"{r[5]:.4f}",
                str(r[6]), str(r[7]), f"{r[8]:.2f}",
                str(r[9]), r[10], r[11]]) + "\n")

    n_flag_low = sum(1 for r in rows if r[11] == "fewer_than_expected")
    n_flag_high = sum(1 for r in rows if r[11] == "more_than_expected")
    sys.stderr.write(
        f"\n{n_clumps} clumps evaluated across {len(files)} files.\n"
        f"  {n_flag_low} flagged fewer_than_expected (obs/exp < {args.low_ratio})\n"
        f"  {n_flag_high} flagged more_than_expected (obs/exp > {args.high_ratio})\n"
        f"Wrote {args.out}\n")
    if n_missing_rsq:
        sys.stderr.write(
            f"NOTE: {n_missing_rsq} clump(s) had members with no parseable RSQ; "
            f"expected counts there rest only on the index.\n")
    sys.stderr.write(
        "CAVEAT: RSQ is from the clumping LD reference panel, not in-sample LD; "
        "a panel/ancestry mismatch widens the expected-vs-observed gap "
        "independently of data quality.\n")


if __name__ == "__main__":
    main()
