#!/usr/bin/env python3
"""
Add per-variant meta-analysis metrics (n_cohorts, dir_concordance, i2) to three
sheets, looking each variant up in the matching <PHENO>_<STRATUM>.tsv.gz file.

Sheets and how each keys to a meta file:
  Clump_index_variants   GWAS column is already PHENO_STRATUM; variant in
                         index_variant. One file per row.
  Ancestry specific loci Phenotype + Ancestry -> <Phenotype>_<Ancestry>.tsv.gz;
                         variant in 'variant'. One file per row.
  Sex specific loci      Phenotype + Variant, spans BOTH sexes -> look up in
                         <Phenotype>_MALE and <Phenotype>_FEMALE, adding
                         male_/female_ prefixed columns.

Position-only matching (chr:pos, alleles ignored), consistent with the pipeline.
Each meta file is scanned once for all the variants that need it.

Usage:
    python3 add_meta_metrics.py \\
        --xlsx SGC_Meta_Analysis_Results.xlsx \\
        --gwas-dir . \\
        --out SGC_Meta_Analysis_Results.xlsx \\
        --suffix .tsv.gz --threads 8
"""

import argparse
import gzip
import io
import os
import re
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from copy import copy

CHR_ALIASES = {"CHR", "CHROMOSOME", "CHROM", "#CHROM"}
POS_ALIASES = {"POS", "POSITION", "GENPOS", "BP", "BASE_PAIR_LOCATION"}
# metric name -> accepted header spellings (case-insensitive)
METRIC_ALIASES = {
    "n_cohorts":       {"N_COHORTS", "N_COHORT", "NCOHORT", "NCOHORTS",
                        "N_STUDIES", "NSTUDY"},
    "dir_concordance": {"DIR_CONCORDANCE", "DIRECTION_CONCORDANCE",
                        "DIRECTIONCONCORDANCE", "CONCORDANCE",
                        "DIR_CONC", "DIRECTION"},
    "i2":              {"I2", "I_2", "ISQ", "I2_HET", "HETISQ", "HETEROGENEITY_I2"},
}
METRICS = ["n_cohorts", "dir_concordance", "i2"]


def norm_chr(c):
    c = re.sub(r"^chr", "", str(c).strip(), flags=re.I).upper()
    return {"23": "X", "24": "Y", "25": "X", "26": "MT", "M": "MT"}.get(c, c)


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


def variant_pos(variant):
    p = str(variant).split(":")
    if len(p) < 2:
        return None
    try:
        return (norm_chr(p[0]), int(p[1]))
    except ValueError:
        return None


def resolve_metric_cols(idx, overrides):
    """idx = {HEADER_UPPER: col}. Return {metric: col_index or None}."""
    out = {}
    for m in METRICS:
        if overrides.get(m):
            out[m] = idx.get(overrides[m].upper())
        else:
            out[m] = next((idx[a] for a in METRIC_ALIASES[m] if a in idx), None)
    return out


def scan_file(job):
    """Read one meta file once; return {(chr,pos): {metric: value}} for wanted."""
    path, want, overrides = job
    fn = os.path.basename(path)
    out = {}
    if not os.path.exists(path):
        return (fn, None, "file not found", {})
    try:
        with open_text(path) as fh:
            first = fh.readline()
            delim = sniff_delim(first)
            header = (first.rstrip("\n").split(delim) if delim
                      else first.split())
            idx = {h.strip().upper(): i for i, h in enumerate(header)}
            ic = next((idx[a] for a in CHR_ALIASES if a in idx), None)
            ip = next((idx[a] for a in POS_ALIASES if a in idx), None)
            mcols = resolve_metric_cols(idx, overrides)
            if ic is None or ip is None:
                return (fn, None, f"no chr/pos column; header={header[:6]}", {})
            found_metrics = {m: c for m, c in mcols.items() if c is not None}
            if not found_metrics:
                return (fn, None,
                        f"none of {METRICS} present; header={header[:10]}", {})
            need = max([ic, ip] + list(found_metrics.values()))
            for line in fh:
                f = line.rstrip("\n").split(delim) if delim else line.split()
                if len(f) <= need:
                    continue
                try:
                    key = (norm_chr(f[ic]), int(f[ip]))
                except ValueError:
                    continue
                if key in want and key not in out:
                    out[key] = {m: f[c].strip() for m, c in found_metrics.items()}
    except OSError as e:
        return (fn, None, f"read error: {e}", {})
    return (fn, out, None, {m: (m in found_metrics) for m in METRICS})


def numify(v):
    if v is None or v == "":
        return None
    try:
        return int(v)
    except ValueError:
        try:
            return float(v)
        except ValueError:
            return v


def style_header(ws, col, name):
    hs = ws.cell(1, 1)
    c = ws.cell(1, col, name)
    c.font = copy(hs.font); c.fill = copy(hs.fill)
    c.alignment = copy(hs.alignment); c.border = copy(hs.border)


def gwas_for_row_clump(ws, r, cols):
    """(gwas_label, variant) for Clump_index_variants."""
    g = ws.cell(r, cols["GWAS"]).value
    v = ws.cell(r, cols["index_variant"]).value
    return (str(g).strip() if g else None, v)


def process_direct_sheet(ws, wb, gwas_dir, suffix, overrides, threads,
                         gwas_col, var_col, sheet_label):
    """Sheets with one (gwas, variant) per row -> add 3 metric columns after
    the variant column. gwas_col/var_col are header names.
    Returns (per_gwas jobs done, warnings)."""
    hdr = [ws.cell(1, c).value for c in range(1, ws.max_column + 1)]
    H = {h: i + 1 for i, h in enumerate(hdr) if h is not None}
    for need in (gwas_col, var_col):
        if need not in H:
            sys.stderr.write(f"  [{sheet_label}] missing {need!r}; skipped\n")
            return
    c_g, c_v = H[gwas_col], H[var_col]

    per_gwas = defaultdict(set)
    row_key = {}
    for r in range(2, ws.max_row + 1):
        g = ws.cell(r, c_g).value
        v = ws.cell(r, c_v).value
        if g is None or v is None:
            continue
        pk = variant_pos(v)
        if pk is None:
            continue
        g = str(g).strip()
        per_gwas[g].add(pk)
        row_key[r] = (g, pk)

    results = scan_all(per_gwas, gwas_dir, suffix, overrides, threads,
                       sheet_label)

    # insert 3 columns after the variant column
    start = c_v + 1
    ws.insert_cols(start, amount=3)
    for i, m in enumerate(METRICS):
        style_header(ws, start + i, m)
    bstyle = ws.cell(2, 1)
    filled = 0
    for r in range(2, ws.max_row + 1):
        rk = row_key.get(r)
        for i, m in enumerate(METRICS):
            cell = ws.cell(r, start + i)
            cell.font = copy(bstyle.font)
        if rk is None:
            continue
        g, pk = rk
        vals = results.get(g, {}).get(pk)
        if vals:
            for i, m in enumerate(METRICS):
                nv = numify(vals.get(m))
                if nv is not None:
                    ws.cell(r, start + i).value = nv
            filled += 1
    sys.stderr.write(f"  [{sheet_label}] filled {filled} rows\n")


def process_ancestry_sheet(ws, gwas_dir, suffix, overrides, threads):
    """Phenotype + Ancestry -> <Phenotype>_<Ancestry>. variant in 'variant'."""
    hdr = [ws.cell(1, c).value for c in range(1, ws.max_column + 1)]
    H = {h: i + 1 for i, h in enumerate(hdr) if h is not None}
    for need in ("Phenotype", "Ancestry", "variant"):
        if need not in H:
            sys.stderr.write(f"  [Ancestry] missing {need!r}; skipped\n")
            return
    c_ph, c_anc, c_v = H["Phenotype"], H["Ancestry"], H["variant"]

    per_gwas = defaultdict(set)
    row_key = {}
    for r in range(2, ws.max_row + 1):
        ph = ws.cell(r, c_ph).value
        anc = ws.cell(r, c_anc).value
        v = ws.cell(r, c_v).value
        if ph is None or anc is None or v is None:
            continue
        pk = variant_pos(v)
        if pk is None:
            continue
        g = f"{str(ph).strip()}_{str(anc).strip()}"
        per_gwas[g].add(pk)
        row_key[r] = (g, pk)

    results = scan_all(per_gwas, gwas_dir, suffix, overrides, threads,
                       "Ancestry")

    start = c_v + 1
    ws.insert_cols(start, amount=3)
    for i, m in enumerate(METRICS):
        style_header(ws, start + i, m)
    bstyle = ws.cell(2, 1)
    filled = 0
    for r in range(2, ws.max_row + 1):
        rk = row_key.get(r)
        for i in range(3):
            ws.cell(r, start + i).font = copy(bstyle.font)
        if rk is None:
            continue
        g, pk = rk
        vals = results.get(g, {}).get(pk)
        if vals:
            for i, m in enumerate(METRICS):
                nv = numify(vals.get(m))
                if nv is not None:
                    ws.cell(r, start + i).value = nv
            filled += 1
    sys.stderr.write(f"  [Ancestry] filled {filled} rows\n")


def process_sex_sheet(ws, gwas_dir, suffix, overrides, threads):
    """Phenotype + Variant spanning both sexes -> look up in _MALE and _FEMALE,
    add male_ and female_ prefixed metric columns."""
    hdr = [ws.cell(1, c).value for c in range(1, ws.max_column + 1)]
    H = {h: i + 1 for i, h in enumerate(hdr) if h is not None}
    for need in ("Phenotype", "Variant"):
        if need not in H:
            sys.stderr.write(f"  [Sex] missing {need!r}; skipped\n")
            return
    c_ph, c_v = H["Phenotype"], H["Variant"]

    per_gwas = defaultdict(set)     # <Phenotype>_MALE / _FEMALE -> positions
    row_key = {}                    # row -> (phenotype, pk)
    for r in range(2, ws.max_row + 1):
        ph = ws.cell(r, c_ph).value
        v = ws.cell(r, c_v).value
        if ph is None or v is None:
            continue
        pk = variant_pos(v)
        if pk is None:
            continue
        ph = str(ph).strip()
        per_gwas[f"{ph}_MALE"].add(pk)
        per_gwas[f"{ph}_FEMALE"].add(pk)
        row_key[r] = (ph, pk)

    results = scan_all(per_gwas, gwas_dir, suffix, overrides, threads, "Sex")

    # insert 6 columns after Variant: male_<m> then female_<m>
    start = c_v + 1
    new_cols = [f"male_{m}" for m in METRICS] + [f"female_{m}" for m in METRICS]
    ws.insert_cols(start, amount=len(new_cols))
    for i, name in enumerate(new_cols):
        style_header(ws, start + i, name)
    bstyle = ws.cell(2, 1)
    filled = 0
    for r in range(2, ws.max_row + 1):
        for i in range(len(new_cols)):
            ws.cell(r, start + i).font = copy(bstyle.font)
        rk = row_key.get(r)
        if rk is None:
            continue
        ph, pk = rk
        got = False
        mvals = results.get(f"{ph}_MALE", {}).get(pk)
        fvals = results.get(f"{ph}_FEMALE", {}).get(pk)
        for i, m in enumerate(METRICS):
            if mvals:
                nv = numify(mvals.get(m))
                if nv is not None:
                    ws.cell(r, start + i).value = nv; got = True
            if fvals:
                nv = numify(fvals.get(m))
                if nv is not None:
                    ws.cell(r, start + 3 + i).value = nv; got = True
        if got:
            filled += 1
    sys.stderr.write(f"  [Sex] filled {filled} rows\n")


def scan_all(per_gwas, gwas_dir, suffix, overrides, threads, label):
    """Scan each needed meta file once (parallel). Returns
    {gwas: {(chr,pos): {metric: value}}}."""
    jobs, order = [], []
    for g, keys in per_gwas.items():
        jobs.append((os.path.join(gwas_dir, f"{g}{suffix}"), keys, overrides))
        order.append(g)
    nt = threads if threads and threads > 0 else (os.cpu_count() or 1)
    nt = min(nt, len(jobs)) if jobs else 1
    sys.stderr.write(f"  [{label}] scanning {len(jobs)} meta files on {nt} "
                     f"process(es)\n")
    results = {}
    missing = []
    metric_seen = {m: False for m in METRICS}
    done = 0
    total = len(jobs)

    def handle(g, res):
        nonlocal done
        fn, out, err, seen = res
        done += 1
        if err:
            sys.stderr.write(f"    WARN [{fn}] {err}\n")
            if "not found" in err:
                missing.append(g)
        else:
            results[g] = out
            for m, s in seen.items():
                metric_seen[m] = metric_seen[m] or s
        if done % 5 == 0 or done == total:
            sys.stderr.write(f"    [{label}] scanned {done}/{total} files\n")
            sys.stderr.flush()

    if nt == 1:
        for g, job in zip(order, jobs):
            handle(g, scan_file(job))
    else:
        with ProcessPoolExecutor(max_workers=nt) as ex:
            for g, res in zip(order, ex.map(scan_file, jobs)):
                handle(g, res)
    if missing:
        u = sorted(set(missing))
        sys.stderr.write(f"    {len(u)} GWAS had no file: " + ", ".join(u[:8])
                         + (" ..." if len(u) > 8 else "") + "\n")
    absent = [m for m in METRICS if not metric_seen[m]]
    if absent:
        sys.stderr.write(f"    NOTE: metric(s) {absent} not found in any file "
                         f"header (check --col overrides)\n")
    return results


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xlsx", required=True)
    ap.add_argument("--gwas-dir", default=".")
    ap.add_argument("--suffix", default=".tsv.gz")
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--n-cohorts-col", default=None,
                    help="override header name for n_cohorts")
    ap.add_argument("--dir-concordance-col", default=None,
                    help="override header name for dir_concordance")
    ap.add_argument("--i2-col", default=None,
                    help="override header name for i2")
    ap.add_argument("--sheets", nargs="*",
                    default=["Clump_index_variants", "Ancestry specific loci",
                             "Sex specific loci"],
                    help="which sheets to annotate")
    args = ap.parse_args()

    overrides = {"n_cohorts": args.n_cohorts_col,
                 "dir_concordance": args.dir_concordance_col,
                 "i2": args.i2_col}

    from openpyxl import load_workbook
    wb = load_workbook(args.xlsx)

    for sheet in args.sheets:
        if sheet not in wb.sheetnames:
            sys.stderr.write(f"  sheet {sheet!r} not found; skipped\n")
            continue
        ws = wb[sheet]
        sys.stderr.write(f"Processing {sheet!r}\n")
        if sheet == "Clump_index_variants":
            process_direct_sheet(ws, wb, args.gwas_dir, args.suffix, overrides,
                                 args.threads, "GWAS", "index_variant", sheet)
        elif sheet == "Ancestry specific loci":
            process_ancestry_sheet(ws, args.gwas_dir, args.suffix, overrides,
                                   args.threads)
        elif sheet == "Sex specific loci":
            process_sex_sheet(ws, args.gwas_dir, args.suffix, overrides,
                              args.threads)
        else:
            sys.stderr.write(f"  no rule for {sheet!r}; skipped\n")

    wb.save(args.out)
    sys.stderr.write(f"\nWrote {args.out}\n")


if __name__ == "__main__":
    main()
