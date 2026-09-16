#!/usr/bin/env python3
"""
Build the full Cross-disease variants sheet (all 101 columns).

Combines three sources into the wide per-locus layout:

  1. LOCI + phenotype/pvalue pairs  -- from cross_disease_verbose.py's wide TSV
     (Variant, n_phenotypes, n_gwsig_in_clump, phenotype_1, pvalue_1, ...).
     Gives: Locus_variant, n_associated_phenotypes, most_gwsig_in_clump,
            phenotype_k / pvalue_k.

  2. NOVELTY ROLLUP  -- each locus's phenotypes matched against the known-loci
     catalogue (Previous Known loci- variants format), per phenotype, with bp
     distance to the nearest known lead. Gives:
       n_phenotypes_total, n_known, n_novel_but_related, n_novel,
       Known phenotypes (dist), Known via related phenotype, Novel phenotypes
     A phenotype is:
       Known               its own phenotype has a catalogue entry near the locus
       Novel but related   a DIFFERENT but mapped-related phenotype does
                           (via --related-map, e.g. PSOR_VULGARIS->PSOR)
       Novel               neither

  3. DIRECTION DIFFERENCE (ALL only)  -- betas at the locus representative from
     each phenotype's _ALL sumstats, aligned to one reference allele. If the
     phenotypes disagree in sign (both a clearly + and a clearly - present),
     Direction_difference = "Yes" and n_positive/n_negative/pos_phenos/
     neg_phenos are filled; otherwise "No" and those four are "NA".
     (Needs --gwas-dir; skipped with "No"/NA if not supplied. Meta
     files have fixed column names -- NO cohort mapping needed.)

rsid and AF are left BLANK (pasted manually), inserted after Locus_variant by
the workbook assembler -- this script does NOT emit them.

Usage:
    python3 build_cross_disease.py \\
        --cross cross_disease.tsv \\
        --known novelty_info.txt \\
        --gwas-dir . \\
        --out cross_disease_full.tsv --threads 8
"""

import argparse
import csv
import glob
import gzip
import io
import math
import os
import re
import sys
from collections import defaultdict, Counter

ANC_MAP = {
    "european": "EUR", "east asian": "EAS", "south asian": "SAS",
    "african american or afro-caribbean": "AFR", "african unspecified": "AFR",
    "hispanic or latin american": "AMR",
}
NEARBY_BP = 500_000     # a known match beyond this distance is tagged "nearby"


def norm_chr(c):
    c = re.sub(r"^chr", "", str(c).strip(), flags=re.I).upper()
    return {"23": "X", "24": "Y", "25": "X", "26": "MT", "M": "MT"}.get(c, c)


def to_float(t):
    try:
        v = float(t); return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


# ---------- known-loci catalogue ----------
def load_known(path):
    """Return by_chr: chr -> list of {pheno, chr, bp, lead, anc}."""
    if path.lower().endswith((".xlsx", ".xlsm")):
        from openpyxl import load_workbook
        wb = load_workbook(path, read_only=True, data_only=True)
        ws = wb["Previous Known loci- variants"] if \
            "Previous Known loci- variants" in wb.sheetnames else wb.active
        data = list(ws.iter_rows(values_only=True))
        hdr = [str(h).strip() if h is not None else "" for h in data[0]]
        recs = [dict(zip(hdr, r)) for r in data[1:]]
    else:
        with open(path, encoding="utf-8-sig", errors="replace") as fh:
            first = fh.readline()
            delim = "\t" if "\t" in first else ("," if "," in first else "\t")
            fh.seek(0)
            recs = list(csv.DictReader(fh, delimiter=delim))
    by_chr = defaultdict(list)
    for r in recs:
        ph, chv, bpv = r.get("Phenotype"), r.get("Chr"), r.get("BP")
        if ph is None or chv is None or bpv is None:
            continue
        try:
            ch = norm_chr(chv); bp = int(float(bpv))
        except (ValueError, TypeError):
            continue
        by_chr[ch].append({
            "pheno": str(ph).strip(), "chr": ch, "bp": bp,
            "lead": r.get("Lead variant") or f"{ch}:{bp}",
            "anc": str(r.get("GWAS") or "").strip().lower()})
    return by_chr


def nearest_known(by_chr, pheno, ch, pos):
    """Nearest catalogue entry for `pheno` on `ch`; return (dist, entry) or None."""
    best = None
    for e in by_chr.get(ch, ()):
        if e["pheno"] != pheno:
            continue
        d = abs(e["bp"] - pos)
        if best is None or d < best[0]:
            best = (d, e)
    return best


# ---------- direction difference from _ALL META betas ----------
def open_text(path):
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8",
                                errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def scan_all_beta(job):
    """Read one <PHENO>_ALL META file for a set of positions; return
    {(chr,pos): (beta, effect_allele, other_allele)}.

    Meta-analysis files have KNOWN, consistent column names, so we read them
    directly -- no cohort column-mapping. Effect allele = alt, non-effect = ref
    (the pipeline's meta convention, same as extract_cross_disease_betas.py).
    Column names are overridable via the tuple passed in `cols`."""
    path, want, cols = job
    c_chr, c_pos, c_ref, c_alt, c_beta = cols
    out = {}
    try:
        with open_text(path) as fh:
            first = fh.readline()
            delim = "\t" if "\t" in first else ("," if "," in first else None)
            header = first.rstrip("\n").split(delim) if delim else first.split()
            idx = {h.strip().lstrip("\ufeff"): i for i, h in enumerate(header)}
            # case-insensitive fallback
            low = {h.lower(): i for h, i in idx.items()}
            def col(name):
                if name in idx: return idx[name]
                return low.get(name.lower())
            ic, ip, ib = col(c_chr), col(c_pos), col(c_beta)
            ialt, iref = col(c_alt), col(c_ref)   # effect=alt, other=ref
            if ic is None or ip is None or ib is None:
                return out
            need = max(x for x in [ic, ip, ib, ialt, iref] if x is not None)
            for line in fh:
                f = line.rstrip("\n").split(delim) if delim else line.split()
                if len(f) <= need:
                    continue
                try:
                    key = (norm_chr(f[ic]), int(f[ip]))
                except ValueError:
                    continue
                if key not in want:
                    continue
                out[key] = (to_float(f[ib]),
                            str(f[ialt]).strip().upper() if ialt is not None else "",
                            str(f[iref]).strip().upper() if iref is not None else "")
    except OSError:
        pass
    return out


def align_sign(beta, ea, oa, ref_ea, ref_oa):
    if beta is None or not ref_ea or not ea:
        return beta
    if ea == ref_ea:
        return beta
    if ea == ref_oa or oa == ref_ea:
        return -beta
    return beta


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cross", required=True,
                    help="cross_disease_verbose.py wide TSV")
    ap.add_argument("--known", required=True, help="known-loci catalogue")
    ap.add_argument("--related-map", default=None,
                    help="optional TSV: child_phenotype<TAB>parent_phenotype "
                         "for 'novel but related' (e.g. PSOR_VULGARIS PSOR)")
    ap.add_argument("--gwas-dir", default=None,
                    help="meta files <PHENO>_ALL.tsv.gz for the direction-"
                         "difference block; if omitted, that block is 'No'/NA")
    ap.add_argument("--suffix", default=".tsv.gz")
    # meta-analysis column names (fixed & consistent; override if yours differ)
    ap.add_argument("--chr-col", default="chromosome")
    ap.add_argument("--pos-col", default="position")
    ap.add_argument("--ref-col", default="ref")        # non-effect allele
    ap.add_argument("--alt-col", default="alt")        # effect allele
    ap.add_argument("--beta-col", default="beta")
    ap.add_argument("--min-beta", type=float, default=0.0,
                    help="ignore |beta| below this when judging sign (default 0)")
    ap.add_argument("--out", default="cross_disease_full.tsv")
    ap.add_argument("--threads", type=int, default=0)
    args = ap.parse_args()

    # ---- read cross-disease wide TSV ----
    with open(args.cross, encoding="utf-8-sig") as fh:
        first = fh.readline()
        delim = "\t" if "\t" in first else ","
        fh.seek(0)
        rows = list(csv.reader(fh, delimiter=delim))
    hdr = rows[0]
    H = {h: i for i, h in enumerate(hdr)}
    v_i = H.get("Variant", 0)
    ng_i = H.get("n_gwsig_in_clump")
    pcol_idx = [i for i, h in enumerate(hdr) if h.startswith("phenotype_")]
    max_np = len(pcol_idx)

    loci = []
    for r in rows[1:]:
        if not r or not r[v_i]:
            continue
        var = r[v_i].strip()
        phenos = []       # (pheno_label, pvalue_str)
        for pi in pcol_idx:
            if pi < len(r) and r[pi] and str(r[pi]).strip():
                pv = r[pi + 1] if pi + 1 < len(r) else ""
                phenos.append((str(r[pi]).strip(), pv))
        most_gw = r[ng_i] if ng_i is not None and ng_i < len(r) else ""
        loci.append({"var": var, "phenos": phenos, "most_gw": most_gw})
    sys.stderr.write(f"{len(loci)} cross-disease loci; up to {max_np} phenotypes\n")

    by_chr = load_known(args.known)

    related = {}
    if args.related_map and os.path.exists(args.related_map):
        with open(args.related_map, encoding="utf-8-sig") as fh:
            for line in fh:
                p = re.split(r"[\t,]", line.strip())
                if len(p) >= 2 and p[0] and p[1]:
                    related[p[0].strip()] = p[1].strip()

    # ---- direction difference: gather ALL-stratum betas at each locus rep ----
    beta_at = {}    # (pheno, chr, pos) -> (beta, ea, oa)
    if args.gwas_dir:
        meta_cols = (args.chr_col, args.pos_col, args.ref_col,
                     args.alt_col, args.beta_col)
        need = defaultdict(set)     # pheno -> {(chr,pos)}
        for l in loci:
            p = l["var"].split(":")
            try:
                key = (norm_chr(p[0]), int(p[1]))
            except (ValueError, IndexError):
                continue
            for lab, _ in l["phenos"]:
                ph = re.sub(r"\s*\(.*", "", lab).strip()
                need[ph].add(key)
        jobs = [(os.path.join(args.gwas_dir, f"{ph}_ALL{args.suffix}"),
                 keys, meta_cols, ph)
                for ph, keys in need.items()
                if os.path.exists(os.path.join(args.gwas_dir,
                                               f"{ph}_ALL{args.suffix}"))]
        sys.stderr.write(f"Direction: scanning {len(jobs)} _ALL meta files "
                         f"for betas\n")
        nt = args.threads if args.threads and args.threads > 0 else (os.cpu_count() or 1)
        nt = min(nt, len(jobs)) if jobs else 1

        results = {}
        if jobs:
            simple = [(path, keys, cols) for (path, keys, cols, _ph) in jobs]
            phs = [j[3] for j in jobs]
            if nt == 1:
                for ph, j in zip(phs, simple):
                    results[ph] = scan_all_beta(j)
            else:
                from concurrent.futures import ProcessPoolExecutor
                with ProcessPoolExecutor(max_workers=nt) as ex:
                    for ph, hits in zip(phs, ex.map(scan_all_beta, simple)):
                        results[ph] = hits
        for ph, hits in results.items():
            for key, (beta, ea, oa) in hits.items():
                beta_at[(ph, key[0], key[1])] = (beta, ea, oa)

    # ---- assemble ----
    base_hdr = ["Locus_variant", "n_associated_phenotypes", "most_gwsig_in_clump"]
    pair_hdr = []
    for k in range(1, max_np + 1):
        pair_hdr += [f"phenotype_{k}", f"pvalue_{k}"]
    tail_hdr = ["n_phenotypes_total", "n_known", "n_novel_but_related", "n_novel",
                "Known phenotypes (dist)", "Known via related phenotype",
                "Novel phenotypes",
                "Direction_difference_across_phenos_(ALL only)",
                "n_positive", "n_negative", "pos_phenos", "neg_phenos"]
    full_hdr = base_hdr + pair_hdr + tail_hdr

    def dist_tag(d):
        return f"{d}bp" + (",nearby" if d > NEARBY_BP else "")

    out_rows = []
    for l in loci:
        var = l["var"]
        p = var.split(":")
        try:
            ch, pos = norm_chr(p[0]), int(p[1])
        except (ValueError, IndexError):
            ch, pos = None, None
        phenos = l["phenos"]
        n_assoc = len(phenos)

        # novelty rollup
        known_list, related_list, novel_list = [], [], []
        for lab, _pv in phenos:
            ph = re.sub(r"\s*\(.*", "", lab).strip()
            hit = nearest_known(by_chr, ph, ch, pos) if ch else None
            if hit:
                known_list.append(f"{ph}({dist_tag(hit[0])})")
                continue
            parent = related.get(ph)
            phit = nearest_known(by_chr, parent, ch, pos) if (parent and ch) else None
            if phit:
                related_list.append(f"{ph}\u2192{parent}({dist_tag(phit[0])})")
            else:
                novel_list.append(ph)

        # direction difference (ALL only)
        dir_diff, n_pos, n_neg, pos_ph, neg_ph = "No", "NA", "NA", "NA", "NA"
        if beta_at and ch is not None:
            signs = []      # (pheno, beta)
            ref_ea = ref_oa = ""
            for lab, _pv in phenos:
                ph = re.sub(r"\s*\(.*", "", lab).strip()
                rec = beta_at.get((ph, ch, pos))
                if not rec:
                    continue
                beta, ea, oa = rec
                if beta is None:
                    continue
                if not ref_ea and ea:
                    ref_ea, ref_oa = ea, oa
                b = align_sign(beta, ea, oa, ref_ea, ref_oa)
                if abs(b) >= args.min_beta:
                    signs.append((ph, b))
            pos_p = [ph for ph, b in signs if b > 0]
            neg_p = [ph for ph, b in signs if b < 0]
            if pos_p and neg_p:
                dir_diff = "Yes"
                n_pos, n_neg = str(len(pos_p)), str(len(neg_p))
                pos_ph = "; ".join(pos_p); neg_ph = "; ".join(neg_p)

        row = [var, str(n_assoc), str(l["most_gw"])]
        for k in range(max_np):
            if k < len(phenos):
                row += [phenos[k][0], phenos[k][1]]
            else:
                row += ["", ""]
        row += [str(n_assoc), str(len(known_list)), str(len(related_list)),
                str(len(novel_list)),
                "; ".join(known_list), "; ".join(related_list),
                "; ".join(novel_list),
                dir_diff, n_pos, n_neg, pos_ph, neg_ph]
        out_rows.append(row)

    with open(args.out, "w") as out:
        out.write("\t".join(full_hdr) + "\n")
        for r in out_rows:
            out.write("\t".join(r) + "\n")

    n_yes = sum(1 for r in out_rows if r[base_len(max_np)] == "Yes") \
        if out_rows else 0
    sys.stderr.write(
        f"Wrote {args.out}: {len(out_rows)} loci, {len(full_hdr)} columns\n"
        f"  direction-difference computed: "
        + ("yes" if beta_at else "no (no --gwas-dir) -> all 'No'/NA")
        + "\n  (rsid / AF left for manual paste; assembler inserts them)\n")


def base_len(max_np):
    # index of Direction_difference column in a row (0-based)
    return 3 + 2 * max_np + 7


if __name__ == "__main__":
    main()
