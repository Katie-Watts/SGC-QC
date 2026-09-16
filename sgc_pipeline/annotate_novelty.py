#!/usr/bin/env python3
"""
Annotate variants with NOVELTY against the known-loci catalogue.

The catalogue has the "Previous Known loci- variants" columns:
    Phenotype, GWAS(=ancestry label), Study, Lead variant, Chr, BP,
    P-value, Clump start, Clump end
(supplied either as novelty_info.txt/tsv or read from that sheet of the xlsx.)

For each query variant (Phenotype, chr, pos[, stratum]) it decides:
    Novelty            "Known" | "Novel but related" | "Novel"
    Matched known variant   the catalogue Lead variant it matched (if any)
    Match details      ancestry|study|dist_to_lead=..bp  (or window note)

Matching rule (position-based, ancestry-aware -- same logic family as
known_loci_recovery.py):
  * a query matches a catalogue entry of the SAME phenotype when its position
    falls in the entry's [Clump start, Clump end] window (or exact BP if no
    window). Exact-position match is reported as dist_to_lead=0.
  * ancestry gating (only when --stratum-col given and the query has a stratum):
      ALL / MALE / FEMALE      -> match catalogue entries of ANY ancestry
      EUR/EAS/SAS/AFR/AMR       -> match only entries whose GWAS ancestry maps
                                   to that stratum, PLUS any-ancestry catalogue
                                   entries (Any Ancestry / Multiancestry)
  * "Known"             matched an entry of the same phenotype
    "Novel but related" no same-phenotype match, but the position matches a
                        catalogue entry of a DIFFERENT phenotype (--related)
    "Novel"             no catalogue match anywhere

Input variant file: TSV with a phenotype column, and either a variant id column
(chr:pos[:a:a]) or separate chr/pos columns. Output = input + 3 novelty columns.

Usage:
    python3 annotate_novelty.py \\
        --variants clump_index.tsv --known novelty_info.txt \\
        --pheno-col GWAS --variant-col index_variant --stratum-from-gwas \\
        --out clump_index_novelty.tsv
"""

import argparse
import re
import sys
from collections import defaultdict

ANC_MAP = {
    "european": "EUR", "east asian": "EAS", "south asian": "SAS",
    "african american or afro-caribbean": "AFR", "african unspecified": "AFR",
    "hispanic or latin american": "AMR",
}
ANY_ANCESTRY_LABELS = {"any ancestry", "multiancestry"}
ANY_ANCESTRY_STRATA = {"ALL", "MALE", "FEMALE"}
KNOWN_STRATA = {"ALL", "EUR", "AFR", "AMR", "EAS", "SAS", "MALE", "FEMALE"}


def norm_chr(c):
    c = re.sub(r"^chr", "", str(c).strip(), flags=re.I).upper()
    return {"23": "X", "24": "Y", "25": "X", "26": "MT", "M": "MT"}.get(c, c)


def strip_stratum(gwas):
    m = re.match(r"^(.*)_([A-Za-z]+)$", str(gwas))
    if m and m.group(2).upper() in KNOWN_STRATA:
        return m.group(1), m.group(2).upper()
    return str(gwas), ""


def load_known(path):
    """Read the catalogue (tsv/csv/xlsx). Return list of dicts and an index
    by (phenotype, chr) -> list of entries."""
    rows = []
    if path.lower().endswith((".xlsx", ".xlsm")):
        from openpyxl import load_workbook
        wb = load_workbook(path, read_only=True, data_only=True)
        ws = wb["Previous Known loci- variants"] if \
            "Previous Known loci- variants" in wb.sheetnames else wb.active
        data = list(ws.iter_rows(values_only=True))
        hdr = [str(h).strip() if h is not None else "" for h in data[0]]
        recs = [dict(zip(hdr, r)) for r in data[1:]]
    else:
        import csv
        with open(path, encoding="utf-8-sig", errors="replace") as fh:
            sample = fh.readline()
            delim = "\t" if "\t" in sample else ("," if "," in sample else "\t")
            fh.seek(0)
            recs = list(csv.DictReader(fh, delimiter=delim))
    for r in recs:
        ph = r.get("Phenotype")
        chro = r.get("Chr"); bp = r.get("BP")
        if ph is None or chro is None or bp is None:
            continue
        try:
            ch = norm_chr(chro); pos = int(float(bp))
        except (ValueError, TypeError):
            continue
        def gi(k):
            v = r.get(k)
            try: return int(float(v))
            except (ValueError, TypeError): return pos
        lo, hi = gi("Clump start"), gi("Clump end")
        anc = str(r.get("GWAS") or "").strip().lower()
        rows.append({"pheno": str(ph).strip(), "chr": ch, "bp": pos,
                     "lo": min(lo, hi), "hi": max(lo, hi),
                     "anc": anc, "lead": r.get("Lead variant") or f"{ch}:{pos}",
                     "study": str(r.get("Study") or "")})
    by_chr = defaultdict(list)
    for e in rows:
        by_chr[e["chr"]].append(e)
    return rows, by_chr


def anc_eligible(entry_anc, stratum):
    """Does a catalogue entry's ancestry apply to this query stratum?"""
    if not stratum:
        return True                       # no stratum given -> any ancestry
    if stratum in ANY_ANCESTRY_STRATA:
        return True
    # ancestry-specific stratum: exact-ancestry OR any-ancestry catalogue entry
    if entry_anc in ANY_ANCESTRY_LABELS:
        return True
    return ANC_MAP.get(entry_anc) == stratum


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variants", required=True)
    ap.add_argument("--known", required=True,
                    help="novelty_info.txt / .tsv / the xlsx with the catalogue")
    ap.add_argument("--out", required=True)
    ap.add_argument("--pheno-col", default="Phenotype")
    ap.add_argument("--variant-col", default=None,
                    help="column holding chr:pos[:a:a]; or use --chr-col/--pos-col")
    ap.add_argument("--chr-col", default=None)
    ap.add_argument("--pos-col", default=None)
    ap.add_argument("--stratum-col", default=None,
                    help="column holding the stratum, for ancestry-aware matching")
    ap.add_argument("--stratum-from-gwas", action="store_true",
                    help="derive phenotype+stratum from the pheno-col (which is "
                         "actually GWAS = PHENO_STRATUM)")
    ap.add_argument("--related", action="store_true",
                    help="tag 'Novel but related' when a different phenotype "
                         "matches the position")
    args = ap.parse_args()

    known, by_chr = load_known(args.known)
    sys.stderr.write(f"Catalogue: {len(known)} known-loci entries\n")

    import csv
    with open(args.variants, encoding="utf-8-sig", errors="replace") as fh:
        first = fh.readline()
        delim = "\t" if "\t" in first else ("," if "," in first else "\t")
        fh.seek(0)
        rdr = csv.DictReader(fh, delimiter=delim)
        cols = rdr.fieldnames
        rows = list(rdr)

    def qpos(row):
        if args.variant_col:
            p = str(row.get(args.variant_col) or "").split(":")
            if len(p) < 2: return None
            try: return norm_chr(p[0]), int(p[1])
            except ValueError: return None
        c = row.get(args.chr_col); p = row.get(args.pos_col)
        try: return norm_chr(c), int(float(p))
        except (ValueError, TypeError): return None

    def qpheno_stratum(row):
        raw = str(row.get(args.pheno_col) or "").strip()
        if args.stratum_from_gwas:
            return strip_stratum(raw)
        st = str(row.get(args.stratum_col) or "").strip().upper() \
            if args.stratum_col else ""
        return raw, st

    out_cols = list(cols) + ["Novelty", "Matched known variant", "Match details"]
    n_known = n_related = n_novel = 0
    with open(args.out, "w") as out:
        out.write("\t".join(out_cols) + "\n")
        for row in rows:
            pos = qpos(row)
            pheno, stratum = qpheno_stratum(row)
            novelty, matched, details = "Novel", "", ""
            if pos is not None:
                ch, bp = pos
                same, other = None, None
                for e in by_chr.get(ch, ()):
                    if not (e["lo"] <= bp <= e["hi"]):
                        continue
                    if not anc_eligible(e["anc"], stratum):
                        continue
                    if e["pheno"] == pheno:
                        # prefer exact / closest lead
                        d = abs(e["bp"] - bp)
                        if same is None or d < same[0]:
                            same = (d, e)
                    elif other is None:
                        other = e
                if same is not None:
                    d, e = same
                    novelty = "Known"
                    matched = e["lead"]
                    details = f"{e['anc']}|{e['study']}|dist_to_lead={d}bp"
                    n_known += 1
                elif args.related and other is not None:
                    novelty = "Novel but related"
                    matched = other["lead"]
                    details = (f"related pheno {other['pheno']}|{other['anc']}"
                               f"|{other['study']}")
                    n_related += 1
                else:
                    n_novel += 1
            else:
                n_novel += 1
            out.write("\t".join(str(row.get(c, "")) for c in cols)
                      + f"\t{novelty}\t{matched}\t{details}\n")

    sys.stderr.write(
        f"Wrote {args.out}: {n_known} Known, "
        + (f"{n_related} Novel-but-related, " if args.related else "")
        + f"{n_novel} Novel\n")


if __name__ == "__main__":
    main()
