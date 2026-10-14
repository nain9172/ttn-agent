#!/usr/bin/env python3

import re
import sys
import json
import time
import pathlib
import requests
import pandas as pd

# ---------------------------------------------------------------- settings --
ROOT = pathlib.Path(__file__).resolve().parents[2]   # repo 根目錄
sys.path.insert(0, str(ROOT))
import config  # noqa: E402  讀 .env 的 PUBMED_EMAIL / PUBMED_API_KEY
SRC = f"{ROOT}/data/pubtator/pubtator_ttn_variants_all_with_clinvar.csv"
OUT = f"{ROOT}/data/pubtator/pubtator_ttn_variants_clinvar_annotated.csv"  # pipeline 直接讀這份
REVIEW = f"{ROOT}/data/pubtator/clinvar_ambiguous_for_review.csv"
CACHE = f"{ROOT}/data/pubtator/clinvar_cache.json"

EMAIL = config.PUBMED_EMAIL        # NCBI 要求附聯絡 email（設定在 .env）
API_KEY = config.PUBMED_API_KEY    # 沒有就留空（速率較慢）
ASSEMBLY = "GRCh38"
GENE = "TTN"      # results are constrained to this gene to cut down false hits

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DELAY = 0.11 if API_KEY else 0.36     # respect 10/s with key, 3/s without

session = requests.Session()


def _params(extra):
    p = {"db": "clinvar", "retmode": "json"}
    if EMAIL:
        p["email"] = EMAIL
    if API_KEY:
        p["api_key"] = API_KEY
    p.update(extra)
    return p


def _get(endpoint, extra, retries=4):
    for attempt in range(retries):
        try:
            r = session.get(f"{BASE}/{endpoint}.fcgi",
                            params=_params(extra), timeout=60)
            if r.status_code == 429:
                time.sleep(2 * (attempt + 1))
                continue
            r.raise_for_status()
            return r.json()
        except Exception as e:
            if attempt == retries - 1:
                print(f"    ! {endpoint} failed: {e}")
                return None
            time.sleep(2 * (attempt + 1))
    return None


# ------------------------------------------------------------ query building --
def clean_rsid(v):
    if not isinstance(v, str):
        return None
    v = v.strip().lower()
    return v if re.fullmatch(r"rs\d+", v) else None


def clean_hgvs(v):
    """Return a ClinVar-searchable query string, or None if hopeless."""
    if not isinstance(v, str):
        return None
    v = v.strip().strip('"')
    if not v:
        return None
    v = re.sub(r",\s*TT?N?$", "", v)          # strip trailing ",TTN" / ",TT"
    v = v.replace(" ", "").strip("[]")
    v = v.replace("c.[", "c.").rstrip("]")

    # "c." followed by a genomic-scale number is a mislabelled GRCh37 coordinate.
    m = re.fullmatch(r"c\.(\d{8,})([ACGT])>([ACGT])", v, re.I)
    if m:
        # NC_000002.11 = chr2 on GRCh37; ClinVar indexes these g. descriptors.
        return f"NC_000002.11:g.{m.group(1)}{m.group(2).upper()}>{m.group(3).upper()}"

    if re.fullmatch(r"g\.\d+[ACGT]>[ACGT]", v, re.I):
        return None                            # no reference sequence, unusable
    return v


def search(term, gene_scoped=True):
    """Return list of ClinVar UIDs for a query term."""
    q = f'{term} AND {GENE}[gene]' if gene_scoped else term
    data = _get("esearch", {"term": q, "retmax": "20"})
    if not data:
        return []
    return data.get("esearchresult", {}).get("idlist", [])


def summarize(uids):
    """uid -> (GRCh38 start, chr, canonical variant name)."""
    out = {}
    for i in range(0, len(uids), 200):
        chunk = uids[i:i + 200]
        data = _get("esummary", {"id": ",".join(chunk)})
        if not data:
            continue
        res = data.get("result", {})
        for uid in res.get("uids", []):
            rec = res[uid]
            for vs in rec.get("variation_set", []):
                for loc in vs.get("variation_loc", []):
                    if loc.get("assembly_name") == ASSEMBLY and loc.get("start"):
                        out[uid] = (loc["start"], loc.get("chr"),
                                    vs.get("variation_name"))
                        break
                if uid in out:
                    break
        time.sleep(DELAY)
    return out


# ---------------------------------------------------------------------- main --
def main():
    df = pd.read_csv(SRC, encoding="latin-1")
    print(f"loaded {len(df)} rows")

    # build one query per row, remembering which field it came from
    queries = []
    for _, row in df.iterrows():
        rs = clean_rsid(row.get("rsid"))
        if rs:
            queries.append((rs, "rsid"))
            continue
        hg = clean_hgvs(row.get("hgvs"))
        queries.append((hg, "hgvs") if hg else (None, "none"))

    uniq = sorted({q for q, _ in queries if q})
    print(f"{len(uniq)} unique query terms to resolve")

    cache = {}
    if pathlib.Path(CACHE).exists():
        cache = json.loads(pathlib.Path(CACHE).read_text())
        print(f"resumed {len(cache)} cached lookups")

    todo = [q for q in uniq if q not in cache]
    for n, q in enumerate(todo, 1):
        uids = search(q)
        if not uids:                       # retry once without the gene filter
            time.sleep(DELAY)
            uids = search(q, gene_scoped=False)
        cache[q] = uids
        time.sleep(DELAY)
        if n % 100 == 0 or n == len(todo):
            print(f"  searched {n}/{len(todo)}")
            pathlib.Path(CACHE).write_text(json.dumps(cache))
    pathlib.Path(CACHE).write_text(json.dumps(cache))

    all_uids = sorted({u for v in cache.values() for u in v})
    print(f"fetching coordinates for {len(all_uids)} ClinVar records")
    coords = summarize(all_uids)

    ids, positions, chroms, names, status = [], [], [], [], []
    ambiguous = []
    for idx, (q, src) in enumerate(queries):
        if not q:
            ids.append(None); positions.append(None); chroms.append(None)
            names.append(None); status.append("no_query_key")
            continue
        uids = cache.get(q, [])
        if len(uids) == 1 and uids[0] in coords:
            start, chrom, name = coords[uids[0]]
            ids.append(uids[0]); positions.append(start); chroms.append(chrom)
            names.append(name); status.append(f"unique_by_{src}")
        elif len(uids) == 0:
            ids.append(None); positions.append(None); chroms.append(None)
            names.append(None); status.append(f"not_found_by_{src}")
        else:
            ids.append(None); positions.append(None); chroms.append(None)
            names.append(None); status.append(f"ambiguous_by_{src}({len(uids)})")
            ambiguous.append({
                "row": idx, "query": q, "source_field": src,
                "n_hits": len(uids), "candidate_ids": ";".join(uids),
                "candidate_names": " | ".join(
                    coords.get(u, (None, None, "?"))[2] for u in uids),
            })

    df["ClinVar_Variation_ID"] = ids
    df["position"] = positions
    df["chromosome"] = chroms
    df["clinvar_variant_name"] = names
    df["clinvar_match_status"] = status

    df.to_csv(OUT, index=False, encoding="utf-8-sig")
    if ambiguous:
        pd.DataFrame(ambiguous).to_csv(REVIEW, index=False, encoding="utf-8-sig")

    print(f"\nwrote {OUT}")
    print(f"filled position for {df['position'].notna().sum()} / {len(df)} rows")
    print("\nbreakdown:")
    print(df["clinvar_match_status"].str.replace(r"\(\d+\)", "", regex=True)
          .value_counts().to_string())
    if ambiguous:
        print(f"\n{len(ambiguous)} ambiguous rows written to {REVIEW}")


if __name__ == "__main__":
    main()
