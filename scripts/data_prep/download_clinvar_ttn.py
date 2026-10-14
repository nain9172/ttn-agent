#!/usr/bin/env python3
"""用 NCBI E-utilities 下載 ClinVar 目前的 TTN 單鹼基變異（pipeline 的輸入清單）。

預設範圍與論文 §3.1.1 的快照（2026-04-03，4,990 顆）相同：germline 分類為
Pathogenic / Likely pathogenic / Benign / Likely benign（含 P/LP、B/LB 合併），且 review
status 為 "criteria provided, multiple submitters, no conflicts"（2 星）。
  --all-review-status  連 1 星（single submitter）等也收（數量約多 3 倍）
  --include-vus        連 VUS、conflicting 等全部收
範圍越大，跑完整 pipeline 的時間越長。

輸出（預設 data/clinvar_ttn.csv）可直接給 batch_run_variants.py 與
score_evo2_for_variants.py 使用：

  input                   chr-pos-ref-alt，ref/alt 為 TTN 轉錄方向（負股）的鹼基，
                          與報告 variant ID 相同（例：genomic A>G 寫成 T>C）
  chr, position, ref, alt 同上，拆開的欄位
  Name, VariationID, Germline classification, Germline review status, Molecular consequence
  clinvar_label           pathogenic / benign / uncertain / conflicting / other

只保留 GRCh38、位在 TTN 參考序列範圍內、且 ref 與 data/sequence.fasta 相符的 SNV。
ClinVar 持續更新，每次下載的數量會不同；檔案裡記下了下載日期。

用法：python scripts/data_prep/download_clinvar_ttn.py [--out data/clinvar_ttn.csv] [--include-vus]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import pandas as pd  # noqa: E402

import config  # noqa: E402

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
COMPLEMENT = str.maketrans("ACGT", "TGCA")
BASE_TERM = 'TTN[gene] AND "single nucleotide variant"[Type of variation]'
PLP_BLB = ("clinsig_pathogenic[prop] OR clinsig_likely_pathogenic[prop] "
           "OR clinsig_benign[prop] OR clinsig_likely_benign[prop]")
KEEP_DEFAULT = {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic",
                "benign", "likely benign", "benign/likely benign"}
BATCH = 500


def label_of(sig: str) -> str:
    s = str(sig or "").lower()
    if "conflicting" in s:
        return "conflicting"
    if "pathogenic" in s and "benign" in s:
        return "conflicting"
    if "pathogenic" in s:
        return "pathogenic"
    if "benign" in s:
        return "benign"
    if "uncertain" in s:
        return "uncertain"
    return "other"


def eutils(endpoint: str, params: dict, retries: int = 5) -> dict:
    params = {**params, "retmode": "json"}
    if config.PUBMED_API_KEY:
        params["api_key"] = config.PUBMED_API_KEY
    if config.PUBMED_EMAIL:
        params["email"] = config.PUBMED_EMAIL
    data = urllib.parse.urlencode(params).encode()
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(f"{EUTILS}/{endpoint}", data=data, timeout=120) as resp:
                return json.load(resp)
        except Exception as e:  # noqa: BLE001
            if attempt == retries - 1:
                raise
            print(f"[clinvar] {endpoint} 失敗（{e}），重試…", flush=True)
            time.sleep(3 * (attempt + 1))
    return {}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=ROOT / "data" / "clinvar_ttn.csv")
    ap.add_argument("--include-vus", action="store_true",
                    help="連 VUS / conflicting 等全部下載（預設只取 P/LP/B/LB）")
    ap.add_argument("--all-review-status", action="store_true",
                    help="不限 review status（預設只取 multiple submitters, no conflicts）")
    args = ap.parse_args()

    term = BASE_TERM if args.include_vus else f"{BASE_TERM} AND ({PLP_BLB})"
    found = eutils("esearch.fcgi", {"db": "clinvar", "term": term, "usehistory": "y", "retmax": 0})
    res = found["esearchresult"]
    total, webenv, qkey = int(res["count"]), res["webenv"], res["querykey"]
    print(f"[clinvar] ClinVar 查到 {total:,} 筆 TTN SNV 紀錄，分 {(total + BATCH - 1) // BATCH} 批下載…",
          flush=True)

    rows = []
    delay = 0.12 if config.PUBMED_API_KEY else 0.4
    for start in range(0, total, BATCH):
        result = eutils("esummary.fcgi", {"db": "clinvar", "WebEnv": webenv, "query_key": qkey,
                                          "retstart": start, "retmax": BATCH})["result"]
        for uid in result.get("uids", []):
            r = result[uid]
            germ = r.get("germline_classification") or {}
            for vs in r.get("variation_set") or []:
                spdi = vs.get("canonical_spdi") or ""
                parts = spdi.split(":")
                if len(parts) != 4 or not parts[0].startswith("NC_000002"):
                    continue
                rows.append({
                    "pos": int(parts[1]) + 1,  # SPDI 是 0-based
                    "ref_plus": parts[2], "alt_plus": parts[3],
                    "Name": r.get("title", ""),
                    "VariationID": uid,
                    "Germline classification": germ.get("description", ""),
                    "Germline review status": germ.get("review_status", ""),
                    "Molecular consequence": "|".join(r.get("molecular_consequence_list") or []),
                })
        print(f"[clinvar]   {min(start + BATCH, total):,} / {total:,}", flush=True)
        time.sleep(delay)

    df = pd.DataFrame(rows)
    if not args.include_vus:
        df = df[df["Germline classification"].str.lower().isin(KEEP_DEFAULT)]
    if not args.all_review_status:
        df = df[df["Germline review status"] == "criteria provided, multiple submitters, no conflicts"]
    df = df[df.ref_plus.str.fullmatch("[ACGT]") & df.alt_plus.str.fullmatch("[ACGT]")]
    start_, end_ = config.TTN_SEQUENCE_START, config.TTN_SEQUENCE_END  # 負股：start > end
    df = df[df.pos.between(end_, start_)]

    # 正股 allele → TTN 轉錄方向，並以參考序列確認
    seq = "".join(l.strip() for l in config.REFERENCE_GENOME_PATH.read_text().splitlines()
                  if not l.startswith(">")).upper()
    df["ref"] = df.ref_plus.str.translate(COMPLEMENT)
    df["alt"] = df.alt_plus.str.translate(COMPLEMENT)
    ok = df.apply(lambda r: seq[start_ - r.pos] == r.ref, axis=1)
    if (~ok).any():
        print(f"[clinvar] 略過 {int((~ok).sum())} 筆 ref 與參考序列不符的紀錄", flush=True)
    df = df[ok]

    out = pd.DataFrame({
        "input": "2-" + df.pos.astype(str) + "-" + df.ref + "-" + df.alt,
        "chr": "2", "position": df.pos, "ref": df.ref, "alt": df.alt,
        "Name": df.Name, "VariationID": df.VariationID,
        "Germline classification": df["Germline classification"],
        "Germline review status": df["Germline review status"],
        "Molecular consequence": df["Molecular consequence"],
        "clinvar_label": df["Germline classification"].map(label_of),
        "download_date": date.today().isoformat(),
    }).sort_values("position").drop_duplicates("input")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"[clinvar] {date.today()} 快照：{len(out):,} 顆 TTN SNV → {args.out}")
    print("[clinvar] " + ", ".join(f"{k} {v:,}" for k, v in out.clinvar_label.value_counts().items()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
