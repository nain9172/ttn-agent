#!/usr/bin/env python3
"""
建立 TTN 轉錄本外顯子模型（HGVS c./p. ↔ genomic 座標轉換用）。

輸出 JSON（預設 data/ttn_transcript_model.json），內容為每個 RefSeq 轉錄本的：
  - strand、外顯子 genomic 區間（轉錄方向排序）、CDS genomic 區間
  - 由外顯子推導的 transcript(r.) 座標與 CDS 偏移（c.1 對應的 r. 座標）

有了這份模型，utils/hgvs_resolver.py 就能完全離線把
`c.98671-2A>G`、`c.29800T>C`、`p.R31056*` 換算成 chr2 genomic 座標。

資料來源（依序嘗試）：
  1. --gff 指定的本機 GFF3（NCBI RefSeq annotation，含 gene=TTN 的行即可）
  2. NCBI E-utilities efetch NC_000002.12 的 GFF3（會下載整條 chr2 的 annotation，
     約 110 MB，只解析 gene=TTN 的行；建議加 --gff 用本機檔案避免重複下載）

用法：
    conda run -n gm4 python build_ttn_transcript_model.py
    conda run -n gm4 python build_ttn_transcript_model.py --gff /tmp/ttn.gff3
    conda run -n gm4 python build_ttn_transcript_model.py --validate   # 用 HGMD 的 c.↔POS 驗證
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional

PROJECT_ROOT = Path(__file__).parent
DEFAULT_OUT = PROJECT_ROOT / "data" / "ttn_transcript_model.json"
CHR2_ACCESSION = "NC_000002.12"
GFF_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    f"?db=nuccore&id={CHR2_ACCESSION}&rettype=gff3&retmode=text"
)
# 只保留 TTN 的 RefSeq 轉錄本（XM_ 預測轉錄本不用於 HGVS 解析）
TRANSCRIPT_RE = re.compile(r"^N[MR]_")


def _attrs(field: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for kv in field.rstrip(";").split(";"):
        if "=" in kv:
            k, _, v = kv.partition("=")
            out[k] = v
    return out


def _iter_gff_lines(gff_path: Optional[Path]) -> Iterable[str]:
    if gff_path:
        opener = gzip.open if gff_path.suffix == ".gz" else open
        with opener(gff_path, "rt", encoding="utf-8", errors="replace") as fh:
            yield from fh
        return
    print(f"下載 chr2 GFF3（{CHR2_ACCESSION}）… 只解析 gene=TTN 的行")
    with urllib.request.urlopen(GFF_URL, timeout=1800) as resp:
        for raw in resp:
            yield raw.decode("utf-8", errors="replace")


def parse_gff(gff_path: Optional[Path]) -> Dict[str, Dict]:
    """從 GFF3 取出 TTN 每個 RefSeq 轉錄本的 exon / CDS genomic 區間。"""
    exons: Dict[str, List[List[int]]] = defaultdict(list)
    cds: Dict[str, List[List[int]]] = defaultdict(list)
    strand: Dict[str, str] = {}
    n_ttn = 0

    for line in _iter_gff_lines(gff_path):
        if line.startswith("#") or "gene=TTN" not in line:
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9 or cols[2] not in ("exon", "CDS"):
            continue
        a = _attrs(cols[8])
        if a.get("gene") != "TTN":
            continue
        tid = a.get("transcript_id") or a.get("Parent", "").replace("rna-", "")
        if not tid or not TRANSCRIPT_RE.match(tid):
            continue
        n_ttn += 1
        start, end = int(cols[3]), int(cols[4])
        strand[tid] = cols[6]
        (exons if cols[2] == "exon" else cds)[tid].append([start, end])

    print(f"解析到 {n_ttn} 行 TTN exon/CDS，轉錄本：{sorted(exons)}")
    return {"exons": exons, "cds": cds, "strand": strand}


def build_model(parsed: Dict[str, Dict]) -> Dict[str, Dict]:
    """把 exon/CDS 區間整理成含 r.（transcript）座標與 CDS 偏移的模型。"""
    model: Dict[str, Dict] = {}
    for tid, ex in parsed["exons"].items():
        st = parsed["strand"].get(tid, "-")
        reverse = st == "-"
        # 依「轉錄方向」排序外顯子：負鏈 = genomic 由大到小
        ex_sorted = sorted(ex, key=lambda p: p[0], reverse=reverse)

        blocks = []
        r = 1
        for g_start, g_end in ex_sorted:
            length = g_end - g_start + 1
            blocks.append(
                {"g_start": g_start, "g_end": g_end, "r_start": r, "r_end": r + length - 1}
            )
            r += length
        tx_len = r - 1

        cds_list = parsed["cds"].get(tid) or []
        cds_offset = None  # c.1 對應的 r. 座標
        cds_end_r = None
        if cds_list:
            c_min = min(p[0] for p in cds_list)
            c_max = max(p[1] for p in cds_list)
            first_g = c_max if reverse else c_min  # CDS 起點（轉錄方向）
            last_g = c_min if reverse else c_max
            cds_offset = _g_to_r(blocks, first_g, reverse)
            cds_end_r = _g_to_r(blocks, last_g, reverse)

        model[tid] = {
            "strand": st,
            "tx_length": tx_len,
            "cds_r_start": cds_offset,
            "cds_r_end": cds_end_r,
            "n_exons": len(blocks),
            "exons": blocks,
        }
    return model


def _g_to_r(blocks: List[Dict], g: int, reverse: bool) -> Optional[int]:
    for b in blocks:
        if b["g_start"] <= g <= b["g_end"]:
            return b["r_start"] + (b["g_end"] - g if reverse else g - b["g_start"])
    return None


def validate_with_hgmd(model_path: Path) -> None:
    """用 HGMD 的 (DNA c.HGVS, POS) 配對驗證 c.→genomic 轉換。"""
    import pandas as pd

    import config
    from utils.hgvs_resolver import TTNTranscriptModel

    tm = TTNTranscriptModel(model_path=model_path)
    df = pd.read_excel(config.HGMD_KB_PATH)
    ok = bad = skip = 0
    bad_examples: List[str] = []
    for _, row in df.iterrows():
        dna = str(row.get("DNA") or "")
        m = re.match(r"(N[MR]_[\d.]+):c\.(\d+)(?:([+-])(\d+))?", dna)
        if not m:
            skip += 1
            continue
        tid, cpos = m.group(1), int(m.group(2))
        offset = 0
        if m.group(3):
            offset = int(m.group(4)) * (1 if m.group(3) == "+" else -1)
        try:
            g = tm.c_to_genomic(tid, cpos, offset)
        except Exception as e:  # noqa: BLE001
            g, err = None, str(e)
        else:
            err = "no-mapping"
        if g is None:
            bad += 1
            if len(bad_examples) < 10:
                bad_examples.append(f"{dna} POS={row['POS']} ({err})")
            continue
        if int(g) == int(row["POS"]):
            ok += 1
        else:
            bad += 1
            if len(bad_examples) < 10:
                bad_examples.append(f"{dna} expect POS={row['POS']} got {g}")

    total = ok + bad
    print(f"\nHGMD c.→genomic 驗證：{ok}/{total} 正確（{ok / total:.1%}），無 c.HGVS 略過 {skip}")
    for e in bad_examples:
        print("  ✗", e)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", type=Path, default=None, help="本機 GFF3（含 gene=TTN 行）")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--validate", action="store_true", help="建好後用 HGMD 驗證")
    args = ap.parse_args()

    if args.gff and not args.gff.exists():
        raise SystemExit(f"找不到 GFF：{args.gff}")

    model = build_model(parse_gff(args.gff))
    if not model:
        raise SystemExit("沒有解析到任何 TTN 轉錄本，請檢查 GFF 來源。")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps({"assembly": "GRCh38", "chrom": "2", "accession": CHR2_ACCESSION,
                    "transcripts": model}, indent=1),
        encoding="utf-8",
    )
    print(f"\n寫入 {args.out}")
    for tid, m in sorted(model.items()):
        print(
            f"  {tid:<20} strand={m['strand']} exons={m['n_exons']:>3} "
            f"tx_len={m['tx_length']:>6} cds_r_start={m['cds_r_start']} "
            f"cds_len={(m['cds_r_end'] - m['cds_r_start'] + 1) if m['cds_r_end'] else None}"
        )

    if args.validate:
        validate_with_hgmd(args.out)


if __name__ == "__main__":
    main()
