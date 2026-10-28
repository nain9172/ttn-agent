"""TTN 變異 ID 的兩種寫法與互轉。

TTN 位於 chr2 負股。這個專案內部（outputs/ 的報告檔名、Evo2 計分、HGMD / PubTator /
PubMind 比對鍵、各種快取）一直用「轉錄方向」的鹼基寫 `2-pos-ref-alt`；但 VCF、ClinVar、
gnomAD 等外部工具的 `chr-pos-ref-alt` 一律是基因組正股。兩者位置相同、鹼基互補：

    內部（轉錄方向） 2-178526348-G-T   ⇔   對外（正股，VCF） 2-178526348-C-A

規則：使用者看到的與內部計算一律用轉錄方向（display_id == internal_id，鹼基與 data/gene.fna、
c. 寫法同一股）；基因組正股寫法只當作對照 VCF / ClinVar 的 vcf_id。
normalize() 兩種都能解析，用參考序列判斷：同一位置的正股鹼基與其互補一定不同，所以 ref 只會
符合其中一種，不會誤判（網站送出分析時只接受轉錄方向，見 spark_job.validate）。

只用標準函式庫（建站、啟動程式的 Spark 端都會 import），參考序列第一次用到時才讀。
"""

from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path
from typing import Dict, Optional

ROOT = Path(__file__).resolve().parent.parent
FASTA = ROOT / "data" / "sequence.fasta"
TTN_SEQUENCE_START = 178807423  # 參考序列第一個鹼基的座標（負股，座標由大到小）
TTN_SEQUENCE_END = 178525989
COMPLEMENT = str.maketrans("ACGT", "TGCA")
ID_RE = re.compile(r"^\s*(?:chr)?(\w+)[-:_](\d+)[-:_]([ACGTacgt])[-:_>]([ACGTacgt])\s*$")


class VariantIdError(ValueError):
    pass


def complement(base: str) -> str:
    return base.upper().translate(COMPLEMENT)


@lru_cache(maxsize=1)
def _sequence() -> str:
    lines = FASTA.read_text().splitlines()
    return "".join(l.strip() for l in lines if not l.startswith(">")).upper()


def transcript_base(pos: int) -> Optional[str]:
    """該位置在轉錄方向（參考序列）的鹼基；不在參考序列範圍內回 None。"""
    i = TTN_SEQUENCE_START - int(pos)
    seq = _sequence()
    return seq[i] if 0 <= i < len(seq) else None


def to_vcf(internal_id: str) -> str:
    """內部 ID（轉錄方向）→ 基因組正股（VCF / ClinVar）寫法。純互補，不需要參考序列。"""
    chrom, pos, ref, alt = internal_id.split("-")
    return f"{chrom}-{pos}-{complement(ref)}-{complement(alt)}"


def to_display(internal_id: str) -> str:
    """對外顯示的 ID：與內部 ID 相同（轉錄方向）。"""
    return internal_id


def normalize(text: str) -> Dict[str, object]:
    """解析使用者輸入的 `chr-pos-ref-alt`（正股或轉錄方向皆可），回傳兩種寫法。

    回傳 dict：chrom, pos, internal_id, display_id, tx_ref, tx_alt（轉錄方向），
    plus_ref, plus_alt（正股），orientation（'plus' 或 'transcript'，代表使用者用的是哪一種）。
    """
    m = ID_RE.match(text or "")
    if not m:
        raise VariantIdError("格式應為 chr-position-ref-alt，例如 2-178527098-G-A")
    chrom, pos, ref, alt = m.group(1).upper(), int(m.group(2)), m.group(3).upper(), m.group(4).upper()
    if chrom != "2":
        raise VariantIdError("TTN 位於 chr2")
    if ref == alt:
        raise VariantIdError("ref 與 alt 相同")
    tx = transcript_base(pos)
    if tx is None:
        # 參考序列外（TTN 兩端少數幾十 bp）無從判斷，依標準寫法視為正股；Evo2 在這裡本來就無法計分
        tx_ref, tx_alt = complement(ref), complement(alt)
        internal = f"2-{pos}-{tx_ref}-{tx_alt}"
        return {
            "chrom": "2", "pos": pos, "internal_id": internal, "display_id": to_display(internal),
            "tx_ref": tx_ref, "tx_alt": tx_alt, "plus_ref": ref, "plus_alt": alt,
            "orientation": "plus-assumed",
        }
    plus = complement(tx)
    if ref == plus:
        orientation, tx_ref, tx_alt = "plus", tx, complement(alt)
    elif ref == tx:
        orientation, tx_ref, tx_alt = "transcript", ref, alt
    else:  # pragma: no cover - 四種鹼基中 ref 必定是 plus、tx 或另外兩種之一
        raise VariantIdError(
            f"chr2:{pos:,} 的參考鹼基是 {plus}（基因組正股）；輸入的 ref {ref} 不符"
        )
    internal = f"2-{pos}-{tx_ref}-{tx_alt}"
    return {
        "chrom": "2", "pos": pos,
        "internal_id": internal, "display_id": to_display(internal),
        "tx_ref": tx_ref, "tx_alt": tx_alt,
        "plus_ref": complement(tx_ref), "plus_alt": complement(tx_alt),
        "orientation": orientation,
    }


def internal_id(text: str) -> str:
    """任一寫法 → 內部 ID；格式不對或不在 TTN 時丟 VariantIdError。"""
    return normalize(text)["internal_id"]


@lru_cache(maxsize=1)
def _transcript_model():
    from utils.hgvs_resolver import TTNTranscriptModel

    return TTNTranscriptModel(ROOT / "data" / "ttn_transcript_model.json", FASTA)


def hgvs_c(pos: int, tx_ref: str, tx_alt: str, transcript: str = "NM_001267550.2") -> Optional[str]:
    """轉錄方向的 c. 寫法（例：NM_001267550.2:c.107743A>C）；查不到時回 None。"""
    try:
        model = _transcript_model()
        c = model.genomic_to_c(transcript, int(pos))
        tx = model.transcripts[model.resolve_transcript_id(transcript) or transcript]
        cds_len = int(tx["cds_r_end"]) - int(tx["cds_r_start"]) + 1
    except Exception:  # noqa: BLE001
        return None
    if not c:
        return None
    # genomic_to_c 把 UTR 也寫成 CDS 座標；HGVS 規定 3′UTR 寫 c.*N、5′UTR 寫 c.-N
    m = re.match(r"c\.(-?\d+)(.*)$", c)
    if m:
        n, rest = int(m.group(1)), m.group(2)
        down = re.fullmatch(r"\+(\d+)", rest)
        if n > cds_len and down:
            # TTN 的 3′UTR 全在最後一個外顯子裡，所以 3′ 端之後的 +N 是基因下游：c.*(UTR 長度 + N)
            c = f"c.*{n - cds_len + int(down.group(1))}"
        elif n > cds_len:
            c = f"c.*{n - cds_len}{rest}"
        elif n <= 0:
            c = f"c.-{1 - n}{rest}"
    return f"{transcript}:{c}{tx_ref}>{tx_alt}"


def add_display_fields(variant_info: Dict) -> Dict:
    """在 pipeline 的 variant_info（內部 ID、轉錄方向鹼基）加上對外顯示用的欄位。"""
    vid = variant_info["variant_id"]
    variant_info.setdefault("display_id", to_display(vid))
    variant_info.setdefault("vcf_id", to_vcf(vid))
    variant_info.setdefault("plus_ref", complement(variant_info["ref"]))
    variant_info.setdefault("plus_alt", complement(variant_info["alt"]))
    if not variant_info.get("hgvs_c_display"):
        variant_info["hgvs_c_display"] = hgvs_c(variant_info["pos"], variant_info["ref"], variant_info["alt"])
    return variant_info
