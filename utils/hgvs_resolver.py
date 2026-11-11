#!/usr/bin/env python3
"""
TTN HGVS → genomic 座標解析器（離線）。

讓整條 pipeline（尤其 Track B 機制推理）除了 `2-178527098-C-T` 之外，
也能直接吃文獻裡最常見的 HGVS 寫法：

    c.98671-2A>G            （intronic / canonical splice site）
    NM_001267550.2:c.29800T>C
    c.104509_104510del      （frameshift）
    p.R31056*  /  p.(Cys13708Trp)  /  p.Arg23029ThrfsX9
    g.178704672A>G  /  NC_000002.12:g.178704672A>G
    2-178527098-C-T

做法：
  1. `data/ttn_transcript_model.json`（由 build_ttn_transcript_model.py 從 NCBI RefSeq
     GFF3 產生）提供每個 TTN RefSeq 轉錄本的外顯子 genomic 區間與 CDS 偏移，
     可把 c./p. 位置換算成 chr2 genomic 座標，intronic 的 ±N 偏移也能處理。
  2. 文獻常省略轉錄本編號，而 TTN 各 isoform 的 c. 編號完全不同。因此會逐一試
     候選轉錄本，並用兩個硬條件驗證：
       - intronic 偏移必須落在真正的外顯子邊界（±1/±2 才是 canonical splice site）
       - 參考鹼基必須與基因組序列一致（data/sequence.fasta，TTN 負鏈已 reverse
         complement，方向與 c./p. 相同）
     p. 則要求該位置轉譯出的胺基酸與 HGVS 的參考胺基酸相符。
  3. 同時由 CDS 轉譯推出蛋白質後果（missense / nonsense / synonymous），
     所以只給 c. 也能定出變異類型，不需要呼叫 VEP。

座標慣例（與 HGMD KB / evo2_predictor / 既有 variant_id 一致）：
  pos 為 chr2 GRCh38 genomic 位置，ref/alt 為「轉錄本方向」鹼基（TTN 在負鏈，
  故與 c.HGVS 直接一致）。plus_ref / plus_alt 另外提供正鏈（VCF）表示。
"""

from __future__ import annotations

import json
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_MODEL_PATH = PROJECT_ROOT / "data" / "ttn_transcript_model.json"
DEFAULT_FASTA_PATH = PROJECT_ROOT / "data" / "sequence.fasta"

# sequence.fasta 的 index 0 對應的 genomic 位置（負鏈，往下遞減）
GENE_SENSE_START = 178807423
GENE_SENSE_END = 178525989
# TTN 區間內 GRCh37 → GRCh38 的位移是常數（hg19 179,390,716 = hg38 178,525,989）。
# 舊文獻常給 hg19 的 g. 座標，換算後仍會用參考鹼基驗證，錯了就不採用。
GRCH37_TO_38_OFFSET = 864_727

# 文獻省略轉錄本時的嘗試順序：meta-transcript（MANE / HGMD 主用）優先
DEFAULT_TRANSCRIPT_PRIORITY = [
    "NM_001267550.2",  # meta transcript（最長，含幾乎所有外顯子）
    "NM_133378.4",     # N2A（骨骼肌）
    "NM_001256850.1",  # N2BA（心肌，長）
    "NM_003319.4",     # N2B（心肌，短）
    "NM_133432.3",
    "NM_133437.4",
    "NM_133379.5",
]

_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}

AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V", "Ter": "*", "Sec": "U", "Xaa": "X",
}
AA1 = set("ACDEFGHIKLMNPQRSTVWY*")

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def complement(base: str) -> str:
    return "".join(_COMPLEMENT.get(b, "N") for b in str(base).upper()[::-1])


class HgvsResolveError(Exception):
    """HGVS 無法解析成 genomic 座標。"""


# ── HGVS 字串解析 ───────────────────────────────────────────────────────────


@dataclass
class ParsedHgvs:
    """解析後的 HGVS（尚未換算 genomic 座標）。"""

    kind: str                  # 'c' | 'p' | 'g' | 'vcf'
    raw: str
    transcript: Optional[str] = None
    # c./g. 用
    start: Optional[int] = None
    start_offset: int = 0
    end: Optional[int] = None
    end_offset: int = 0
    edit: Optional[str] = None  # 'sub' | 'del' | 'dup' | 'ins' | 'delins' | 'unknown'
    ref: Optional[str] = None
    alt: Optional[str] = None
    utr: Optional[str] = None   # '5' (c.-N) | '3' (c.*N)
    # p. 用
    aa_pos: Optional[int] = None
    aa_ref: Optional[str] = None
    aa_alt: Optional[str] = None
    p_kind: Optional[str] = None  # 'missense' | 'nonsense' | 'frameshift' | 'synonymous' | 'other'
    # vcf 用
    pos: Optional[int] = None


def _norm_text(text: str) -> str:
    """去掉空白與常見雜訊，統一大於符號。"""
    s = str(text or "").strip()
    s = s.replace("→", ">").replace("＞", ">").replace("∗", "*")
    s = re.sub(r"\(Expired\)", "", s, flags=re.IGNORECASE)
    s = re.sub(r"\s+", "", s)
    # ClinVar 標題寫法：NM_001267550.2(TTN):c.107909C>T(p.Thr35970Ile) → NM_001267550.2:c.107909C>T
    s = re.sub(r"^([A-Z]{2}_\d+(?:\.\d+)?)\([A-Za-z0-9-]+\)(?=:)", r"\1", s)
    s = re.sub(r"(?<=.)\(p\.[^)]*\)$", "", s)
    # 文獻常見的千分位逗號：c.104,515C>T → c.104515C>T
    s = re.sub(r"(?<=\d),(?=\d{3}(?!\d))", "", s)
    return s


def _aa_token_to_1(token: Optional[str]) -> Optional[str]:
    if not token:
        return None
    t = token.strip()
    if t in ("*", "X", "Ter", "ter"):
        return "*"
    if len(t) == 1 and t.upper() in AA1:
        return t.upper()
    key = t[:1].upper() + t[1:3].lower()
    return AA3_TO_1.get(key)


_C_SUB_RE = re.compile(
    r"^c\.(?P<pos>[*-]?\d+)(?P<off>[+-]\d+)?(?P<ref>[ACGT]+)>(?P<alt>[ACGT]+)$",
    re.IGNORECASE,
)
_C_INDEL_RE = re.compile(
    r"^c\.(?P<pos>[*-]?\d+)(?P<off>[+-]\d+)?"
    r"(?:_(?P<pos2>[*-]?\d+)(?P<off2>[+-]\d+)?)?"
    r"(?P<edit>delins|del|dup|ins)(?P<seq>[ACGT]*)",
    re.IGNORECASE,
)
_G_RE = re.compile(
    r"^g\.(?P<pos>\d+)(?:(?P<ref>[ACGT]+)>(?P<alt>[ACGT]+)|(?P<edit>delins|del|dup|ins)(?P<seq>[ACGT]*))?",
    re.IGNORECASE,
)
_P_RE = re.compile(
    r"^p\.\(?(?P<ref>[A-Z][a-z]{2}|[A-Z])(?P<pos>\d+)(?P<rest>.*?)\)?$"
)
_VCF_RE = re.compile(r"^(?:chr)?(?P<chrom>\d+|[XYM]T?)[-_:](?P<pos>\d+)[-_:](?P<ref>[ACGT]+)[-_:](?P<alt>[ACGT]+)$",
                     re.IGNORECASE)


def parse_hgvs(text: str) -> ParsedHgvs:
    """把各種變異寫法解析成 ParsedHgvs（不查座標）。"""
    raw = str(text or "").strip()
    s = _norm_text(raw)
    if not s:
        raise HgvsResolveError("空的變異字串")

    transcript = None
    # NM_001267550.2:c.123A>G / NC_000002.12(NM_...):c.123A>G / NP_...:p.xxx
    m = re.match(r"^(?P<acc>[A-Z]{2}_[\d.]+)(?:\((?P<sel>[A-Z]{2}_[\d.]+)\))?:(?P<rest>.+)$", s)
    if m:
        acc, sel, s = m.group("acc"), m.group("sel"), m.group("rest")
        transcript = sel or (acc if acc.startswith(("NM_", "NR_", "XM_")) else None)

    low = s.lower()

    if low.startswith("c.") or low.startswith("n."):
        s = "c." + s[2:]
        sub = _C_SUB_RE.match(s)
        if sub:
            pos, utr = _split_utr(sub.group("pos"))
            return ParsedHgvs(
                kind="c", raw=raw, transcript=transcript, start=pos, utr=utr,
                start_offset=int(sub.group("off") or 0), edit="sub",
                ref=sub.group("ref").upper(), alt=sub.group("alt").upper(),
            )
        indel = _C_INDEL_RE.match(s)
        if indel:
            pos, utr = _split_utr(indel.group("pos"))
            pos2, _ = _split_utr(indel.group("pos2")) if indel.group("pos2") else (None, None)
            edit = indel.group("edit").lower()
            seq = (indel.group("seq") or "").upper()
            return ParsedHgvs(
                kind="c", raw=raw, transcript=transcript, start=pos, utr=utr,
                start_offset=int(indel.group("off") or 0), end=pos2,
                end_offset=int(indel.group("off2") or 0), edit=edit,
                ref=seq if edit in ("del", "dup") else None,
                alt=seq if edit in ("ins", "delins") else None,
            )
        # 只有位置（例如 c.12345）→ 仍可定位
        only = re.match(r"^c\.([*-]?\d+)([+-]\d+)?$", s)
        if only:
            pos, utr = _split_utr(only.group(1))
            return ParsedHgvs(kind="c", raw=raw, transcript=transcript, start=pos,
                              utr=utr, start_offset=int(only.group(2) or 0), edit="unknown")
        raise HgvsResolveError(f"無法解析 c.HGVS：{raw}")

    if low.startswith("g.") or low.startswith("m."):
        g = _G_RE.match("g." + s[2:])
        if not g:
            raise HgvsResolveError(f"無法解析 g.HGVS：{raw}")
        return ParsedHgvs(
            kind="g", raw=raw, transcript=transcript, pos=int(g.group("pos")),
            edit="sub" if g.group("ref") else (g.group("edit") or "unknown"),
            ref=(g.group("ref") or "").upper() or None,
            alt=(g.group("alt") or "").upper() or None,
        )

    if low.startswith("p."):
        p = _P_RE.match(s)
        if not p:
            raise HgvsResolveError(f"無法解析 p.HGVS：{raw}")
        aa_ref = _aa_token_to_1(p.group("ref"))
        rest = p.group("rest") or ""
        p_kind, aa_alt = _classify_protein_change(rest)
        return ParsedHgvs(
            kind="p", raw=raw, transcript=transcript, aa_pos=int(p.group("pos")),
            aa_ref=aa_ref, aa_alt=aa_alt, p_kind=p_kind,
        )

    v = _VCF_RE.match(s)
    if v:
        return ParsedHgvs(
            kind="vcf", raw=raw, pos=int(v.group("pos")), edit="sub",
            ref=v.group("ref").upper(), alt=v.group("alt").upper(),
        )

    raise HgvsResolveError(f"不認得的變異寫法：{raw}")


def _split_utr(token: Optional[str]) -> Tuple[Optional[int], Optional[str]]:
    if token is None:
        return None, None
    t = str(token)
    if t.startswith("*"):
        return int(t[1:]), "3"
    if t.startswith("-"):
        return int(t[1:]), "5"
    return int(t), None


def _classify_protein_change(rest: str) -> Tuple[str, Optional[str]]:
    """由 p. 的尾段判斷類型：p.R31056* / p.L34837EfsX12 / p.I35986T / p.Q123=。"""
    r = (rest or "").strip()
    if not r:
        return "other", None
    if re.search(r"fs", r, re.IGNORECASE):
        return "frameshift", None
    if r in ("=", "="):
        return "synonymous", None
    if r in ("*", "X") or r.lower().startswith("ter"):
        return "nonsense", "*"
    if re.match(r"^(?:del|dup|ins)", r, re.IGNORECASE):
        return "other", None
    alt = _aa_token_to_1(r[:3] if len(r) >= 3 and r[:3].isalpha() else r[:1])
    if alt == "*":
        return "nonsense", "*"
    if alt:
        return "missense", alt
    return "other", None


# ── 轉錄本模型 ──────────────────────────────────────────────────────────────


class TTNTranscriptModel:
    """TTN RefSeq 轉錄本的 c./p. ↔ genomic 座標換算（離線）。"""

    def __init__(
        self,
        model_path: Optional[Path] = None,
        fasta_path: Optional[Path] = None,
    ):
        self.model_path = Path(model_path or DEFAULT_MODEL_PATH)
        self.fasta_path = Path(fasta_path or DEFAULT_FASTA_PATH)
        self._model: Optional[Dict[str, Any]] = None
        self._seq: Optional[str] = None
        self._cds_cache: Dict[str, str] = {}
        self._protein_cache: Dict[str, str] = {}

    # -- 載入 ---------------------------------------------------------------

    @property
    def transcripts(self) -> Dict[str, Any]:
        if self._model is None:
            if not self.model_path.exists():
                raise HgvsResolveError(
                    f"找不到轉錄本模型 {self.model_path}；"
                    "請先執行 build_ttn_transcript_model.py 建立。"
                )
            self._model = json.loads(self.model_path.read_text(encoding="utf-8"))
        return self._model["transcripts"]

    def has_transcript(self, tid: str) -> bool:
        return tid in self.transcripts

    def resolve_transcript_id(self, tid: Optional[str]) -> Optional[str]:
        """允許不帶版本號（NM_001267550 → NM_001267550.2）。"""
        if not tid:
            return None
        if tid in self.transcripts:
            return tid
        base = tid.split(".")[0]
        for known in self.transcripts:
            if known.split(".")[0] == base:
                return known
        return None

    @property
    def sequence(self) -> Optional[str]:
        """TTN 基因方向（負鏈已 reverse complement）的參考序列。"""
        if self._seq is None:
            if not self.fasta_path.exists():
                logger.warning(f"參考序列不存在：{self.fasta_path}（將略過鹼基驗證）")
                self._seq = ""
            else:
                lines = []
                with open(self.fasta_path, "rt", encoding="utf-8") as fh:
                    for line in fh:
                        if line.startswith(">"):
                            continue
                        lines.append(line.strip())
                self._seq = "".join(lines).upper()
        return self._seq or None

    def base_at(self, pos: int) -> Optional[str]:
        """genomic 位置的「轉錄本方向」鹼基（TTN 負鏈 → 與 c.HGVS 同向）。"""
        seq = self.sequence
        if not seq:
            return None
        idx = GENE_SENSE_START - int(pos)
        if idx < 0 or idx >= len(seq):
            return None
        return seq[idx]

    def bases_at(self, pos: int, length: int) -> Optional[str]:
        """從 pos 起、沿轉錄方向取 length 個鹼基。"""
        seq = self.sequence
        if not seq:
            return None
        idx = GENE_SENSE_START - int(pos)
        if idx < 0 or idx + length > len(seq):
            return None
        return seq[idx: idx + length]

    # -- 座標換算 -----------------------------------------------------------

    def _tx(self, tid: str) -> Dict[str, Any]:
        t = self.transcripts.get(tid)
        if not t:
            raise HgvsResolveError(f"模型中沒有轉錄本 {tid}")
        return t

    def c_to_r(self, tid: str, cpos: int, utr: Optional[str] = None) -> int:
        """c. 位置 → r.（transcript）位置。"""
        t = self._tx(tid)
        cds_start, cds_end = t["cds_r_start"], t["cds_r_end"]
        if cds_start is None:
            raise HgvsResolveError(f"{tid} 沒有 CDS 資訊")
        if utr == "5":
            return cds_start - int(cpos)
        if utr == "3":
            return cds_end + int(cpos)
        return cds_start + int(cpos) - 1

    def _block_for_r(self, tid: str, r: int) -> Optional[Tuple[int, Dict[str, Any]]]:
        for i, b in enumerate(self._tx(tid)["exons"]):
            if b["r_start"] <= r <= b["r_end"]:
                return i, b
        return None

    def c_to_genomic(
        self,
        tid: str,
        cpos: int,
        offset: int = 0,
        utr: Optional[str] = None,
        strict_boundary: bool = True,
    ) -> Optional[int]:
        """
        c. 位置（可帶 intronic 偏移）→ chr2 genomic 位置。

        strict_boundary=True 時，帶偏移的 c.X-N / c.X+N 必須真的落在外顯子邊界
        （這是判斷「文獻用的是哪個轉錄本編號」最有力的線索）。
        """
        tid = self.resolve_transcript_id(tid) or tid
        try:
            r = self.c_to_r(tid, cpos, utr)
        except HgvsResolveError:
            return None
        found = self._block_for_r(tid, r)
        if not found:
            return None
        idx, b = found
        minus = self._tx(tid)["strand"] == "-"

        if not offset:
            return b["g_end"] - (r - b["r_start"]) if minus else b["g_start"] + (r - b["r_start"])

        exons = self._tx(tid)["exons"]
        if offset < 0:
            # intron 在該外顯子之前（轉錄方向）
            if strict_boundary and (r != b["r_start"] or idx == 0):
                return None
            return b["g_end"] - offset if minus else b["g_start"] + offset
        if strict_boundary and (r != b["r_end"] or idx == len(exons) - 1):
            return None
        return b["g_start"] - offset if minus else b["g_end"] + offset

    def genomic_to_c(self, tid: str, pos: int) -> Optional[str]:
        """genomic → c. 表示（外顯子內才給精確值，intron 給最近邊界 ±N）。"""
        tid = self.resolve_transcript_id(tid) or tid
        t = self._tx(tid)
        minus = t["strand"] == "-"
        cds_start = t["cds_r_start"]
        for i, b in enumerate(t["exons"]):
            if b["g_start"] <= pos <= b["g_end"]:
                r = b["r_start"] + (b["g_end"] - pos if minus else pos - b["g_start"])
                return f"c.{r - cds_start + 1}"
        # 落在 intron：找最近的外顯子邊界
        best = None
        for b in t["exons"]:
            for g_edge, is_start in ((b["g_start"], not minus), (b["g_end"], minus)):
                d = abs(pos - g_edge)
                if best is None or d < best[0]:
                    r = b["r_start"] if is_start else b["r_end"]
                    sign = "-" if is_start else "+"
                    best = (d, f"c.{r - cds_start + 1}{sign}{d}")
        return best[1] if best else None

    # -- 序列 / 蛋白質 ------------------------------------------------------

    def cds_sequence(self, tid: str) -> Optional[str]:
        """轉錄本 CDS 序列（轉錄方向，A of ATG 起）。"""
        tid = self.resolve_transcript_id(tid) or tid
        if tid in self._cds_cache:
            return self._cds_cache[tid]
        seq = self.sequence
        if not seq:
            return None
        t = self._tx(tid)
        cds_start, cds_end = t["cds_r_start"], t["cds_r_end"]
        if cds_start is None or cds_end is None:
            return None
        minus = t["strand"] == "-"
        chunks: List[str] = []
        for b in t["exons"]:
            lo, hi = max(b["r_start"], cds_start), min(b["r_end"], cds_end)
            if lo > hi:
                continue
            g_first = (b["g_end"] - (lo - b["r_start"])) if minus else (b["g_start"] + (lo - b["r_start"]))
            piece = self.bases_at(g_first, hi - lo + 1) if minus else None
            if piece is None:
                return None
            chunks.append(piece)
        cds = "".join(chunks)
        self._cds_cache[tid] = cds
        return cds

    def protein_sequence(self, tid: str) -> Optional[str]:
        tid = self.resolve_transcript_id(tid) or tid
        if tid in self._protein_cache:
            return self._protein_cache[tid]
        cds = self.cds_sequence(tid)
        if not cds:
            return None
        aa = "".join(
            CODON_TABLE.get(cds[i: i + 3], "X") for i in range(0, len(cds) - 2, 3)
        )
        self._protein_cache[tid] = aa
        return aa

    def aa_at(self, tid: str, aa_pos: int) -> Optional[str]:
        prot = self.protein_sequence(tid)
        if not prot or aa_pos < 1 or aa_pos > len(prot):
            return None
        return prot[aa_pos - 1]

    def predict_substitution_consequence(
        self, tid: str, cpos: int, ref: str, alt: str
    ) -> Dict[str, Any]:
        """CDS 內單鹼基置換 → 蛋白質後果（離線轉譯，不需 VEP）。"""
        out: Dict[str, Any] = {"vclass": "unknown", "hgvs_p": None}
        cds = self.cds_sequence(tid)
        if not cds or cpos < 1 or cpos > len(cds):
            return out
        if cds[cpos - 1].upper() != str(ref).upper():
            out["ref_mismatch"] = True
            return out
        aa_pos = (cpos - 1) // 3 + 1
        codon_start = (aa_pos - 1) * 3
        codon = cds[codon_start: codon_start + 3]
        if len(codon) < 3:
            return out
        i = (cpos - 1) % 3
        new_codon = codon[:i] + str(alt).upper() + codon[i + 1:]
        aa_ref = CODON_TABLE.get(codon, "X")
        aa_alt = CODON_TABLE.get(new_codon, "X")
        out.update({"aa_pos": aa_pos, "aa_ref": aa_ref, "aa_alt": aa_alt})
        if aa_ref == aa_alt:
            out["vclass"] = "synonymous"
            out["hgvs_p"] = f"p.{aa_ref}{aa_pos}="
        elif aa_alt == "*":
            out["vclass"] = "nonsense"
            out["hgvs_p"] = f"p.{aa_ref}{aa_pos}*"
        elif aa_ref == "*":
            out["vclass"] = "stop_loss"
            out["hgvs_p"] = f"p.*{aa_pos}{aa_alt}"
        else:
            out["vclass"] = "missense"
            out["hgvs_p"] = f"p.{aa_ref}{aa_pos}{aa_alt}"
        return out

    def aa_to_genomic(self, tid: str, aa_pos: int) -> Optional[int]:
        """胺基酸位置 → 該 codon 第一個鹼基的 genomic 位置。"""
        return self.c_to_genomic(tid, (int(aa_pos) - 1) * 3 + 1)


# ── 解析器 ──────────────────────────────────────────────────────────────────


@dataclass
class ResolvedVariant:
    """HGVS 解析結果（座標慣例與 HGMD KB / pipeline 一致）。"""

    variant_id: str
    pos: int
    ref: Optional[str]
    alt: Optional[str]
    transcript: Optional[str]
    hgvs_c: Optional[str] = None
    hgvs_p: Optional[str] = None
    vclass: str = "unknown"
    vgroup: str = "unknown"
    resolved_by: str = ""
    confidence: str = "low"        # high | medium | low
    warnings: List[str] = field(default_factory=list)
    chrom: str = "2"
    input: str = ""

    @property
    def plus_ref(self) -> Optional[str]:
        return complement(self.ref) if self.ref else None

    @property
    def plus_alt(self) -> Optional[str]:
        return complement(self.alt) if self.alt else None

    def to_variant_info(self) -> Dict[str, Any]:
        """轉成 pipeline 慣用的 variant_info dict。"""
        return {
            "variant_id": self.variant_id,
            "chrom": self.chrom,
            "pos": self.pos,
            "ref": self.ref,
            "alt": self.alt,
            "transcript": self.transcript,
            "hgvs_c": self.hgvs_c,
            "hgvs_p": self.hgvs_p,
            "hgvs_vclass": self.vclass,
            "hgvs_vgroup": self.vgroup,
            "hgvs_resolved_by": self.resolved_by,
            "hgvs_confidence": self.confidence,
            "hgvs_warnings": list(self.warnings),
            "input": self.input,
        }


_TRUNCATING = {"nonsense", "frameshift", "splice"}


def _vgroup_for(vclass: str) -> str:
    if vclass in _TRUNCATING:
        return "truncating"
    if vclass == "missense":
        return "missense"
    if vclass in ("synonymous", "inframe", "stop_loss", "intron", "utr", "splice_region"):
        return "other"
    return "unknown"


class HgvsResolver:
    """把各種 TTN 變異寫法解析成 genomic 座標 + 變異類型。"""

    def __init__(
        self,
        model_path: Optional[Path] = None,
        fasta_path: Optional[Path] = None,
        transcript_priority: Optional[List[str]] = None,
        hgmd_index: Optional[Dict[str, Tuple[int, str, str]]] = None,
        allow_network: bool = False,
        cache_path: Optional[Path] = None,
    ):
        self.model = TTNTranscriptModel(model_path=model_path, fasta_path=fasta_path)
        self.transcript_priority = list(transcript_priority or DEFAULT_TRANSCRIPT_PRIORITY)
        # 可選：HGMD 的 {正規化 HGVS -> (pos, ref, alt)}，命中即直接採用
        self.hgmd_index = hgmd_index or {}
        # rsID 沒有座標資訊，只能查 dbSNP（唯一需要連網的路徑，結果會快取）
        self.allow_network = bool(allow_network)
        self.cache_path = Path(cache_path) if cache_path else (
            PROJECT_ROOT / "outputs" / "hgvs_rsid_cache.json"
        )
        self._rsid_cache: Optional[Dict[str, Any]] = None

    # -- 對外入口 -----------------------------------------------------------

    def resolve(self, query: str, hgvs_p: Optional[str] = None) -> ResolvedVariant:
        """
        解析單一變異。query 可為 c./p./g. HGVS、`2-178527098-C-T` 或 rsID
        （rsID 需 allow_network=True）。
        另外給 hgvs_p 時，會用它交叉驗證 c. 該用哪個轉錄本編號，並補上蛋白質類型。
        """
        if re.match(r"^rs\d+$", _norm_text(query), re.IGNORECASE):
            return self._resolve_rsid(_norm_text(query))
        parsed = parse_hgvs(query)
        p_extra = None
        if hgvs_p:
            try:
                p_extra = parse_hgvs(hgvs_p)
                if p_extra.kind != "p":
                    p_extra = None
            except HgvsResolveError:
                p_extra = None

        if parsed.kind == "vcf":
            return self._resolve_vcf(parsed)
        if parsed.kind == "g":
            return self._resolve_genomic(parsed)
        if parsed.kind == "c":
            return self._resolve_cdna(parsed, p_extra)
        if parsed.kind == "p":
            return self._resolve_protein(parsed)
        raise HgvsResolveError(f"不支援的 HGVS 類型：{parsed.kind}")

    def try_resolve(self, query: str, hgvs_p: Optional[str] = None) -> Optional[ResolvedVariant]:
        try:
            return self.resolve(query, hgvs_p=hgvs_p)
        except HgvsResolveError as e:
            logger.debug(f"HGVS 解析失敗 {query}: {e}")
            return None

    def resolve_best(self, candidates: List[Optional[str]]) -> Optional[ResolvedVariant]:
        """依序嘗試多個寫法（例如 c. → p. → g. → rsid 之外的別名），取第一個成功者。"""
        cands = [c for c in candidates if c and str(c).strip().lower() != "nan"]
        prot = next((c for c in cands if _norm_text(c).lower().startswith("p.")), None)
        best: Optional[ResolvedVariant] = None
        for c in cands:
            r = self.try_resolve(c, hgvs_p=prot if not _norm_text(c).lower().startswith("p.") else None)
            if r is None:
                continue
            if r.confidence == "high":
                return r
            if best is None:
                best = r
        return best

    # -- 各種輸入 -----------------------------------------------------------

    def _finish(
        self,
        pos: int,
        ref: Optional[str],
        alt: Optional[str],
        parsed: ParsedHgvs,
        transcript: Optional[str],
        vclass: str,
        hgvs_p: Optional[str],
        resolved_by: str,
        confidence: str,
        warnings: List[str],
    ) -> ResolvedVariant:
        vid = (
            f"2-{pos}-{ref}-{alt}" if ref and alt and len(ref) == 1 and len(alt) == 1
            else f"2-{pos}-{(ref or '.')}-{(alt or '.')}"
        )
        hgvs_c = None
        if parsed.kind == "c":
            hgvs_c = _norm_text(parsed.raw).split(":")[-1]
        return ResolvedVariant(
            variant_id=vid, pos=int(pos), ref=ref, alt=alt, transcript=transcript,
            hgvs_c=hgvs_c, hgvs_p=hgvs_p, vclass=vclass, vgroup=_vgroup_for(vclass),
            resolved_by=resolved_by, confidence=confidence, warnings=warnings,
            input=parsed.raw,
        )

    def _resolve_vcf(self, parsed: ParsedHgvs) -> ResolvedVariant:
        pos, ref, alt = parsed.pos, parsed.ref, parsed.alt
        warnings: List[str] = []
        base = self.model.base_at(pos)
        # 既有 pipeline 的 variant_id 是轉錄本方向（負鏈）；若對不上就試正鏈
        if base and ref and base != ref:
            if base == complement(ref):
                ref, alt = complement(ref), complement(alt or "")
                warnings.append("輸入的 ref/alt 看起來是正鏈，已轉成轉錄本方向")
            else:
                warnings.append(f"參考鹼基不符：基因組為 {base}，輸入為 {ref}")
        vclass, hgvs_p = self._consequence_from_genomic(pos, ref, alt)
        return self._finish(
            pos, ref, alt, parsed, None, vclass, hgvs_p, "vcf_coordinates",
            "high" if base and base == ref else "low", warnings,
        )

    def _resolve_genomic(self, parsed: ParsedHgvs) -> ResolvedVariant:
        pos = int(parsed.pos or 0)
        warnings: List[str] = []
        ref, alt = parsed.ref, parsed.alt
        if not (GENE_SENSE_END <= pos <= GENE_SENSE_START):
            # 舊文獻常給 GRCh37/hg19 座標（TTN 在 chr2:179.3–179.7 Mb）→ 常數位移換算，
            # 換算後必須通得過參考鹼基驗證才採用
            lifted = pos - GRCH37_TO_38_OFFSET
            base37 = self.model.base_at(lifted)
            ok = base37 is not None and (
                not ref or base37 == ref or base37 == complement(ref)
            )
            if not ok:
                hint = ("（看起來是 GRCh37/hg19 座標，但換算後參考鹼基仍不符）"
                        if GENE_SENSE_END <= lifted <= GENE_SENSE_START else "")
                raise HgvsResolveError(f"g.{pos} 不在 GRCh38 的 TTN 區間內{hint}")
            warnings.append(f"g.{pos} 判定為 GRCh37 座標，已換算成 GRCh38 g.{lifted}")
            pos = lifted
        base = self.model.base_at(pos)
        if ref and base and base != ref:
            if base == complement(ref):
                ref, alt = complement(ref), complement(alt or "")
                warnings.append("g.HGVS 為正鏈，已轉成轉錄本方向")
            else:
                warnings.append(f"參考鹼基不符：基因組為 {base}，輸入為 {ref}")
        vclass, hgvs_p = self._consequence_from_genomic(pos, ref, alt)
        return self._finish(
            pos, ref, alt, parsed, None, vclass, hgvs_p, "genomic_hgvs",
            "high" if ref and base == ref else "medium", warnings,
        )

    def _consequence_from_genomic(
        self, pos: int, ref: Optional[str], alt: Optional[str]
    ) -> Tuple[str, Optional[str]]:
        """已知 genomic 座標時，用主要轉錄本推類型（外顯子內轉譯 / 否則 intron）。"""
        if not ref or not alt or len(ref) != 1 or len(alt) != 1:
            return "unknown", None
        for tid in self.transcript_priority:
            if not self.model.has_transcript(tid):
                continue
            c_desc = self.model.genomic_to_c(tid, pos)
            if not c_desc:
                continue
            m = re.match(r"^c\.(\d+)$", c_desc)
            if not m:
                # intronic：±1/±2 視為 canonical splice site
                m2 = re.match(r"^c\.(\d+)([+-])(\d+)$", c_desc)
                if m2 and int(m2.group(3)) <= 2:
                    return "splice", None
                continue
            cons = self.model.predict_substitution_consequence(tid, int(m.group(1)), ref, alt)
            if cons["vclass"] != "unknown":
                return cons["vclass"], cons["hgvs_p"]
        return "unknown", None

    def _resolve_cdna(
        self, parsed: ParsedHgvs, p_extra: Optional[ParsedHgvs] = None
    ) -> ResolvedVariant:
        """c.HGVS：逐一試候選轉錄本，用外顯子邊界 + 參考鹼基（+ p. 交叉驗證）挑最合理者。"""
        key = _norm_text(parsed.raw).split(":")[-1].upper()
        if key in self.hgmd_index:
            pos, ref, alt = self.hgmd_index[key]
            vclass = self._vclass_from_parsed(parsed)
            return self._finish(pos, ref, alt, parsed, parsed.transcript, vclass, None,
                                "hgmd_kb_exact", "high", [])

        explicit = self.model.resolve_transcript_id(parsed.transcript)
        candidates = [explicit] if explicit else [
            t for t in self.transcript_priority if self.model.has_transcript(t)
        ]

        scored: List[Tuple[int, Dict[str, Any]]] = []
        for tid in candidates:
            pos = self.model.c_to_genomic(
                tid, int(parsed.start or 0), parsed.start_offset, parsed.utr,
                strict_boundary=True,
            )
            if pos is None:
                continue
            score = 0
            warnings: List[str] = []
            if explicit:
                score += 4
            if parsed.start_offset:
                score += 3  # 通得過外顯子邊界檢查本身就是強證據

            # 參考鹼基驗證
            ref_ok = None
            if parsed.ref and len(parsed.ref) == 1:
                base = self.model.base_at(pos)
                ref_ok = base is not None and base == parsed.ref
                score += 4 if ref_ok else -6

            # 蛋白質後果（外顯子內置換）
            vclass, hgvs_p = "unknown", None
            if parsed.edit == "sub" and not parsed.start_offset and parsed.utr is None:
                cons = self.model.predict_substitution_consequence(
                    tid, int(parsed.start or 0), parsed.ref or "", parsed.alt or ""
                )
                vclass, hgvs_p = cons["vclass"], cons["hgvs_p"]
                # 與文獻給的 p. 交叉驗證（同一顆變異、不同轉錄本 aa 編號不同）
                if p_extra and cons.get("aa_pos"):
                    same_ref = not p_extra.aa_ref or p_extra.aa_ref == cons.get("aa_ref")
                    same_alt = not p_extra.aa_alt or p_extra.aa_alt == cons.get("aa_alt")
                    if p_extra.aa_pos == cons["aa_pos"] and same_ref:
                        score += 5 + (1 if same_alt else 0)
                    else:
                        score -= 2
            if vclass == "unknown":
                vclass = self._vclass_from_parsed(parsed)
            if parsed.ref and ref_ok is False:
                warnings.append(
                    f"{tid}: c.{parsed.start} 的參考鹼基不符"
                    f"（基因組 {self.model.base_at(pos)} vs HGVS {parsed.ref}）"
                )
            scored.append((score, {
                "tid": tid, "pos": pos, "vclass": vclass, "hgvs_p": hgvs_p,
                "warnings": warnings, "ref_ok": ref_ok,
            }))

        if not scored:
            raise HgvsResolveError(
                f"{parsed.raw}：所有候選轉錄本都無法定位"
                f"{'（intronic 位置不在任何外顯子邊界）' if parsed.start_offset else ''}"
            )

        scored.sort(key=lambda x: -x[0])
        best_score, best = scored[0]
        warnings = list(best["warnings"])
        if best["ref_ok"] is False and not explicit:
            raise HgvsResolveError(
                f"{parsed.raw}：找不到參考鹼基相符的轉錄本（試過 {len(scored)} 個）"
            )
        # 多個轉錄本同分 → 位置可能不唯一
        tied = [d["pos"] for s, d in scored if s == best_score]
        if len(set(tied)) > 1:
            warnings.append(
                f"有 {len(set(tied))} 個轉錄本同樣相符，採用 {best['tid']}（可能有座標歧義）"
            )
        confidence = "high" if (best["ref_ok"] or parsed.start_offset) else "medium"
        if len(set(tied)) > 1:
            confidence = "medium"
        return self._finish(
            best["pos"], parsed.ref, parsed.alt, parsed, best["tid"], best["vclass"],
            best["hgvs_p"], f"cdna_hgvs:{best['tid']}", confidence, warnings,
        )

    def _resolve_protein(self, parsed: ParsedHgvs) -> ResolvedVariant:
        """p.HGVS：用轉譯出的胺基酸驗證是哪個轉錄本編號，回傳 codon 起點座標。"""
        explicit = self.model.resolve_transcript_id(parsed.transcript)
        candidates = [explicit] if explicit else [
            t for t in self.transcript_priority if self.model.has_transcript(t)
        ]
        matches: List[Tuple[str, int]] = []
        for tid in candidates:
            aa = self.model.aa_at(tid, int(parsed.aa_pos or 0))
            if aa is None:
                continue
            if parsed.aa_ref and aa != parsed.aa_ref:
                continue
            pos = self.model.aa_to_genomic(tid, int(parsed.aa_pos or 0))
            if pos is not None:
                matches.append((tid, pos))
        if not matches:
            raise HgvsResolveError(
                f"{parsed.raw}：沒有轉錄本在 aa {parsed.aa_pos} 有 {parsed.aa_ref}"
            )
        tid, pos = matches[0]
        warnings: List[str] = []
        if len({p for _, p in matches}) > 1:
            warnings.append(
                f"有 {len(matches)} 個轉錄本的 aa {parsed.aa_pos} 都是 {parsed.aa_ref}，"
                f"採用 {tid}（可能有座標歧義）"
            )
        vclass = {
            "missense": "missense", "nonsense": "nonsense",
            "frameshift": "frameshift", "synonymous": "synonymous",
        }.get(parsed.p_kind or "", "unknown")
        return self._finish(
            pos, None, None, parsed, tid, vclass, _norm_text(parsed.raw).split(":")[-1],
            f"protein_hgvs:{tid}", "medium" if not warnings else "low", warnings,
        )

    # -- rsID（唯一需要連網的路徑）-------------------------------------------

    def _load_rsid_cache(self) -> Dict[str, Any]:
        if self._rsid_cache is None:
            try:
                self._rsid_cache = json.loads(self.cache_path.read_text(encoding="utf-8"))
            except Exception:  # noqa: BLE001
                self._rsid_cache = {}
        return self._rsid_cache

    def _save_rsid_cache(self) -> None:
        try:
            self.cache_path.parent.mkdir(parents=True, exist_ok=True)
            self.cache_path.write_text(
                json.dumps(self._load_rsid_cache(), ensure_ascii=False), encoding="utf-8"
            )
        except Exception as e:  # noqa: BLE001
            logger.debug(f"rsID 快取寫入失敗：{e}")

    def _resolve_rsid(self, rsid: str) -> ResolvedVariant:
        """用 dbSNP（NCBI Variation Services）取 GRCh38 座標。"""
        rsid = rsid.lower()
        cache = self._load_rsid_cache()
        hit = cache.get(rsid, "missing")
        if hit == "missing":
            if not self.allow_network:
                raise HgvsResolveError(f"{rsid}：rsID 需要連網查 dbSNP（allow_network=False）")
            hit = self._fetch_rsid(rsid)
            cache[rsid] = hit
            self._save_rsid_cache()
        if not hit:
            raise HgvsResolveError(f"{rsid}：dbSNP 查不到 GRCh38 chr2 座標")
        pos, plus_ref, plus_alt = int(hit[0]), str(hit[1]), str(hit[2])
        ref, alt = complement(plus_ref), complement(plus_alt)
        warnings: List[str] = []
        base = self.model.base_at(pos)
        if base and len(ref) == 1 and base != ref:
            warnings.append(f"參考鹼基不符：基因組為 {base}，dbSNP 為 {ref}")
        vclass, hgvs_p = self._consequence_from_genomic(pos, ref, alt)
        parsed = ParsedHgvs(kind="vcf", raw=rsid, pos=pos, ref=ref, alt=alt, edit="sub")
        return self._finish(
            pos, ref, alt, parsed, None, vclass, hgvs_p, f"dbsnp:{rsid}",
            "high" if base == ref else "low", warnings,
        )

    @staticmethod
    def _fetch_rsid(rsid: str) -> Optional[List[Any]]:
        """回傳第一個 alt allele 的 (pos, plus_ref, plus_alt)；查不到回 None。"""
        alleles = fetch_rsid_alleles(rsid)
        return list(alleles[0]) if alleles else None

    # -- 類型（無法轉譯時的規則） -------------------------------------------

    @staticmethod
    def _vclass_from_parsed(parsed: ParsedHgvs) -> str:
        if parsed.kind == "c":
            if parsed.start_offset and abs(parsed.start_offset) <= 2:
                return "splice"
            if parsed.start_offset:
                return "intron"
            if parsed.utr:
                return "utr"
            edit = parsed.edit or "unknown"
            if edit in ("del", "dup", "ins", "delins"):
                length = HgvsResolver._indel_length(parsed)
                if length is None:
                    return "frameshift"  # TTN 文獻的 del/dup 絕大多數是 frameshift
                return "inframe" if length % 3 == 0 else "frameshift"
            return "unknown"
        return "unknown"

    @staticmethod
    def _indel_length(parsed: ParsedHgvs) -> Optional[int]:
        if parsed.ref:
            return len(parsed.ref)
        if parsed.alt and parsed.edit in ("ins", "delins"):
            return len(parsed.alt)
        if parsed.start and parsed.end:
            return abs(int(parsed.end) - int(parsed.start)) + 1
        if parsed.start and not parsed.end and parsed.edit in ("del", "dup"):
            return 1
        return None


# ── HGMD 索引（可選的加速 / 校準來源）───────────────────────────────────────


def fetch_rsid_alleles(rsid: str) -> Optional[List[Tuple[int, str, str]]]:
    """dbSNP（NCBI Variation Services）在 GRCh38 chr2（NC_000002.12）上的所有 alt allele。

    回傳 [(pos, plus_ref, plus_alt), ...]，pos 為 1-based、鹼基為正股；一個 rsID 可能有多個 alt。
    查詢失敗（例如沒有網路）回 None，查得到但不在 GRCh38 chr2 上回 []。
    """
    import urllib.error
    import urllib.request

    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{str(rsid).lower().removeprefix('rs')}"
    try:
        with urllib.request.urlopen(url, timeout=60) as resp:
            data = json.loads(resp.read())
    except urllib.error.HTTPError as e:
        if e.code in (404, 410):  # dbSNP 沒有這個（或已撤銷的）rsID
            return []
        logger.debug(f"dbSNP 查詢失敗 {rsid}: {e}")
        return None
    except Exception as e:  # noqa: BLE001
        logger.debug(f"dbSNP 查詢失敗 {rsid}: {e}")
        return None
    out: List[Tuple[int, str, str]] = []
    snap = data.get("primary_snapshot_data") or {}
    for placement in snap.get("placements_with_allele", []):
        # 只用 GRCh38；NC_000002.11 是 GRCh37，座標不同
        if str(placement.get("seq_id") or "") != "NC_000002.12":
            continue
        for allele in placement.get("alleles", []):
            spdi = allele.get("allele", {}).get("spdi") or {}
            deleted, inserted = spdi.get("deleted_sequence"), spdi.get("inserted_sequence")
            if deleted and inserted and deleted != inserted:
                # SPDI position 為 0-based → +1 轉成 1-based genomic
                out.append((int(spdi["position"]) + 1, deleted, inserted))
    return out


def build_hgmd_hgvs_index(kb_rows: List[Dict[str, Any]]) -> Dict[str, Tuple[int, str, str]]:
    """由 MechanisticInference._load_kb() 的結果建立 {正規化 c.HGVS -> (pos, ref, alt)}。"""
    index: Dict[str, Tuple[int, str, str]] = {}
    for r in kb_rows:
        dna = str(r.get("dna") or "")
        if ":" in dna:
            dna = dna.split(":", 1)[1]
        key = _norm_text(dna).upper()
        if key.startswith("C.") and r.get("pos"):
            index.setdefault(key, (int(r["pos"]), str(r.get("ref") or ""), str(r.get("alt") or "")))
    return index


__all__ = [
    "HgvsResolveError",
    "HgvsResolver",
    "ParsedHgvs",
    "ResolvedVariant",
    "TTNTranscriptModel",
    "build_hgmd_hgvs_index",
    "parse_hgvs",
]
