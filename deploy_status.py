"""docker 服務的 volume 檢查：列出 pipeline 需要的資料、模型與結果目前在不在。

web_server.py 的設定頁（/setup）與 /api/status 用它；也可以直接執行看文字版：

    python deploy_status.py

每一項標明來源：
  auto     可以由網頁上的「下載公開資料」或「從頭重算」自動取得（對應 setup_and_run.sh 的步驟）
  manual   需要授權或沒有公開下載點，必須手動放進 volume（見 docker_readme/DEPLOY_DOCKER.md）
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Dict, List

import config

ROOT = Path(__file__).resolve().parent
DATA = config.DATA_DIR
OUT = config.OUTPUT_DIR
SITE = Path(os.environ.get("TTN_SITE_DIR") or ROOT / "site")
HF_HOME = Path(os.environ.get("HF_HOME") or Path.home() / ".cache" / "huggingface")


def _hf_model(repo: str) -> Path:
    return HF_HOME / "hub" / ("models--" + repo.replace("/", "--"))


def _model_ok(repo: str) -> bool:
    snap = _hf_model(repo) / "snapshots"
    return snap.is_dir() and any(f.suffix == ".safetensors" or f.suffix == ".pt"
                                 for f in snap.rglob("*") if f.is_file() or f.is_symlink())


# (key, 名稱, 路徑, 來源, 必要程度, 說明)
ITEMS = [
    ("sequence", "TTN reference sequence", DATA / "sequence.fasta", "auto:data", "required",
     "GRCh38 chr2:178,525,989-178,807,423 (minus strand) from NCBI."),
    ("transcript_model", "TTN transcript model", DATA / "ttn_transcript_model.json", "auto:data", "required",
     "RefSeq exon model for HGVS ↔ genomic conversion, built from the NCBI annotation."),
    ("transcript_interval", "Isoform intervals", DATA / "transcript_interval.xlsx", "auto:image", "required",
     "Ships with the image and is copied into the volume on first start."),
    ("clinvar", "ClinVar TTN variant list", DATA / "clinvar_ttn.csv", "auto:data", "recommended",
     "Input list for the full recompute (about 5,000 variants)."),
    ("pubmind", "PubMind TTN records", DATA / "pubmind_ttn.json", "auto:data", "recommended",
     "Public version downloads automatically; the full licensed version replaces it if you place it here."),
    ("pubtator", "PubTator TTN variant corpus", DATA / "pubtator" / "pubtator_ttn_variants_clinvar_annotated.csv",
     "auto:pubtator", "recommended",
     "Text-mined TTN variants from PubMed; the first harvest takes many hours and needs the LLM."),
    ("hgmd_csv", "HGMD export (licensed)", DATA / "HGMD_latest.csv", "manual", "recommended",
     "HGMD Professional export for TTN. Needed for database lookup and analogy inference."),
    ("hgmd_xlsx", "HGMD processed table (licensed)", DATA / "HGMD_TTN_processed.xlsx", "manual", "recommended",
     "Provides REF/ALT so HGMD can be matched allele-exactly."),
    ("cardiodb", "cardiodb exon table (licensed)", DATA / "cardiodb.xlsx", "manual", "recommended",
     "Exon coordinates, bands and measured heart / skeletal PSI used by the inference."),
    ("exon_psi", "Exon PSI table", DATA / "ttn_exon_psi.csv", "manual", "optional",
     "Measured PSI for 323 exons; without it cardiodb PSI is used."),
]

MODELS = [
    ("evo2_model", "Evo2-7B weights", "arcinstitute/" + config.EVO2_MODEL, "required",
     "About 14 GB from Hugging Face."),
    ("llm_model", "Gemma 4 weights", config.LOCAL_LLM_MODEL, "required",
     "About 31 GB from Hugging Face; needs HF_TOKEN if the model is gated."),
]


def status() -> Dict:
    items: List[Dict] = []
    for key, name, path, source, level, note in ITEMS:
        ok = path.exists() and path.stat().st_size > 0
        items.append({"key": key, "name": name, "path": _rel(path), "source": source, "level": level,
                      "note": note, "ok": ok})
    for key, name, repo, level, note in MODELS:
        items.append({"key": key, "name": name, "path": _rel(_hf_model(repo)), "source": "auto:models",
                      "level": level, "note": note, "ok": _model_ok(repo)})

    reports = sum(1 for _ in OUT.glob("2-*.html")) if OUT.exists() else 0
    site_ok = (SITE / "index.html").exists()
    evo2_table = (OUT / "evo2_table" / "evo2_table.js").exists()
    missing_required = [i["name"] for i in items if i["level"] == "required" and not i["ok"]]
    return {
        "items": items,
        "reports": reports,
        "site_built": site_ok,
        "evo2_table": evo2_table,
        "missing_required": missing_required,
        "ready_for_analysis": not missing_required,
        "hf_token": bool(os.environ.get("HF_TOKEN") or os.environ.get("HUGGING_FACE_HUB_TOKEN")),
        "ncbi_email": bool(getattr(config, "PUBMED_EMAIL", "")),
    }


def _rel(path: Path) -> str:
    vol = os.environ.get("TTN_VOLUME")
    real = str(path.resolve())  # 不存在的檔案也會解開 /app/data → /vol/data 的連結
    if vol and real.startswith(vol.rstrip("/") + "/"):
        # 容器內 /vol/... → 主機上的資料夾路徑（docker-compose 傳入 TTN_HOST_DATA_DIR）
        host = (os.environ.get("TTN_HOST_DATA_DIR") or "").rstrip("/")
        return (host or "<volume>") + "/" + real[len(vol.rstrip("/")) + 1:]
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


if __name__ == "__main__":
    st = status()
    for i in st["items"]:
        print(f"{'OK ' if i['ok'] else '-- '} {i['level']:<11} {i['source']:<13} {i['name']}  ({i['path']})")
    print(json.dumps({k: v for k, v in st.items() if k != "items"}, indent=1, ensure_ascii=False))
