#!/usr/bin/env python3
"""Spark 端的任務佇列：接收使用者從 site 送來的變異、跑完整 pipeline、打包成 site 補丁。

使用者電腦上的 TTN-Site 啟動程式（launcher/）透過 SSH 呼叫本腳本，所有指令的輸出
都是單一 JSON 物件（stdout），方便對方解析：

    python spark_job.py submit            # 從 stdin 讀 {"variants": [...]}，排入佇列
    python spark_job.py status [JOB_ID…]  # 列出任務狀態（不指定就列出全部）
    python spark_job.py check             # 環境自我檢查（連線測試用）
    python spark_job.py worker            # 依序處理佇列（submit 會自動在背景啟動）

每個任務在 outputs/site_jobs/<job_id>/ 下：
    job.json       狀態（queued / running / done / failed）
    input.csv      變異清單（batch_run_variants.py 的輸入）
    result.csv     batch_run_variants.py 的輸出；report_file 有值代表該顆完成
    run.log        pipeline 的完整輸出
    patch/         完成後的 site 補丁：reports/、assets/、records.json（見 build_site.build_patch）

GPU 一次只能跑一個 pipeline，所以 worker 用檔案鎖保證全機只有一個在執行。
pipeline 走 batch_run_variants.py 的子程序模式（每顆變異一個 main.py，含 Evo2 與 LLM），
已在 outputs/ 有報告的變異會直接沿用。
"""

from __future__ import annotations

import csv
import fcntl
import json
import os
import re
import subprocess
import sys
import time
import uuid
from pathlib import Path
from typing import Dict, List, Optional

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import config  # noqa: E402

JOBS_DIR = config.OUTPUT_DIR / "site_jobs"
LOCK_PATH = JOBS_DIR / "worker.lock"
VARIANT_RE = re.compile(r"^(?:chr)?2-(\d+)-([ACGT])-([ACGT])$", re.IGNORECASE)
JOB_ID_RE = re.compile(r"^\d{8}-\d{6}-[0-9a-f]{6}$")
MAX_VARIANTS_PER_JOB = 200
# 開始前要空出多少記憶體：Gemma 4 由 vLLM 一次預留「總量 × LLM_GPU_MEMORY_UTILIZATION」，
# 不夠時 main.py 會在載入模型時被系統 OOM killer 直接殺掉。這台 Spark 也有別人在跑 GPU 工作，
# 所以先等記憶體空出來再開始，最多等 MEMORY_WAIT_HOURS 小時。
MEMORY_WAIT_HOURS = 12
MEMORY_POLL_SECONDS = 60
_PROBE = (
    "import json, torch; f, t = torch.cuda.mem_get_info(0); "
    "print(json.dumps([f / 2**30, t / 2**30]))"
)


def emit(obj: Dict) -> None:
    print(json.dumps(obj, ensure_ascii=False))


def now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S")


def load_sequence() -> str:
    lines = config.REFERENCE_GENOME_PATH.read_text().splitlines()
    return "".join(l.strip() for l in lines if not l.startswith(">")).upper()


def validate(raw: str, seq: str) -> tuple[Optional[str], Optional[str]]:
    """回傳 (正規化後的 variant_id, 錯誤原因)。目前只支援 TTN 上的單鹼基變異。"""
    text = raw.strip()
    m = VARIANT_RE.match(text)
    if not m:
        return None, "format must be 2-<position>-<ref>-<alt> (single-nucleotide variant on chr2)"
    pos, ref, alt = int(m.group(1)), m.group(2).upper(), m.group(3).upper()
    if ref == alt:
        return None, "ref and alt are the same base"
    if not (config.TTN_SEQUENCE_END <= pos <= config.TTN_SEQUENCE_START):
        return None, (f"position is outside TTN "
                      f"(chr2:{config.TTN_SEQUENCE_END:,}-{config.TTN_SEQUENCE_START:,})")
    # sequence.fasta 是負股（轉錄方向），與報告 variant ID 的 ref/alt 同一個方向
    actual = seq[config.TTN_SEQUENCE_START - pos]
    if actual != ref:
        return None, f"reference base at {pos:,} is {actual}, not {ref} (alleles are on the TTN coding strand)"
    return f"2-{pos}-{ref}-{alt}", None


def read_job(job_dir: Path) -> Dict:
    return json.loads((job_dir / "job.json").read_text(encoding="utf-8"))


def write_job(job_dir: Path, job: Dict) -> None:
    tmp = job_dir / "job.json.tmp"
    tmp.write_text(json.dumps(job, ensure_ascii=False, indent=2), encoding="utf-8")
    tmp.replace(job_dir / "job.json")


def completed_variants(job_dir: Path) -> List[str]:
    result = job_dir / "result.csv"
    if not result.exists():
        return []
    try:
        with result.open(newline="", encoding="utf-8") as fh:
            return [row["input"] for row in csv.DictReader(fh) if (row.get("report_file") or "").strip()]
    except (OSError, KeyError, csv.Error):
        return []


def summarize(job_dir: Path) -> Dict:
    job = read_job(job_dir)
    job["completed"] = completed_variants(job_dir)
    patch = job_dir / "patch"
    if job.get("state") == "done" and patch.exists():
        job["patch_dir"] = str(patch)
        job["patch_files"] = sorted(
            str(p.relative_to(patch)).replace(os.sep, "/") for p in patch.rglob("*") if p.is_file()
        )
    return job


def worker_running() -> bool:
    JOBS_DIR.mkdir(parents=True, exist_ok=True)
    with LOCK_PATH.open("a") as fh:
        try:
            fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return True
        fcntl.flock(fh, fcntl.LOCK_UN)
    return False


def start_worker() -> None:
    if worker_running():
        return
    log = (JOBS_DIR / "worker.log").open("a")
    subprocess.Popen(
        [sys.executable, str(Path(__file__).resolve()), "worker"],
        cwd=str(ROOT), stdin=subprocess.DEVNULL, stdout=log, stderr=log,
        start_new_session=True,  # SSH 斷線後繼續跑
    )


# ── 指令 ─────────────────────────────────────────────────────────────────────


def cmd_check() -> int:
    issues = []
    if not config.REFERENCE_GENOME_PATH.exists():
        issues.append(f"missing {config.REFERENCE_GENOME_PATH}")
    import importlib.util

    for module in ("pandas", "torch", "evo2", "vllm"):
        if importlib.util.find_spec(module) is None:
            issues.append(f"Python package {module} is not installed in {sys.executable}")
    emit({
        "ok": not issues,
        "issues": issues,
        "python": sys.executable,
        "repo": str(ROOT),
        "worker_running": worker_running(),
        "queued": sum(1 for d in JOBS_DIR.glob("*/job.json")
                      if json.loads(d.read_text()).get("state") in ("queued", "running"))
        if JOBS_DIR.exists() else 0,
    })
    return 0


def cmd_submit() -> int:
    try:
        payload = json.loads(sys.stdin.read() or "{}")
    except json.JSONDecodeError as e:
        emit({"ok": False, "error": f"invalid JSON on stdin: {e}"})
        return 1
    raw = payload.get("variants") or []
    if not isinstance(raw, list) or not raw:
        emit({"ok": False, "error": "no variants given"})
        return 1
    if len(raw) > MAX_VARIANTS_PER_JOB:
        emit({"ok": False, "error": f"at most {MAX_VARIANTS_PER_JOB} variants per job"})
        return 1

    seq = load_sequence()
    accepted, rejected = [], []
    for item in raw:
        vid, err = validate(str(item), seq)
        if err:
            rejected.append({"variant": str(item), "reason": err})
        elif vid not in accepted:
            accepted.append(vid)
    if not accepted:
        emit({"ok": False, "error": "no valid variants", "rejected": rejected})
        return 1

    job_id = time.strftime("%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:6]
    job_dir = JOBS_DIR / job_id
    job_dir.mkdir(parents=True)
    with (job_dir / "input.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["input"])
        writer.writerows([v] for v in accepted)
    write_job(job_dir, {
        "job_id": job_id, "state": "queued", "variants": accepted, "rejected": rejected,
        # rerun：即使 outputs/ 已有報告也重跑整個 pipeline（使用者在網頁勾選重新分析）
        "rerun": bool(payload.get("rerun")),
        "client": str(payload.get("client") or ""), "created": now(),
    })
    start_worker()
    emit({"ok": True, "job": summarize(job_dir)})
    return 0


def cmd_status(job_ids: List[str]) -> int:
    jobs = []
    if JOBS_DIR.exists():
        job_ids = [j for j in job_ids if JOB_ID_RE.match(j)]
        dirs = [JOBS_DIR / j for j in job_ids] if job_ids else sorted(JOBS_DIR.iterdir())
        for d in dirs:
            if (d / "job.json").exists():
                jobs.append(summarize(d))
    emit({"ok": True, "jobs": jobs, "worker_running": worker_running()})
    return 0


def gpu_occupants() -> str:
    """nvidia-smi 列出的 GPU 行程與擁有者，例：lilyhsiao VLLM::EngineCore 100.7 GiB。"""
    try:
        out = subprocess.run(["nvidia-smi", "--query-compute-apps=pid,used_memory,process_name",
                              "--format=csv,noheader,nounits"],
                             capture_output=True, text=True, timeout=15).stdout
    except (OSError, subprocess.SubprocessError):
        return ""
    rows = []
    for line in out.strip().splitlines():
        parts = [x.strip() for x in line.split(",")]
        if len(parts) < 3:
            continue
        user = subprocess.run(["ps", "-o", "user=", "-p", parts[0]],
                              capture_output=True, text=True).stdout.strip() or "?"
        try:
            mem = f"{float(parts[1]) / 1024:.1f} GiB"
        except ValueError:
            mem = parts[1]
        rows.append(f"{user} {parts[2].split('/')[-1]} {mem}")
    return "; ".join(rows)


def memory_headroom() -> Optional[tuple]:
    """(可用 GiB, 需要 GiB)；先釋放本機可釋放的檔案快取。查不到就回 None（不擋）。"""
    subprocess.run([sys.executable, str(ROOT / "scripts" / "release_page_cache.py")],
                   capture_output=True, timeout=600)
    try:
        out = subprocess.run([sys.executable, "-c", _PROBE], capture_output=True, text=True,
                             timeout=120).stdout.strip().splitlines()[-1]
        free, total = json.loads(out)
    except Exception:  # noqa: BLE001
        return None
    return free, total * float(getattr(config, "LLM_GPU_MEMORY_UTILIZATION", 0.85))


def wait_for_memory(job_dir: Path) -> bool:
    """記憶體不夠就把任務標成 waiting 並定期重查；夠了回 True，等太久回 False。"""
    deadline = time.time() + MEMORY_WAIT_HOURS * 3600
    while True:
        head = memory_headroom()
        if head is None or head[0] >= head[1]:
            return True
        free, need = head
        who = gpu_occupants()
        job = read_job(job_dir)
        job.update(state="waiting",
                   note=(f"Waiting for memory on the Spark: {free:.0f} GiB free, {need:.0f} GiB needed"
                         + (f". In use by: {who}" if who else "")))
        write_job(job_dir, job)
        if time.time() >= deadline:
            return False
        time.sleep(MEMORY_POLL_SECONDS)


def run_job(job_dir: Path) -> None:
    if not wait_for_memory(job_dir):
        job = read_job(job_dir)
        job.update(state="failed", finished=now(),
                   error=f"gave up after waiting {MEMORY_WAIT_HOURS} h for memory. " + job.get("note", ""))
        write_job(job_dir, job)
        return

    job = read_job(job_dir)
    job.pop("note", None)
    job.update(state="running", started=now())
    write_job(job_dir, job)
    rerun = bool(job.get("rerun"))
    run_started = time.time()

    env = dict(os.environ)
    # 經 SSH 執行時 conda 沒有 activate；vLLM 需要環境 bin 裡的 ninja
    env["PATH"] = str(Path(sys.executable).parent) + os.pathsep + env.get("PATH", "")
    cmd = [
        sys.executable, str(ROOT / "batch_run_variants.py"),
        "--input", str(job_dir / "input.csv"),
        "--out", str(job_dir / "result.csv"),
    ]
    if rerun:
        cmd.append("--no-reuse")
    with (job_dir / "run.log").open("a", encoding="utf-8") as log:
        log.write(f"[{now()}] {' '.join(cmd)}\n")
        log.flush()
        code = subprocess.call(cmd, cwd=str(ROOT), env=env, stdout=log, stderr=subprocess.STDOUT)

    job = read_job(job_dir)
    done = completed_variants(job_dir)
    reports = []
    if rerun:
        # 重跑時只收這次新產生的報告。batch_run_variants 會依模型 tag 偏好挑報告，
        # 若這次失敗或產生不同 tag，可能挑到舊報告，被誤當成新結果送回去。
        done = []
        for vid in job["variants"]:
            fresh = [f for f in config.OUTPUT_DIR.glob(f"{vid}_*.html")
                     if f.stat().st_mtime >= run_started]
            if fresh:
                reports.append(max(fresh, key=lambda f: f.stat().st_mtime))
                done.append(vid)
    else:
        result = job_dir / "result.csv"
        if result.exists():
            with result.open(newline="", encoding="utf-8") as fh:
                for row in csv.DictReader(fh):
                    name = (row.get("report_file") or "").strip()
                    if name and (config.OUTPUT_DIR / name).exists():
                        reports.append(config.OUTPUT_DIR / name)

    if reports:
        import build_site

        build_site.build_patch(reports, job_dir / "patch", getattr(config, "PUBMIND_JSON_PATH", None))
    missing = [v for v in job["variants"] if v not in done]
    job.update(
        state="done" if reports else "failed",
        finished=now(),
        exit_code=code,
        failed_variants=missing,
        error="" if reports else failure_reason(job_dir, code),
    )
    write_job(job_dir, job)


def failure_reason(job_dir: Path, code: int) -> str:
    """把 run.log 裡最常見的失敗原因翻成使用者看得懂的一句話。"""
    log = (job_dir / "run.log").read_text(encoding="utf-8", errors="replace") if (job_dir / "run.log").exists() else ""
    if "SIGKILL" in log:
        who = gpu_occupants()
        return ("the Spark ran out of memory and the system stopped the pipeline"
                + (f" (GPU in use by: {who})" if who else "") + "; try again when the Spark is free")
    if "insufficient free GPU memory" in log or "Not enough free GPU memory" in log:
        return "not enough free GPU memory on the Spark to load the language model; try again later"
    return f"pipeline produced no report (exit code {code}); see run.log on the Spark"


def cmd_worker() -> int:
    JOBS_DIR.mkdir(parents=True, exist_ok=True)
    with LOCK_PATH.open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return 0  # 已經有 worker 在跑
        # 上次被中斷而停在 running 的任務重新排隊
        for d in JOBS_DIR.iterdir():
            if (d / "job.json").exists() and read_job(d).get("state") in ("running", "waiting"):
                job = read_job(d)
                job["state"] = "queued"
                write_job(d, job)
        while True:
            queued = sorted(d for d in JOBS_DIR.iterdir()
                            if (d / "job.json").exists() and read_job(d).get("state") == "queued")
            if not queued:
                return 0
            try:
                run_job(queued[0])
            except Exception as e:  # noqa: BLE001
                job = read_job(queued[0])
                job.update(state="failed", finished=now(), error=f"{type(e).__name__}: {e}")
                write_job(queued[0], job)


def main(argv: List[str]) -> int:
    if not argv:
        print(__doc__)
        return 2
    command, rest = argv[0], argv[1:]
    if command == "check":
        return cmd_check()
    if command == "submit":
        return cmd_submit()
    if command == "status":
        return cmd_status(rest)
    if command == "worker":
        return cmd_worker()
    emit({"ok": False, "error": f"unknown command {command!r}"})
    return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
