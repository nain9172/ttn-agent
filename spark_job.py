#!/usr/bin/env python3
"""Spark 端的任務佇列：接收使用者從 site 送來的變異、跑完整 pipeline、打包成 site 補丁。

兩種用法共用同一個佇列：
  - docker 服務（web_server.py + 本腳本的 worker --forever）：網頁直接呼叫 submit_job()，
    worker 跑完後把結果直接併進 TTN_SITE_DIR 指定的網站（見 merge_patch_into_site）。
  - 舊版 Windows 啟動程式（launcher/）透過 SSH 呼叫本腳本；所有指令的輸出都是單一
    JSON 物件（stdout），方便對方解析：

    python spark_job.py submit            # 從 stdin 讀 {"variants": [...]}，排入佇列
    python spark_job.py status [JOB_ID…]  # 列出任務狀態（不指定就列出全部）
    python spark_job.py check             # 環境自我檢查（連線測試用）
    python spark_job.py worker            # 依序處理佇列（submit 會自動在背景啟動）
    python spark_job.py worker --forever  # 常駐 worker（docker），佇列空了也不結束

除了變異分析，佇列也接受「setup」任務（kind=setup）：在空的 volume 上下載資料、
從頭重算或重建網站，實際執行 setup_and_run.sh 的指定步驟。

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
from utils import variant_ids  # noqa: E402

JOBS_DIR = config.OUTPUT_DIR / "site_jobs"
# docker 服務設定 TTN_SITE_DIR：任務完成後直接把報告併進這個網站
SITE_DIR = Path(os.environ["TTN_SITE_DIR"]) if os.environ.get("TTN_SITE_DIR") else None
LOCK_PATH = JOBS_DIR / "worker.lock"
VARIANT_RE = re.compile(r"^(?:chr)?2-(\d+)-([ACGT])-([ACGT])$", re.IGNORECASE)
JOB_ID_RE = re.compile(r"^\d{8}-\d{6}-[0-9a-f]{6}$")
MAX_VARIANTS_PER_JOB = 200
# 開始前要空出多少記憶體：Gemma 4 由 vLLM 一次預留「總量 × LLM_GPU_MEMORY_UTILIZATION」，
# 不夠時 main.py 會在載入模型時被系統 OOM killer 直接殺掉。這台 Spark 也有別人在跑 GPU 工作，
# 所以先等記憶體空出來再開始，最多等 MEMORY_WAIT_HOURS 小時。
MEMORY_WAIT_HOURS = 12
MEMORY_POLL_SECONDS = 60
POLL_SECONDS = 5  # 常駐 worker 檢查新任務的間隔
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
    """回傳 (內部 variant_id, 錯誤原因)。目前只支援 TTN 上的單鹼基變異。

    輸入可用基因組正股（VCF / ClinVar 慣例）或 TTN 轉錄方向的鹼基，依參考序列自動判斷；
    回傳的是內部（轉錄方向）ID，因為 outputs/ 的報告檔名用它。對外顯示時再轉回正股（summarize）。
    """
    text = raw.strip()
    m = VARIANT_RE.match(text)
    if not m:
        return None, "format must be 2-<position>-<ref>-<alt> (single-nucleotide variant on chr2, GRCh38)"
    pos, ref, alt = int(m.group(1)), m.group(2).upper(), m.group(3).upper()
    if ref == alt:
        return None, "ref and alt are the same base"
    if not (config.TTN_SEQUENCE_END <= pos <= config.TTN_SEQUENCE_START):
        return None, (f"position is outside TTN "
                      f"(chr2:{config.TTN_SEQUENCE_END:,}-{config.TTN_SEQUENCE_START:,})")
    # sequence.fasta 是負股（轉錄方向）
    tx = seq[config.TTN_SEQUENCE_START - pos]
    if ref == variant_ids.complement(tx):  # 正股寫法（標準）
        return f"2-{pos}-{tx}-{variant_ids.complement(alt)}", None
    if ref == tx:  # 轉錄方向寫法
        return f"2-{pos}-{ref}-{alt}", None
    return None, (f"reference base at chr2:{pos:,} is {variant_ids.complement(tx)} "
                  f"(GRCh38 + strand), not {ref}")


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
    if job.get("kind") == "setup":
        return job
    job["completed"] = completed_variants(job_dir)
    # job.json 存內部（轉錄方向）ID；回給網頁的一律是正股寫法
    for key in ("variants", "completed", "failed_variants"):
        if job.get(key):
            job[key] = [variant_ids.to_display(v) for v in job[key]]
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


def submit_job(raw: List, rerun: bool = False, client: str = "") -> Dict:
    """驗證變異並排入佇列；回傳 {"ok", "job" | "error", "rejected"}。不負責啟動 worker。"""
    if not isinstance(raw, list) or not raw:
        return {"ok": False, "error": "no variants given"}
    if len(raw) > MAX_VARIANTS_PER_JOB:
        return {"ok": False, "error": f"at most {MAX_VARIANTS_PER_JOB} variants per job"}
    if not config.REFERENCE_GENOME_PATH.exists():
        return {"ok": False, "error": "the TTN reference sequence is not installed yet; finish the data setup first"}

    seq = load_sequence()
    accepted, rejected = [], []
    for item in raw:
        vid, err = validate(str(item), seq)
        if err:
            rejected.append({"variant": str(item), "reason": err})
        elif vid not in accepted:
            accepted.append(vid)
    if not accepted:
        return {"ok": False, "error": "no valid variants", "rejected": rejected}

    job_dir = new_job_dir()
    with (job_dir / "input.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["input"])
        writer.writerows([v] for v in accepted)
    write_job(job_dir, {
        "job_id": job_dir.name, "kind": "variants", "state": "queued", "variants": accepted, "rejected": rejected,
        # rerun：即使 outputs/ 已有報告也重跑整個 pipeline（使用者在網頁勾選重新分析）
        "rerun": bool(rerun),
        "client": str(client or ""), "created": now(),
    })
    return {"ok": True, "job": summarize(job_dir)}


# setup 任務可以執行的 setup_and_run.sh 步驟（env / config 在 docker 映像裡已經處理好）
SETUP_STAGES = ("models", "data", "pubtator", "evo2", "agent", "backfill", "site")


def submit_setup_job(stages: List[str], client: str = "", extra: Optional[List[str]] = None) -> Dict:
    """排入一個 setup 任務（下載資料、從頭重算、重建網站）。同時只允許一個未完成的 setup 任務。"""
    stages = [s for s in stages if s in SETUP_STAGES]
    if not stages:
        return {"ok": False, "error": "no valid setup steps given"}
    for job in list_jobs():
        if job.get("kind") == "setup" and job.get("state") in ("queued", "running", "waiting"):
            return {"ok": False, "error": f"setup job {job['job_id']} is still {job['state']}"}
    job_dir = new_job_dir()
    write_job(job_dir, {
        "job_id": job_dir.name, "kind": "setup", "state": "queued", "stages": stages,
        "extra_args": list(extra or []), "client": str(client or ""), "created": now(),
    })
    return {"ok": True, "job": summarize(job_dir)}


def new_job_dir() -> Path:
    job_id = time.strftime("%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:6]
    job_dir = JOBS_DIR / job_id
    job_dir.mkdir(parents=True)
    return job_dir


def list_jobs(job_ids: Optional[List[str]] = None) -> List[Dict]:
    jobs = []
    if JOBS_DIR.exists():
        ids = [j for j in (job_ids or []) if JOB_ID_RE.match(j)]
        dirs = [JOBS_DIR / j for j in ids] if ids else sorted(JOBS_DIR.iterdir())
        for d in dirs:
            if (d / "job.json").exists():
                try:
                    jobs.append(summarize(d))
                except (OSError, ValueError):
                    continue  # 正在寫入的 job.json
    return jobs


def cancel_job(job_id: str) -> Dict:
    """排隊中的直接取消；執行中的設旗標，由 worker 結束子程序。"""
    if not JOB_ID_RE.match(job_id or "") or not (JOBS_DIR / job_id / "job.json").exists():
        return {"ok": False, "error": "no such job"}
    job_dir = JOBS_DIR / job_id
    job = read_job(job_dir)
    if job.get("state") == "queued":
        job.update(state="cancelled", finished=now())
    elif job.get("state") in ("running", "waiting"):
        job["cancel_requested"] = True
    else:
        return {"ok": False, "error": f"job is already {job.get('state')}"}
    write_job(job_dir, job)
    return {"ok": True, "job": summarize(job_dir)}


def cancel_requested(job_dir: Path) -> bool:
    try:
        return bool(read_job(job_dir).get("cancel_requested"))
    except (OSError, ValueError):
        return False


def run_logged(job_dir: Path, cmd: List[str], env: Dict) -> int:
    """執行子程序（輸出寫入 run.log），期間每幾秒檢查取消旗標。被取消回傳 -1。"""
    import signal

    with (job_dir / "run.log").open("a", encoding="utf-8") as log:
        log.write(f"[{now()}] {' '.join(cmd)}\n")
        log.flush()
        proc = subprocess.Popen(cmd, cwd=str(ROOT), env=env, stdout=log, stderr=subprocess.STDOUT,
                                start_new_session=True)
        while True:
            try:
                return proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                pass
            write_heartbeat("busy", job_dir.name)
            if cancel_requested(job_dir):
                log.write(f"[{now()}] cancelled by user\n")
                for sig in (signal.SIGTERM, signal.SIGKILL):
                    try:
                        os.killpg(proc.pid, sig)
                    except ProcessLookupError:
                        break
                    try:
                        proc.wait(timeout=60)
                        break
                    except subprocess.TimeoutExpired:
                        continue
                return -1


def cmd_submit() -> int:
    try:
        payload = json.loads(sys.stdin.read() or "{}")
    except json.JSONDecodeError as e:
        emit({"ok": False, "error": f"invalid JSON on stdin: {e}"})
        return 1
    res = submit_job(payload.get("variants") or [], bool(payload.get("rerun")), str(payload.get("client") or ""))
    if res["ok"]:
        start_worker()
    emit(res)
    return 0 if res["ok"] else 1


def cmd_status(job_ids: List[str]) -> int:
    emit({"ok": True, "jobs": list_jobs(job_ids), "worker_running": worker_running()})
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
        if cancel_requested(job_dir):
            return False
        free, need = head
        who = gpu_occupants()
        job = read_job(job_dir)
        job.update(state="waiting",
                   note=(f"Waiting for memory on the Spark: {free:.0f} GiB free, {need:.0f} GiB needed"
                         + (f". In use by: {who}" if who else "")))
        write_job(job_dir, job)
        if time.time() >= deadline:
            return False
        write_heartbeat("waiting", job_dir.name)
        time.sleep(MEMORY_POLL_SECONDS)


def run_job(job_dir: Path) -> None:
    if read_job(job_dir).get("kind") == "setup":
        run_setup_job(job_dir)  # setup_and_run.sh 自己會在需要 GPU 的步驟前等記憶體
        return
    if not wait_for_memory(job_dir):
        job = read_job(job_dir)
        if job.get("cancel_requested"):
            job.update(state="cancelled", finished=now())
            write_job(job_dir, job)
            return
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

    env = pipeline_env()
    cmd = [
        sys.executable, str(ROOT / "batch_run_variants.py"),
        "--input", str(job_dir / "input.csv"),
        "--out", str(job_dir / "result.csv"),
    ]
    if rerun:
        cmd.append("--no-reuse")
    code = run_logged(job_dir, cmd, env)
    if code == -1:
        job = read_job(job_dir)
        job.update(state="cancelled", finished=now())
        write_job(job_dir, job)
        return

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

    merged = {}
    if reports:
        import build_site

        build_site.build_patch(reports, job_dir / "patch", getattr(config, "PUBMIND_JSON_PATH", None))
        if SITE_DIR is not None:
            # docker 服務：直接併進 volume 上的網站，所有使用者重新整理就看得到
            try:
                merged = {"merged": True, "reports": merge_patch_into_site(job_dir / "patch", SITE_DIR)}
            except Exception as e:  # noqa: BLE001
                merged = {"merge_error": f"could not add the reports to the site: {e}"}
    missing = [v for v in job["variants"] if v not in done]
    job.update(
        state="done" if reports else "failed",
        finished=now(),
        exit_code=code,
        failed_variants=missing,
        error="" if reports else failure_reason(job_dir, code),
        **merged,
    )
    write_job(job_dir, job)


def pipeline_env() -> Dict[str, str]:
    env = dict(os.environ)
    # 經 SSH 執行時 conda 沒有 activate；vLLM 需要環境 bin 裡的 ninja
    env["PATH"] = str(Path(sys.executable).parent) + os.pathsep + env.get("PATH", "")
    return env


def merge_patch_into_site(patch: Path, site: Path) -> List[str]:
    """把 build_site.build_patch 的補丁併進網站：報告、圖檔照路徑放，索引列寫進 data/local-reports.js。

    與舊版 Windows 啟動程式（launcher/app.go mergePatch）做的事相同。回傳加入的報告檔名。
    """
    import shutil

    for f in patch.rglob("*"):
        if f.is_file() and f.name != "records.json":
            dest = site / f.relative_to(patch)
            dest.parent.mkdir(parents=True, exist_ok=True)
            tmp = dest.with_name(dest.name + ".tmp")
            shutil.copyfile(f, tmp)
            os.replace(tmp, dest)
    records = json.loads((patch / "records.json").read_text(encoding="utf-8"))
    data_dir = site / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    local_json = data_dir / "local-reports.json"
    with (data_dir / "local-reports.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        existing = json.loads(local_json.read_text(encoding="utf-8")) if local_json.exists() else []
        for r in records:
            ids = {r.get("variant_id"), r.get("internal_id")} - {None}
            existing = [e for e in existing if e.get("variant_id") not in ids]
            existing.append(r)
        payload = json.dumps(existing, ensure_ascii=False, separators=(",", ":"))
        for path, text in ((local_json, payload),
                           (data_dir / "local-reports.js",
                            "window.__TTN_LOCAL_REPORTS__ = " + payload.replace("</", "<\\/") + ";\n")):
            tmp = path.with_name(path.name + ".tmp")
            tmp.write_text(text, encoding="utf-8")
            os.replace(tmp, path)
    return [r["file"] for r in records]


def run_setup_job(job_dir: Path) -> None:
    """執行 setup_and_run.sh 的指定步驟（下載資料、從頭重算、重建網站）。"""
    job = read_job(job_dir)
    job.update(state="running", started=now())
    write_job(job_dir, job)
    cmd = ["bash", str(ROOT / "setup_and_run.sh"), "--yes"]
    for stage in job.get("stages") or []:
        cmd += ["--only", stage]
    cmd += [str(a) for a in job.get("extra_args") or []]
    code = run_logged(job_dir, cmd, pipeline_env())
    job = read_job(job_dir)
    if code == -1:
        job.update(state="cancelled", finished=now())
    else:
        job.update(state="done" if code == 0 else "failed", finished=now(), exit_code=code,
                   error="" if code == 0 else f"setup_and_run.sh exited with code {code}; see the log")
    write_job(job_dir, job)


def write_heartbeat(state: str, job_id: str = "") -> None:
    """worker 每次輪詢都寫一次，網頁據此判斷 worker 是否還活著。"""
    path = JOBS_DIR / "worker_heartbeat.json"
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps({"time": time.time(), "state": state, "job_id": job_id, "pid": os.getpid()}))
    os.replace(tmp, path)


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


def cmd_worker(forever: bool = False) -> int:
    """依序處理佇列。forever（docker 常駐 worker）：佇列空了繼續等，並定期寫 heartbeat。"""
    JOBS_DIR.mkdir(parents=True, exist_ok=True)
    with LOCK_PATH.open("a") as lock:
        while True:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                if not forever:
                    return 0  # 已經有 worker 在跑
                time.sleep(30)
        # 上次被中斷而停在 running 的任務重新排隊（使用者已按取消的就直接結束）
        for d in JOBS_DIR.iterdir():
            if (d / "job.json").exists() and read_job(d).get("state") in ("running", "waiting"):
                job = read_job(d)
                job["state"] = "cancelled" if job.get("cancel_requested") else "queued"
                write_job(d, job)
        while True:
            queued = sorted(d for d in JOBS_DIR.iterdir()
                            if (d / "job.json").exists() and read_job(d).get("state") == "queued")
            if not queued:
                if not forever:
                    return 0
                write_heartbeat("idle")
                time.sleep(POLL_SECONDS)
                continue
            write_heartbeat("busy", queued[0].name)
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
        return cmd_worker(forever="--forever" in rest)
    emit({"ok": False, "error": f"unknown command {command!r}"})
    return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
