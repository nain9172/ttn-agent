#!/usr/bin/env python3
"""TTN 報告網站的 HTTP 服務（docker 的 web 容器執行它）。

  - 提供 site/ 的靜態網頁（報告、索引、Information 頁）
  - /api/*：送出變異分析、查任務、取消任務、volume 狀態與 setup 任務
  - /setup：volume 設定頁（缺哪些資料、下載公開資料、從頭重算、重建網站）

實際分析由 worker 容器（python spark_job.py worker --forever）執行；兩者透過 volume 上的
outputs/site_jobs/ 溝通，所以 web 不需要 GPU，GPU 被佔用或 pipeline 當掉時網站照樣能開。

只用標準函式庫。存取控制：只接受 TTN_ALLOWED_NETWORKS 內的來源 IP（預設私有網段與本機）。

環境變數：
  TTN_SITE_DIR          網站目錄（預設 ./site）
  TTN_PORT / TTN_BIND   監聽埠與位址（預設 8080 / 0.0.0.0）
  TTN_ALLOWED_NETWORKS  允許的來源網段，逗號分隔，例 "10.0.0.0/8,140.112.0.0/16"；"*" 代表不限制
"""

from __future__ import annotations

import ipaddress
import json
import os
import sys
import time
from http import HTTPStatus
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Dict, List, Optional
from urllib.parse import parse_qs, urlparse

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
os.environ.setdefault("TTN_SITE_DIR", str(ROOT / "site"))

import deploy_status  # noqa: E402
import spark_job  # noqa: E402

SITE = Path(os.environ["TTN_SITE_DIR"])
DEFAULT_NETWORKS = "127.0.0.0/8,::1/128,10.0.0.0/8,172.16.0.0/12,192.168.0.0/16,fc00::/7"
MAX_BODY = 1 << 20
HEARTBEAT_STALE_SECONDS = 180

# 設定頁按鈕 → setup_and_run.sh 的步驟
SETUP_ACTIONS = {
    "download": (["models", "data"], [], "Download public data and models"),
    "trial": (["models", "data", "evo2", "agent", "backfill", "site"], ["--limit", "5"],
              "Trial run on 5 ClinVar variants"),
    "recompute": (["models", "data", "pubtator", "evo2", "agent", "backfill", "site"], [],
                  "Recompute everything from scratch"),
    "site": (["site"], [], "Rebuild the website from existing results"),
}


def parse_networks(text: str) -> Optional[List]:
    if text.strip() == "*":
        return None
    nets = []
    for part in text.split(","):
        part = part.strip()
        if part:
            nets.append(ipaddress.ip_network(part, strict=False))
    return nets


NETWORKS = parse_networks(os.environ.get("TTN_ALLOWED_NETWORKS") or DEFAULT_NETWORKS)


def client_allowed(addr: str) -> bool:
    if NETWORKS is None:
        return True
    try:
        ip = ipaddress.ip_address(addr)
    except ValueError:
        return False
    if ip.version == 6 and ip.ipv4_mapped:
        ip = ip.ipv4_mapped
    return any(ip in n for n in NETWORKS)


def worker_info() -> Dict:
    hb = spark_job.JOBS_DIR / "worker_heartbeat.json"
    try:
        data = json.loads(hb.read_text())
    except (OSError, ValueError):
        return {"alive": False, "state": "never seen"}
    age = time.time() - float(data.get("time", 0))
    return {"alive": age < HEARTBEAT_STALE_SECONDS, "state": data.get("state"), "job_id": data.get("job_id"),
            "seconds_since_heartbeat": round(age)}


def job_log_tail(job_id: str, lines: int = 200) -> str:
    if not spark_job.JOB_ID_RE.match(job_id or ""):
        return ""
    log = spark_job.JOBS_DIR / job_id / "run.log"
    if not log.exists():
        return ""
    with log.open("rb") as fh:
        fh.seek(0, os.SEEK_END)
        size = fh.tell()
        fh.seek(max(0, size - 256 * 1024))
        text = fh.read().decode("utf-8", errors="replace")
    # tqdm 進度列用 \r 覆寫同一行，只留最後一段
    rows = [r.split("\r")[-1] for r in text.splitlines()]
    return "\n".join(rows[-lines:])


class Handler(SimpleHTTPRequestHandler):
    server_version = "TTNSite/1.0"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(SITE), **kwargs)

    # ── 共用 ──
    def _deny(self) -> bool:
        if client_allowed(self.client_address[0]):
            return False
        self.send_error(HTTPStatus.FORBIDDEN, "This service only accepts connections from the internal network")
        return True

    def _json(self, obj: Dict, status: int = 200) -> None:
        body = json.dumps(obj, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Cache-Control", "no-store")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _html(self, text: str) -> None:
        body = text.encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Cache-Control", "no-store")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def end_headers(self) -> None:
        # 索引與本機報告清單會隨分析結果更新，不要被瀏覽器快取
        path = urlparse(self.path).path
        if path in ("/", "/index.html") or path.startswith("/data/"):
            self.send_header("Cache-Control", "no-cache")
        self.send_header("X-Content-Type-Options", "nosniff")
        super().end_headers()

    def log_message(self, fmt: str, *args) -> None:  # noqa: D401
        sys.stderr.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} {self.client_address[0]} {fmt % args}\n")

    # ── GET ──
    def do_GET(self) -> None:  # noqa: N802
        if self._deny():
            return
        url = urlparse(self.path)
        path = url.path
        if path.startswith("/api/"):
            return self._api_get(path, parse_qs(url.query))
        if path in ("/setup", "/setup/"):
            return self._html(SETUP_HTML)
        if path in ("/", "/index.html") and not (SITE / "index.html").exists():
            self.send_response(HTTPStatus.FOUND)
            self.send_header("Location", "/setup")
            self.end_headers()
            return
        if any(part.startswith(".") for part in path.split("/")):
            self.send_error(HTTPStatus.NOT_FOUND)
            return
        super().do_GET()

    def do_HEAD(self) -> None:  # noqa: N802
        if self._deny():
            return
        super().do_HEAD()

    def _api_get(self, path: str, query: Dict) -> None:
        if path == "/api/ping":
            return self._json({"ok": True, "mode": "server"})
        if path == "/api/status":
            st = deploy_status.status()
            return self._json({"ok": True, "setup": st, "worker": worker_info(),
                               "setup_jobs": [j for j in spark_job.list_jobs() if j.get("kind") == "setup"][-5:]})
        if path == "/api/jobs":
            jobs = [dict(j, host=j.get("client", "")) for j in spark_job.list_jobs() if j.get("kind") != "setup"]
            return self._json({"ok": True, "jobs": jobs[-50:], "worker": worker_info()})
        if path == "/api/local-reports":
            local = SITE / "data" / "local-reports.json"
            try:
                reports = json.loads(local.read_text(encoding="utf-8")) if local.exists() else []
            except ValueError:
                reports = []
            return self._json({"ok": True, "reports": reports})
        if path == "/api/job-log":
            job_id = (query.get("job") or [""])[0]
            return self._json({"ok": True, "log": job_log_tail(job_id)})
        self._json({"ok": False, "error": "not found"}, 404)

    # ── POST ──
    def do_POST(self) -> None:  # noqa: N802
        if self._deny():
            return
        # 自訂標頭擋掉跨站表單送出（瀏覽器跨來源加自訂標頭要先 preflight，而這裡不回 CORS）
        if self.headers.get("X-TTN-Site") != "1":
            return self._json({"ok": False, "error": "missing X-TTN-Site header"}, 403)
        try:
            length = int(self.headers.get("Content-Length") or 0)
        except ValueError:
            length = 0
        if length > MAX_BODY:
            return self._json({"ok": False, "error": "request too large"}, 413)
        try:
            body = json.loads(self.rfile.read(length) or b"{}")
        except ValueError:
            return self._json({"ok": False, "error": "invalid JSON"}, 400)
        path = urlparse(self.path).path
        client = self.client_address[0]
        if path == "/api/submit":
            return self._json(spark_job.submit_job(body.get("variants") or [], bool(body.get("rerun")), client))
        if path == "/api/cancel":
            return self._json(spark_job.cancel_job(str(body.get("job_id") or "")))
        if path == "/api/setup":
            action = SETUP_ACTIONS.get(str(body.get("action") or ""))
            if not action:
                return self._json({"ok": False, "error": "unknown action"}, 400)
            stages, extra, _ = action
            return self._json(spark_job.submit_setup_job(stages, client, extra))
        self._json({"ok": False, "error": "not found"}, 404)


SETUP_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>TTN Site Setup</title>
<style>
* { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: 'Segoe UI', 'Noto Sans', 'Noto Sans TC', 'Microsoft JhengHei', -apple-system, Roboto, Arial, sans-serif;
       color: #24303f; background: #f4f6fa; line-height: 1.55; }
header { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: #fff; padding: 28px 24px 24px; }
header .inner, main { max-width: 1100px; margin: 0 auto; }
header h1 { font-size: 1.7em; }
header p { opacity: 0.9; margin-top: 4px; }
header a { color: #fff; }
main { padding: 22px 24px 60px; }
section { background: #fff; border-radius: 12px; padding: 20px 22px; margin-bottom: 18px; box-shadow: 0 6px 20px rgba(31,45,68,0.07); }
h2 { font-size: 1.2em; margin-bottom: 12px; color: #2c3e50; }
table { width: 100%; border-collapse: collapse; font-size: 0.92em; }
th { text-align: left; color: #6b7a90; font-size: 0.8em; text-transform: uppercase; letter-spacing: 0.04em; padding: 6px 8px; border-bottom: 1px solid #e4e9f0; }
td { padding: 8px; border-bottom: 1px solid #f1f4f8; vertical-align: top; }
td.path { font-family: Consolas, Menlo, monospace; font-size: 0.85em; color: #55637a; word-break: break-all; }
.ok { color: #1e8449; font-weight: 700; }
.miss { color: #c0392b; font-weight: 700; }
.miss.optional, .miss.recommended { color: #b26a00; }
.pill { display: inline-block; padding: 1px 8px; border-radius: 999px; font-size: 0.78em; font-weight: 700; background: #eef1f6; color: #55637a; }
.pill.manual { background: #fff3cd; color: #8a5a00; }
.note { color: #6b7a90; font-size: 0.88em; }
.actions { display: flex; flex-wrap: wrap; gap: 10px; margin-top: 6px; }
button { border: 0; border-radius: 8px; padding: 9px 16px; font-size: 0.95em; cursor: pointer; background: #667eea; color: #fff; font-weight: 600; }
button.secondary { background: #eef1f6; color: #3d4d63; }
button:disabled { opacity: 0.5; cursor: default; }
.msg { margin-top: 10px; font-size: 0.9em; }
.msg.err { color: #c0392b; }
pre { background: #1e2633; color: #d6e2f0; border-radius: 8px; padding: 12px 14px; font-size: 0.8em; max-height: 420px; overflow: auto; white-space: pre-wrap; margin-top: 10px; }
.state { font-weight: 700; }
.banner { padding: 12px 16px; border-radius: 8px; margin-bottom: 14px; }
.banner.ready { background: #e8f8ef; color: #1e6b34; }
.banner.notready { background: #fff4e0; color: #6b4500; }
</style>
</head>
<body>
<header><div class="inner">
  <h1>TTN Site Setup</h1>
  <p>Data, models and results stored in this server's data folder. <a href="/">&larr; Back to the reports</a></p>
</div></header>
<main>
  <div id="banner"></div>
  <section>
    <h2>Volume contents</h2>
    <table><thead><tr><th></th><th>Item</th><th>Source</th><th>Location</th><th>Notes</th></tr></thead>
    <tbody id="items"></tbody></table>
    <p class="note" style="margin-top:10px">Licensed files (<span class="pill manual">manual</span>) cannot be downloaded
    automatically. Copy them to the location shown above (see <b>docker_readme/DEPLOY_DOCKER.md</b>, section &ldquo;授權資料&rdquo;), then reload this page.</p>
    <p class="note" id="summary" style="margin-top:6px"></p>
  </section>
  <section>
    <h2>Actions</h2>
    <p class="note">Each action runs on the GPU worker in the background; you can close this page. Steps that are already
    done are skipped. A full recompute takes several days (about 5,000 ClinVar variants).</p>
    <div class="actions" id="actions"></div>
    <div class="msg" id="msg"></div>
  </section>
  <section>
    <h2>Setup jobs <span class="note" id="worker"></span></h2>
    <table><thead><tr><th>Created</th><th>Steps</th><th>State</th><th></th></tr></thead><tbody id="jobs"></tbody></table>
    <pre id="log" hidden></pre>
  </section>
</main>
<script>
(function () {
  var ACTIONS = __ACTIONS__;
  var $ = function (id) { return document.getElementById(id); };
  var logJob = null;
  function esc(s) { return String(s == null ? '' : s).replace(/[&<>"']/g, function (c) {
    return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]; }); }
  function api(method, path, body) {
    return fetch(path, { method: method, headers: { 'X-TTN-Site': '1', 'Content-Type': 'application/json' },
                         body: body ? JSON.stringify(body) : undefined }).then(function (r) { return r.json(); });
  }
  Object.keys(ACTIONS).forEach(function (k) {
    var b = document.createElement('button');
    b.textContent = ACTIONS[k];
    if (k === 'site' || k === 'trial') b.className = 'secondary';
    b.addEventListener('click', function () {
      if (k === 'recompute' && !confirm('Recompute everything from scratch? This runs for several days and uses the GPU the whole time.')) return;
      $('msg').textContent = 'Submitting…'; $('msg').className = 'msg';
      api('POST', '/api/setup', { action: k }).then(function (res) {
        $('msg').textContent = res.ok ? 'Queued: ' + ACTIONS[k] : (res.error || 'Failed');
        $('msg').className = 'msg' + (res.ok ? '' : ' err');
        if (res.ok) { logJob = res.job.job_id; refresh(); }
      });
    });
    $('actions').appendChild(b);
  });
  function refresh() {
    api('GET', '/api/status').then(function (res) {
      var st = res.setup;
      $('items').innerHTML = st.items.map(function (i) {
        var src = i.source.indexOf('manual') === 0 ? '<span class="pill manual">manual</span>'
                                                   : '<span class="pill">auto · ' + esc(i.source.split(':')[1]) + '</span>';
        return '<tr><td class="' + (i.ok ? 'ok' : 'miss ' + i.level) + '">' + (i.ok ? '✓' : '✗') + '</td>' +
          '<td><b>' + esc(i.name) + '</b><br><span class="note">' + esc(i.level) + '</span></td><td>' + src + '</td>' +
          '<td class="path">' + esc(i.path) + '</td><td class="note">' + esc(i.note) + '</td></tr>';
      }).join('');
      $('summary').innerHTML = st.reports + ' report(s) in the data folder · website ' + (st.site_built ? 'built' : '<b>not built yet</b>') +
        ' · Hugging Face token ' + (st.hf_token ? 'set' : 'not set') + ' · NCBI email ' + (st.ncbi_email ? 'set' : '<b>not set</b>');
      var w = res.worker;
      $('worker').textContent = w.alive ? '· GPU worker ' + w.state : '· GPU worker is not running (check the worker container)';
      $('banner').innerHTML = st.ready_for_analysis
        ? '<div class="banner ready">Required data and models are in place; new variants can be analyzed.' +
          (st.site_built ? '' : ' The website has not been built yet: run &ldquo;Rebuild the website&rdquo; or a recompute.') + '</div>'
        : '<div class="banner notready">Missing before variants can be analyzed: ' + esc(st.missing_required.join(', ')) +
          '. Start with &ldquo;Download public data and models&rdquo;.</div>';
      var jobs = (res.setup_jobs || []).slice().reverse();
      if (!logJob && jobs.length) logJob = jobs[0].job_id;
      $('jobs').innerHTML = jobs.map(function (j) {
        var active = j.state === 'queued' || j.state === 'running';
        return '<tr><td>' + esc((j.created || '').replace('T', ' ')) + '</td><td>' + esc((j.stages || []).join(' → ')) +
          ((j.extra_args || []).length ? ' <span class="note">' + esc(j.extra_args.join(' ')) + '</span>' : '') + '</td>' +
          '<td class="state">' + esc(j.state) + (j.error ? '<br><span class="note">' + esc(j.error) + '</span>' : '') + '</td>' +
          '<td><button class="secondary" data-log="' + esc(j.job_id) + '">Log</button>' +
          (active ? ' <button class="secondary" data-cancel="' + esc(j.job_id) + '">Cancel</button>' : '') + '</td></tr>';
      }).join('');
      if (logJob) api('GET', '/api/job-log?job=' + encodeURIComponent(logJob)).then(function (r) {
        var pre = $('log'); pre.hidden = !r.log; var bottom = pre.scrollTop + pre.clientHeight >= pre.scrollHeight - 20;
        pre.textContent = r.log; if (bottom) pre.scrollTop = pre.scrollHeight;
      });
    });
  }
  $('jobs').addEventListener('click', function (e) {
    var t = e.target;
    if (t.dataset.log) { logJob = t.dataset.log; refresh(); }
    if (t.dataset.cancel && confirm('Cancel this setup job?')) api('POST', '/api/cancel', { job_id: t.dataset.cancel }).then(refresh);
  });
  refresh();
  setInterval(refresh, 5000);
})();
</script>
</body>
</html>
""".replace("__ACTIONS__", json.dumps({k: v[2] for k, v in SETUP_ACTIONS.items()}))


def main() -> int:
    port = int(os.environ.get("TTN_PORT") or 8080)
    bind = os.environ.get("TTN_BIND") or "0.0.0.0"
    SITE.mkdir(parents=True, exist_ok=True)
    spark_job.JOBS_DIR.mkdir(parents=True, exist_ok=True)
    nets = "any address" if NETWORKS is None else ", ".join(str(n) for n in NETWORKS)
    print(f"TTN site on http://{bind}:{port}/ serving {SITE}; allowed clients: {nets}", flush=True)
    ThreadingHTTPServer.daemon_threads = True
    with ThreadingHTTPServer((bind, port), Handler) as httpd:
        httpd.serve_forever()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
