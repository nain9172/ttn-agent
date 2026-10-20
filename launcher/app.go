package main

import (
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"log"
	"net/http"
	"os"
	"path"
	"path/filepath"
	"regexp"
	"strings"
	"sync"
	"time"

	"github.com/gen2brain/beeep"
	"github.com/pkg/sftp"
	"golang.org/x/crypto/ssh"
)

const pollInterval = 30 * time.Second

var variantRe = regexp.MustCompile(`^(?i)(?:chr)?2-\d+-[ACGT]-[ACGT]$`)

type app struct {
	site  string
	port  int
	store *store
	auth  *sshAuth
	poll  sync.Mutex // 同一時間只做一輪輪詢／合併
	wake  chan struct{}

	errMu     sync.Mutex
	pollError string
}

func newApp(site string) (*app, error) {
	st, err := newStore(site)
	if err != nil {
		return nil, err
	}
	return &app{site: site, store: st, auth: &sshAuth{}, wake: make(chan struct{}, 1)}, nil
}

// ── HTTP ─────────────────────────────────────────────────────────────────────

func (a *app) routes() http.Handler {
	mux := http.NewServeMux()
	mux.Handle("/", http.FileServer(http.Dir(a.site)))
	api := map[string]http.HandlerFunc{
		"GET /api/ping":          a.handlePing,
		"GET /api/config":        a.handleGetConfig,
		"POST /api/config":       a.handleSetConfig,
		"POST /api/test":         a.handleTest,
		"POST /api/trust":        a.handleTrust,
		"POST /api/passphrase":   a.handlePassphrase,
		"POST /api/submit":       a.handleSubmit,
		"GET /api/jobs":          a.handleJobs,
		"GET /api/local-reports": a.handleLocalReports,
		"POST /api/setup-key":    a.handleSetupKey,
	}
	for pattern, h := range api {
		mux.Handle(pattern, a.guard(h))
	}
	return noCache(mux)
}

// guard 擋掉不是來自本機這個網頁的請求：Host 必須是 127.0.0.1/localhost（防 DNS rebinding），
// 而且要帶自訂標頭 X-TTN-Site（其他網站的跨站請求帶不了這個標頭，因為我們不回應 CORS 預檢）。
func (a *app) guard(h http.HandlerFunc) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		host := r.Host
		okHost := host == fmt.Sprintf("127.0.0.1:%d", a.port) || host == fmt.Sprintf("localhost:%d", a.port)
		if !okHost || r.Header.Get("X-TTN-Site") != "1" {
			http.Error(w, "forbidden", http.StatusForbidden)
			return
		}
		h(w, r)
	})
}

func noCache(h http.Handler) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		// 報告與索引會被合併更新，不能讓瀏覽器拿到舊版
		w.Header().Set("Cache-Control", "no-store")
		h.ServeHTTP(w, r)
	})
}

func writeJSONResponse(w http.ResponseWriter, v any) {
	w.Header().Set("Content-Type", "application/json; charset=utf-8")
	_ = json.NewEncoder(w).Encode(v)
}

func errorResponse(w http.ResponseWriter, err error) {
	resp := map[string]any{"ok": false, "error": err.Error()}
	var unknown errUnknownHost
	switch {
	case errors.As(err, &unknown):
		resp["need"] = "trust"
		resp["fingerprint"] = unknown.Fingerprint
	case errors.Is(err, errNeedPassphrase):
		resp["need"] = "passphrase"
	case errors.Is(err, errNoKey), strings.Contains(err.Error(), "unable to authenticate"):
		// 沒有金鑰，或 Spark 不認得這把金鑰：請使用者用密碼做一次金鑰設定
		resp["need"] = "setup"
	}
	writeJSONResponse(w, resp)
}

func decode(r *http.Request, v any) error {
	return json.NewDecoder(io.LimitReader(r.Body, 1<<20)).Decode(v)
}

func (a *app) handlePing(w http.ResponseWriter, _ *http.Request) {
	writeJSONResponse(w, map[string]any{"ok": true, "version": version})
}

func (a *app) handleGetConfig(w http.ResponseWriter, _ *http.Request) {
	writeJSONResponse(w, map[string]any{"ok": true, "config": a.store.config(), "known_hosts": knownHostsPath()})
}

func (a *app) handleSetConfig(w http.ResponseWriter, r *http.Request) {
	var c Config
	if err := decode(r, &c); err != nil {
		errorResponse(w, err)
		return
	}
	c.Host, c.User = strings.TrimSpace(c.Host), strings.TrimSpace(c.User)
	if err := a.store.setConfig(c); err != nil {
		errorResponse(w, err)
		return
	}
	writeJSONResponse(w, map[string]any{"ok": true})
}

func (a *app) handleTest(w http.ResponseWriter, _ *http.Request) {
	a.handleTestWith(w, nil)
}

// handleTestWith 連線並在 Spark 上跑 spark_job.py check；extra 會併進回應。
func (a *app) handleTestWith(w http.ResponseWriter, extra map[string]any) {
	cfg := a.store.config()
	client, err := a.auth.dial(cfg)
	if err != nil {
		errorResponse(w, err)
		return
	}
	defer client.Close()
	res, err := sparkJob(client, cfg, nil, "check")
	if err != nil {
		errorResponse(w, err)
		return
	}
	resp := map[string]any{"ok": res["ok"] == true, "check": res}
	for k, v := range extra {
		resp[k] = v
	}
	writeJSONResponse(w, resp)
}

func (a *app) handleTrust(w http.ResponseWriter, r *http.Request) {
	var body struct{ Fingerprint string }
	if err := decode(r, &body); err != nil {
		errorResponse(w, err)
		return
	}
	a.auth.mu.Lock()
	ok := body.Fingerprint != "" && body.Fingerprint == a.auth.lastHostKey
	if ok {
		a.auth.trusted = body.Fingerprint
	}
	a.auth.mu.Unlock()
	if !ok {
		errorResponse(w, errors.New("fingerprint does not match the key the Spark presented; test the connection again"))
		return
	}
	a.handleTest(w, r) // 立即連線一次，把指紋寫進 known_hosts
}

func (a *app) handlePassphrase(w http.ResponseWriter, r *http.Request) {
	var body struct{ Passphrase string }
	if err := decode(r, &body); err != nil {
		errorResponse(w, err)
		return
	}
	a.auth.mu.Lock()
	a.auth.passphrase = []byte(body.Passphrase)
	a.auth.mu.Unlock()
	a.handleTest(w, r)
}

func (a *app) handleSubmit(w http.ResponseWriter, r *http.Request) {
	var body struct {
		Variants []string
		Rerun    bool // 已在 site 裡的變異也要重新跑（Spark 端不沿用舊報告）
	}
	if err := decode(r, &body); err != nil {
		errorResponse(w, err)
		return
	}
	var variants []string
	var rejected []map[string]any
	for _, v := range body.Variants {
		v = strings.TrimSpace(v)
		if v == "" {
			continue
		}
		if !variantRe.MatchString(v) {
			rejected = append(rejected, map[string]any{"variant": v, "reason": "format must be 2-<position>-<ref>-<alt>"})
			continue
		}
		variants = append(variants, strings.ToUpper(strings.TrimPrefix(strings.ToLower(v), "chr")))
	}
	if len(variants) == 0 {
		writeJSONResponse(w, map[string]any{"ok": false, "error": "no valid variants", "rejected": rejected})
		return
	}

	cfg := a.store.config()
	client, err := a.auth.dial(cfg)
	if err != nil {
		errorResponse(w, err)
		return
	}
	defer client.Close()
	host, _ := os.Hostname()
	res, err := sparkJob(client, cfg, map[string]any{"variants": variants, "client": host, "rerun": body.Rerun}, "submit")
	if err != nil {
		errorResponse(w, err)
		return
	}
	if res["ok"] != true {
		if rj, ok := res["rejected"].([]any); ok {
			for _, x := range rj {
				if m, ok := x.(map[string]any); ok {
					rejected = append(rejected, m)
				}
			}
		}
		writeJSONResponse(w, map[string]any{"ok": false, "error": res["error"], "rejected": rejected})
		return
	}
	job := jobFromRemote(res["job"])
	job.Host = cfg.User + "@" + cfg.Host
	job.Rejected = append(rejected, job.Rejected...)
	if err := a.store.addJob(job); err != nil {
		errorResponse(w, err)
		return
	}
	a.kick()
	writeJSONResponse(w, map[string]any{"ok": true, "job": job})
}

func (a *app) handleJobs(w http.ResponseWriter, _ *http.Request) {
	a.errMu.Lock()
	pollErr := a.pollError
	a.errMu.Unlock()
	writeJSONResponse(w, map[string]any{"ok": true, "jobs": a.store.snapshot(), "poll_error": pollErr})
}

func (a *app) handleLocalReports(w http.ResponseWriter, _ *http.Request) {
	records, err := a.loadLocalRecords()
	if err != nil {
		errorResponse(w, err)
		return
	}
	writeJSONResponse(w, map[string]any{"ok": true, "reports": records})
}

// ── 輪詢與合併 ───────────────────────────────────────────────────────────────

func (a *app) kick() {
	select {
	case a.wake <- struct{}{}:
	default:
	}
}

func (a *app) pollLoop() {
	ticker := time.NewTicker(pollInterval)
	defer ticker.Stop()
	for {
		a.pollOnce()
		select {
		case <-ticker.C:
		case <-a.wake:
		}
	}
}

func (a *app) setPollError(err error) {
	a.errMu.Lock()
	defer a.errMu.Unlock()
	if err == nil {
		a.pollError = ""
	} else {
		a.pollError = err.Error()
	}
}

func (a *app) pollOnce() {
	a.poll.Lock()
	defer a.poll.Unlock()

	var pending []string
	for _, j := range a.store.snapshot() {
		if !j.finished() {
			pending = append(pending, j.JobID)
		}
	}
	if len(pending) == 0 {
		a.setPollError(nil)
		return
	}
	cfg := a.store.config()
	client, err := a.auth.dial(cfg)
	if err != nil {
		a.setPollError(err)
		return
	}
	defer client.Close()
	res, err := sparkJob(client, cfg, nil, append([]string{"status"}, pending...)...)
	if err != nil {
		a.setPollError(err)
		return
	}
	a.setPollError(nil)

	remote, _ := res["jobs"].([]any)
	for _, raw := range remote {
		m, ok := raw.(map[string]any)
		if !ok {
			continue
		}
		r := jobFromRemote(m)
		patchDir, _ := m["patch_dir"].(string)
		patchFiles := toStrings(m["patch_files"])
		var justFinished bool
		_ = a.store.update(r.JobID, func(j *Job) {
			wasDone := j.State == "done" || j.State == "failed"
			j.State, j.Completed, j.FailedVariants = r.State, r.Completed, r.FailedVariants
			j.Started, j.Finished, j.Error, j.Note = r.Started, r.Finished, r.Error, r.Note
			justFinished = !wasDone && (r.State == "done" || r.State == "failed")
		})
		if r.State == "done" {
			job := a.findJob(r.JobID)
			if job != nil && !job.Merged {
				reports, err := a.mergePatch(client, patchDir, patchFiles, r.JobID)
				_ = a.store.update(r.JobID, func(j *Job) {
					if err != nil {
						j.MergeError = err.Error()
					} else {
						j.Merged, j.MergeError, j.Reports = true, "", reports
					}
				})
				if err == nil {
					notify("TTN analysis finished",
						fmt.Sprintf("%d of %d variant(s) were added to the site.", len(reports), len(job.Variants)))
				} else {
					log.Printf("merge %s: %v", r.JobID, err)
				}
			}
		} else if r.State == "failed" && justFinished {
			notify("TTN analysis failed", r.Error)
		}
	}
}

func (a *app) findJob(id string) *Job {
	for _, j := range a.store.snapshot() {
		if j.JobID == id {
			jj := j
			return &jj
		}
	}
	return nil
}

// mergePatch 把 Spark 上的補丁下載進 site：報告與圖檔照路徑放，索引列併進 data/local-reports.js。
func (a *app) mergePatch(client *ssh.Client, patchDir string, files []string, jobID string) ([]string, error) {
	if patchDir == "" || len(files) == 0 {
		return nil, errors.New("the Spark reported no patch files")
	}
	sc, err := sftp.NewClient(client)
	if err != nil {
		return nil, fmt.Errorf("cannot open SFTP: %w", err)
	}
	defer sc.Close()

	var records []map[string]any
	var reports []string
	for _, rel := range files {
		clean := path.Clean(rel)
		allowed := clean == "records.json" || strings.HasPrefix(clean, "reports/") || strings.HasPrefix(clean, "assets/")
		if !allowed || strings.Contains(clean, "..") || path.IsAbs(clean) {
			return nil, fmt.Errorf("refusing unexpected patch path %q", rel)
		}
		data, err := readRemote(sc, path.Join(patchDir, clean))
		if err != nil {
			return nil, err
		}
		if clean == "records.json" {
			if err := json.Unmarshal(data, &records); err != nil {
				return nil, fmt.Errorf("bad records.json: %w", err)
			}
			continue
		}
		dest := filepath.Join(a.site, filepath.FromSlash(clean))
		if err := os.MkdirAll(filepath.Dir(dest), 0o755); err != nil {
			return nil, err
		}
		if err := os.WriteFile(dest+".tmp", data, 0o644); err != nil {
			return nil, err
		}
		if err := os.Rename(dest+".tmp", dest); err != nil {
			return nil, err
		}
		if strings.HasPrefix(clean, "reports/") {
			reports = append(reports, strings.TrimPrefix(clean, "reports/"))
		}
	}
	for _, rec := range records {
		rec["local_job"] = jobID
	}
	if err := a.mergeLocalRecords(records); err != nil {
		return nil, err
	}
	return reports, nil
}

func readRemote(sc *sftp.Client, p string) ([]byte, error) {
	f, err := sc.Open(p)
	if err != nil {
		return nil, fmt.Errorf("cannot download %s: %w", p, err)
	}
	defer f.Close()
	return io.ReadAll(io.LimitReader(f, 200<<20))
}

func (a *app) localRecordsPath() string {
	return filepath.Join(a.site, "data", "local-reports.json")
}

func (a *app) loadLocalRecords() ([]map[string]any, error) {
	var records []map[string]any
	err := readJSON(a.localRecordsPath(), &records)
	if errors.Is(err, os.ErrNotExist) {
		return []map[string]any{}, nil
	}
	return records, err
}

// mergeLocalRecords 依 variant_id 取代舊列，再同時寫出 .json（給啟動程式自己讀）與
// .js（index.html 用 <script> 載入，直接開檔案時也看得到新報告）。
func (a *app) mergeLocalRecords(add []map[string]any) error {
	existing, err := a.loadLocalRecords()
	if err != nil {
		return err
	}
	byID := map[string]int{}
	for i, r := range existing {
		if id, ok := r["variant_id"].(string); ok {
			byID[id] = i
		}
	}
	for _, r := range add {
		id, _ := r["variant_id"].(string)
		i, ok := byID[id]
		if !ok {
			// 舊版啟動程式寫入的列以轉錄方向 ID 為 variant_id；新列的 internal_id 就是它
			if internal, _ := r["internal_id"].(string); internal != "" {
				i, ok = byID[internal]
				if ok {
					byID[id] = i
				}
			}
		}
		if ok {
			existing[i] = r
		} else {
			byID[id] = len(existing)
			existing = append(existing, r)
		}
	}
	if err := writeJSON(a.localRecordsPath(), existing); err != nil {
		return err
	}
	data, err := json.Marshal(existing)
	if err != nil {
		return err
	}
	js := "window.__TTN_LOCAL_REPORTS__ = " + strings.ReplaceAll(string(data), "</", `<\/`) + ";\n"
	jsPath := filepath.Join(a.site, "data", "local-reports.js")
	if err := os.WriteFile(jsPath+".tmp", []byte(js), 0o644); err != nil {
		return err
	}
	return os.Rename(jsPath+".tmp", jsPath)
}

func jobFromRemote(raw any) *Job {
	m, _ := raw.(map[string]any)
	j := &Job{}
	if m == nil {
		return j
	}
	j.JobID, _ = m["job_id"].(string)
	j.State, _ = m["state"].(string)
	j.Created, _ = m["created"].(string)
	j.Started, _ = m["started"].(string)
	j.Finished, _ = m["finished"].(string)
	j.Error, _ = m["error"].(string)
	j.Note, _ = m["note"].(string)
	j.Rerun, _ = m["rerun"].(bool)
	j.Variants = toStrings(m["variants"])
	j.Completed = toStrings(m["completed"])
	j.FailedVariants = toStrings(m["failed_variants"])
	if rj, ok := m["rejected"].([]any); ok {
		for _, x := range rj {
			if mm, ok := x.(map[string]any); ok {
				j.Rejected = append(j.Rejected, mm)
			}
		}
	}
	return j
}

func toStrings(v any) []string {
	list, _ := v.([]any)
	out := make([]string, 0, len(list))
	for _, x := range list {
		if s, ok := x.(string); ok {
			out = append(out, s)
		}
	}
	return out
}

func notify(title, message string) {
	if err := beeep.Notify(title, message, ""); err != nil {
		log.Printf("notification: %v", err)
	}
}
