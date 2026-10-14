package main

import (
	"encoding/json"
	"errors"
	"os"
	"path/filepath"
	"sync"
)

// Config 是 Spark 連線設定，存在使用者的設定資料夾（%AppData%\TTN-Site\config.json），
// 換一份新下載的 site 也不用重填。
type Config struct {
	Host    string `json:"host"`
	Port    int    `json:"port"`
	User    string `json:"user"`
	KeyPath string `json:"key_path"` // 空白 = 依序嘗試 ~/.ssh/id_ed25519、id_ecdsa、id_rsa
	Repo    string `json:"repo"`     // Spark 上 ttn-agent 的路徑
	Python  string `json:"python"`   // Spark 上 gm4 環境的 python
}

func defaultConfig() Config {
	return Config{
		Port:   22,
		Repo:   "~/ttn-agent",
		Python: "~/miniconda3/envs/gm4/bin/python",
	}
}

// Job 是本機記錄的任務。State/Completed 等欄位來自 Spark 的 spark_job.py status；
// Merged 表示補丁已經下載並合併進這份 site。
type Job struct {
	JobID          string           `json:"job_id"`
	State          string           `json:"state"`
	Variants       []string         `json:"variants"`
	Rerun          bool             `json:"rerun,omitempty"`
	Rejected       []map[string]any `json:"rejected,omitempty"`
	Completed      []string         `json:"completed"`
	FailedVariants []string         `json:"failed_variants,omitempty"`
	Created        string           `json:"created,omitempty"`
	Started        string           `json:"started,omitempty"`
	Finished       string           `json:"finished,omitempty"`
	Error          string           `json:"error,omitempty"`
	Note           string           `json:"note,omitempty"` // 例：等待 Spark 記憶體空出來（誰在佔用）
	Host           string           `json:"host"`
	Merged         bool             `json:"merged"`
	MergeError     string           `json:"merge_error,omitempty"`
	Reports        []string         `json:"reports,omitempty"` // 合併後的報告檔名（reports/ 底下）
}

func (j *Job) finished() bool {
	return (j.State == "done" && j.Merged) || j.State == "failed"
}

type store struct {
	mu         sync.Mutex
	configPath string
	jobsPath   string
	cfg        Config
	jobs       []*Job
}

func newStore(site string) (*store, error) {
	dir, err := os.UserConfigDir()
	if err != nil {
		dir = site
	}
	dir = filepath.Join(dir, "TTN-Site")
	if err := os.MkdirAll(dir, 0o700); err != nil {
		return nil, err
	}
	s := &store{
		configPath: filepath.Join(dir, "config.json"),
		// 任務紀錄跟著 site 走：合併進哪一份 site，就記在哪一份裡
		jobsPath: filepath.Join(site, "data", "local-jobs.json"),
		cfg:      defaultConfig(),
	}
	if err := readJSON(s.configPath, &s.cfg); err != nil && !errors.Is(err, os.ErrNotExist) {
		return nil, err
	}
	if err := readJSON(s.jobsPath, &s.jobs); err != nil && !errors.Is(err, os.ErrNotExist) {
		return nil, err
	}
	return s, nil
}

func (s *store) config() Config {
	s.mu.Lock()
	defer s.mu.Unlock()
	return s.cfg
}

func (s *store) setConfig(c Config) error {
	s.mu.Lock()
	defer s.mu.Unlock()
	if c.Port == 0 {
		c.Port = 22
	}
	s.cfg = c
	return writeJSON(s.configPath, c)
}

func (s *store) snapshot() []Job {
	s.mu.Lock()
	defer s.mu.Unlock()
	out := make([]Job, len(s.jobs))
	for i, j := range s.jobs {
		out[i] = *j
	}
	return out
}

func (s *store) addJob(j *Job) error {
	s.mu.Lock()
	defer s.mu.Unlock()
	s.jobs = append(s.jobs, j)
	return writeJSON(s.jobsPath, s.jobs)
}

// update 在鎖內修改指定任務並存檔。
func (s *store) update(id string, fn func(*Job)) error {
	s.mu.Lock()
	defer s.mu.Unlock()
	for _, j := range s.jobs {
		if j.JobID == id {
			fn(j)
		}
	}
	return writeJSON(s.jobsPath, s.jobs)
}

func readJSON(path string, v any) error {
	data, err := os.ReadFile(path)
	if err != nil {
		return err
	}
	return json.Unmarshal(data, v)
}

// writeJSON 先寫暫存檔再改名，避免寫到一半被讀到。
func writeJSON(path string, v any) error {
	data, err := json.MarshalIndent(v, "", "  ")
	if err != nil {
		return err
	}
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	tmp := path + ".tmp"
	if err := os.WriteFile(tmp, data, 0o644); err != nil {
		return err
	}
	return os.Rename(tmp, path)
}
