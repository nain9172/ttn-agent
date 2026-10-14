package main

import (
	"bytes"
	"encoding/json"
	"errors"
	"fmt"
	"net"
	"os"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"
	"sync"
	"time"

	"golang.org/x/crypto/ssh"
	"golang.org/x/crypto/ssh/knownhosts"
)

// errUnknownHost：第一次連到這台 Spark，要使用者在網頁上確認主機指紋。
type errUnknownHost struct{ Fingerprint string }

func (e errUnknownHost) Error() string {
	return "this Spark is not in known_hosts yet; confirm its fingerprint " + e.Fingerprint
}

// errNeedPassphrase：私鑰有密碼，而且 ssh-agent 裡沒有可用的金鑰。
var errNeedPassphrase = errors.New("the SSH private key is protected by a passphrase")

// errNoKey：這台電腦還沒有 SSH 金鑰，可以在網頁上用 Spark 密碼做一次金鑰設定。
var errNoKey = errors.New("no SSH key on this computer yet")

type sshAuth struct {
	mu          sync.Mutex
	passphrase  []byte // 只留在記憶體，不寫檔
	trusted     string // 使用者在網頁上確認過、但還沒寫入 known_hosts 的指紋
	lastHostKey string
}

func sshDir() string {
	home, err := os.UserHomeDir()
	if err != nil {
		return ".ssh"
	}
	return filepath.Join(home, ".ssh")
}

func knownHostsPath() string { return filepath.Join(sshDir(), "known_hosts") }

func (a *sshAuth) hostKeyCallback() (ssh.HostKeyCallback, error) {
	path := knownHostsPath()
	if _, err := os.Stat(path); errors.Is(err, os.ErrNotExist) {
		if err := os.MkdirAll(filepath.Dir(path), 0o700); err != nil {
			return nil, err
		}
		if err := os.WriteFile(path, nil, 0o600); err != nil {
			return nil, err
		}
	}
	check, err := knownhosts.New(path)
	if err != nil {
		return nil, fmt.Errorf("cannot read %s: %w", path, err)
	}
	return func(hostname string, remote net.Addr, key ssh.PublicKey) error {
		err := check(hostname, remote, key)
		var keyErr *knownhosts.KeyError
		if err == nil || !errors.As(err, &keyErr) {
			return err
		}
		if len(keyErr.Want) > 0 {
			// 指紋跟以前記錄的不同：可能是重灌，也可能是中間人攻擊，一律不自動接受
			return fmt.Errorf("HOST KEY CHANGED for %s: known_hosts has a different key. "+
				"If the Spark was reinstalled, remove its line from %s and try again", hostname, path)
		}
		fp := ssh.FingerprintSHA256(key)
		a.mu.Lock()
		trusted := a.trusted == fp
		a.lastHostKey = fp
		a.mu.Unlock()
		if !trusted {
			return errUnknownHost{Fingerprint: fp}
		}
		line := knownhosts.Line([]string{knownhosts.Normalize(hostname)}, key)
		f, err := os.OpenFile(path, os.O_APPEND|os.O_WRONLY, 0o600)
		if err != nil {
			return err
		}
		defer f.Close()
		_, err = f.WriteString(line + "\n")
		return err
	}, nil
}

func (a *sshAuth) methods(cfg Config) ([]ssh.AuthMethod, error) {
	var signers []ssh.Signer
	if ag := dialAgent(); ag != nil {
		if s, err := ag.Signers(); err == nil {
			signers = append(signers, s...)
		}
	}

	paths := []string{cfg.KeyPath}
	if cfg.KeyPath == "" {
		paths = []string{
			filepath.Join(sshDir(), "id_ed25519"),
			filepath.Join(sshDir(), "id_ecdsa"),
			filepath.Join(sshDir(), "id_rsa"),
		}
	}
	needPass := false
	for _, p := range paths {
		data, err := os.ReadFile(expandLocal(p))
		if err != nil {
			continue
		}
		signer, err := ssh.ParsePrivateKey(data)
		var missing *ssh.PassphraseMissingError
		if errors.As(err, &missing) {
			a.mu.Lock()
			pass := a.passphrase
			a.mu.Unlock()
			if pass == nil {
				needPass = true
				continue
			}
			signer, err = ssh.ParsePrivateKeyWithPassphrase(data, pass)
		}
		if err != nil {
			continue
		}
		signers = append(signers, signer)
	}
	if len(signers) == 0 {
		if needPass {
			return nil, errNeedPassphrase
		}
		return nil, fmt.Errorf("%w (looked for %s)", errNoKey, strings.Join(paths, ", "))
	}
	return []ssh.AuthMethod{ssh.PublicKeys(signers...)}, nil
}

func expandLocal(p string) string {
	if strings.HasPrefix(p, "~/") || strings.HasPrefix(p, `~\`) {
		if home, err := os.UserHomeDir(); err == nil {
			return filepath.Join(home, p[2:])
		}
	}
	return p
}

func (a *sshAuth) dial(cfg Config) (*ssh.Client, error) {
	if cfg.Host == "" || cfg.User == "" {
		return nil, errors.New("set the Spark host and user first")
	}
	hk, err := a.hostKeyCallback()
	if err != nil {
		return nil, err
	}
	auth, err := a.methods(cfg)
	if err != nil {
		return nil, err
	}
	addr := net.JoinHostPort(cfg.Host, strconv.Itoa(cfg.Port))
	client, err := ssh.Dial("tcp", addr, &ssh.ClientConfig{
		User:            cfg.User,
		Auth:            auth,
		HostKeyCallback: hk,
		Timeout:         15 * time.Second,
	})
	if err != nil {
		var unknown errUnknownHost
		if errors.As(err, &unknown) {
			return nil, unknown
		}
		return nil, fmt.Errorf("cannot connect to %s: %w", addr, err)
	}
	return client, nil
}

// shellPath 把設定裡的路徑轉成安全的 shell 片段；開頭的 ~/ 改成 "$HOME"/，其餘一律單引號。
func shellPath(p string) string {
	if p == "~" {
		return `"$HOME"`
	}
	if strings.HasPrefix(p, "~/") {
		return `"$HOME"/` + shellQuote(p[2:])
	}
	return shellQuote(p)
}

func shellQuote(s string) string {
	return "'" + strings.ReplaceAll(s, "'", `'"'"'`) + "'"
}

var jobIDRe = regexp.MustCompile(`^\d{8}-\d{6}-[0-9a-f]{6}$`)

// sparkJob 在 Spark 上執行 spark_job.py <args>，stdin 傳 JSON，回傳解析後的 JSON。
// 參數只會是固定的子指令或驗證過格式的 job id；變異清單一律走 stdin，不進命令列。
func sparkJob(client *ssh.Client, cfg Config, stdin any, args ...string) (map[string]any, error) {
	for _, a := range args[1:] {
		if !jobIDRe.MatchString(a) {
			return nil, fmt.Errorf("invalid job id %q", a)
		}
	}
	cmd := "cd " + shellPath(cfg.Repo) + " && " + shellPath(cfg.Python) + " spark_job.py " + strings.Join(args, " ")
	session, err := client.NewSession()
	if err != nil {
		return nil, err
	}
	defer session.Close()
	var out, errOut bytes.Buffer
	session.Stdout, session.Stderr = &out, &errOut
	if stdin != nil {
		data, _ := json.Marshal(stdin)
		session.Stdin = bytes.NewReader(data)
	}
	runErr := session.Run(cmd)
	// spark_job.py 只印一行 JSON；前面若有其他輸出（例如 .bashrc 的訊息）就略過
	lines := strings.Split(strings.TrimSpace(out.String()), "\n")
	last := strings.TrimSpace(lines[len(lines)-1])
	var result map[string]any
	if err := json.Unmarshal([]byte(last), &result); err != nil {
		msg := strings.TrimSpace(errOut.String())
		if msg == "" {
			msg = strings.TrimSpace(out.String())
		}
		if runErr != nil {
			return nil, fmt.Errorf("spark_job.py failed on the Spark: %v: %s", runErr, tail(msg, 600))
		}
		return nil, fmt.Errorf("unexpected output from spark_job.py: %s", tail(msg, 600))
	}
	return result, nil
}

func tail(s string, n int) string {
	if len(s) <= n {
		return s
	}
	return "…" + s[len(s)-n:]
}

// signerSource 是 ssh-agent 用到的最小介面。
type signerSource interface {
	Signers() ([]ssh.Signer, error)
}
