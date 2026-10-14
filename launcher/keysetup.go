package main

import (
	"crypto/ed25519"
	"crypto/rand"
	"encoding/pem"
	"errors"
	"fmt"
	"net"
	"net/http"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"time"

	"golang.org/x/crypto/ssh"
)

// 讓使用者不必自己開 PowerShell 跑 ssh-keygen / ssh-copy-id：
//  1. 本機沒有金鑰就產生 ~/.ssh/id_ed25519（與 Windows 內建 OpenSSH 相容的格式）
//  2. 用使用者輸入的 Spark 密碼登入一次，把公鑰加進 Spark 的 ~/.ssh/authorized_keys
// 密碼只用於這一次連線，不寫檔、不保留在記憶體。

func defaultKeyPath() string { return filepath.Join(sshDir(), "id_ed25519") }

// ensureKey 回傳私鑰路徑與公鑰（authorized_keys 格式）；沒有就產生一把。
func ensureKey(cfg Config) (string, string, bool, error) {
	keyPath := expandLocal(cfg.KeyPath)
	if keyPath == "" {
		keyPath = defaultKeyPath()
	}
	if data, err := os.ReadFile(keyPath + ".pub"); err == nil {
		return keyPath, strings.TrimSpace(string(data)), false, nil
	}
	if _, err := os.Stat(keyPath); err == nil {
		return "", "", false, fmt.Errorf("%s exists but %s.pub is missing; "+
			"recreate the .pub file or leave the key field blank to use a new key", keyPath, keyPath)
	}

	pub, priv, err := ed25519.GenerateKey(rand.Reader)
	if err != nil {
		return "", "", false, err
	}
	host, _ := os.Hostname()
	comment := "ttn-site@" + host
	block, err := ssh.MarshalPrivateKey(priv, comment)
	if err != nil {
		return "", "", false, err
	}
	sshPub, err := ssh.NewPublicKey(pub)
	if err != nil {
		return "", "", false, err
	}
	pubLine := strings.TrimSpace(string(ssh.MarshalAuthorizedKey(sshPub))) + " " + comment
	if err := os.MkdirAll(filepath.Dir(keyPath), 0o700); err != nil {
		return "", "", false, err
	}
	if err := os.WriteFile(keyPath, pem.EncodeToMemory(block), 0o600); err != nil {
		return "", "", false, err
	}
	if err := os.WriteFile(keyPath+".pub", []byte(pubLine+"\n"), 0o644); err != nil {
		return "", "", false, err
	}
	return keyPath, pubLine, true, nil
}

func (a *app) handleSetupKey(w http.ResponseWriter, r *http.Request) {
	var body struct{ Password string }
	if err := decode(r, &body); err != nil {
		errorResponse(w, err)
		return
	}
	cfg := a.store.config()
	if cfg.Host == "" || cfg.User == "" {
		errorResponse(w, errors.New("fill in the Spark host and user name first"))
		return
	}
	if body.Password == "" {
		errorResponse(w, errors.New("enter the password you use to log in to the Spark"))
		return
	}
	keyPath, pubLine, created, err := ensureKey(cfg)
	if err != nil {
		errorResponse(w, err)
		return
	}

	hk, err := a.auth.hostKeyCallback()
	if err != nil {
		errorResponse(w, err)
		return
	}
	password := body.Password
	addr := net.JoinHostPort(cfg.Host, strconv.Itoa(cfg.Port))
	client, err := ssh.Dial("tcp", addr, &ssh.ClientConfig{
		User: cfg.User,
		Auth: []ssh.AuthMethod{
			ssh.Password(password),
			// 有些 sshd 只開 keyboard-interactive（PAM），題目通常就是一行密碼
			ssh.KeyboardInteractive(func(_, _ string, questions []string, _ []bool) ([]string, error) {
				answers := make([]string, len(questions))
				for i := range answers {
					answers[i] = password
				}
				return answers, nil
			}),
		},
		HostKeyCallback: hk,
		Timeout:         15 * time.Second,
	})
	if err != nil {
		var unknown errUnknownHost
		if errors.As(err, &unknown) {
			errorResponse(w, unknown)
			return
		}
		errorResponse(w, fmt.Errorf("password login to %s failed: %w "+
			"(check the user name and password; if the Spark does not allow password login, "+
			"ask its admin to add this public key: %s)", addr, err, pubLine))
		return
	}
	defer client.Close()

	session, err := client.NewSession()
	if err != nil {
		errorResponse(w, err)
		return
	}
	defer session.Close()
	// 冪等：同一把公鑰已經在 authorized_keys 裡就不重複加
	session.Stdin = strings.NewReader(pubLine + "\n")
	cmd := `umask 077; mkdir -p "$HOME/.ssh" && touch "$HOME/.ssh/authorized_keys" && ` +
		`read -r key && (grep -qxF "$key" "$HOME/.ssh/authorized_keys" || printf '%s\n' "$key" >> "$HOME/.ssh/authorized_keys")`
	if out, err := session.CombinedOutput(cmd); err != nil {
		errorResponse(w, fmt.Errorf("could not install the key on the Spark: %v %s", err, strings.TrimSpace(string(out))))
		return
	}

	// 用新的金鑰設定，接著跑一次正式的連線測試
	if cfg.KeyPath == "" && keyPath != defaultKeyPath() {
		cfg.KeyPath = keyPath
		_ = a.store.setConfig(cfg)
	}
	a.handleTestWith(w, map[string]any{"key_created": created, "key_path": keyPath})
}
