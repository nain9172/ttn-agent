package main

import (
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"golang.org/x/crypto/ssh"
)

// 產生的金鑰要能被 OpenSSH（Windows 內建的 ssh.exe 也是 OpenSSH）直接讀取。
func TestEnsureKeyIsOpenSSHCompatible(t *testing.T) {
	dir := t.TempDir()
	keyPath := filepath.Join(dir, "id_ed25519")
	path, pub, created, err := ensureKey(Config{KeyPath: keyPath})
	if err != nil || !created || path != keyPath {
		t.Fatalf("ensureKey: %v created=%v path=%s", err, created, path)
	}
	data, _ := os.ReadFile(keyPath)
	signer, err := ssh.ParsePrivateKey(data)
	if err != nil {
		t.Fatalf("parse private key: %v", err)
	}
	if !strings.HasPrefix(pub, strings.TrimSpace(string(ssh.MarshalAuthorizedKey(signer.PublicKey())))) {
		t.Fatalf("public key does not match private key")
	}
	if keygen, err := exec.LookPath("ssh-keygen"); err == nil {
		out, err := exec.Command(keygen, "-y", "-f", keyPath).Output()
		if err != nil || !strings.HasPrefix(pub, strings.TrimSpace(string(out))) {
			t.Fatalf("ssh-keygen cannot read the key: %v %s", err, out)
		}
	}
	// 第二次呼叫沿用同一把，不會覆寫
	_, pub2, created2, err := ensureKey(Config{KeyPath: keyPath})
	if err != nil || created2 || pub2 != pub {
		t.Fatalf("second ensureKey should reuse the key: %v created=%v", err, created2)
	}
}
