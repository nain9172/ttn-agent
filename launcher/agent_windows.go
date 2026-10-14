//go:build windows

package main

import (
	"time"

	"github.com/Microsoft/go-winio"
	"golang.org/x/crypto/ssh/agent"
)

// dialAgent 連到 Windows 內建 OpenSSH 的 ssh-agent（服務沒開就回傳 nil）。
func dialAgent() signerSource {
	timeout := 2 * time.Second
	conn, err := winio.DialPipe(`\\.\pipe\openssh-ssh-agent`, &timeout)
	if err != nil {
		return nil
	}
	return agent.NewClient(conn)
}
