//go:build !windows

package main

import (
	"net"
	"os"

	"golang.org/x/crypto/ssh/agent"
)

// dialAgent 連到 SSH_AUTH_SOCK 指定的 ssh-agent（沒有就回傳 nil）。
func dialAgent() signerSource {
	sock := os.Getenv("SSH_AUTH_SOCK")
	if sock == "" {
		return nil
	}
	conn, err := net.Dial("unix", sock)
	if err != nil {
		return nil
	}
	return agent.NewClient(conn)
}
