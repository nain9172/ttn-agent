#!/usr/bin/env bash
# 編譯 Windows 版 TTN-Site 啟動程式到 launcher/dist/；build_site.py 會自動把它們放進 site/。
# 需要 Go（https://go.dev/dl/）；在 Linux / DGX Spark 上就能交叉編譯 Windows 版。
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p dist
go test ./...
GOOS=windows GOARCH=amd64 go build -trimpath -ldflags "-s -w" -o dist/TTN-Site.exe .
GOOS=windows GOARCH=arm64 go build -trimpath -ldflags "-s -w" -o dist/TTN-Site-arm64.exe .
ls -l dist
