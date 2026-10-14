// TTN-Site launcher：放在下載的 site 資料夾裡，雙擊執行。
//
// 它在本機 127.0.0.1 提供 site 網頁，並代替瀏覽器做網頁做不到的事：
// 用 SSH 把使用者輸入的變異送到 DGX Spark、追蹤進度、把跑完的報告下載回
// site 資料夾，最後跳出 Windows 通知。網頁介面就是 site 首頁上的
// 「Analyze new variants」面板。
package main

import (
	"errors"
	"flag"
	"fmt"
	"log"
	"net"
	"net/http"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"time"
)

const version = "0.1.0"

func main() {
	siteFlag := flag.String("site", "", "site folder (default: the folder this program is in)")
	portFlag := flag.Int("port", 8765, "first local port to try")
	noBrowser := flag.Bool("no-browser", false, "do not open the browser")
	flag.Parse()

	site, err := findSite(*siteFlag)
	if err != nil {
		fail(err)
	}
	app, err := newApp(site)
	if err != nil {
		fail(err)
	}

	ln, port, err := listen(*portFlag)
	if err != nil {
		fail(err)
	}
	app.port = port
	url := fmt.Sprintf("http://127.0.0.1:%d/", port)

	fmt.Printf("TTN-Site launcher %s\n", version)
	fmt.Printf("  site   : %s\n", site)
	fmt.Printf("  open   : %s\n", url)
	fmt.Println("  Keep this window open while you use the site. Close it to stop.")
	fmt.Println()

	go app.pollLoop()
	if !*noBrowser {
		go func() {
			time.Sleep(400 * time.Millisecond)
			openBrowser(url)
		}()
	}
	log.Fatal(http.Serve(ln, app.routes()))
}

// findSite 找出 site 資料夾：參數 > 執行檔所在資料夾 > 目前目錄。
func findSite(flagValue string) (string, error) {
	var candidates []string
	if flagValue != "" {
		candidates = append(candidates, flagValue)
	} else {
		if exe, err := os.Executable(); err == nil {
			candidates = append(candidates, filepath.Dir(exe))
		}
		if wd, err := os.Getwd(); err == nil {
			candidates = append(candidates, wd)
		}
	}
	for _, dir := range candidates {
		abs, err := filepath.Abs(dir)
		if err != nil {
			continue
		}
		if _, err := os.Stat(filepath.Join(abs, "index.html")); err == nil {
			if _, err := os.Stat(filepath.Join(abs, "reports")); err == nil {
				return abs, nil
			}
		}
	}
	return "", errors.New("cannot find the site folder (index.html + reports/). " +
		"Put TTN-Site.exe inside the downloaded site folder, or pass -site <folder>")
}

func listen(first int) (net.Listener, int, error) {
	for port := first; port < first+20; port++ {
		ln, err := net.Listen("tcp", fmt.Sprintf("127.0.0.1:%d", port))
		if err == nil {
			return ln, port, nil
		}
	}
	return nil, 0, fmt.Errorf("no free local port between %d and %d", first, first+19)
}

func openBrowser(url string) {
	var cmd *exec.Cmd
	switch runtime.GOOS {
	case "windows":
		cmd = exec.Command("rundll32", "url.dll,FileProtocolHandler", url)
	case "darwin":
		cmd = exec.Command("open", url)
	default:
		cmd = exec.Command("xdg-open", url)
	}
	_ = cmd.Start()
}

// fail 印出錯誤並等使用者按 Enter，避免雙擊執行時視窗一閃就關掉。
func fail(err error) {
	fmt.Fprintln(os.Stderr, "Error:", err)
	fmt.Fprintln(os.Stderr, "Press Enter to close.")
	fmt.Scanln()
	os.Exit(1)
}
