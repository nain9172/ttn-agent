#!/usr/bin/env bash
# 初始化：預先計算整條 TTN 上所有 SNV 的 Evo2 分數（site 的「未收錄變異」查表用）。
#
# 第一次 clone 後執行一次即可。計算在背景進行，GB10 上約 31 天；
# 會先算外顯子 + 剪接區（約 14 天），可隨時停止、重跑本腳本就會續算。
#
# 用法（先 conda activate gm4）：
#   ./init_evo2_table.sh               準備資料並在背景開始／續算
#   ./init_evo2_table.sh --foreground  在前景執行（Ctrl-C 會存檔後離開）
#   ./init_evo2_table.sh --status      查看進度
#   ./init_evo2_table.sh --stop        請背景計算在目前 batch 結束後存檔停止
#   ./init_evo2_table.sh --export      把目前結果匯出成 site 用的 evo2_table.js
#
# 環境變數：PYTHON（預設 python）、OUT_DIR（預設 outputs/evo2_table）
set -euo pipefail
cd "$(dirname "$0")"

PYTHON=${PYTHON:-python}
OUT_DIR=${OUT_DIR:-outputs/evo2_table}
FASTA=data/sequence.fasta
TX_MODEL=data/ttn_transcript_model.json
PID_FILE="$OUT_DIR/run.pid"
TTN_LENGTH=281435
# GRCh38 chr2:178,525,989-178,807,423，負股（strand=2），與 config.TTN_SEQUENCE_START/END 一致
FASTA_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000002.12&seq_start=178525989&seq_stop=178807423&strand=2&rettype=fasta&retmode=text"

say() { echo "[init] $*"; }
die() { echo "[init] 錯誤：$*" >&2; exit 1; }

running_pid() {
    if [[ -f "$PID_FILE" ]]; then
        local pid
        pid=$(cat "$PID_FILE")
        if kill -0 "$pid" 2>/dev/null; then
            echo "$pid"
            return 0
        fi
    fi
    return 1
}

precompute() { "$PYTHON" precompute_evo2_table.py --out-dir "$OUT_DIR" "$@"; }

mode=${1:-start}
case "$mode" in
    --status)
        if pid=$(running_pid); then say "背景計算執行中（PID $pid）"; else say "目前沒有背景計算"; fi
        precompute --status
        [[ -f "$OUT_DIR/progress.log" ]] && { echo; tail -n 3 "$OUT_DIR/progress.log"; }
        exit 0 ;;
    --stop)
        if pid=$(running_pid); then
            touch "$OUT_DIR/STOP"
            say "已通知 PID $pid 停止；目前 batch 結束後會存檔離開（通常一分鐘內）"
        else
            say "目前沒有背景計算"
        fi
        exit 0 ;;
    --export)
        precompute --export
        exit 0 ;;
    start|--foreground) ;;
    *) die "未知參數 $mode（可用：--foreground --status --stop --export）" ;;
esac

# ── 1. 環境檢查 ───────────────────────────────────────────────────────────────
say "檢查 Python 環境（$PYTHON）"
"$PYTHON" - <<'EOF' || die "環境不完整，請先 conda activate gm4（需要 evo2、torch 與可用的 CUDA GPU）"
import sys
import evo2, numpy, torch
if not torch.cuda.is_available():
    sys.exit("找不到 CUDA GPU")
print(f"[init]   evo2 OK，GPU：{torch.cuda.get_device_name(0)}")
EOF

if pid=$(running_pid); then
    die "已經有計算在跑（PID $pid）；用 --status 看進度，或 --stop 停止"
fi

# ── 2. TTN 參考序列（data/ 不在 git 裡，第一次要下載）───────────────────────
mkdir -p data
if [[ ! -s "$FASTA" ]]; then
    say "下載 TTN 參考序列（NCBI efetch，約 290 KB）"
    if command -v curl >/dev/null 2>&1; then
        curl -fsSL "$FASTA_URL" -o "$FASTA.tmp"
    else
        wget -qO "$FASTA.tmp" "$FASTA_URL"
    fi
    mv "$FASTA.tmp" "$FASTA"
fi
"$PYTHON" - "$FASTA" "$TTN_LENGTH" <<'EOF' || die "$FASTA 內容不符（應為 GRCh38 chr2:178525989-178807423 負股）；刪掉後重跑會重新下載"
import sys
path, expected = sys.argv[1], int(sys.argv[2])
lines = open(path).read().splitlines()
header = lines[0] if lines else ""
seq = "".join(l.strip() for l in lines[1:])
assert header.startswith(">NC_000002.12:c178807423-178525989"), header
assert len(seq) == expected, f"length {len(seq)} != {expected}"
assert set(seq.upper()) <= set("ACGTN"), "unexpected characters"
print(f"[init]   參考序列 OK（{len(seq):,} bp）")
EOF

# ── 3. 轉錄本模型（只用來決定先算外顯子；失敗也能照序列順序算）────────────
if [[ ! -s "$TX_MODEL" ]]; then
    say "建立 TTN 轉錄本模型（會從 NCBI 下載 chr2 GFF3，約 110 MB）"
    if ! "$PYTHON" build_ttn_transcript_model.py; then
        say "警告：轉錄本模型建立失敗，改為照序列順序計算（結果相同，只是外顯子不會先算完）"
    fi
fi

# ── 4. 開始／續算 ─────────────────────────────────────────────────────────────
mkdir -p "$OUT_DIR"
rm -f "$OUT_DIR/STOP"

if [[ "$mode" == "--foreground" ]]; then
    precompute
    exit $?
fi

nohup "$PYTHON" precompute_evo2_table.py --out-dir "$OUT_DIR" >> "$OUT_DIR/run.out" 2>&1 &
echo $! > "$PID_FILE"
sleep 2
kill -0 "$(cat "$PID_FILE")" 2>/dev/null || die "背景程序啟動失敗，請看 $OUT_DIR/run.out"

say "已在背景開始計算（PID $(cat "$PID_FILE")），完整約需 31 天（GB10 實測）"
say "  進度：./init_evo2_table.sh --status"
say "  停止：./init_evo2_table.sh --stop   （之後重跑本腳本即續算）"
say "  記錄：$OUT_DIR/progress.log"
say "每次停止或完成時會自動更新 $OUT_DIR/evo2_table.js；計算期間 GPU 會被佔用"
