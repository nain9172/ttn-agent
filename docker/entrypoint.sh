#!/usr/bin/env bash
# docker 容器的進入點：把 volume 接到程式預期的路徑，再執行指定的服務。
#
#   web        網站與 API（web_server.py），不需要 GPU
#   worker     GPU 分析佇列（spark_job.py worker --forever）
#   setup ...  直接執行 setup_and_run.sh（例：setup --only data）
#   status     列出 volume 裡的資料與模型狀態
#   其他       當成指令執行（例：bash）
set -euo pipefail

VOL=${TTN_VOLUME:-/vol}
if [[ ! -w $VOL ]]; then
    echo "✗ 資料目錄 $VOL 沒有寫入權限（容器以 uid $(id -u) 執行）。" >&2
    echo "  使用主機資料夾（TTN_DATA_DIR）時，請先用執行帳號建立它：mkdir -p <資料夾>，" >&2
    echo "  並在 .env 設定 TTN_UID=\$(id -u)、TTN_GID=\$(id -g)。見 docker_readme/DEPLOY_DOCKER.md。" >&2
    exit 1
fi
mkdir -p "$VOL"/data "$VOL"/outputs "$VOL"/site "$VOL"/models "$VOL"/cache "$VOL"/llm_logs

# /app/data、/app/outputs、/app/site、/app/llm_logs 在映像裡已是指向 $VOL 的連結（見 Dockerfile）
# 映像內附、可公開的資料（目前只有 transcript_interval.xlsx）：資料目錄沒有才複製
if [[ -d /app/seed/data ]]; then cp -r --update=none --no-preserve=ownership /app/seed/data/. "$VOL/data/"; fi

# 以一般帳號（TTN_UID）執行時，容器的 /etc/passwd 裡沒有這個 uid：
#   - 查不到使用者名稱，PyTorch 載入時（getpass.getuser）會失敗，vLLM 因此無法載入 → 補上 USER / LOGNAME
#   - 沒有可寫的家目錄 → 各種工具的設定與快取都放在資料目錄
if ! id -un >/dev/null 2>&1; then
    export USER=ttn LOGNAME=ttn
fi
if [[ ! -w ${HOME:-/} ]]; then
    export HOME=$VOL/cache/home
    mkdir -p "$HOME"
fi

# 模型權重與各種編譯快取也放在 volume，重啟容器不必重新下載或編譯
export HF_HOME=${HF_HOME:-$VOL/models/huggingface}
export XDG_CACHE_HOME=${XDG_CACHE_HOME:-$VOL/cache}
export TRITON_CACHE_DIR=${TRITON_CACHE_DIR:-$VOL/cache/triton}
export VLLM_CACHE_ROOT=${VLLM_CACHE_ROOT:-$VOL/cache/vllm}
export TTN_SITE_DIR=${TTN_SITE_DIR:-$VOL/site}
export TTN_SERVER_MODE=1

cd /app
case "${1:-web}" in
    web)    exec python web_server.py ;;
    worker)
        # 先用實際執行的帳號載入一次 pipeline 需要的套件；失敗就停下來並說明，
        # 免得每個使用者送出的分析都失敗（網頁會顯示 worker 沒在執行）
        if ! python -c "import cv2, torch, transformers, evo2, vllm; from vllm import LLM" >/tmp/import-check.log 2>&1; then
            echo "✗ worker 無法載入 pipeline 需要的套件（uid $(id -u)）：" >&2
            tail -5 /tmp/import-check.log >&2
            exit 1
        fi
        exec python spark_job.py worker --forever ;;
    setup)  shift; exec bash setup_and_run.sh --yes "$@" ;;
    status) exec python deploy_status.py ;;
    *)      exec "$@" ;;
esac
