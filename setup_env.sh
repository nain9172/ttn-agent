#!/usr/bin/env bash
# 只建立（或修補）TTN agent 的 conda 環境，不下載資料、不跑 pipeline。
#
#   1. conda env create -f environment.yml（環境已存在就沿用；缺套件時用 conda env update 補齊）
#   2. 安裝 flash-attn（Evo2 需要；PyPI 沒有 aarch64 版），依序嘗試：
#        FLASH_ATTN_WHEEL（檔案路徑或網址）→ wheels/ 裡 Python 版本相符的 wheel → 從原始碼編譯
#   3. 移除 evo2 依賴帶進來、只有空殼的 transformer-engine（有它 Evo2 反而無法載入）
#   4. 檢查：主要套件都在、PyTorch 看得到 GPU、evo2 載入得了
#
# 可以重複執行：已經裝好的部分會跳過，只補缺的。setup_and_run.sh 的 env 步驟就是呼叫這支。
#
# 用法：
#   ./setup_env.sh                       # 環境名稱用 environment.yml 裡的（gm4）
#   ./setup_env.sh --env gm4test         # 指定名稱
#   FLASH_ATTN_WHEEL=/path/or/url.whl ./setup_env.sh
#   MAX_JOBS=2 ./setup_env.sh            # 必須從原始碼編譯 flash-attn 時限制平行數
set -uo pipefail
cd "$(dirname "$0")"

ENV_NAME=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --env) ENV_NAME="$2"; shift 2 ;;
        -h|--help) sed -n '2,17p' "$0"; exit 0 ;;
        *) echo "未知選項：$1（--help 看用法）" >&2; exit 2 ;;
    esac
done
[[ -z $ENV_NAME ]] && ENV_NAME=$(awk '/^name:/{print $2; exit}' environment.yml)

say()  { printf '\n\033[1;34m[%s] %s\033[0m\n' "$(date +%H:%M:%S)" "$*"; }
info() { printf '    %s\n' "$*"; }
warn() { printf '\033[1;33m    ⚠ %s\033[0m\n' "$*"; }
die()  { printf '\033[1;31m    ✗ %s\033[0m\n' "$*"; exit 1; }

CONDA_BASE=$(conda info --base 2>/dev/null || true)
if [[ -z $CONDA_BASE ]]; then
    for d in "$HOME/miniforge3" "$HOME/miniconda3" "$HOME/anaconda3"; do
        [[ -f $d/etc/profile.d/conda.sh ]] && CONDA_BASE=$d && break
    done
fi
[[ -n $CONDA_BASE ]] || die "找不到 conda。請先安裝 Miniforge：https://github.com/conda-forge/miniforge"
# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"
activate_env() { conda activate "$ENV_NAME" || die "無法啟用 conda 環境 $ENV_NAME"; }

say "建立 conda 環境 $ENV_NAME"

# ── 1. conda 環境 ──────────────────────────────────────────────────────────
if [[ -d $CONDA_BASE/envs/$ENV_NAME ]]; then
    info "環境已存在"
else
    info "conda env create -f environment.yml（約 20–40 分鐘）"
    conda env create -n "$ENV_NAME" -f environment.yml \
        || warn "conda env create 沒有完全成功，下面改用 conda env update 補齊"
fi
activate_env
# 上次建到一半（例如 pip 某個套件失敗）時主要套件會缺，用 update 補齊。
# 只看套件在不在、不 import：evo2 一 import 就會去載 flash-attn，它還沒裝時會誤判
if ! python -c "import importlib.util as u, sys; sys.exit(any(u.find_spec(m) is None for m in ('torch', 'vllm', 'evo2', 'pandas', 'bs4')))" 2>/dev/null; then
    info "環境不完整，conda env update -f environment.yml 補齊缺少的套件"
    conda env update -n "$ENV_NAME" -f environment.yml || die "環境補齊失敗，見上方訊息"
    activate_env
fi

# ── 2. flash-attn ─────────────────────────────────────────────────────────
if python -c "import flash_attn" 2>/dev/null; then
    info "flash-attn 已安裝（$(python -c 'import flash_attn; print(flash_attn.__version__)')）"
else
    pytag=$(python -c 'import sys; print(f"cp{sys.version_info.major}{sys.version_info.minor}")')
    wheel=${FLASH_ATTN_WHEEL:-}
    if [[ -z $wheel ]]; then
        wheel=$(ls wheels/flash_attn-2.8.3-"$pytag"-"$pytag"-linux_"$(uname -m)".whl 2>/dev/null | head -1)
    fi
    if [[ -n $wheel ]]; then
        info "安裝預先編好的 flash-attn：$wheel"
        pip install "$wheel" || warn "這個 wheel 裝不起來，改從原始碼編譯"
    fi
fi
if ! python -c "import flash_attn" 2>/dev/null; then
    info "沒有相符的預編 wheel（需要 ${pytag:-cp3xx} / $(uname -m)），從原始碼編譯 flash-attn 2.8.3（GB10 上要一小時以上）"
    [[ -z ${CUDA_HOME:-} && -x /usr/local/cuda/bin/nvcc ]] && export CUDA_HOME=/usr/local/cuda
    # 每個編譯工作（內含多個 nvcc 執行緒）實測超過 16 GB：118 GB 可用、7 個平行仍被 OOM killer 砍掉。
    # 依目前可用記憶體決定平行數（flash-attn 官方建議記憶體少於 96 GB 時設 MAX_JOBS=4）
    if [[ -z ${MAX_JOBS:-} ]]; then
        avail_gb=$(awk '/MemAvailable/ {print int($2 / 1048576)}' /proc/meminfo)
        MAX_JOBS=$(( avail_gb / 32 )); (( MAX_JOBS < 1 )) && MAX_JOBS=1; (( MAX_JOBS > $(nproc) )) && MAX_JOBS=$(nproc)
        info "可用記憶體 ${avail_gb} GiB → 同時編譯 ${MAX_JOBS} 個（可用 MAX_JOBS 覆寫）"
    fi
    MAX_JOBS=$MAX_JOBS pip install flash-attn==2.8.3 --no-build-isolation \
        || die "flash-attn 編譯失敗。log 裡有 Killed 代表記憶體不足（等其他工作結束或設小 MAX_JOBS 後重跑）；找不到 nvcc 則設 CUDA_HOME 指向 CUDA 目錄"
fi

# ── 3. transformer-engine 空殼 ────────────────────────────────────────────
# evo2 的依賴會裝進只有空殼的 transformer-engine（沒有 PyTorch 擴充），Evo2 一偵測到就報錯。
# gm4 不裝 transformer-engine（Evo2 走非 FP8 路徑，分數與 site 使用的 gm4 一致），這裡比照移除。
if pip show transformer_engine >/dev/null 2>&1 && ! python -c "import transformer_engine.pytorch" 2>/dev/null; then
    info "移除不完整的 transformer-engine（Evo2 不需要；有它反而無法載入）"
    pip uninstall -y transformer_engine transformer-engine-cu12 transformer-engine-cu13 >/dev/null 2>&1 || true
fi

# ── 4. 檢查 ────────────────────────────────────────────────────────────────
activate_env
python - <<'PY' || die "環境檢查失敗"
import importlib.util, sys
missing = [m for m in ("torch", "vllm", "evo2", "flash_attn", "pandas", "bs4") if importlib.util.find_spec(m) is None]
if missing:
    sys.exit(f"缺少套件：{missing}")
import torch
if not torch.cuda.is_available():
    sys.exit("PyTorch 看不到 GPU")
try:
    import evo2  # noqa: F401  連 flash-attn / transformer-engine 的問題一起在這裡抓出來
except Exception as e:  # noqa: BLE001
    sys.exit(f"evo2 無法載入：{e}")
print(f"    GPU：{torch.cuda.get_device_name(0)}，torch {torch.__version__}")
PY
info "環境 $ENV_NAME 可以使用：conda activate $ENV_NAME"
