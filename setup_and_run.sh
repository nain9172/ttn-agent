#!/usr/bin/env bash
# TTN agent：從頭到尾的一鍵流程（在 DGX Spark / Linux + NVIDIA GPU 上執行）。
#
#   1 config    設定 NCBI email / API key（寫入 .env）與 Hugging Face 登入
#   2 env       建立 conda 環境（呼叫 setup_env.sh：environment.yml → gm4、flash-attn、檢查 GPU）
#   3 models    預先下載 Evo2-7B 與 Gemma 4 權重
#   4 data      下載 TTN 參考序列、轉錄本模型、ClinVar TTN 變異、PubMind TTN 資料；
#               檢查需授權、無法自動下載的檔案（HGMD、cardiodb）
#   5 pubtator  PubTator 有新的 TTN 文獻才重跑 harvest → ClinVar ID → 座標標註
#   6 evo2      對所有 ClinVar 變異算 Evo2 分數
#   7 agent     完整 pipeline（文獻、Gemma 擷取、資料庫比對、Track B）；LLM 只載入一次
#   8 backfill  把 Evo2 分數與最新的資料庫比對寫回報告
#   9 site      打包成網頁 site/（有 Go 的話一併編譯 Windows 版 TTN-Site.exe）
#
# 每一步都可以中斷後重跑，已完成的部分會沿用。
#
# 用法：
#   ./setup_and_run.sh                    # 全部（第一次會問 email / API key / HF token）
#   ./setup_and_run.sh --limit 5          # 先用前 5 顆 ClinVar 變異試跑整個流程
#   ./setup_and_run.sh --from evo2        # 從某一步開始
#   ./setup_and_run.sh --only site        # 只跑某一步（可重複指定）
#   ./setup_and_run.sh --yes              # 不問問題，沿用 .env 與現有登入（適合放背景跑）
#
# 其他選項：
#   --env NAME            conda 環境名稱（預設 environment.yml 裡的 gm4）
#   --include-vus         ClinVar 連 VUS / conflicting 一起跑（預設只跑 2 星 P/LP/B/LB，約 5,600 顆）
#   --all-review-status   ClinVar 不限 review status
#   --refresh-data        重新下載 ClinVar / PubMind（預設已存在就沿用）
#   --force-pubtator      不管有沒有新文獻都重跑 PubTator harvest
#   --variant ID          只跑指定的變異（可重複；例 --variant 2-178527200-T-C），取代 ClinVar 清單
#   --skip STEP           跳過某一步（可重複；例 --skip pubtator）
#   --max-wait-hours N    GPU 記憶體被別人佔用時最多等幾小時（預設 12）
#
# 全部跑完需要數天（約 5,600 顆變異：Evo2 每顆約 6 秒；Gemma 擷取視文獻數，每顆數分鐘），
# 建議在 tmux 裡執行；記錄寫在 outputs/run_logs/。
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$(pwd)

STAGES=(config env models data pubtator evo2 agent backfill site)
FROM=""; ONLY=(); SKIP=(); VARIANTS=(); YES=0; LIMIT=0; ENV_NAME=""; CLINVAR_FLAGS=(); REFRESH=0; FORCE_PT=0
MAX_WAIT_HOURS=12
while [[ $# -gt 0 ]]; do
    case "$1" in
        --from) FROM="$2"; shift 2 ;;
        --only) ONLY+=("$2"); shift 2 ;;
        --skip) SKIP+=("$2"); shift 2 ;;
        --variant) VARIANTS+=("$2"); shift 2 ;;
        --max-wait-hours) MAX_WAIT_HOURS="$2"; shift 2 ;;
        --yes|-y) YES=1; shift ;;
        --limit) LIMIT="$2"; shift 2 ;;
        --env) ENV_NAME="$2"; shift 2 ;;
        --include-vus) CLINVAR_FLAGS+=(--include-vus); shift ;;
        --all-review-status) CLINVAR_FLAGS+=(--all-review-status); shift ;;
        --refresh-data) REFRESH=1; shift ;;
        --force-pubtator) FORCE_PT=1; shift ;;
        -h|--help) sed -n '2,33p' "$0"; exit 0 ;;
        *) echo "未知選項：$1（--help 看用法）" >&2; exit 2 ;;
    esac
done
[[ -z $ENV_NAME ]] && ENV_NAME=$(awk '/^name:/{print $2; exit}' environment.yml)

mkdir -p outputs/run_logs
LOG="outputs/run_logs/setup_and_run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

say()  { printf '\n\033[1;34m[%s] %s\033[0m\n' "$(date +%H:%M:%S)" "$*"; }
info() { printf '    %s\n' "$*"; }
warn() { printf '\033[1;33m    ⚠ %s\033[0m\n' "$*"; }
die()  { printf '\033[1;31m    ✗ %s\033[0m\n' "$*"; exit 1; }

want() {  # 這一步要不要跑
    local s=$1
    [[ " ${SKIP[*]} " == *" $s "* ]] && return 1
    if ((${#ONLY[@]})); then [[ " ${ONLY[*]} " == *" $s "* ]]; return; fi
    if [[ -n $FROM ]]; then
        local hit=0
        for x in "${STAGES[@]}"; do [[ $x == "$FROM" ]] && hit=1; [[ $x == "$s" ]] && return $((1 - hit)); done
    fi
    return 0
}
for s in "${ONLY[@]}" "${SKIP[@]}" ${FROM:+"$FROM"}; do
    [[ " ${STAGES[*]} " == *" $s "* ]] || die "沒有這個步驟：$s（可用：${STAGES[*]}）"
done

# ── conda ─────────────────────────────────────────────────────────────────
CONDA_BASE=$(conda info --base 2>/dev/null || true)
if [[ -z $CONDA_BASE ]]; then
    for d in "$HOME/miniforge3" "$HOME/miniconda3" "$HOME/anaconda3"; do
        [[ -f $d/etc/profile.d/conda.sh ]] && CONDA_BASE=$d && break
    done
fi
activate_env() {
    [[ -n $CONDA_BASE ]] || die "找不到 conda。請先安裝 Miniforge：https://github.com/conda-forge/miniforge"
    # shellcheck disable=SC1091
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    conda activate "$ENV_NAME" || die "無法啟用 conda 環境 $ENV_NAME（先跑 --only env）"
}
env_exists() { [[ -n $CONDA_BASE ]] && [[ -d $CONDA_BASE/envs/$ENV_NAME ]]; }
# GB10 的 CPU / GPU 共用記憶體，檔案快取會被算成「已使用」而擋下 vLLM；每個 GPU 步驟前先釋放
release_cache() { python scripts/release_page_cache.py || true; }
# Gemma 4 要求「總記憶體 × LLM_GPU_MEMORY_UTILIZATION」都是空的；不夠就先停下來說明，
# 免得 pipeline 每顆都失敗、一直自動重啟
llm_headroom_or_die() {
    # 記憶體被別的工作（可能是其他使用者）佔住時先等，每 10 分鐘回報一次佔用者。
    # 參數：需要的 GiB（省略 = Gemma 需要的「總量 × LLM_GPU_MEMORY_UTILIZATION」）
    local NEED_GIB=${1:-}
    export NEED_GIB
    local deadline=$(( $(date +%s) + MAX_WAIT_HOURS * 3600 )) last=0
    while true; do
        release_cache >/dev/null
        if python - <<'PY'
import sys
import torch
import config
free, total = torch.cuda.mem_get_info(0)
import os
need = float(os.environ["NEED_GIB"]) * 2**30 if os.environ.get("NEED_GIB") else total * config.LLM_GPU_MEMORY_UTILIZATION
print(f"    GPU 可用 {free / 2**30:.1f} GiB，需要 {need / 2**30:.1f} GiB", end="")
sys.exit(0 if free >= need else 1)
PY
        then echo "：OK"; return 0; fi
        echo
        if (( $(date +%s) - last >= 600 )); then
            info "記憶體不足，等待其他工作結束。目前佔用 GPU 的行程：$(python -c 'import spark_job; print(spark_job.gpu_occupants() or "未知")')"
            last=$(date +%s)
        fi
        if (( $(date +%s) >= deadline )); then
            info "已等待 ${MAX_WAIT_HOURS} 小時。可以："
            info "  1. 清除整機檔案快取：sudo sh -c 'sync; echo 3 > /proc/sys/vm/drop_caches'"
            info "  2. 結束其他佔用 GPU / 記憶體的程式（nvidia-smi、free -g 查看）"
            info "  3. 在 config.py 調低 LLM_GPU_MEMORY_UTILIZATION（例如 0.75；KV cache 會變小、批次變慢）"
            die "記憶體不足，Gemma 4 無法載入"
        fi
        sleep 60
    done
}

INPUT=data/clinvar_ttn.csv
if ((LIMIT > 0)); then INPUT="data/clinvar_ttn_first${LIMIT}.csv"; fi
if ((${#VARIANTS[@]})); then
    INPUT=data/custom_variants.csv
    mkdir -p data
    { echo input; printf '%s\n' "${VARIANTS[@]}"; } > "$INPUT"
fi

say "TTN agent 一鍵流程｜環境 $ENV_NAME｜輸入 $INPUT｜記錄 $LOG"

# ── 1. 設定 ────────────────────────────────────────────────────────────────
if want config; then
    say "1/9 設定 NCBI 與 Hugging Face"
    touch .env && chmod 600 .env
    get_env() { grep -E "^$1=" .env | tail -1 | cut -d= -f2-; }
    set_env() {
        local tmp; tmp=$(mktemp)
        grep -vE "^$1=" .env > "$tmp" || true
        echo "$1=$2" >> "$tmp"; mv "$tmp" .env; chmod 600 .env
    }
    email=$(get_env PUBMED_EMAIL); key=$(get_env PUBMED_API_KEY)
    if ((YES == 0)); then
        info "NCBI 要求每個 API 請求附上聯絡 email；API key 可把速率從每秒 3 次提高到 10 次"
        info "（https://www.ncbi.nlm.nih.gov/account/settings/ → API Key Management）。直接按 Enter 保留目前的值。"
        read -r -p "    NCBI email [${email:-未設定}]: " v; [[ -n $v ]] && email=$v
        read -r -s -p "    NCBI API key [${key:+已設定（隱藏）}${key:-未設定}]: " v; echo; [[ -n $v ]] && key=$v
    fi
    [[ -n $email ]] || die "必須設定 NCBI email（PUBMED_EMAIL）"
    set_env PUBMED_EMAIL "$email"
    [[ -n $key ]] && set_env PUBMED_API_KEY "$key"
    info "已寫入 .env（權限 600，不會進 git）"

    # Hugging Face：目前用到的模型都不需要授權，但登入可避免下載被限速，
    # 之後改用需要同意授權的模型（例如 gated 的 Gemma 版本）也必須登入。
    if env_exists; then
        activate_env
        if hf auth whoami >/dev/null 2>&1; then
            info "Hugging Face：已登入（$(hf auth whoami 2>/dev/null | head -1)）"
        elif ((YES == 0)); then
            info "Hugging Face token：https://huggingface.co/settings/tokens（Read 權限即可）；不登入就直接按 Enter"
            read -r -s -p "    HF token: " tok; echo
            if [[ -n $tok ]]; then
                hf auth login --token "$tok" >/dev/null && info "Hugging Face 登入成功" || warn "Hugging Face 登入失敗，稍後可再跑 --only config"
            else
                info "略過 Hugging Face 登入"
            fi
        fi
    else
        info "conda 環境還沒建立；Hugging Face 登入會在建立環境後進行（--only config 再跑一次即可）"
        HF_LOGIN_PENDING=1
    fi
fi

# ── 2. 環境 ────────────────────────────────────────────────────────────────
if want env; then
    say "2/9 建立 conda 環境 $ENV_NAME"
    # 建環境（含 flash-attn、transformer-engine 處理與檢查）都在 setup_env.sh，也可以單獨執行它
    ./setup_env.sh --env "$ENV_NAME" || die "環境建立失敗，見上方訊息（修正後可單獨重跑 ./setup_env.sh --env $ENV_NAME）"
    activate_env
    if [[ ${HF_LOGIN_PENDING:-0} == 1 && $YES == 0 ]] && ! hf auth whoami >/dev/null 2>&1; then
        read -r -s -p "    Hugging Face token（不登入就按 Enter）: " tok; echo
        [[ -n $tok ]] && { hf auth login --token "$tok" >/dev/null && info "Hugging Face 登入成功"; }
    fi
fi

# ── 3. 模型 ────────────────────────────────────────────────────────────────
if want models; then
    say "3/9 下載模型權重（Evo2-7B 約 14 GB、Gemma 4 31B NVFP4 約 31 GB；已下載的會跳過）"
    activate_env
    read -r LLM EVO2 < <(python -c "import config; print(config.LOCAL_LLM_MODEL, 'arcinstitute/' + config.EVO2_MODEL)")
    for repo in "$EVO2" "$LLM"; do
        info "hf download $repo"
        hf download "$repo" >/dev/null || die "下載 $repo 失敗（網路或 Hugging Face 權限）"
    done
fi

# ── 4. 資料 ────────────────────────────────────────────────────────────────
if want data; then
    say "4/9 下載資料"
    activate_env
    mkdir -p data/pubtator
    if [[ ! -s data/sequence.fasta ]]; then
        info "TTN 參考序列（NCBI efetch，chr2:178,525,989-178,807,423 負股）"
        curl -fsS --retry 3 -o data/sequence.fasta \
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000002.12&seq_start=178525989&seq_stop=178807423&strand=2&rettype=fasta&retmode=text" \
            || die "參考序列下載失敗"
    fi
    n=$(grep -v '^>' data/sequence.fasta | tr -d '\n' | wc -c)
    [[ $n -eq 281435 ]] || die "參考序列長度 $n，應為 281,435"
    info "參考序列 OK（281,435 bp）"

    if [[ ! -s data/ttn_transcript_model.json ]]; then
        info "建立 TTN RefSeq 轉錄本模型（會下載約 110 MB 的 chr2 註解）"
        python build_ttn_transcript_model.py || die "轉錄本模型建立失敗"
    fi
    info "轉錄本模型 OK"

    if [[ ! -s data/clinvar_ttn.csv || $REFRESH == 1 ]]; then
        python scripts/data_prep/download_clinvar_ttn.py --out data/clinvar_ttn.csv "${CLINVAR_FLAGS[@]}" \
            || die "ClinVar 下載失敗"
    else
        info "ClinVar：沿用 data/clinvar_ttn.csv（$(($(wc -l < data/clinvar_ttn.csv) - 1)) 顆，--refresh-data 可重新下載）"
    fi
    if ((LIMIT > 0)); then
        head -n $((LIMIT + 1)) data/clinvar_ttn.csv > "$INPUT"
        info "試跑：只用前 $LIMIT 顆 → $INPUT"
    fi

    # PubMind：完整版（含逐篇疾病與致病性）需向 PubMind 申請；這裡下載公開版。
    # 已經是完整版（有 disease 欄位）就絕不覆蓋；公開版在 --refresh-data 時更新。
    pm_full=$( [[ -s data/pubmind_ttn.json ]] && python -c "import json; print(int('disease' in json.load(open('data/pubmind_ttn.json'))[0]))" 2>/dev/null || echo 0)
    if [[ $pm_full == 1 ]]; then
        info "PubMind：data/pubmind_ttn.json 是完整版，保留不動"
    elif [[ ! -s data/pubmind_ttn.json || $REFRESH == 1 ]]; then
        info "PubMind TTN 資料（公開版 /download_json）"
        curl -fsS --retry 3 -o data/pubmind_ttn.json.part \
            "https://pubmind.wglab.org/download_json?query=TTN&field=gene" \
            && python -c "import json,sys; d=json.load(open('data/pubmind_ttn.json.part')); assert isinstance(d,list) and d" \
            && mv data/pubmind_ttn.json.part data/pubmind_ttn.json \
            || warn "PubMind 下載失敗，報告與網頁會略過 PubMind 比對"
    fi
    [[ -s data/pubmind_ttn.json ]] && info "PubMind OK（$(python -c "import json; print(len(json.load(open('data/pubmind_ttn.json'))))") 筆）"

    # 需要授權、無法自動下載的檔案
    for f in "data/HGMD_latest.csv:HGMD 完整匯出（資料庫比對與 Track B 類比庫）" \
             "data/HGMD_TTN_processed.xlsx:HGMD 已處理表（提供 REF/ALT，讓 HGMD 能精確比對 allele）" \
             "data/cardiodb.xlsx:cardiodb 外顯子表（Track B 的外顯子 PSI）"; do
        path=${f%%:*}; what=${f#*:}
        if [[ -s $path ]]; then info "$path OK"
        else warn "缺 $path — $what。需自行取得授權後放入；沒有的話報告仍會產生，但缺這部分。"; fi
    done
fi

# ── 5. PubTator ────────────────────────────────────────────────────────────
if want pubtator; then
    say "5/9 PubTator（有新的 TTN 文獻才重跑）"
    activate_env
    PT=data/pubtator/pubtator_ttn_variants_clinvar_annotated.csv
    run_pt=0
    if [[ ! -s $PT ]]; then
        info "還沒有 $PT，需要完整 harvest（第一次會跑很久：逐篇抓 PubTator 標註並用 Gemma 擷取表型）"
        run_pt=1
    elif ((FORCE_PT)); then
        run_pt=1
    else
        python scripts/data_prep/pubtator_check_new.py; code=$?
        if [[ $code -eq 0 ]]; then run_pt=1
        elif [[ $code -ne 10 ]]; then warn "無法查詢 PubTator（網路？），沿用現有資料"; fi
    fi
    if ((run_pt)); then
        # harvest 逐篇快取、LLM 結果有檢查點，重跑只會處理新文獻
        llm_headroom_or_die
        python scripts/data_prep/pubtator_harvest.py --max-articles 100000 --out data/pubtator/pubtator_ttn_variants_all.csv \
            && python scripts/data_prep/clinvar_variation_id.py --pubtator \
            && python scripts/data_prep/clinvar_annotate.py \
            && python scripts/data_prep/pubtator_check_new.py --record \
            || die "PubTator 更新失敗（可重跑 --only pubtator 續跑）"
        info "PubTator 資料已更新：$PT"
    fi
fi

# ── 6. Evo2 ───────────────────────────────────────────────────────────────
if want evo2; then
    say "6/9 Evo2 計分（每顆約 6 秒；已算過或報告已有分數的會跳過）"
    activate_env
    [[ -s $INPUT ]] || die "找不到 $INPUT，先跑 --only data"
    out="outputs/evo2_scores_$(basename "${INPUT%.csv}").csv"
    llm_headroom_or_die 40   # Evo2-7B 載入時約需 30 GiB
    python score_evo2_for_variants.py --input "$INPUT" --out "$out" || die "Evo2 計分失敗（重跑會從斷點續算）"
fi

# ── 7. 完整 pipeline ───────────────────────────────────────────────────────
if want agent; then
    say "7/9 完整 pipeline（Gemma 4 只載入一次；中斷後重跑會續跑）"
    activate_env
    [[ -s $INPUT ]] || die "找不到 $INPUT，先跑 --only data"
    # Evo2 已在上一步獨立算好（第 8 步寫回報告），這裡跳過以免每顆重複載入模型
    llm_headroom_or_die
    ./run_batch_inprocess.sh --input "$INPUT" --skip-evo2 --require-llm --restart-every 200 \
        || die "pipeline 中止（重跑 --from agent 會續跑）"
fi

# ── 8. 寫回 ────────────────────────────────────────────────────────────────
if want backfill; then
    say "8/9 寫回 Evo2 分數與最新的資料庫比對"
    activate_env
    python backfill_evo2_scores.py || die "Evo2 寫回失敗"
    python backfill_database_annotations.py || die "資料庫比對寫回失敗"
fi

# ── 9. 網頁 ────────────────────────────────────────────────────────────────
if want site; then
    say "9/9 打包網頁"
    activate_env
    if [[ ${TTN_SERVER_MODE:-} == 1 ]]; then
        :  # docker 服務：網站由 web 容器提供，不需要 Windows 啟動程式
    elif command -v go >/dev/null 2>&1 || [[ -x $HOME/.local/go/bin/go ]]; then
        PATH="$HOME/.local/go/bin:$PATH" launcher/build.sh >/dev/null && info "已編譯 TTN-Site.exe"
    else
        info "沒有 Go，略過 Windows 版 TTN-Site.exe（安裝 https://go.dev/dl/ 後重跑 --only site）"
    fi
    python build_site.py || die "網頁打包失敗"
    if [[ ${TTN_SERVER_MODE:-} == 1 ]]; then info "完成：網站已更新（重新整理瀏覽器即可）"
    else info "完成：$ROOT/site/index.html（Windows 使用者雙擊 site/TTN-Site.exe）"; fi
fi

say "全部完成。記錄：$LOG"
