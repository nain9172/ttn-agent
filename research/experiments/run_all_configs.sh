#!/usr/bin/env bash
# 依序跑完實驗二／三的所有設定（每個設定一個獨立行程，各自完整載入一次模型）。
#
# 用法（先 conda activate gm4）：
#   ./research/experiments/run_all_configs.sh              # 等 GPU 空出來後依序跑完
#   PILOT=3 ./research/experiments/run_all_configs.sh      # 每個設定只跑 3 顆（先估時間用）
#   WAIT=0  ./research/experiments/run_all_configs.sh      # 不等 GPU，立刻開始
#
# 每個設定的輸出是 research/experiments/runs/<label>.jsonl，逐顆變異 append，中斷可續跑，
# 重跑這個腳本會從斷點接續（已完成的變異會被跳過）。
set -u
cd "$(dirname "$0")/../.."

# SET 決定跑哪一輪的評估集：空 = 第一輪（隨機 60 顆），hub = 第二輪（hub 文獻的 57 顆）
SET=${SET:-}
if [[ -n $SET ]]; then
    export TTN_EXP_SET="$SET"
    VARIANTS=research/experiments/llm_eval_set_${SET}.csv
    RUNS_DIR=research/experiments/runs_${SET}
else
    VARIANTS=research/experiments/round1_random60/llm_eval_set.csv
    RUNS_DIR=research/experiments/round1_random60/runs
fi
set_args=(--variants "$VARIANTS" --out-dir "$RUNS_DIR")
# hub 系列只餵「提及最多變異的那篇」文獻，其他文獻不進 LLM；
# 評估集沒有 hub_pmids 欄（例如 Set A）時就餵完整文獻集。
if head -1 "$VARIANTS" | grep -q "hub_pmids"; then
    set_args+=(--pmid-column "${PMID_COLUMN:-hub_pmids}")
fi

PILOT=${PILOT:-0}                 # >0 時每個設定只跑前 N 顆
WAIT=${WAIT:-1}                   # 1 = 等到 GPU free 記憶體足夠再開始
GPU_UTIL=${GPU_UTIL:-0.78}        # vLLM 會佔「總記憶體 × 這個比例」；所有設定用同一個值
NEED_GIB=${NEED_GIB:-96}          # 開跑前要求的 free 記憶體（略大於 GPU_UTIL × 122 GiB）
limit_args=()
((PILOT > 0)) && limit_args=(--limit "$PILOT")

# label 與旗標。all-guards 用 config.py 的預設模型（NVFP4 的 Gemma 4 31B），
# 它同時是實驗三的基準與實驗二的「現行部署」參考點，所以只跑一次。
configs=(
    "all-guards|"
    "medgemma-27b|--model google/medgemma-27b-it"
    "gemma4-31b-bf16|--model google/gemma-4-31B-it"
    "no-alias-gate|--no-alias-gate"
    "no-evidence|--no-verify-evidence"
    "no-guided-json|--no-guided-json"
    "sc-n1|--self-consistency 1"
    "sc-n3|--self-consistency 3"
    "none|--no-alias-gate --no-guided-json --no-verify-evidence --self-consistency 1"
)

free_gib() {
    python - <<'PY'
import torch
free, total = torch.cuda.mem_get_info()
print(int(free / 1024**3))
PY
}

if ((WAIT)); then
    echo "[queue] 等 GPU 釋放（需要 ${NEED_GIB} GiB free）…"
    while :; do
        f=$(free_gib 2>/dev/null || echo 0)
        [[ ${f:-0} -ge $NEED_GIB ]] && break
        echo "[queue] $(date +%H:%M) 目前 free ${f} GiB，60s 後再看" >&2
        sleep 60
    done
fi

for entry in "${configs[@]}"; do
    label=${entry%%|*}
    flags=${entry#*|}
    echo
    echo "=========================================================="
    echo "[queue] $(date +%H:%M:%S) 開始 $label $flags"
    echo "=========================================================="
    # shellcheck disable=SC2086  # flags 需要分詞成多個參數
    python research/experiments/run_extraction_config.py --label "$label" ${flags} \
        --gpu-memory-utilization "$GPU_UTIL" "${set_args[@]}" "${limit_args[@]}"
    code=$?
    echo "[queue] $label 結束，exit $code"
    if [[ $code -ne 0 ]]; then
        echo "[queue] 非 0 結束碼，停在這裡；修好後重跑本腳本會從斷點續跑" >&2
        exit $code
    fi
    sleep 20     # 讓上一個行程的顯存完全歸還
done

# 實驗二的 NVFP4 參考點與實驗三的基準是同一份結果，做個連結讓兩邊的 notebook 都讀得到
ln -sf all-guards.jsonl "$RUNS_DIR/gemma4-31b-nvfp4.jsonl"

# 抽取跑完後直接把兩個 notebook 執行一遍，圖表與摘要表就都生好了（這段不需要 GPU）。
# SKIP_ANALYSIS=1 可以只跑抽取、不碰 notebook。
if [[ ${SKIP_ANALYSIS:-0} -eq 0 && $PILOT -eq 0 ]]; then
    echo
    echo "[queue] $(date +%H:%M:%S) 執行 notebook 分析…"
    mkdir -p research/experiments/executed
    round_tag=${SET:-round1_random60}
    for nb in research/exp2_medgemma_vs_gemma4.ipynb research/exp3_hallucination_guards.ipynb; do
        # 執行後的副本另存，原始 notebook 保持乾淨、兩輪結果都留著
        out="$(basename "${nb%.ipynb}")_${round_tag}.ipynb"
        echo "[queue] nbconvert $nb -> research/experiments/executed/$out"
        jupyter nbconvert --to notebook --execute \
            --output-dir research/experiments/executed --output "$out" \
            --ExecutePreprocessor.timeout=1800 "$nb" \
            || echo "[queue] $nb 執行失敗，手動開 notebook 逐格跑" >&2
    done
fi

echo
echo "[queue] 全部完成。輸出："
echo "  抽取結果  research/experiments/runs/*.jsonl"
echo "  摘要表    research/experiments/exp2_model_summary.csv, research/experiments/exp3_guard_summary.csv"
echo "  圖        research/experiments/figures/*.png"
echo "  notebook  exp2_medgemma_vs_gemma4.ipynb, exp3_hallucination_guards.ipynb（已含輸出）"
