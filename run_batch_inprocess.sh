#!/usr/bin/env bash
# 以 in-process 模式跑 batch_run_variants.py，遇到 exit 3（vLLM 引擎壞掉／連續失敗）
# 或 exit 4（--restart-every 定期重啟）時自動重啟續跑。
#
# 用法（先 conda activate gm4）：
#   ./run_batch_inprocess.sh --input data/clinvar_ttn.csv --skip-evo2 --require-llm [--force] [--restart-every 200]
#
# 環境變數：MAX_RESTARTS（預設 100）、RESTART_DELAY 秒（預設 30）
set -u
cd "$(dirname "$0")"

if ! command -v ninja >/dev/null 2>&1; then
    echo "[wrapper] 找不到 ninja：請先 conda activate gm4，否則 vLLM 無法編譯 kernel" >&2
    exit 1
fi

args=("$@")
max_restarts=${MAX_RESTARTS:-100}
delay=${RESTART_DELAY:-30}

for ((i = 0; i <= max_restarts; i++)); do
    python -u batch_run_variants.py --in-process "${args[@]}"
    code=$?
    if [[ $code -ne 3 && $code -ne 4 ]]; then
        echo "[wrapper] batch_run_variants.py 結束，exit code $code"
        exit $code
    fi

    # 續跑不能再帶 --force（會從原始表格重來、覆蓋已完成的列）。
    # 原本有 --force 代表不要沿用舊報告，所以改加 --no-reuse。
    resumed=()
    had_force=0
    for a in "${args[@]}"; do
        if [[ $a == --force ]]; then had_force=1; else resumed+=("$a"); fi
    done
    ((had_force)) && resumed+=(--no-reuse)
    args=("${resumed[@]}")

    echo "[wrapper] exit $code，${delay}s 後續跑（第 $((i + 1)) 次重啟）：${args[*]}"
    sleep "$delay"
done

echo "[wrapper] 已重啟 $max_restarts 次，停止" >&2
exit 1
