#!/usr/bin/env bash
# 重現 Evo2 分數：git clone 之後執行 ./reproduce_evo2.sh 就好。
# 適用 DGX Spark（aarch64、CUDA 13）；GPU 需要約 30 GiB 可用記憶體。
#
# 步驟（每一步都可以重跑，已完成的會跳過；計分中斷後重跑會從斷點續算）：
#   1 建 conda 環境 gm4（setup_env.sh，約 20–40 分鐘；flash-attn 用 repo 裡 wheels/ 的預編 wheel）
#   2 下載 TTN 參考序列（NCBI）
#   3 下載 ClinVar 的 TTN 變異 → data/clinvar_ttn.csv（約 5 分鐘），條件：
#       Germline classification：Benign、Likely benign、Likely pathogenic、Pathogenic
#       Variation type：Single nucleotide；Review status：multiple submitters, no conflicts
#     ClinVar 持續更新，每次下載的數量會不同；已經有這個檔案就沿用，要重抓先刪掉它
#   4 下載 Evo2-7B 權重（約 14 GB）
#   5 對全部 ClinVar 變異計分 → outputs/evo2_scores_clinvar_ttn.csv（每顆約 6 秒，約 5,700 顆共 9–10 小時，
#     建議在 tmux 裡跑）
#   6 用 ClinVar 標籤評估：AUROC，以及 cutoff（delta < -0.000714 判為 pathogenic）的敏感度、特異度
set -eo pipefail
cd "$(dirname "$0")"

ENV_NAME=gm4
CLINVAR=data/clinvar_ttn.csv
OUT=outputs/evo2_scores_clinvar_ttn.csv
say() { printf '\n\033[1;34m[%s] %s\033[0m\n' "$(date +%H:%M:%S)" "$*"; }

say "1/6 conda 環境 $ENV_NAME"
./setup_env.sh --env "$ENV_NAME"
CONDA_BASE=$(conda info --base 2>/dev/null || true)
for d in "$CONDA_BASE" "$HOME/miniforge3" "$HOME/miniconda3" "$HOME/anaconda3"; do
    [[ -n $d && -f $d/etc/profile.d/conda.sh ]] && { source "$d/etc/profile.d/conda.sh"; break; }
done
conda activate "$ENV_NAME"

say "2/6 TTN 參考序列（chr2:178,525,989-178,807,423 負股）"
mkdir -p data outputs
if [[ ! -s data/sequence.fasta ]]; then
    curl -fsS --retry 3 -o data/sequence.fasta \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000002.12&seq_start=178525989&seq_stop=178807423&strand=2&rettype=fasta&retmode=text"
fi
n=$(grep -v '^>' data/sequence.fasta | tr -d '\n' | wc -c)
[[ $n -eq 281435 ]] || { echo "參考序列長度 $n，應為 281,435。刪掉 data/sequence.fasta 後重跑"; exit 1; }
echo "    OK（281,435 bp）"

say "3/6 ClinVar TTN 變異"
if [[ ! -s $CLINVAR ]]; then
    python scripts/data_prep/download_clinvar_ttn.py --out "$CLINVAR"
else
    echo "    沿用 $CLINVAR（$(($(wc -l < "$CLINVAR") - 1)) 顆；要重新下載先刪掉它）"
fi

say "4/6 Evo2-7B 權重（第一次約 14 GB）"
hf download arcinstitute/evo2_7b_base >/dev/null
echo "    OK"

say "5/6 Evo2 計分（每顆約 6 秒，中斷後重跑會續算）"
python score_evo2_for_variants.py --input "$CLINVAR" --out "$OUT" --all

say "6/6 以 ClinVar 標籤評估"
python - "$CLINVAR" "$OUT" <<'PY'
import sys
import pandas as pd
from sklearn.metrics import roc_auc_score
from config import PATHOGENIC_THRESHOLD as THR

cv = pd.read_csv(sys.argv[1])
cv["vid"] = "2-" + cv["position"].astype(str) + "-" + cv["tx_ref"] + "-" + cv["tx_alt"]
df = pd.read_csv(sys.argv[2]).merge(cv[["vid", "clinvar_label"]], on="vid")
df = df[df.clinvar_label.isin(["pathogenic", "benign"])]
y = df.clinvar_label.eq("pathogenic")
call = df.delta_score < THR
print(f"    變異數：{len(df):,}（pathogenic {int(y.sum()):,}、benign {int((~y).sum()):,}）")
print(f"    AUROC：{roc_auc_score(y, -df.delta_score):.4f}")
print(f"    cutoff {THR}：敏感度 {(call & y).sum() / y.sum():.3f}、特異度 {(~call & ~y).sum() / (~y).sum():.3f}、"
      f"準確率 {(call == y).mean():.3f}")
PY
say "完成：分數在 $OUT"
