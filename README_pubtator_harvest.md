# `pubtator_harvest.py` 使用說明

PubTator3 → TTN 致病變異彙整器（獨立測試版）。

與主 pipeline（variant-centric：先有變異 → 找文獻）相反，本程式是 **literature-centric**：

> 從 PubTator 上 TTN 基因的 journal article 文獻出發 → 抽出文中提及的 TTN 變異
> → 用本地 LLM 判定致病性並擷取 disease / 影響組織 / 發病年齡 / 遺傳模式
> → 把所有文獻的結果彙整成單一 CSV。

---

## 1. 環境需求

- **conda 環境 `gm4`**（已裝好 `vllm`、`transformers`、`requests`、`pandas` 等）。
- **GPU**：完整流程（含 LLM）需要 GPU；只抓變異（`--no-llm`）則不需要。
- **`config.py`**（與主 pipeline 共用）提供以下關鍵設定：

| 設定 | 說明 | 目前值 |
|---|---|---|
| `LOCAL_LLM_MODEL` | 本地 LLM 模型 | `google/gemma-4-31B-it`（BF16, ~58GB） |
| `LOCAL_LLM_MAX_MODEL_LEN` | vLLM 單 prompt 上限 (tokens) | `32768` |
| `LOCAL_LLM_MAX_CONTEXT_LENGTH` | 餵給 LLM 的文本上限 (字元) | `30000` |
| `LLM_SELF_CONSISTENCY_N` | 每個 prompt sample 幾次後投票 | `5` |
| `LLM_USE_GUIDED_JSON` | 強制輸出符合 JSON schema | `True` |
| `LLM_MAX_NUM_SEQS` | vLLM 同時解碼序列上限（控 KV cache，防 OOM） | `16`（沒設則預設 16） |
| `TRITON_PTXAS_PATH` | vLLM 編譯需要 | `/usr/local/cuda/bin/ptxas` |
| `PUBMED_API_KEY` / `PUBMED_EMAIL` | NCBI esummary 用；設了 API key 可從 3 req/s 提速到 10 req/s | — |

---

## 2. 快速開始

```bash
# 只跑 Stage 1–2（不需 GPU），先驗證變異抓取是否正常
conda run -n gm4 python pubtator_harvest.py --max-articles 5 --no-llm

# 完整流程（需 GPU / vLLM）：處理前 100 篇 journal article
conda run -n gm4 python pubtator_harvest.py --max-articles 100
```

背景跑大批量並即時看輸出（log 也會自動寫檔，見第 6 節）：

```bash
PYTHONUNBUFFERED=1 conda run --no-capture-output -n gm4 \
  python -u pubtator_harvest.py --max-articles 2000 \
  --out outputs/pubtator_ttn_variants.csv
```

---

## 3. 命令列參數

| 參數 | 預設 | 說明 |
|---|---|---|
| `--max-articles` | `100` | 要處理的 journal article 篇數。設成 `15000`（>全庫 14,199）= 全取。 |
| `--out` | `outputs/pubtator_ttn_variants.csv` | 輸出 CSV 路徑。 |
| `--no-llm` | off | 跳過 LLM，只抓變異清單（不需 GPU）。 |
| `--no-fulltext` | off | 不抓全文（full=true / EuropePMC），變異與 context 只用 abstract。 |
| `--include-reviews` | off | 一併納入 Review 類型文獻（預設只收 Journal Article）。 |
| `--refresh` | off | 忽略快取，重新下載 biocjson 全文/abstract。 |
| `--llm-model` | config 值 | 覆寫 LLM 模型。 |
| `--llm-tokenizer` | None | 覆寫 tokenizer。 |
| `--llm-batch-size` | `32` | 每批送進 vLLM 的（文獻×變異）數量。越小越省 KV cache、越大吞吐越高。 |
| `--refresh-llm` | off | 清掉 LLM 檢查點，重新推論所有變異。 |

---

## 4. 處理流程（Pipeline）

| 階段 | 動作 |
|---|---|
| **Stage 1** | 用基因實體查詢 `@GENE_TTN`（非 free-text，避免雜訊），依 relevance 取候選 PMID。緩衝量為 `max_articles × 3`。 |
| **Stage 1b** | 用 NCBI `esummary` 取 pubtype/年份，過濾出 Journal Article（PubTator 內建 filter 實測無效，故改此法）。 |
| **Stage 2** | 批次抓 `full=true` 全文 biocjson；缺全文者退回 abstract 層級；再缺則試 EuropePMC 補全文（僅供 LLM context）。解析 Variant 標註中 `CorrespondingGene:7273`（TTN）者，含 HGVS / rsID。 |
| **Stage 3+4** | 本地 LLM（vLLM）判定致病性並擷取表型（disease / 影響組織 / 發病年齡 / 遺傳模式 / 病患數 / 證據句）。**分批 + 檢查點**，可續跑。 |
| **Stage 5** | 彙整輸出單一 CSV（每列 = 一個「文獻 × 變異」對）。 |

---

## 5. 輸出 CSV 欄位

`outputs/pubtator_ttn_variants.csv`（UTF-8-SIG，每列為一組「文獻×變異」）：

`pmid, pmcid, year, journal, title, variant_text, hgvs, rsid, gene, mention_count,`
`pathogenicity, disease, tissue_affected, age_onset, inheritance, patient_count,`
`evidence_sentence, reasoning, text_source, pubmed_link`

- `pathogenicity`：Pathogenic / Likely pathogenic / Uncertain significance / Likely benign / Benign / Not specified
- `tissue_affected`：Cardiac / Skeletal / Both / Not specified
- `text_source`：`pubtator_full`（全文）/ `abstract` / `europepmc`

---

## 6. 快取、檢查點與 Log（可續跑）

| 路徑 | 內容 |
|---|---|
| `outputs/pubtator_cache/fulltext/<pmid>.json` | full=true 全文 biocjson 快取 |
| `outputs/pubtator_cache/abstract/<pmid>.json` | abstract 層級 biocjson 快取 |
| `outputs/pubtator_cache/eupmc/<pmid>.txt` | EuropePMC 全文快取 |
| `outputs/pubtator_cache/llm_results.jsonl` | **LLM 檢查點**：每筆 `{key, info}`，中斷後重跑會自動略過已完成的 |
| `outputs/pubtator_harvest.log` | **自動寫檔的執行 log**（append 模式，每次執行前有 `=== 新一次執行開始 ===` 分隔線） |

- **續跑**：直接用相同指令再執行即可，已下載的全文與已完成的 LLM 結果都會跳過。
- **重來**：`--refresh` 重抓全文；`--refresh-llm` 清掉 LLM 檢查點重新推論。

---

## 7. 跑「全部文獻」

PubTator `@GENE_TTN` 目前命中 **14,199 篇**（會在過濾 Journal Article 後略少）。建議分兩階段：

```bash
# 1) 先全量收割變異（不需 GPU，約 1–2 小時，含全文下載快取）
conda run -n gm4 python pubtator_harvest.py --max-articles 15000 --no-llm \
  --out outputs/pubtator_ttn_variants_all.csv

# 2) 再跑 LLM（全文已快取、檢查點可續跑，可在背景分段進行）
PYTHONUNBUFFERED=1 conda run --no-capture-output -n gm4 \
  python -u pubtator_harvest.py --max-articles 15000 \
  --out outputs/pubtator_ttn_variants_all.csv
```

> ⚠️ 全量 LLM 推論非常耗時（外推約 7 萬個 job × n=5，需數天 GPU）。因有檢查點，可隨時中斷再續。

---

## 8. 效能 / OOM 調校

- **吞吐量**主要取決於同時解碼的序列數：把 `--llm-batch-size` 加大、或調高 `config.LLM_MAX_NUM_SEQS` 可提速，但 KV cache 峰值也會上升。
- 若出現 `Cannot get free blocks`（KV cache OOM）：調小 `LLM_MAX_NUM_SEQS` 或 `--llm-batch-size`，或縮短 `LOCAL_LLM_MAX_CONTEXT_LENGTH`。
- 想降成本：把 `config.LLM_SELF_CONSISTENCY_N` 從 5 降到 3，序列數直接少 ~40%（代價是投票穩定度略降）。

---

## 9. 常見問題

- **抓不到變異**：很多 TTN 文獻是綜述/功能研究，內文沒有可抽取的明確變異，屬正常；只有部分文獻會產出變異列。
- **沒有 API key 很慢**：在 `config.py` 設 `PUBMED_API_KEY`（NCBI 帳號可申請）可把限速從 3 req/s 提到 10 req/s，加速 Stage 1/1b。
- **全文 doc id 對不上**：`full=true` 回傳的 `doc['id']` 其實是 PMCID 且批次會重排，程式用 passage 的 `article-id_pmid` 對應回 PMID（已處理）。
