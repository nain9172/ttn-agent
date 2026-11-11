# TTN Variant AI Agent

## 網站服務（Docker）

在 DGX Spark 上用 Docker 架設報告網站：內網使用者用瀏覽器開 `http://<伺服器 IP>:8080` 就能瀏覽報告、送出新的變異分析，
分析在同一台 Spark 的 GPU 上執行，結果直接加進網站。資料、模型與結果都放在主機上的一個資料夾（`.env` 的 `TTN_DATA_DIR`）。

```bash
cp deploy.env.example .env        # 填入 NCBI email、允許連線的內網網段、資料夾位置與帳號
mkdir -p /data/ttn                # 資料夾（與 .env 的 TTN_DATA_DIR 相同）
docker compose build
docker compose up -d
```

資料夾是空的時候，網站會導向設定頁 `/setup`：可以下載公開資料與模型、試跑、從頭重算，並列出授權資料要放在哪裡。
完整步驟（從 NAS 套件準備資料、內網限制、授權資料、日常操作）見 [docker_readme/DEPLOY_DOCKER.md](docker_readme/DEPLOY_DOCKER.md)。

## 安裝環境

```bash
./setup_env.sh                           # 建立或修補 gm4 conda 環境（包含 flash-attn）
conda activate gm4
```

`setup_env.sh` 只負責建立環境，不會下載資料或執行任何分析，而且可以重複執行：已經裝好的部分會跳過，只補上缺少的。
用 `--env NAME` 可以指定其他環境名稱。`setup_and_run.sh` 的 `env` 步驟就是呼叫這個腳本。

在 repo 根目錄建立 `.env` 檔，填入你的 NCBI 聯絡資訊。這個檔案不會進 git。NCBI 規定要附上 email；
API key 可以把請求速率從每秒 3 次提高到 10 次。

```
PUBMED_EMAIL=you@example.com
PUBMED_API_KEY=your_ncbi_api_key
```




## 網頁介面（舊版：Windows 啟動程式 + SSH）

> 已由上面的 Docker 網站服務取代；以下保留給還在用 `TTN-Site.exe` 的人。

報告網站可以把還沒有報告的變異送到 DGX Spark，由 Spark 跑完整 pipeline，新報告會自動加進網站。

**在 Spark 上（只需一次）：** clone 這個 repo，建立好 `gm4` 環境並準備好模型，確認 `python main.py <variant>` 能正常執行，再檢查任務腳本：

```bash
~/miniconda3/envs/gm4/bin/python ~/ttn-agent/spark_job.py check
```

**在 Windows 電腦上：** 把site資料夾下載並打開 `TTN-Site.exe`（ARM 筆電請用 `TTN-Site-arm64.exe`）。瀏覽器會打開網站，最上方有 **Analyze variants on your DGX Spark** 面板。

1. 填入 Spark 的 IP(140.xxx.xxx.xxx) 和你的使用者名稱，按 **Save & test connection**。
2. 第一次連線時，把畫面上顯示的指紋和在 Spark 上執行 `ssh-keygen -lf /etc/ssh/ssh_host_ed25519_key.pub` 的結果比對，一致後按 **Trust this Spark**。
3. 如果還沒有 SSH 金鑰，按 **Set up SSH key…** 並輸入一次 Spark 密碼。TTN-Site 會建立 `%USERPROFILE%\.ssh\id_ed25519` 並安裝到 Spark 上，密碼不會被儲存。
4. 以 `2-<position>-<ref>-<alt>` 格式輸入變異（GRCh38 基因組正股，與 VCF / ClinVar 相同，例如 ClinVar 332662 輸入 `2-178525997-T-C`；TTN 轉錄方向（c. 寫法的鹼基）`2-178525997-A-G` 會被拒絕，並提示對應的正股寫法）。也可以用 HGVS c. / g.（例如 `NM_001267550.2:c.*1015A>G`）或 dbSNP rsID（例如 `rs72629798`），送出後由伺服器換算；p. 寫法不接受。各種寫法的範例見網頁上的 **Download template**，一行一個，也可以載入 CSV 或文字檔，然後按 **Run on Spark**。
   已經有報告的變異只會顯示連結；勾選 **Re-run variants already in this site** 可以對它們重新跑完整 pipeline。
   新報告會取代舊報告；如果重跑失敗，舊報告會保留。

請讓 TTN-Site 的視窗保持開著。任務完成時 Windows 會跳出通知，新報告會出現在表格中並帶有 **Spark** 標籤。
這些報告存在網站資料夾的 `data/local-reports.js`，所以直接打開 `index.html` 也看得到

==================================主要用以上這些==================================
## 快速開始：一個腳本跑完全部(包含下載ClinVar資料、reference sequence，執行PubTator文獻挖掘，把生成報告整理成網站)

```bash
./setup_and_run.sh --limit 5      # 先用前 5 顆 ClinVar 變異試跑
./setup_and_run.sh                # 完整執行（建議在 tmux 裡跑，需要好幾天）
```

## 使用方式

### 單一變異

```bash
python main.py 2-178527121-A-G                    # chr-pos-ref-alt（GRCh38 基因組正股，與 VCF / ClinVar 相同；轉錄方向寫法也能解析）
python main.py "NM_001267550.2:c.29800T>C"         # 也接受 HGVS c./p./g.
python main.py 2-178527121-A-G --skip-evo2         # 跳過 Evo2（較快）
python main.py 2-178527121-A-G --skip-clinical-extraction   # 不用 LLM；Track B 改用規則式推論
```

其他選項：`--output PATH`、`--llm-backend {vllm,transformers,ollama}`、`--llm-model NAME`，
以及 `--require-llm`（LLM 無法執行時直接失敗，而不是產生沒有 LLM 分析的報告）。

### 批次執行

```bash
python batch_run_variants.py --input data/clinvar_ttn.csv
./run_batch_inprocess.sh --input data/clinvar_ttn.csv --skip-evo2 --require-llm   # LLM 只載入一次，出錯會自動重啟
```

輸入檔需要有 `input` 欄（`chr-pos-ref-alt`），或是 `chr`/`position`/`ref`/`alt` 這幾欄。每跑完一顆就會寫回結果，
中斷後重跑會從斷點繼續。`outputs/` 裡已有報告的變異會直接沿用；修改過分析引擎後要加 `--force` 才會重跑。

## 一次性初始化：Evo2 分數表

報告網站可以查詢還沒有報告的變異的 Evo2 分數，這需要先算好 TTN 上所有 SNV 的分數表
（281,435 個位置 × 3 種替代鹼基，8,192 bp 窗口，計分方式與 `main.py` 相同）。clone 之後執行一次：

```bash
conda activate gm4
./init_evo2_table.sh            # 缺 TTN 序列時會先下載，然後在背景開始計算
./init_evo2_table.sh --status   # 查看進度
./init_evo2_table.sh --stop     # 目前這批算完後停止；重新執行第一行指令即可續算
```

在 GB10 上全部算完約需 31 天，期間會一直佔用 GPU，所以不要同時跑 vLLM 批次。外顯子與剪接區（±10 bp）會先算，約需 14 天。
每次停止或算完都會寫出網站用的 `outputs/evo2_table/evo2_table.js`；也可以執行 `./init_evo2_table.sh --export` 手動更新。

分數和執行環境有關：gm3（FP8）和 gm4 的分數約差 1e-4，足以讓接近門檻的變異判定翻轉。請使用與 `main.py` 相同的環境。
所用的套件版本記錄在 `outputs/evo2_table/meta.json`。](https://www.ncbi.nlm.nih.gov/research/litvar2/)

## 專案結構

| 路徑 | 說明 |
|---|---|
| `main.py` | 對單一變異執行完整 pipeline，輸出 `outputs/<variant>_<model>.html` |
| `batch_run_variants.py`、`run_batch_inprocess.sh` | 從 CSV/XLSX 批次執行，可中斷續跑；每份報告的結果會寫回表格 |
| `build_site.py`、`site_info.py` | 把 `outputs/` 的報告打包成 `site/` 網站，並產生說明頁 `info.html`（band、isoform、分數與推論說明） |
| `web_server.py`、`deploy_status.py` | Docker 網站服務：網頁、API、volume 設定頁 `/setup` |
| `Dockerfile`、`docker-compose.yml`、`docker/`、`deploy.env.example` | Docker 映像與服務設定（見 `docker_readme/DEPLOY_DOCKER.md`） |
| `spark_job.py` | 分析任務佇列（Docker 的 worker 容器；也供舊版啟動程式透過 SSH 呼叫） |
| `launcher/` | 舊版 Windows 啟動程式 `TTN-Site.exe`（Go，透過 SSH 連 Spark） |
| `setup_and_run.sh` | 一鍵完成安裝設定與完整執行（見「快速開始」） |
| `setup_env.sh` | 只建立 conda 環境（`environment.yml`、flash-attn wheel、環境檢查） |
| `score_evo2_for_variants.py`、`backfill_evo2_scores.py` | 批次計算 Evo2 分數，並寫回以 `--skip-evo2` 產生的報告 |
| `reproduce_evo2.sh` | 重現 Evo2 分數：建環境、下載參考序列與 ClinVar TTN 變異、對全部變異計分，並以 ClinVar 標籤評估（AUROC、敏感度、特異度） |
| `backfill_database_annotations.py` | 不重跑 pipeline，直接重新產生既有報告的資料庫比對區塊 |
| `backfill_report_ui.py` | 不重跑 pipeline，把報告改版（正股 ID、Evo2 說明、isoform 示意圖、推論標示、遺傳模式、側邊目錄）套到既有報告 |
| `init_evo2_table.sh`、`precompute_evo2_table.py` | 一次性預先計算 TTN 上所有 SNV 的 Evo2 分數 |
| `build_ttn_transcript_model.py` | 建立 `data/ttn_transcript_model.json`（HGVS ↔ 基因組座標用的 RefSeq 外顯子模型） |
| `config.py` | 所有設定：模型、門檻、資料路徑、功能開關 |
| `utils/` | pipeline 的各個模組（`variant_ids.py`：正股（VCF / ClinVar）與轉錄方向 ID 互轉；`inheritance_inference.py`：AD / AR 推論） |
| `scripts/data_prep/` | 輸入資料：下載 ClinVar TTN、PubTator harvest 與新文獻檢查、ClinVar 標註、HGMD ↔ ClinVar ID |
| `research/` | 評估、實驗、論文與簡報圖表、notebook；pipeline 不會用到 |
