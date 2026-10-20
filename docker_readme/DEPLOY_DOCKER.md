# TTN 報告網站：Docker 部署說明

這份文件說明如何在 DGX Spark 上用 Docker 架設 TTN 報告網站。架好之後，內網的使用者用瀏覽器開
`http://<伺服器 IP>:8080` 就能瀏覽報告、送出新的變異分析，不需要在自己電腦上安裝任何東西。

## 1. 架構

```
瀏覽器 ──HTTP──▶ web 容器（網站與 API，不用 GPU）
                     │  任務寫進資料夾
                     ▼
               主機上的資料夾 /data/ttn（資料、模型、報告、網站、任務佇列）
                     ▲
                     │  取任務、寫結果
                 worker 容器（GPU：Evo2、Gemma 4、完整 pipeline）
```

- **web**：提供網站與 API（`web_server.py`）。
- **worker**：一次處理一個任務（`spark_job.py worker --forever`），完成的報告直接加進網站，所有人重新整理就看得到。
- 兩個容器共用主機上的一個資料夾（`.env` 的 `TTN_DATA_DIR`，以下以 `/data/ttn` 為例），在容器裡看到的是 `/vol`。
  GPU 忙碌或 pipeline 出錯時，網站仍可以正常瀏覽。
- 容器以主機帳號的身分執行（`TTN_UID`／`TTN_GID`），資料夾裡的檔案都屬於這個帳號，可以直接用 `cp`、`ls` 管理，不需要 sudo。
- 映像檔只包含程式與執行環境；資料、模型權重、分析結果都放在資料夾。**授權資料（HGMD、cardiodb）不會進映像檔。**

## 2. 需求

| 項目 | 說明 |
|---|---|
| 機器 | DGX Spark（aarch64、GB10）。映像與 flash-attn 都是 aarch64 版，不能直接用在 x86 主機 |
| Docker | Docker Engine 與 Docker Compose v2（`docker compose version` ≥ 2.24） |
| GPU | NVIDIA Container Toolkit（DGX OS 已內建）；`docker run --rm --gpus all nvidia/cuda:13.0.0-base-ubuntu24.04 nvidia-smi` 要能看到 GPU |
| 權限 | 執行的帳號要在 `docker` 群組：`sudo usermod -aG docker $USER`，然後重新登入 |
| 磁碟 | 映像約 20 GB；資料夾約 70 GB 起（模型權重約 43 GB、報告與網站約 15 GB），建議預留 100 GB |
| 網路 | 第一次建置與下載資料時需要連到 GitHub、conda-forge、PyPI、Hugging Face、NCBI |
| 記憶體 | Gemma 4 需要約 104 GiB 可用的統一記憶體；同一台機器若有其他人在跑 GPU 工作，任務會排隊等待 |

## 3. 安裝

```bash
git clone https://github.com/nain9172/ttn-agent.git
cd ttn-agent
cp deploy.env.example .env        # 編輯 .env，見下方說明
mkdir -p /data/ttn                # 資料夾，用執行 docker 的帳號建立（路徑要與 .env 的 TTN_DATA_DIR 相同）
docker compose build              # 第一次約 30–60 分鐘（建立 conda 環境）
docker compose up -d
docker compose ps                 # web、worker 都應該是 running
```

如果拿到的是匯出好的映像檔（`ttn-agent-image.tar`），用 `docker load -i ttn-agent-image.tar` 取代 `docker compose build`，可以省下建置時間，也不需要連 conda-forge／PyPI。

`.env` 要填的欄位：

| 變數 | 說明 |
|---|---|
| `PUBMED_EMAIL` | 必填。NCBI 規定每個請求都要附聯絡 email |
| `PUBMED_API_KEY` | 建議填。速率從每秒 3 次提高到 10 次（NCBI 帳號 → API Key Management） |
| `HF_TOKEN` | Hugging Face token（Read 權限）。下載需要同意授權的模型時必須 |
| `TTN_ALLOWED_NETWORKS` | 允許連線的內網網段，見第 4 節 |
| `TTN_PORT`、`TTN_BIND_ADDRESS` | 網站埠（預設 8080）與要開放的網路介面 |
| `TTN_DATA_DIR` | 資料夾的完整路徑，例如 `/data/ttn`。必須先用執行帳號 `mkdir` 建立 |
| `TTN_UID`、`TTN_GID` | 容器以哪個帳號執行。填 `id -u`、`id -g` 的結果，資料夾裡的檔案就屬於這個帳號 |

接著用瀏覽器開 `http://<伺服器 IP>:8080`：

- 資料夾還是空的：會自動導向設定頁 `/setup`，照第 5 節建立資料。
- 已經有網站：直接顯示報告列表。

## 4. 限制只有內網可以連線（請院內資訊人員設定）

網站只接受 `TTN_ALLOWED_NETWORKS` 裡的來源 IP，其他來源一律回 403。

- **預設值**：只允許私有網段（`10.0.0.0/8`、`172.16.0.0/12`、`192.168.0.0/16`）與本機。
- **院內使用公開 IP 網段時**：如果醫院內網的電腦使用公開 IP（例如 `140.112.x.x`），預設值會把他們擋掉，必須把該網段加進去，例如：

  ```
  TTN_ALLOWED_NETWORKS=127.0.0.0/8,10.0.0.0/8,172.16.0.0/12,192.168.0.0/16,140.112.0.0/16
  ```

  實際網段請向院內資訊室確認。改完執行 `docker compose up -d` 套用。

- **只綁定內網網卡**：伺服器有多張網卡（例如同時接內網與外網）時，把 `TTN_BIND_ADDRESS` 設成內網那張的 IP，網站就只在內網介面上開放。

- **注意防火牆**：Docker 對外開放的埠會**繞過 ufw**（Docker 自己寫 iptables 規則），所以只靠 `ufw deny` 擋不住。要在主機防火牆層限制，請把規則加在 iptables 的 `DOCKER-USER` chain，或交給院內的網路防火牆處理。`TTN_ALLOWED_NETWORKS` 是應用程式層的第二道防線。

> 報告中含有 HGMD（授權資料庫）的紀錄，網站只能開放給授權範圍內的使用者。

## 5. 準備資料夾裡的資料

資料夾（`TTN_DATA_DIR`，例如 `/data/ttn`）的結構如下。容器裡看到的路徑是 `/vol/...`，對應到主機的 `/data/ttn/...`：

```
/data/ttn/
├── data/                    pipeline 輸入資料（含授權檔）
│   ├── sequence.fasta, ttn_transcript_model.json, transcript_interval.xlsx
│   ├── clinvar_ttn.csv, pubmind_ttn.json, pubtator/
│   ├── HGMD_latest.csv, HGMD_TTN_processed.xlsx, cardiodb.xlsx   ← 授權資料，手動放
│   └── ttn_exon_psi.csv                                           ← 選用，手動放
├── outputs/                 報告（*.html）、示意圖、Evo2 分數、任務佇列（site_jobs/）、執行紀錄
├── site/                    網站（由 outputs/ 打包而成）
├── models/huggingface/      Evo2、Gemma 4 權重
├── cache/                   Triton / vLLM 編譯快取
└── llm_logs/                LLM 推論紀錄
```

設定頁 `/setup` 會逐項列出哪些檔案已經有、哪些還缺，並顯示每個檔案在主機上應該放的完整路徑。
有三種方式建立資料，**有 NAS 套件的話用方式 A 最快**。

### 方式 A：使用提供的完整套件（解開就能用）

完整套件 `ttn-full-<日期>.tar`（約 50 GB）就是一個準備好的資料夾，內部結構與 `TTN_DATA_DIR` 相同：

```
ttn-full-<日期>.tar
├── data/        參考序列、轉錄本模型、ClinVar、PubMind、PubTator、HGMD、cardiodb、外顯子 PSI
├── outputs/     7,792 份 HTML 報告、Evo2 分數、PubTator 快取
├── models/      Evo2-7B 與 Gemma 4 權重（約 43 GB）
└── site/        已打包好的網站
```

```bash
sha256sum -c ttn-full-<日期>.tar.sha256          # 確認傳輸沒有損壞（套件與 .sha256 放在同一個目錄）
mkdir -p /data/ttn
tar -xf ttn-full-<日期>.tar -C /data/ttn       # 用執行 docker 的帳號解壓縮，檔案就屬於這個帳號
```

解壓縮後照第 3 節啟動服務（`docker compose build`、`docker compose up -d`），用瀏覽器開 `http://<伺服器 IP>:8080`
就會看到全部報告，也可以直接送出新的變異分析，不需要再下載或重建任何東西。

> 只有資料與報告、沒有模型與網站的套件（`ttn-data-<日期>.tar`）也是同樣解壓縮到 `/data/ttn`；
> 之後在 `/setup` 依序按 **Download public data and models**（下載模型，約 1 小時）與
> **Rebuild the website from existing results**（約 10–20 分鐘）。

### 方式 B：從零開始（網頁操作）

1. 開 `/setup`，按 **Download public data and models**。會下載：
   - 參考序列、轉錄本模型、ClinVar TTN 變異清單、PubMind 公開版
   - Evo2 與 Gemma 4 權重（約 45 GB）
2. 把授權資料放進資料夾（見第 6 節）。
3. 按 **Trial run on 5 ClinVar variants**，用 5 顆變異試跑整個流程（約 30–60 分鐘），確認環境沒問題。
4. 按 **Recompute everything from scratch**。依序執行 PubTator 文獻挖掘、Evo2 計分、完整 pipeline、寫回與建站，**需要好幾天**，期間 GPU 會一直被佔用。可以關掉網頁，進度與紀錄在 `/setup` 下方的 Setup jobs。worker 一次只跑一個任務，所以重算期間使用者送出的變異分析會排在後面，等重算結束才會開始。

任務中斷（例如重開機）後，worker 會自動把它重新排隊，已完成的部分會沿用。

### 方式 C：指令列

```bash
docker compose run --rm worker status                       # 列出資料夾裡的資料狀態
docker compose run --rm worker setup --only models --only data
docker compose run --rm worker setup --limit 5              # 試跑 5 顆
docker compose run --rm worker setup --only site            # 重建網站
```

`setup` 之後的參數與 `setup_and_run.sh` 相同（`--help` 可查看）。

## 6. 授權資料

以下檔案受授權限制，無法自動下載，需要由有授權的單位提供後放進資料夾的 `data/`：

| 檔案 | 放置位置 | 用途 | 沒有的影響 |
|---|---|---|---|
| `HGMD_latest.csv` | `/data/ttn/data/HGMD_latest.csv` | HGMD Professional 的 TTN 匯出：資料庫比對、類比推論 | 報告沒有 HGMD 比對，推論只能靠其他證據 |
| `HGMD_TTN_processed.xlsx` | `/data/ttn/data/HGMD_TTN_processed.xlsx` | 提供 REF/ALT，讓 HGMD 能精確比對 allele | HGMD 只能比對到位置 |
| `cardiodb.xlsx` | `/data/ttn/data/cardiodb.xlsx` | 364 個外顯子的座標、band、心肌／骨骼肌 PSI | 推論改用轉錄本收錄度代理，準確度較低 |
| `pubmind_ttn.json`（完整版） | `/data/ttn/data/pubmind_ttn.json` | 逐篇疾病與致病性；放了就不會被公開版覆蓋 | 使用公開版 |
| `ttn_exon_psi.csv`（選用） | `/data/ttn/data/ttn_exon_psi.csv` | 323 個外顯子的實測 PSI | 改用 cardiodb 的 PSI |

直接複製進去即可：

```bash
cp HGMD_latest.csv HGMD_TTN_processed.xlsx cardiodb.xlsx /data/ttn/data/
```

放好後重新整理 `/setup`，確認這幾項都是 ✓。已經產生的報告不會自動更新；要套用到全部報告，請在 `/setup` 執行 Recompute，或用指令列 `setup --only backfill --only site` 只更新資料庫比對並重建網站。

## 7. 日常操作

| 要做的事 | 指令 |
|---|---|
| 看狀態 | `docker compose ps` |
| 看記錄 | `docker compose logs -f web`、`docker compose logs -f worker` |
| 停止／啟動 | `docker compose stop`、`docker compose start` |
| 更新程式 | `git pull && docker compose build && docker compose up -d`（資料夾不受影響） |
| 備份 | `tar -C /data/ttn -cf /backup/ttn-$(date +%F).tar data outputs`（模型與快取可以重新下載，不必備份） |

- 送出分析：在網站首頁的 **Analyze variants** 區塊輸入變異，一行一個，GRCh38 基因組正股寫法，例如 `2-178527245-T-G`；TTN 轉錄方向的寫法也接受。
- 取消任務：任務列表裡的 **Cancel**。

## 8. 疑難排解

| 狀況 | 原因與處理 |
|---|---|
| 瀏覽器顯示 403 | 來源 IP 不在 `TTN_ALLOWED_NETWORKS`，見第 4 節 |
| 首頁顯示「analysis worker is not running」 | worker 容器沒在跑或拿不到 GPU。執行 `docker compose logs worker`，並在容器內測試 `docker compose exec worker nvidia-smi` |
| 任務一直是 waiting for memory | Gemma 4 需要約 104 GiB 可用的統一記憶體。worker 會先自動清出檔案快取（GB10 上檔案快取會被算成 GPU 已使用）；仍不夠代表被其他工作佔用，任務會自動等待最多 12 小時，網頁會顯示原因。常態不夠時，可在 `.env` 設 `LLM_GPU_MEMORY_UTILIZATION=0.80` 後執行 `docker compose up -d` |
| 下載 Gemma 失敗（401／403） | 在 `.env` 設定 `HF_TOKEN`，並到該模型的 Hugging Face 頁面同意授權條款 |
| 建置映像失敗 | 確認機器是 aarch64（`uname -m`）、能連到 conda-forge 與 PyPI；重跑 `docker compose build` 會從失敗的那一層繼續 |
| 報告缺 HGMD 比對 | 授權檔沒放或位置不對，到 `/setup` 檢查 |
| 容器一直重啟，記錄顯示「資料目錄沒有寫入權限」 | 資料夾不存在（Docker 會自動用 root 建立）或擁有者不是 `TTN_UID` 的帳號。用執行帳號重新 `mkdir`，或 `sudo chown -R $(id -u):$(id -g) /data/ttn` |


## 附錄：改用 docker volume

`.env` 不設定 `TTN_DATA_DIR`（並刪掉 `TTN_UID`、`TTN_GID`）時，資料改放在 Docker 管理的 volume `ttn-data`，
實際位置在 `/var/lib/docker/volumes/ttn-data/_data`，只有 root 看得到，存取要透過 Docker，例如：

```bash
docker compose cp HGMD_latest.csv worker:/vol/data/HGMD_latest.csv
docker compose exec worker ls /vol/data
```
