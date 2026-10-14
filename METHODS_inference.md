# Track B：類比式機制推理（Analogy-Based Mechanistic Inference）

> 論文方法章節草稿。對應程式：`utils/mechanistic_inference.py`、
> 評估：`research/benchmark_inference.py`、圖表：`research/inference_benchmark_figures.ipynb`。

---

## 1. 動機與問題設定

本系統的主要流程（Track A）是 evidence-grounded 的文獻擷取：由 ClinVar 取得該變異被
引用的文獻，再由 LLM 從全文擷取表型描述。此流程有一個結構性的失效條件——**當
ClinVar 未引用任何文獻，或所引文獻並未提及該變異時，Track A 無法產生任何有依據的
判斷**，只能回報「無證據」。

在 TTN 上這並非邊緣案例。TTN 是人類最大的基因（363 個編碼外顯子、約 36,000 個
胺基酸），變異數量龐大而文獻覆蓋稀疏；本研究的 PubTator 評估集中有 **91.5%** 的變異
在 HGMD 中查不到同位置的收錄紀錄。

Track B 針對此缺口，回答一個較弱但仍有臨床意義的問題：

> 在完全沒有該變異的直接文獻證據時，能否依據其**結構位置、變異類型、所影響的
> 轉錄本**，以及**基因內已知致病變異的類比**，推論出最可能的疾病與受累組織？

輸出明確標記為 `evidence_strength = inferred_no_literature`，與 Track A 的
evidence-grounded 結果在報告上分軌呈現，不混為一談。此設計原則與 VarChat
（Bioinformatics 2024）明確標示「答案來自文獻或來自模型知識」的做法一致。

方法論上，這屬於**類比推理／案例式推理（case-based reasoning）**，在臨床遺傳學中
對應 paralogue annotation 的邏輯（一顆變異在某位置已知致病，可用以註釋臨床資料為零
的等價位置），在 ACMG/AMP 框架中則對應 PS1（同一胺基酸改變）、PM5（同一密碼子的
不同改變）與 PM1（突變熱點／關鍵功能域）三條規則。

---

## 2. 系統概觀

給定一顆無文獻證據的變異，Track B 依序執行：

```
輸入變異 (HGVS 或 chr-pos-ref-alt)
   │
   ├─ 階段一：特徵建構
   │    ├─ 定位：結構 band、meta 外顯子編號
   │    ├─ 變異類型：truncating vs missense
   │    ├─ 外顯子使用率 (PSI) → 組織 likelihood ratio
   │    ├─ HGMD 同位置查詢
   │    ├─ 外顯子熱點比對 (ACMG PM1)
   │    └─ 類比檢索：同 band + 同類型的已知致病變異
   │
   └─ 階段二：決策
        ├─ 第 1 層　HGMD 完全命中 → 直接採用
        ├─ 第 2 層　外顯子熱點 → 硬規則覆寫
        └─ 第 3 層　距離加權投票 → 疾病 + 組織
```

決策可由**規則式**（rule-based）完成，亦可將階段一的特徵組成 prompt 交由本地 LLM
（MedGemma / Gemma，經 vLLM 推論）撰寫推論說明。本章的量化評估**全部基於規則式路徑**：
規則式是可重現、可歸因的那一層，LLM 會把 prompt 變動的效果混入，不利於歸因分析。

---

## 3. 資料來源

| 資料 | 內容 | 用途 |
|---|---|---|
| HGMD TTN 匯出 | 1,094 顆已知致病變異（truncating 792、missense 264、類型未定 38；心肌 743、骨骼肌 322、兩者 20、未定 9） | 類比案例庫 |
| `data/ttn_transcript_model.json` | 7 條 curated RefSeq 轉錄本的逐外顯子基因體座標 | 外顯子定位、isoform 收錄度 |
| `data/ttn_exon_psi.csv` | 323 個外顯子的心肌／骨骼肌實測 PSI、metatranscript-only 旗標 | 組織 likelihood ratio |
| `TTN_BAND_EDGES` | 四個結構 band 的基因體邊界 | 粗定位、統計先驗 |

**外顯子 PSI 來源**：Circ Genom Precis Med (2025) 之補充表 S2。該表使用
`NM_001267550.2` 的外顯子編號，與本研究一致。解析結果通過兩項驗證：(1) 可重現論文
的兩個關鍵數字——心肌高 PSI（≥90%）外顯子 **187** 個、「骨骼肌高 PSI 但心肌低 PSI」
**106** 個且全部位於 I-band；(2) 323 個外顯子中有 **320** 個的 band 判定與
`TTN_BAND_EDGES` 一致，3 個不一致（exon 29／252／358）全部落在 band 邊界，屬座標
切點的一個外顯子誤差，非編號錯位。

---

## 4. 特徵建構

### 4.1 定位

**結構 band**：以 `TTN_BAND_EDGES` 對基因體座標做二分查找，分入 M-band／A-band／
I-band／Z-disk 四類。

**外顯子編號**：以 meta-transcript `NM_001267550.2` 的外顯子區間做二分查找。模型含
363 個編碼外顯子，而文獻與 RefSeq 使用 364 個外顯子的編號（含非編碼的 exon 1），
兩者相差 1；此偏移以三顆已知變異校正並固定（HMERF 的 `p.C31712R`、`p.N31786K` →
文獻 exon 344；TMD 的 FINmaj `c.107780_107790` → 文獻 exon 364）。落在內含子者
回傳 `None`。

### 4.2 變異類型分類

由蛋白質／核苷酸 HGVS 以正則式離線判定，不呼叫 VEP：

- **truncating**：nonsense（`p.X###*`／`Ter`）、frameshift（`fs`）、canonical splice
  （`c.####±1`、`c.####±2`）
- **missense**：`p.Xxx###Yyy` 或 `p.X###Y` 且非終止

此二分法對應 TTN 的機制分水嶺：TTNtv 經 nonsense-mediated decay 導致全長 titin 單倍
劑量不足，致病機率較 missense 高約一個數量級；missense 在族群中普遍且多為良性或
意義未明。

### 4.3 外顯子使用率與組織 likelihood ratio

這是本方法相對於單純位置類比的主要擴充。

**動機**：TTN 的心肌與骨骼肌表型差異來自 isoform 的組織特異性表現。心肌主要表現
N2B／N2BA、骨骼肌表現 N2A；A-band 與 M-band 外顯子為 constitutive（所有 isoform 皆
保留），替代性剪接集中於 I-band 的 spring 區。因此**組織判別力理論上應集中在 I-band**。

**為何不能用二元 isoform 歸屬**：本研究實測發現，363 個 meta 外顯子中僅 **1 個**
（exon 48，N2B 專屬元件）為心肌專屬，**0 個**為骨骼肌專屬——因為長心肌 isoform
N2BA 涵蓋了 N2A 的全部 122 個外顯子。二元歸屬在整條基因上幾乎恆為真，不具判別力。
判別訊號必須是**量化的 PSI**，而非轉錄本的有無。

**分桶**：依論文自身的定義（PSI ≥ 90% 為 high PSI）將外顯子分為

| 桶 | 條件 | 外顯子數 |
|---|---|---|
| `both-high` | 心肌與骨骼肌 PSI 皆 ≥ 90% | 185 |
| `skeletal-high` | 骨骼肌 ≥ 90%、心肌 < 90% | 106（全部 I-band） |
| `meta-only` | 補充表的 metatranscript-only 旗標 | 17 |
| `both-low` | 兩者皆 < 90% | 13 |
| `cardiac-high` | 心肌 ≥ 90%、骨骼肌 < 90% | 2 |
| `absent` | 不在補充表涵蓋範圍者 | 40 |

`absent` 為補充表刻意排除的外顯子（exon 45–46、147、149、164、172–201），皆為
I-band 替代性剪接區；此類併為單一桶處理，若依轉錄本收錄度再細分會使各桶樣本量
低於下述門檻而全數棄權，反而損失訊號。

**Likelihood ratio 估計**：對桶 $b$ 與組織 $t \in \{\text{cardiac}, \text{skeletal}, \text{both}\}$，

$$
\mathrm{LR}_t(e) = \left[\frac{\left(n^{(-e)}_{b,t}+1\right)\big/\left(N^{(-e)}_{b}+3\right)}
{\left(n^{(-e)}_{t}+1\right)\big/\left(N^{(-e)}+3\right)}\right]^{\alpha}
$$

其中 $n_{b,t}$ 為桶 $b$ 內組織為 $t$ 的 HGMD 變異數，$n_t$ 為全庫對應數，分子分母
皆採 Laplace 平滑（+1／+3）。指數 $\alpha = 2$，結果夾在 $[1/4,\,4]$。

兩項防護措施：

1. **防循環**：上標 $(-e)$ 表示估計時**扣除查詢變異所在外顯子的全部 HGMD 變異**。
   查詢變異本身與其外顯子鄰居皆不參與自身先驗的估計。此設計回應變異效應預測評估
   中已知的 circularity 問題（同一基因／區域的變異同時出現在訓練與評估兩側會系統性
   高估效能）。
2. **樣本量門檻**：桶內有效樣本數 $N^{(-e)}_{b} < 15$ 時回傳 $\mathrm{LR}=1$（不表態），
   避免以極少樣本推出極端 LR。實務上 `meta-only`（n=5）、`both-low`（n=4）、
   `cardiac-high`（n=3）三桶恆為棄權。

**實際估計值**（扣除自身外顯子後）：

| 桶 | n | 心肌比例 | LR(心肌) | LR(骨骼肌) |
|---|---|---|---|---|
| `both-high` | 726 | 74% | 1.16 | 0.67 |
| `skeletal-high` | 97 | 44% | 0.41 | 3.17 |
| `absent` | 27 | 22% | 0.25 | 4.00 |
| 全庫基準 | 862 | 69% | — | — |

以 LR 而非 High／Low 標籤輸出，形式上對齊 ACMG 證據強度的量化校準慣例（PM5 的
positive LR 經校準為 8.4–16.3；同源位置類比為 13.0）。

### 4.4 外顯子熱點（ACMG PM1）

兩個文獻已確立、表型專一性遠高於所在 band 先驗的外顯子，以硬規則處理：

| 外顯子 | 結構域 | 適用類型 | 疾病 | 組織 |
|---|---|---|---|---|
| 344 | FN3-119 (A150) | **僅 missense** | HMERF | 骨骼肌 |
| 364 | Mex6（末端外顯子） | 不限類型 | TMD／LGMD2J | 骨骼肌 |

exon 344 限定 missense 係依原始文獻的條件——HMERF 是 FN3-119 域的 missense 專屬
表型，該外顯子的 truncating 變異走一般 A-band TTNtv 機制。exon 364 為末端外顯子，
其等位基因逃脫 NMD，機制與 A-band TTNtv 不同；異型合子造成 TMD、同型合子造成
LGMD2J，二者為等位疾病故合併報告。

### 4.5 類比檢索

自 HGMD 檢索 **同一 band 且同一變異類型** 的已知致病變異，依基因體距離
$d = |\mathrm{POS}_{\text{analog}} - \mathrm{POS}_{\text{query}}|$ 由近至遠排序，取前
$k = 8$ 顆。若同類型無任何結果，才放寬為 band-only 並標記 `match_tier = "band"`。

系統另實作了外顯子分層檢索（`exon+type` → `band+type` → `exon` → `band`），但實測
其疾病命中率並未優於 band 檢索（見 §7）。原因為「同外顯子但不同變異類型」的類比
命中率僅 47.6%，顯著低於對照；**truncating／missense 的機制分野比共處同一外顯子
更能決定表型**，band 檢索因強制同類型而避開此問題。

---

## 5. 決策程序

三層優先序，上層命中即回傳：

### 第 1 層　HGMD 完全命中

若 HGMD 收錄 `pos + ref + alt` 完全相同的變異，直接採用其疾病與組織標註，
confidence 設為 High。此類變異本身即為已知致病變異，僅因 ClinVar 未引用文獻而
進入本流程；此層為查表而非推論，故在評估中必須嚴格排除（見 §8.2）。

### 第 2 層　外顯子熱點

落入 §4.4 之熱點且變異類型相符者，以熱點標註覆寫，confidence 依 `strength`
設為 High／Moderate。此層屬 ACMG PM1 等級的區域級定位證據。

### 第 3 層　距離加權投票

絕大多數變異由此層決定。每顆類比案例 $a$ 的票重為

$$
w(a) = \underbrace{0.5^{\,d_a/\lambda}}_{\text{距離}}
\times \underbrace{\tau(\text{tier}_a)}_{\text{定位層級}}
\times \underbrace{\mathrm{LR}_{t(a)}(e)}_{\text{PSI 組織}}
\times \underbrace{\rho(a)}_{\text{RANKSCORE}}
\times \underbrace{\nu(a)}_{\text{籠統病名}}
$$

各項定義：

- **距離權重**：$\lambda = 150$ bp，即每 150 bp 減半。此設定使投票實質上由**最近鄰
  主導**——3 kb 外的類比權重已衰減至 $2^{-20} \approx 10^{-6}$。此性質與空間鄰近度
  方法的既有結論一致：具預測力的是與致病變異**群集**的鄰近度，而非個別遠鄰。
- **定位層級**：$\tau = \{$`exon+type`: 1.0, `band+type`: 0.6, `exon`: 0.4, `band`: 0.2$\}$。
- **PSI 組織 LR**：以**查詢變異所在外顯子**的 LR，對**該類比案例的組織標註**加權。
  例如查詢變異落於 `skeletal-high` 外顯子時，骨骼肌類比乘 3.17、心肌類比乘 0.41。
- **RANKSCORE 折扣**：$\rho = (1-s) + s \cdot \mathrm{rankscore}$，$s$ 為強度參數；
  缺值者給中位數 0.75（不獎不罰）。$s=0$ 時停用。
- **籠統病名折扣**：類比之疾病標註為上位詞（如 "Cardiomyopathy"、"Titinopathy"）
  時乘以 $\nu < 1$。$\nu = 1$ 時停用。

**疾病**取票重加總最大者，$\hat{D} = \arg\max_D \sum_{a:\,\mathrm{phen}(a)=D} w(a)$；
差異診斷取第 2、3 名。**組織**再投一次票，僅由「疾病等於 $\hat{D}$」的類比參與；
若 isoform 覆蓋明確指向單一組織則直接採用，否則依序退回全體組織票數、band 統計先驗。

---

## 6. LLM 推論路徑

當提供 LLM 時，階段一的特徵組成結構化 prompt，內容包含：TTN 機制背景（TTNtv 與
haploinsufficiency、isoform-組織對應、splice 變異的 Evo2 分數偏低警示）、目標變異的
band／外顯子／類型／PSI 意涵／Evo2 分數與信心桶、類比案例表（含距離與疾病）、以及
band 層級的統計先驗。

Prompt 明確要求模型**不得杜撰引用**，並說明其輸出為 inference 而非 citation。
輸出以 guided decoding 約束於固定 JSON schema（`inferred_disease`、`inferred_tissue`
∈ {Cardiac, Skeletal, Both, Uncertain}、`confidence`、`differential_diagnoses`、
`rationale`、`analog_variants_used`），sampling 參數為 temperature 0.2、top-p 0.95、
max_tokens 1536。

後處理有兩道：HGMD 完全命中時以資料庫標註覆寫模型輸出；模型回答 Uncertain 但確有
類比案例時，以規則式多數決回填。此二者皆記錄於 `source` 欄位以利歸因。

> **方法學說明**：強制模型必須給出單一疾病、不得回答 Uncertain，是刻意的設計選擇，
> 目的是產生可評估的預測；其代價是移除了模型的棄權（abstention）能力，校準責任因此
> 完全落在下游。此為已知限制，見 §9。

---

## 7. 評估設計

### 7.1 兩套獨立的 ground truth

| | 來源 | n | 與類比庫的關係 |
|---|---|---|---|
| **HGMD** | 每顆變異 leave-one-out，答案為其自身 PHEN 標註 | 1,072 | 同源，會高估 |
| **PubTator** | 文獻標註，僅採信原文判為 Pathogenic／Likely pathogenic 者 | 503 | 完全獨立 |

PubTator 的 P/LP 過濾並非可有可無的調校：其 disease 欄位由 LLM 自全文擷取，大量
標註實為「收案世代的疾病」而非變異表型（最明顯者為 HCM——TTN 並非公認的 HCM 致病
基因，卻佔未過濾集的 17%，其中 57% 原文分級為 VUS 或良性）。不過濾不僅使分數低估
約 10 個百分點，更會使結論反轉。

**任何改動須在兩套 ground truth 上同時成立方予採納。**

### 7.2 Leave-one-out 的範圍

HGMD 評估採 `loo_scope = "position"`：預測前將**查詢變異所在位置的全部變異**自類比
庫移除。此範圍的選擇至關重要——僅移除變異本身（`"exact"`）並不足夠，因為同位置的
其他變異會以 $d = 0$ 的類比進入投票，而 $0.5^{0} = 1.0$ 正是距離權重的最大值，答案
仍會由側門洩漏。

### 7.3 指標

- **hit**（疾病命中率）：預測與答案的 canonical label 相同（`exact`），或為上位詞
  ↔ 特定病名的相容關係（`compatible`）。
- **specific**：hit 中再排除「預測本身即為籠統上位詞」者。此指標與 hit 的差值即為
  「以上位詞作答」的比例。
- **macro**：先算每個答案類別各自的 hit 率再平均。TTN 的答案分佈極度偏斜（DCM 佔
  大宗），hit 會被多數類別灌水，macro 方能反映少數表型。
- **tissue**：組織相容率（exact + partial），僅計有組織標註者。

---

## 8. 結果

### 8.1 主要結果

| 組態 | HGMD hit | HGMD specific | HGMD tissue | PubTator hit | PubTator specific |
|---|---|---|---|---|---|
| baseline（band 檢索、無 PSI） | 70.1% | 50.7% | 73.1% | 56.9% | 38.8% |
| 逐項累加至 `+psi` | 69.4% | **54.2%** | **74.7%** | 55.5% | **41.9%** |
| **`band-hit`（命中率最佳）** | **71.5%** | 51.7% | 74.3% | **58.4%** | 40.0% |

`band-hit` 之設定為：band 檢索、外顯子熱點開啟、RANKSCORE 折扣停用（$s=0$）、
籠統病名折扣停用（$\nu=1$）、PSI 組織 LR 強度 $\alpha=2$。

**須明確交代的取捨**：`band-hit` 的命中率優勢有相當部分來自「允許以上位詞作答」。
實測其預測使用籠統病名的比例為 **27.4%**，而 `+psi` 為 **21.5%**；由於 hit 的定義
將「上位詞 ↔ 特定病名」計為命中，放寬此項等同一個穩賺 hit 的安全策略，`specific`
同步下降 2.5 個百分點即為其代價。若研究目標強調「能給出具體病名」，應採 `+psi`。

### 8.2 消融一：PSI 組織 LR 的貢獻

在其餘三個槓桿（檢索模式 × RANKSCORE × 籠統病名折扣）構成的 **8 種組合**上，分別
比較停用與啟用 PSI 組織 LR。結果為 **8/8 全部改善**（HGMD 命中率 +0.8 至 +2.1 個
百分點），且增益集中於 **I-band**（60.4% → 63.7%，n=328）——正是文獻預測判別外顯子
所在的位置。此一致性排除了單點噪音的可能。

### 8.3 消融二：移除資料庫查表

為確認報告的準確率衡量的是推論而非查表，設計三個層次的消融：

| 條件 | HGMD hit | PubTator hit |
|---|---|---|
| 現行 | 71.5% | 58.4% |
| A. 停用 `hgmd_lookup` | 71.5% | 58.4% |
| B. 同位置變異移出類比庫 | 71.5% | 58.1% |
| C. A ＋ B | 71.5% | 58.1% |

三項發現：

1. **HGMD 組本就不含查表**。三條件命中率完全相同，因 `loo_scope="position"` 早已
   移除同位置變異，`hgmd_lookup` 無從命中。此結果同時驗證了 leave-one-out 的實作。
2. **僅停用 `hgmd_lookup` 是無效的消融**（PubTator 58.4% → 58.4%）。同一顆變異會以
   $d=0$ 的類比回到投票，權重為最大值 1.0。有效的消融必須自類比庫移除。
3. **最嚴格條件下命中率僅降 0.3 個百分點**（58.4% → 58.1%）。PubTator 集中僅 5.4%
   的變異在 HGMD 有完全命中、3.2% 有同位置命中，91.5% 完全查不到。

**結論：報告的準確率衡量的是類比推理，而非資料庫查表。** 此點對 Track B 的價值
主張（無文獻證據時仍能給出有依據的推測）具決定性。

須一併交代的界線：上述消融中 `exon_hotspot_rules` 仍為開啟。該規則屬文獻確立的
**區域級**知識（ACMG PM1 等級），而非查詢該變異本身的答案，性質與 `hgmd_lookup`
不同，但應於論文中分開說明。

---

## 9. 限制

1. **HGMD leave-one-out 會系統性高估**。查詢與類比庫同源，且 HGMD 的收錄密度與
   研究強度相關——密集區本身即為容易的區域。故所有結論均同時報告獨立的 PubTator
   數字。

2. **距離尺度未經正規化**。$\lambda = 150$ bp 為固定值，未依局部變異密度調整，可能
   將「資料密度」誤認為「證據強度」。文獻建議以未標註資料集約束最小搜尋半徑。

3. **band 的解析度過粗**。四個 band 跨度達 28 萬 bp；領域規範（ACMG 的 DCM 分類
   條件）雖仍使用 A-band 座標，但其作者已自承此舉「排除了其他可能相關的區域」，
   而研究文獻早已轉向 exon-level PSI。本方法在組織軸已改用 PSI，但疾病軸的統計
   先驗仍為 band 層級。

4. **`absent` 桶樣本稀薄**。補充表排除的 40 個 I-band 外顯子僅有 27 顆 HGMD 變異
   支撐，其心肌比例雖低至 22%（訊號強），但估計不穩。

5. **組織標籤本身可能低估 "Both"**。近期研究指出 DCM 患者所帶的單套 TTNtv 亦可
   造成輕度骨骼肌病變，HGMD 的 phenotype 欄記錄的是主要就診表型而非完整的組織
   受累，故「預測 Cardiac 而真值為 Both」的案例可能被低估為錯誤。

6. **強制作答移除了棄權能力**（見 §6）。較合適的校準形式是分層報告 likelihood
   ratio 而非 High／Low 標籤；PSI 組織 LR 已採此形式，但整體 confidence 尚未。

---

## 附：可重現性

```bash
# 主要結果與消融
conda run -n gm4 python research/benchmark_inference.py --truth hgmd
conda run -n gm4 python research/benchmark_inference.py --truth pubtator

# 全部圖表（含 §8.2、§8.3 消融）
conda run -n gm4 jupyter nbconvert --to notebook --execute \
    --inplace research/inference_benchmark_figures.ipynb
```

| 參數 | 值 | 位置 |
|---|---|---|
| 類比數 $k$ | 8 | `config.INFERENCE_TOP_K_ANALOGS` |
| 距離半衰期 $\lambda$ | 150 bp | `distance_scale` |
| LR 指數 $\alpha$ | 2.0 | `inclusion_lr_strength` |
| LR 上下限 | [0.25, 4.0] | `inclusion_lr_cap` |
| LR 樣本量門檻 | 15 | `inclusion_min_n` |
| high PSI 門檻 | 90% | `PSI_HIGH` |
