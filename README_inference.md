# Track B：無文獻證據時的 TTN 變異表型推論

> 簡報用說明文件。對應程式 `utils/mechanistic_inference.py`、
> 評估 `eval_inference_data.py`、圖表 `inference_benchmark_figures.ipynb`。

---

## 1. 要解決的問題

主流程（Track A）從 ClinVar 引用的文獻中擷取表型。它有一個結構性的失效條件：
**ClinVar 沒引文獻、或引的文獻沒提到這顆變異時，Track A 什麼都給不出來。**

TTN 上這是常態而非例外——本研究的 PubTator 評估集中，**49.2%** 的變異在 HGMD 有
同位置紀錄，另外一半完全查不到。

Track B 回答一個較弱但仍有臨床意義的問題：

> 沒有這顆變異的直接證據時，能否從它的**外顯子位置**、**該外顯子的組織使用率**、
> 以及**同基因內已知致病變異的統計**，推論最可能的疾病與受累組織？

輸出一律標記 `evidence_strength = inferred_no_literature`，與 Track A 的結果分軌呈現。

---

## 2. 系統概觀

三層決策，上層命中即回傳：

```
輸入變異 (chr2:pos)
   │
   ├─ 第 1 層  HGMD 收錄同一顆變異？  → 直接採用其標註（查表，非推論）
   │
   ├─ 第 2 層  落在 exon 344 / 364？  → 套用文獻確立的熱點表型
   │
   └─ 第 3 層  階層式收縮             → 疾病 + 組織          ← 主要方法
```

第 3 層另有一個備用實作（距離加權類比投票，§6），以 `predict_mode` 切換。

---

## 3. 資料來源

| 角色 | 檔案 | 規模 |
|---|---|---|
| 計數來源（類比庫） | `TTN_all_mutaion_20260703.csv`（HGMD） | 2,172 顆，其中 1,926 顆可定位到外顯子 |
| 外顯子座標 / band / PSI | `cardiodb.xlsx` | 364 個外顯子 |
| Ground truth | `pubtator_ttn_variants_clinvar_annotated.csv` | 419 顆（有座標且原文判 P/LP） |

**PSI**（percent spliced-in）是該外顯子在某組織的轉錄本中被保留的比例。
cardiodb 對 36 個外顯子沒有骨骼肌 PSI（exon 148、150、172–205，全在 I-band），
缺值視為 0。

---

## 4. 特徵：外顯子使用率分類

依 cardiodb 實測 PSI，以 90% 為高／低門檻，把 364 個外顯子分成四類：

| 類別 | 條件 | 外顯子數 | HGMD 中心肌表型比例 |
|---|---|---|---|
| `both-high` | 心肌、骨骼肌 PSI 皆 ≥ 90% | 186 | 79% |
| `skeletal-high` | 骨骼肌 ≥ 90%、心肌 < 90% | 106（全在 I-band） | 44% |
| `both-low` | 兩者皆 < 90% | 69 | 25% |
| `cardiac-high` | 心肌 ≥ 90%、骨骼肌 < 90% | 3 | 89% |
| | | | 全庫基準 **74%** |

`skeletal-high` 這 106 個外顯子正是文獻報告的「骨骼肌高 PSI、心肌低 PSI」集合，
且全部落在 I-band——替代性剪接集中的彈簧區。它們的心肌表型比例只有 44%，
遠低於 74% 的基準，是組織判別的主要訊號來源。

---

## 5. 第 3 層（主要方法）：階層式收縮

### 5.1 核心問題

最直接的做法是查「這個外顯子的已知變異都造成什麼病」。但 **1,926 顆變異分到
364 個外顯子，平均每個只有 5.3 顆**——直接數次數會極度不穩定，而且沒有變異的
外顯子完全無法估計。

解法是**借力**（borrowing strength）：從粗的層級開始估計，逐層往細的層級收縮。

### 5.2 記號

| 符號 | 意義 |
|---|---|
| $v$ | 查詢變異，位於基因體座標 $p$ |
| $\mathcal{K}$ | HGMD 計數來源（2,172 顆） |
| $\mathcal{K}^{(-p)}$ | 移除座標 $p$ 上全部變異後的 $\mathcal{K}$（防洩漏，見 §7.1） |
| $e(p),\ b(p)$ | 座標 $p$ 所屬的外顯子與結構 band |
| $\psi_H(e),\ \psi_S(e)$ | 外顯子 $e$ 的心肌／骨骼肌 PSI |
| $\mathcal{D}$ | 疾病標籤集合 |
| $\mathcal{T}$ | 組織標籤集合 $\{\text{cardiac},\ \text{skeletal},\ \text{both}\}$ |
| $\mathcal{H}$ | 熱點外顯子 $\{344\ (\text{僅 missense}),\ 364\}$ |

### 5.3 PSI 類別函數

$$c(e) = \begin{cases}
\texttt{both-high} & \psi_H(e) \ge 90 \ \wedge\ \psi_S(e) \ge 90\\[2pt]
\texttt{skeletal-high} & \psi_S(e) \ge 90 > \psi_H(e)\\[2pt]
\texttt{cardiac-high} & \psi_H(e) \ge 90 > \psi_S(e)\\[2pt]
\texttt{both-low} & \text{otherwise}
\end{cases}$$

### 5.4 四個由粗到細的分組

$$\begin{aligned}
L_0 &= \mathcal{K}^{(-p)} && \text{全庫}\\
L_1 &= \{r \in \mathcal{K}^{(-p)} : \mathrm{band}(r) = b(p)\} && \text{同 band}\\
L_2 &= \{r \in \mathcal{K}^{(-p)} : c(\mathrm{exon}(r)) = c(e(p))\} && \text{同 PSI 類別}\\
L_3 &= \{r \in \mathcal{K}^{(-p)} : \mathrm{exon}(r) = e(p)\} && \text{同外顯子}
\end{aligned}$$

計數 $n_\ell(d) = \bigl|\{r \in L_\ell : \mathrm{phen}(r) = d\}\bigr|$，$n_\ell = |L_\ell|$。

> **注意**：這四層是逐漸specific的**分組**，不是嚴格的巢狀集合——$L_2 \not\subseteq L_1$，
> 因為同一個 PSI 類別的外顯子可以跨越不同 band。遞迴只是把每一層當作下一層的先驗。

### 5.5 收縮遞迴

$$P_0(d) = \frac{n_0(d)}{n_0}, \qquad
P_\ell(d) = \frac{n_\ell(d) + m_\ell \, P_{\ell-1}(d)}{n_\ell + m_\ell}, \quad \ell = 1, 2, 3$$

其中 $m_1 = m_{\text{band}} = 8$、$m_2 = m_{\text{bucket}} = 4$、$m_3 = m_{\text{exon}} = 4$。

**展開式**：

$$P_3(d) = \cfrac{n_3(d) + m_3 \cdot
\cfrac{n_2(d) + m_2 \cdot
\cfrac{n_1(d) + m_1 \cdot \frac{n_0(d)}{n_0}}{n_1 + m_1}}
{n_2 + m_2}}
{n_3 + m_3}$$

組織用完全相同的遞迴，只把 $n_\ell(d)$ 換成組織計數 $n_\ell(t)$，記作 $Q_3(t)$。

> 標籤集合固定取自 $L_0$，所以在細層從未出現的疾病仍能透過先驗分到機率質量——
> 這正是收縮能處理「該外顯子完全沒有變異」的原因。

### 5.6 完整決策函數

三層依序判斷，上層命中即回傳：

$$\hat{y}(v) = \begin{cases}
\bigl(\mathrm{phen}(r^*),\ \mathrm{tis}(r^*)\bigr)
& \exists\, r^* \in \mathcal{K}:\ (\mathrm{pos},\mathrm{ref},\mathrm{alt})(r^*) = (p,\mathrm{ref},\mathrm{alt})(v)\\[6pt]
\bigl(D_h,\ T_h\bigr)
& e(p) \in \mathcal{H}\ \wedge\ \text{變異類型符合該熱點}\\[6pt]
\Bigl(\arg\max\limits_{d \in \mathcal{D}} P_3(d),\ \ \arg\max\limits_{t \in \mathcal{T}} Q_3(t)\Bigr)
& \text{otherwise}
\end{cases}$$

- **層 1**：HGMD 完全命中，$(\mathrm{phen}, \mathrm{tis})$ 為該筆紀錄的標註。屬查表而非推論，
  評估時必須排除。
- **層 2**：$(D_h, T_h)$ 為熱點的文獻表型（exon 344 → HMERF；exon 364 → TMD／LGMD2J）。
- **層 3**：階層式收縮，絕大多數變異由此決定。

### 5.7 $m$ 的意義

$m_\ell$ 是**等效樣本數**：該層要累積到 $m_\ell$ 顆變異，其自身計數才與上一層的
先驗等重。

- $m_\ell \to 0$：完全相信該層自己的計數（樣本少時會過擬合）
- $m_\ell \to \infty$：完全忽略該層，直接沿用上一層的估計

預設 $m_{\text{exon}} = 4$、$m_{\text{bucket}} = 4$、$m_{\text{band}} = 8$。

### 5.8 信心與差異診斷

$$\mathrm{conf}(v) = \begin{cases}
\text{High} & P_3(\hat{d}) \ge 0.6 \\
\text{Moderate} & 0.4 \le P_3(\hat{d}) < 0.6 \\
\text{Low} & P_3(\hat{d}) < 0.4
\end{cases}
\qquad
\mathrm{diff}(v) = \operatorname*{top\text{-}2}_{d \ne \hat{d}} P_3(d)$$

---

## 6. 備用方法：距離加權類比投票

保留作為對照與 ensemble 的一半。從 HGMD 取同 band、同變異類型的 $k = 8$ 顆最近變異，
加權投票。類比案例 $a$ 的票重為五項連乘：

$$w(a) = \underbrace{0.5^{\,d_a/\lambda}}_{\text{距離}}
\cdot \underbrace{\tau(\text{tier}_a)}_{\text{定位層級}}
\cdot \underbrace{\mathrm{LR}_{t(a)}(e)}_{\text{PSI 組織}}
\cdot \underbrace{\rho(a)}_{\text{RANKSCORE}}
\cdot \underbrace{\nu(a)}_{\text{籠統病名}}$$

其中 $d_a = |\mathrm{POS}_a - \mathrm{POS}_{\text{query}}|$，半衰期 $\lambda = 150$ bp；
$\tau = \{$`exon+type`: 1.0, `band+type`: 0.6, `exon`: 0.4, `band`: 0.2$\}$。

**PSI 組織 likelihood ratio**（查詢變異所在外顯子 $e$、其所屬類別 $b(e)$）：

$$\mathrm{LR}_t(e) = \mathrm{clip}\left(
\left[\frac{\bigl(n^{(-e)}_{b,t}+1\bigr)\big/\bigl(N^{(-e)}_{b}+3\bigr)}
{\bigl(n^{(-e)}_{t}+1\bigr)\big/\bigl(N^{(-e)}+3\bigr)}\right]^{\alpha},
\ \tfrac{1}{c},\ c \right)$$

- 分子分母皆採 Laplace 平滑（+1／+3，3 為組織類別數），避免某類為 0 時 $\mathrm{LR}=0$
  把該類的類比完全消音
- 上標 $(-e)$：估計時扣除查詢變異所在外顯子的全部變異（防循環）
- $\alpha = 2$（訊號放大）、$c = 4$（夾住極端值）
- 類別內樣本數 $N^{(-e)}_b < 15$ 時回傳 $\mathrm{LR} = 1$（不表態）

疾病取 $\arg\max_D \sum_{a:\,\mathrm{phen}(a)=D} w(a)$。

---

## 7. 評估設計

### 7.1 洩漏處理

PubTator 的 419 顆變異中有 **206 顆（49.2%）在 HGMD 有同位置紀錄**。若不處理，
它們會以 $d = 0$ 的類比進入投票（距離權重 $0.5^0 = 1$，最大值），或直接被計入該
外顯子的統計——等於把答案讀出來。

**所有評估都移除查詢位置的全部 HGMD 變異。** 影響實測如下：

| | hit | specific |
|---|---|---|
| 保留同位置（會洩漏） | 58.7% | 44.6% |
| 移除同位置（採用） | 55.1% | 38.2% |

### 7.2 四個指標

疾病名不是受控詞彙——同一個病可能寫成 `Dilated cardiomyopathy`、`Cardiomyopathy, dilated`
或只寫 `Cardiomyopathy`，而 `Titinopathy` 這種上位詞又涵蓋好幾個具體診斷。因此
預測與答案都先映射成 **canonical label** 再比對。

比對結果分成幾個層級（`utils/disease_normalizer.compare_disease`）：

| level | 意義 |
|---|---|
| `exact` | canonical label 完全相同（`DCM` vs `Cardiomyopathy, dilated`） |
| `compatible` | 上位詞 ↔ 特定病名（`Cardiomyopathy` vs `DCM`、`Titinopathy` vs `TMD`） |
| `differential` | 主答案沒中，但列在差異診斷裡 |
| `same_tissue` | 疾病不同但組織相同（`DCM` vs `HCM`） |
| `miss` | 組織也不同 |

**hit** 與 **specific** 只採計 `exact` 與 `compatible`。

---

#### hit —— 疾病命中率

$$\text{hit} = \frac{\bigl|\{v : \mathrm{level}(v) \in \{\texttt{exact},\ \texttt{compatible}\}\}\bigr|}{|V|}$$

**刻意寬鬆**：答「Cardiomyopathy」而真值是 DCM 也算對。理由是這在臨床上仍有用——
把變異放進正確的大類就已經縮小了範圍。

但這也代表系統可以靠**答得含糊**來灌高分。

#### specific —— 命中且具體

$$\text{specific} = \frac{\bigl|\{v : \mathrm{level}(v) \in \{\texttt{exact},\ \texttt{compatible}\}\ \wedge\ \neg\,\mathrm{vague}(\hat{d}_v)\}\bigr|}{|V|}$$

其中 $\mathrm{vague}(\cdot)$ 判斷預測本身是否為籠統上位詞，目前的上位詞集合是
`{TITINOPATHY, NEUROMUSCULAR, CARDIAC_OTHER, MYOPATHY, CARDIOMYOPATHY}`。

**這個指標的存在就是為了抓 hit 的漏洞。** `hit − specific` 即「靠上位詞作答」的比例。
兩者差距越大，代表系統越常用模糊答案換分數。

> 實例：距離加權類比投票的 hit 55.1%、specific 38.9%，差距 16.2 個百分點——
> 表示它有相當比例的「命中」是靠答 Cardiomyopathy 這類詞得來的。階層式收縮的
> 差距只有 6.9 個百分點（55.6% / 48.7%），答得具體得多。

#### macro —— 逐類別平均

$$\text{macro} = \frac{1}{|\mathcal{C}|} \sum_{c \in \mathcal{C}}
\frac{\bigl|\{v \in c : \mathrm{level}(v) \in \{\texttt{exact},\ \texttt{compatible}\}\}\bigr|}{|c|},
\qquad \mathcal{C} = \{c : |c| \ge 3\}$$

先算每個答案類別各自的命中率，再對類別取平均（樣本數 < 3 的類別排除）。

**為什麼必要**：TTN 的表型分佈極度偏斜，DCM 在 PubTator 集中佔 184/419（44%）。
一個永遠答 DCM 的預測器 hit 可以到 50.1%，但對骨骼肌表型完全無用。macro 讓每個
類別等權，是唯一能反映「罕見表型有沒有被抓到」的指標。

> 多數類別基準的 macro 只有 21.7%，遠低於它 50.1% 的 hit——差距正好量化了
> 「靠多數類別得分」這件事。

#### tissue —— 組織正確率

$$\text{tissue} = \frac{\bigl|\{v \in V_T : \mathrm{tmatch}(\hat{t}_v,\ t_v) \in \{\texttt{exact},\ \texttt{partial}\}\}\bigr|}{|V_T|}$$

$V_T$ 是有組織標註的子集。比對採**部分給分**：

| 預測 \ 真值 | Cardiac | Skeletal | Both |
|---|---|---|---|
| **Cardiac** | exact | miss | partial |
| **Skeletal** | miss | exact | partial |
| **Both** | partial | partial | exact |

`Both` 與 `Cardiac`／`Skeletal` 互為 partial，兩者都算對。

**為什麼要分開算**：疾病層與組織層的錯誤性質不同。把 DCM 誤判成 LVNC 是心肌軸
內部的錯誤；把骨骼肌病誤判成心臟病則是完全不同的臨床後果。組織軸也是 PSI 特徵
最直接影響的地方，分開量才看得出 PSI 有沒有用。

> **但要注意**：組織標籤是由疾病名推導的（`tissue_of_labels(canonical_labels(d))`），
> 所以 `tissue` 與 `hit` 有相關性，不是完全獨立的評估軸。見 §9。

---

## 8. 結果（PubTator，n = 419）

| 方法 | hit | specific | macro | tissue |
|---|---|---|---|---|
| 多數類別基準（一律答 DCM／Cardiac） | 50.1% | **50.1%** | 21.7% | 70.9% |
| 距離加權類比投票 | 55.1% | 38.9% | **35.9%** | **76.2%** |
| **階層式收縮** | **55.6%** | 48.7% | 32.3% | 74.7% |

### 8.1 逐層貢獻

| 層級 | hit | specific | macro | tissue |
|---|---|---|---|---|
| 只用全庫 | 52.0% | 52.0% | 27.1% | 72.8% |
| ＋band | 52.0% | 52.0% | 27.1% | 72.8% |
| ＋PSI 類別 | 55.1% | 49.6% | 30.4% | 74.0% |
| ＋外顯子（完整） | 55.6% | 48.7% | 32.3% | 74.7% |

**band 層完全沒有貢獻**——四個指標一個都沒動。增益全部來自 PSI 類別（+3.1 hit、
+3.3 macro）與外顯子（+0.5 hit、+1.9 macro）。原因是 TTN 的 band 分佈極度偏斜
（A-band 佔 61%），band 層的疾病分佈幾乎等同全庫分佈。

### 8.2 參數敏感度

$m_{\text{exon}}$ 與 $m_{\text{bucket}}$ 在 1～16 之間掃描，PubTator 命中率**全距
僅 1.0%**（54.7%–55.6%）。結果不依賴調參。

---

## 9. 已知限制

1. **`specific` 仍略低於多數類別基準**（48.7% vs 50.1%）。階層式已比類比投票的
   38.9% 好很多，但尚未跨過「一律答 DCM」這條線。

2. **macro 與 tissue 輸給類比投票**。階層式預測的是該外顯子的眾數疾病，罕見表型
   抓得較差；類比投票靠最近鄰反而保留了少數表型的訊號。兩者互補，合併（以階層式
   機率為先驗、類比票為似然）尚未實作。

3. **band 層冗餘**，可以移除以簡化模型（見 §8.1）。

4. **HGMD 收錄密度不均**：研究得多的區域變異密集，估計自然較準。這是資料本身的
   性質，非方法可解。

5. **組織標籤由疾病名推導**：`tissue` 指標因此與 `hit` 有相關性，不是完全獨立的
   評估軸。

---

## 10. 參數總表與重現

| 參數 | 值 | 位置 |
|---|---|---|
| $m_{\text{exon}}$ | 4.0 | `m_exon` |
| $m_{\text{bucket}}$ | 4.0 | `m_bucket` |
| $m_{\text{band}}$ | 8.0 | `m_band` |
| high PSI 門檻 | 90% | `PSI_HIGH` |
| 類比數 $k$ | 8 | `INFERENCE_TOP_K_ANALOGS` |
| 距離半衰期 $\lambda$ | 150 bp | `distance_scale` |
| LR 指數 $\alpha$ | 2.0 | `inclusion_lr_strength` |
| LR 上下限 $c$ | 4.0 | `inclusion_lr_cap` |
| LR 最小樣本 | 15 | `inclusion_min_n` |
| 熱點外顯子 | 344（僅 missense）、364 | `EXON_HOTSPOTS` |

```bash
# 評估
conda run -n gm4 python eval_inference_data.py

# 方法比較
conda run -n gm4 python compare_exon_psi_methods.py

# 全部圖表
conda run -n gm4 jupyter nbconvert --to notebook --execute \
    --inplace inference_benchmark_figures.ipynb
```

```python
# 直接呼叫
from utils.mechanistic_inference import MechanisticInference
mi = MechanisticInference(..., predict_mode="hierarchical",
                          m_exon=4.0, m_bucket=4.0, m_band=8.0)
result = mi.infer(extractor=None, variant_info="c.98671-2A>G")
```
