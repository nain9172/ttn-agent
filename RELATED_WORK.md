# 相關文獻：距離加權類比推理與變異表型預測

對照對象（本專案 Track B 推理）：從 HGMD 取同 band、同變異類型的已知致病變異 → 依 genomic 距離加權投票 → 預測疾病與組織，並以 PubTator 文獻標註當 ground truth 驗證。

檢索日期 2026-08-24；§5–§9 為 2026-08-25 補充（限 2020 年後文獻）。實驗數據對應 `test_inference_strategies.ipynb`。

---

## 最接近的一篇：PathProx

**Comprehensive Analysis of Constraint on the Spatial Distribution of Missense Variants in Human Protein Structures** — [10.1101/109652](https://doi.org/10.1101/109652)

PathProx 定義了一個「致病鄰近度」指標：用 NeighborWeight 把蛋白質 **3D 歐氏距離**轉成權重（下界 8 Å、上界 24 Å），把 VUS 到已知致病變異與中性變異的加權平均距離相減當分數，並用 leave-one-out 交叉驗證評估，中位 ROC AUC 0.75。這和本專案 `utils/mechanistic_inference.py` 的 `_distance_weight`（每 150 bp 減半）是同一個構造，差別只在距離定義於蛋白質 3D 空間而非 genomic bp。

關鍵是它報告的失敗條件：

> 不先篩掉「致病變異沒有群聚」的基因就直接套用，中位 ROC AUC 從 0.75 掉到 0.55（N = 398）。作者結論是**有預測力的是「與致病變異群集的鄰近度」，不是「與個別致病變異的鄰近度」**。

本專案剔除 500 bp 內類比後準確率不掉，量到的是同一個現象——這是這類方法的已知性質，不是實作缺陷。

| 指標 | 值 |
| --- | --- |
| PathProx LOO 中位 ROC AUC（有群聚的基因） | 0.75 |
| 不篩群聚時 | 0.55（≈ 隨機） |
| 距離權重上下界 | 8–24 Å（本專案類比：150 bp 半衰期） |

---

## 方法元件對照

| 元件 | 本專案實作 | 文獻上的名稱／對應概念 |
| --- | --- | --- |
| 檢索條件 | 同 band（4 段）＋ 同變異類型 | PM1「mutational hotspot／critical domain」；TTN 領域已改用 exon-level PSI |
| 相似度 | genomic bp 距離，指數衰減半衰期 150 bp | PathProx 的 NeighborWeight（3D Å）、SPDV 的 K 近鄰平均距離、MutScore 的位置群聚特徵；另一族是序列比對而非距離：HomolVar 的 Pfam 對齊、Paralogue Annotation 的同源等價位置（見 §5） |
| 決策 | 距離加權多數決 → 疾病 ＋ 組織 | kNN 加權投票（P-KNN、VUS.Life）；但它們幾乎都只輸出「是否致病」 |
| 驗證 | PubTator 文獻標註當答案 ＋ HGMD leave-one-out | LitVar／PubTator 抽取的 variant–disease 關聯；PM3-Bench、VDAL 這類文獻證據基準 |

---

## 1. 距離加權近鄰投票（方法核心）

這條線很成熟，但幾乎全部都在預測「是否致病」這個二元標籤，而不是「是哪一種疾病」。

| 論文／工具 | 做法 | 和本專案的差異 |
| --- | --- | --- |
| **PathProx**<br>[10.1101/109652](https://doi.org/10.1101/109652) | 3D 距離的 NeighborWeight 加權，致病 − 中性鄰近度差值，LOO 驗證（AUC 0.75） | 距離在蛋白質 3D 空間；二元輸出。已明確報告「群集 > 個別鄰居」的結論 |
| **P-KNN**<br>[10.1101/2025.09.24.678417](https://doi.org/10.1101/2025.09.24.678417)<br>[GitHub](https://github.com/Brandes-Lab/P-KNN) | 在多個預測工具分數構成的空間裡做加權 kNN，輸出可用於 ClinGen 框架的後驗機率 | 距離定義在「預測分數空間」而非基因座標。值得參考的設計：用一份未標註資料集正規化最小搜尋半徑，避免鄰域過窄 |
| **RyR1 SPDV**<br>[10.1155/humu/1834898](https://doi.org/10.1155/humu/1834898) | 某殘基到 K 個最近「已知致病殘基」的平均 3D 距離，當單一特徵（SPDV_4 最佳，AUC 0.762） | K 是超參數、掃法與本專案相同；同樣是單基因、二元輸出。作者自陳這類做法在罕病臨床框架裡很少被使用 |
| **MutScore**<br>[PMC8948164](https://pmc.ncbi.nlm.nih.gov/articles/PMC8948164/) | 把基因內位置群聚轉成特徵餵給 random forest；840 個 ClinVar 基因中 46% 有顯著非隨機分布 | 提供基準線：位置群聚是普遍現象，但它把位置當特徵而非直接投票 |
| **3DNT（3D 鄰域檢定）**<br>[10.64898/2026.05.29.26354366](https://doi.org/10.64898/2026.05.29.26354366) | 全基因體掃描 case/control 變異在 3D 結構上的局部富集，Mendelian 疾病找到 872 個基因有顯著群集 | 目標是「找出群集區域」而非預測單顆變異；可用於取代 band 當檢索條件 |

---

## 2. ACMG PM1／PM5：同一套邏輯的臨床規則版

PM1（變異落在突變熱點／關鍵功能域）與 PM5（同一殘基已有其他致病 missense）就是「同位置／鄰近位置類比」的規則化版本。這條線最有用的地方是它已經把證據強度量化成 likelihood ratio。

| 來源 | 重點 | 為什麼直接相關 |
| --- | --- | --- |
| Navigating the nuances of clinical sequence variant interpretation in Mendelian disease<br>[10.1038/s41436-018-0100-y](https://doi.org/10.1038/s41436-018-0100-y) | PM1、PS1、PM5 都在主張「變異改動了蛋白質的重要區域」，同時套用會人為放大證據（double-counting） | 本專案量到「exact-match 查表只比強迫推理好 2.8 個百分點、macro 完全相同」，正是這種非正交性的實測 |
| PM5 證據強度校準（ICR）<br>[全文 PDF](https://repository.icr.ac.uk/server/api/core/bitstreams/6a7c6fa5-96d2-4721-9bd6-54dcdeb118b2/content) | 把 PM5 量成 pLR = 8.4（95% CI 7.2–…），加上「比同位置已知變異更 damaging」條件後升到 16.3（10.6–24.9）；並建議 PM5 與 PM1 不要併用 | 這是 confidence 校準該走的形式：報 likelihood ratio，而不是直接標 High／Low |
| ClinGen CSpec 規格登錄<br>[GN008](https://cspec.genome.network/cspec/ui/svi/doc/GN008) | VCEP 只在**逐一列舉的殘基**上把 PM1 升成 PM1_strong，並規定用了 PM1 就不得再用 PM5／PS1 | 領域專家對「熱點有多大」的答案是人工逐殘基認定，不是一個全域的 bp 半衰期——支持「單一 150 bp 尺度是錯的參數化」這個結論 |

---

## 3. TTN 專屬：這個領域已經從 band 走到 exon-level PSI

把 genomic bp 換成外顯子／domain 層級，不只合理，而是 TTN 領域十年前就走過的路，且已寫進 ACMG 的 DCM 判讀規則。

| 論文 | 結論 | 意義 |
| --- | --- | --- |
| Roberts et al. 2015, Sci Transl Med<br>[10.1126/scitranslmed.3010134](https://doi.org/10.1126/scitranslmed.3010134) | TTNtv 的致病性由「該外顯子是否被成人心臟轉錄本收錄」決定，用 PSI 量化；PSI ＋ 位置才能區分致病與良性 | PSI 就是需要的那個「更好的相似度定義」，且有現成的 cardiodb.org/titin 註釋 |
| Titin Cardiomyopathy, Emerging Evidence（2026 回顧）<br>[10.1007/s11886-025-02309-5](https://doi.org/10.1007/s11886-025-02309-5) | 過去認為越靠 A-band 越嚴重，但目前最大規模的研究並未發現不同 TTN domain 之間預後有顯著差異 | 領域自己的結論就是 band 位置預測力很弱——與「band 過濾形同虛設」一致 |
| 2025 預印本<br>[10.1101/2025.11.02.25338804](https://doi.org/10.1101/2025.11.02.25338804) | 長讀定序的 PSI-LR 中 PSI ≥ 0.75 與進展到重度心衰相關（⚠️ 未能取得全文核實，medRxiv 擋抓取；該 DOI 目前檢索到的標題為「Exon Utilization Improves Risk Stratification for Advanced Heart Failure in Titin Cardiomyopathy」，與此處原記的標題不符，引用前需自行確認版本） | 直接可用的替代條件；也提醒 A-band 有很多高 PSI 外顯子其實沒有 DCM 富集 |
| Savarese et al. 2020, J Neuromuscul Dis 7(2):153–166<br>*Improved Criteria for the Classification of Titin Variants in Inherited Skeletal Myopathies*<br>[10.3233/JND-190423](https://journals.sagepub.com/doi/full/10.3233/JND-190423) | TTN 骨骼肌表型的 PM1 只認兩個外顯子：exon 344（HMERF，FN3-119／A150 domain）、exon 364（TMD／Udd myopathy 及遲發型隱性 titinopathy），強度 Moderate。**限定僅適用於 missense 變異**，TTNtv 走另一套規則（已於 2026-08-25 核實原文） | 「組織」預測軸在領域裡的做法是外顯子級、手工認定的兩個熱點，粒度比 band 細得多。那個 missense-only 限定本專案目前沒做：`retrieve_analogs` 有依 `vgroup` 過濾，但 `band_phenotype_prior` 把 truncating 與 missense 混在一起算比例 |
| A comprehensive framework for the interpretation of TTN missense variants（Genome Med 2026）<br>[10.1186/s13073-026-01605-1](https://doi.org/10.1186/s13073-026-01605-1) | 最新的 TTN missense 判讀框架；分析時特意排除 exon 344／364，因為那兩個外顯子的表型已確立 | 目前最該對照的領域基準線，也是 missense 類別的現成比較對象 |

### ⚠️ 更正（2026-08-25）：規範仍用 A-band，改用 PSI 的是回報政策

此節原本寫「ACMG 的 DCM 條件已從必須在 A-band 改成高 PSI 外顯子」，並掛在 2025 預印本下。核實後這是把兩份不同文件混為一談：

| 文件 | 實際內容 |
| --- | --- |
| **Morales et al. 2020, Circ Genom Precis Med**<br>*Variant Interpretation for Dilated Cardiomyopathy: Refinement of the ACMG/ClinGen Guidelines for the DCM Precision Medicine Study*<br>[PMC8070981](https://pmc.ncbi.nlm.nih.gov/articles/PMC8070981/) | 這才是 ACMG/ClinGen 的 DCM **分類**規範。用 **A-band 座標**（chr2:179,483,218–179,400,709，GRCh37／NM_001267550.1），強度 **PVS1_Strong**；**完全未提 PSI 或 constitutive exon**。TTN missense 另設 `BP1_Stand_Alone`（單獨即可判為良性）。原文自承限制：restricting pathogenicity to the A-band「excludes potentially relevant regions that can be identified with more granular analyses」 |
| **ACMG SF v3.2**, Miller et al., Genet Med 2023<br>[PMC10524344](https://pmc.ncbi.nlm.nih.gov/articles/PMC10524344/) | TTN 納入次要發現清單時，明訂只有 frameshift／nonsense 與影響剪接的變異、**且落在高 PSI 外顯子**者才需評估並回報。這是**回報政策**，不是分類條件 |

**這個更正對本專案的論證是加分而非扣分**：研究文獻十年前就證明 PSI 優於 band 位置、ACMG 的回報政策 2023 年也改用 PSI，但 DCM 的**分類規範至今仍停在 A-band**，且規範作者自己承認這是限制。本專案量到的「band 過濾形同虛設」正好落在這個缺口上——是對一份仍在使用的規範提出量化證據，而不是重述已知結論。

---

## 4. 用文獻挖掘當 ground truth

| 論文／工具 | 做法 | 關聯 |
| --- | --- | --- |
| LitVar 2.0, Nat Genet 2023<br>[10.1038/s41588-023-01414-x](https://doi.org/10.1038/s41588-023-01414-x) | tmVar 3.0 ＋ PubTator，從 PubMed／PMC 全文與補充材料抽出約 1,400 萬個變異提及並正規化 | 本專案 ground truth 的上游基礎設施；補充材料常是變異資訊最密集的地方，值得納入 |
| PubMind, Nat Commun 2026<br>[全文](https://www.nature.com/articles/s41467-026-76834-4) | LLM 從 4,100 萬篇摘要抽 variant–function–disease，建成 130 萬變異的資料庫；只有約 10% 與 ClinVar 重疊，其中 >80% 致病性標註一致 | 本專案的 5.6% HGMD∩PubTator 重疊率與 63.9% 病名相符率，在更大規模上是同一個現象——文獻與 curation 互補而非重複 |
| AutoPM3, Bioinformatics 2025<br>[10.1093/bioinformatics/btaf382](https://doi.org/10.1093/bioinformatics/btaf382)<br>[GitHub](https://github.com/HKU-BAL/AutoPM3) | 開源 LLM ＋ RAG ＋ Text2SQL 從全文（含表格）抽 PM3 證據；PM3-Bench 1,027 對，變異命中率 86.1%、in trans recall 72.5% | LLM 抽取步驟最直接的對照，也示範了該怎麼建基準與做 ablation |
| VDAL／VDANet, ACM-BCB 2023<br>[10.1145/3584371.3612995](https://doi.org/10.1145/3584371.3612995) | 用 DisGeNET 的 variant–disease 標籤 ＋ PubTator 標註建成文件級關聯抽取基準 | 若要把「文獻標註當答案」寫進論文，這是可引用的先例 |

### LLM ＋ 檢索（Track B 的 LLM 路徑）

| 系統 | 做法 | 差異 |
| --- | --- | --- |
| ChatTogoVar, Hum Genome Var 2026<br>[10.1038/s41439-026-00344-4](https://doi.org/10.1038/s41439-026-00344-4) | 查 TogoVar 取得該變異的結構化紀錄，塞進 prompt，要求 LLM 引用資料庫證據 | 檢索的是**同一顆變異**的紀錄；本專案檢索的是**類比變異**——這是實質差異 |
| VarChat, Bioinformatics 2024<br>[10.1093/bioinformatics/btae183](https://doi.org/10.1093/bioinformatics/btae183) | 生成式 AI 助理，為人類基因變異產生文獻摘要 | 同上，變異層級檢索 |

---

## 5. 同源位置類比轉移：把「距離」從 bp 換成 domain 對齊

這條線與本專案是同一個想法，但相似度定義在**蛋白質家族的序列比對**上：一顆變異在某位置已知致病，就用來註釋臨床資料為零的等價位置。它是 `retrieve_analogs` 目前最有力的替代方案。

| 論文 | 做法與數字 | 對本專案的意義 |
| --- | --- | --- |
| **Conserved missense variant pathogenicity and correlated phenotypes across paralogous genes**<br>Genome Biology 2025<br>[10.1186/s13059-025-03663-x](https://pmc.ncbi.nlm.nih.gov/articles/PMC12235991/) | 59 萬顆致病 ＋ 190 萬顆族群變異映到 9,928 基因、2,054 個 paralogous family，可判讀殘基數增加 **5.1 倍**；「同源基因等價位置有致病變異」的 **pLR = 13.0**；在 10 個電位門控鈉離子通道基因、22 個專家審訂疾病上，依 3D 空間位置找出 **cross-paralog correlated phenotypes**，且把表型納入類比選擇會改善分類 | 三重意義：(1) 5.1 倍＝Track B「擴大可判讀範圍」的量化版；(2) pLR 13.0 比 PM5 的 8.4–16.3 更適合當校準模板，因為量的就是**類比位置**而非同一位置；(3) **它已經佔掉「用位置鄰近度預測是哪一種疾病」這一格**——見下方修正 |
| **HomolVar**<br>Bioinformatics 2025<br>[btaf305](https://academic.oup.com/bioinformatics/article/41/5/btaf305/8131454) | 用 **Pfam family alignment** 找等價位置，分 strict pair（同一置換）與 similar pair（BLOSUM62 相近）。292 個 Pfam domain、2,976 致病 ＋ 17,555 非致病同源變異，同源位置致病性一致率 **95%**（序列相似度 >30% 時 98%），對比 27 個傳統預測工具的 0.34–0.79 | 可直接照抄的機制，且**不需要 3D 結構**。TTN 的 Ig domain 屬 Pfam I-set（PF07679）、FN3 屬 PF00041，全基因約 280 個重複 domain——把 HGMD 致病變異投影到 Pfam HMM match state，相似度就從 genomic bp 變成 domain 內等價殘基位置 |
| para-PS1／PM5<br>medRxiv 2023<br>[10.1101/2023.08.22.23294353](https://www.medrxiv.org/content/10.1101/2023.08.22.23294353v1) | 把 ClinVar／HGMD／gnomAD 的 250 萬顆變異依基因家族對齊，延伸後的 para-PS1／PM5 讓可判讀胺基酸增加 **3.6 倍** | 證明這套邏輯已被寫成 ACMG 條件的延伸形式 |
| Paralogue Annotation<br>bioRxiv 2023<br>[10.1101/2023.08.07.552236](https://www.biorxiv.org/content/10.1101/2023.08.07.552236v1.full) | 跨蛋白質家族轉移致病註釋，精確度優於當時的 in silico 工具，且廣泛適用於多個蛋白質家族 | 這條線的方法論總述 |

**為什麼這對 TTN 特別重要**：本專案量到 bp 距離幾乎沒有加值。一個具體解釋是——同一個 Ig domain 折疊上的等價殘基，在基因座標上可能相隔數萬 bp，所以 bp 距離量的是「HGMD 收錄密度」而不是「結構等價性」。Pfam 對齊把這兩者分開。

---

## 6. 表型（而非致病性）預測

| 論文 | 做法 | 差異 |
| --- | --- | --- |
| **PheMART**<br>Nature Biomedical Engineering 2026<br>[10.1038/s41551-026-01636-4](https://www.nature.com/articles/s41551-026-01636-4)<br>[GitHub](https://github.com/celehs/PheMART) | 對 missense 變異預測 **4,179 個表型**；contrastive learning 把變異與表型投影到同一低維度量空間，特徵含 protein domain、protein language model、PPI、醫學知識圖譜、EHR | 動機論述與本專案一字不差：二元致病／良性標籤會掩蓋同一基因內的表型異質性。粒度遠細於 V2P 的 23 個 HPO 頂層表型，**應取代 V2P 當主要對照**。同樣沒有使用「與已知致病變異的鄰近度」 |
| V2P<br>[10.1038/s41467-025-66607-w](https://doi.org/10.1038/s41467-025-66607-w) | 多標籤輸出 23 個 HPO 頂層表型機率，訓練資料為 HGMD 252,125 顆致病變異 | 基因／蛋白質／網路層級特徵，無鄰近度 |

---

## 7. Evo2：推論輸入的來源與已知弱點

`evo2_confidence_bucket()` 與 prompt 中的 Evo2 分數目前無出處，需補：

| 論文 | 重點 | 用途 |
| --- | --- | --- |
| **Evo 2**, Nature 2026<br>[10.1038/s41586-026-10176-5](https://www.nature.com/articles/s41586-026-10176-5) | 9 兆 bp 訓練的基因體基礎模型，以 reference／variant 對數概度差做 zero-shot 變異效應預測 | 分數來源的正式引用 |
| **Biological Blind Spots in Evo2 Variant-Effect Prediction**<br>bioRxiv 2026<br>[連結](https://www.biorxiv.org/content/10.64898/2026.03.10.710786v1.full) | 直接質疑 Evo2 zero-shot 致病性預測的臨床成熟度 | 支持本專案把 Evo2 當「弱訊號、須與類比證據併用」，也是 prompt 裡「splice 變異 Evo2 分數異常低、不得因此降權」那條 caveat 的依據 |

---

## 8. 基準偏差：為何 HGMD 自我驗證會高估

本專案用 HGMD leave-one-out 驗證，查詢與知識庫兩側落在同一基因、甚至同一 band，會系統性高估近鄰法的表現；而 HGMD 收錄密度本身與研究強度相關。這是已知的方法學問題，該在論文裡先講明。

| 論文 | 重點 |
| --- | --- |
| **Updated benchmarking of variant effect predictors using deep mutational scanning**<br>Mol Syst Biol 2023<br>[10.15252/msb.202211474](https://www.embopress.org/doi/full/10.15252/msb.202211474) | 用 DMS 而非臨床標註 benchmark 55 個 VEP，以最小化 data circularity；核心論點是「拿臨床觀察當基準本身就會引入偏差」 |
| **Benchmarking VEPs by their ability to infer human traits**<br>Genome Biology 2024<br>[10.1186/s13059-024-03314-7](https://link.springer.com/article/10.1186/s13059-024-03314-7) | 提供不依賴臨床標註的替代評估設計 |

這與 §1 的 P-KNN 半徑正規化建議是同一個問題的兩面：固定 bp 尺度會把「資料密度」誤當成「證據強度」。

---

## 9. LLM 端：analogical prompting 與強制 commit 的代價

`build_prompt()` 要求模型必須給出單一疾病、禁止回答 Uncertain，`_backfill_if_uncertain()` 再用多數決填掉殘餘的 Uncertain。

| 論文／系統 | 重點 | 關聯 |
| --- | --- | --- |
| **Large Language Models as Analogical Reasoners**<br>ICLR 2024<br>[arXiv:2310.01714](https://arxiv.org/abs/2310.01714) | analogical prompting：讓模型自生類比範例再解題，優於 0-shot CoT 與人工 few-shot CoT | 本專案是它的 retrieval 版本——**類比案例來自 HGMD 而非模型自生**，這正是避免捏造的關鍵設計差異，該寫成貢獻 |
| **VarChat**, Bioinformatics 2024<br>[10.1093/bioinformatics/btae183](https://academic.oup.com/bioinformatics/article/40/4/btae183/7641533) | 作者自陳在無文獻變異上會 hallucinate；明確標示答案來自文獻或來自模型知識 | 與本專案 `evidence_strength: inferred_no_literature` ＋ `source` 欄位的分軌呈現是同一設計原則，應寫成明確貢獻而非實作細節 |
| **Benchmarking GPT-4o／Llama 3.1／Qwen 2.5 for cancer genetic variant classification**<br>npj Precision Oncology 2025<br>[10.1038/s41698-025-00935-4](https://www.nature.com/articles/s41698-025-00935-4) | LLM 變異分類的量化基準 | 「強制 commit 是否造成過度自信」這個必然被問到的問題，需要這類基準當對照 |

移除 abstention 選項等於把校準責任完全交給下游——正確的解法是 §2 的 likelihood ratio，而不是保留 Uncertain。

---

## 尚未被做過的部分（2026-08-25 修正）

原本的宣稱是「多類別表型預測 × 鄰近度」這一格全空。**這已不再成立**：Genome Biology 2025（§5）在 10 個鈉離子通道基因、22 個疾病上，依 3D 空間位置做出了 cross-paralog correlated phenotypes，等於已經用位置鄰近度預測「是哪一種疾病」。

仍然空著的是更窄的一格：**基因內部（intra-gene）、跨重複 domain 的同源位置類比，用於多類別表型預測**。既有工作的分佈是——

- PathProx／SPDV：用了鄰近度，但只做二元致病／良性判斷。
- PheMART／V2P：做多類別表型，但特徵在基因／蛋白質／網路層級，沒有鄰近度。
- Genome Biology 2025／HomolVar：用了同源位置類比，但相似度定義在**跨基因**的 paralog 對齊上。

TTN 的特殊之處在於它單一基因內就有約 280 個重複的 Ig／FN3 domain，構成一個基因內部的「蛋白質家族」——這使得 paralogue annotation 的邏輯可以完全在基因內部運作，同時保有多類別（DCM／HCM／骨骼肌病）的預測目標。這是目前的空白。

反面證據也要記住：LMNA 的研究把致病取代位點全部映到結構上，結論是[看不到明顯的基因型—表型對應](https://pmc.ncbi.nlm.nih.gov/articles/PMC8536974/)（三個 domain 的變異都能造成同一種病，同一個位置的不同取代又能造成不同病）。本專案在 TTN 上量到「基礎機率主導、近鄰無加值」，與 LMNA 這個結論屬於同一類現象——位置本身在多數情況下不足以決定是哪一種疾病。

反面證據也要記住：LMNA 的研究把致病取代位點全部映到結構上，結論是[看不到明顯的基因型—表型對應](https://pmc.ncbi.nlm.nih.gov/articles/PMC8536974/)（三個 domain 的變異都能造成同一種病，同一個位置的不同取代又能造成不同病）。本專案在 TTN 上量到「基礎機率主導、近鄰無加值」，與 LMNA 這個結論屬於同一類現象——位置本身在多數情況下不足以決定是哪一種疾病。

---

## 文獻直接支持的四個下一步

1. **把 bp 距離換成 Pfam domain 對齊（優先序最高，唯一能改變結果的）。** HomolVar 用 Pfam family alignment 定義等價位置，同源位置致病性一致率 95%，且不需要 3D 結構。TTN 的 Ig（PF07679 I-set）與 FN3（PF00041）合計約 280 個重複 domain，把 HGMD 致病變異投影到 Pfam HMM 的 match state，就得到一個與 genomic 距離正交的相似度。這同時提供了「為何 bp 距離無加值」的機制性解釋：折疊上的等價殘基在座標上可相隔數萬 bp。
2. **把 band 換成 PSI／外顯子。** Roberts 2015 的 PSI 是現成的量化指標，cardiodb.org/titin 有現成註釋；ACMG SF v3.2 的 TTN **回報**政策已改用高 PSI 外顯子（分類規範 Morales 2020 則仍用 A-band，見 §3 更正）。骨骼肌側是 exon 344／364 兩個手工熱點，且**只適用於 missense**——`band_phenotype_prior` 應比照拆開 truncating／missense 分別計算。這比任何參數調整都有領域依據。
3. **confidence 改報 likelihood ratio。** PM5 校準論文把「同位置類比」量成 pLR 8.4–16.3、Genome Biology 2025 把「同源位置類比」量成 pLR 13.0，可直接搬用：拿 `nearest_analog_bp`（或 domain 對齊後的等價位置）分層估每一層的 LR，取代目前的 High／Low 標籤。
4. **鄰域半徑要正規化，不要用固定半衰期。** P-KNN 用一份未標註資料集約束最小搜尋半徑，正好對應本專案發現的問題——HGMD 密集區本來就是容易的區域，固定 bp 尺度會把「資料密度」誤當成「證據強度」。
