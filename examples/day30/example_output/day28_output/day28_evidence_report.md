# Day 28｜合成資料的限制與證據追蹤報告

> **研究教學示例（OFFLINE / TOY_NOT_GRCH38）**。不是 AlphaGenome API 回應、病人資料或臨床結論。

## 一、來源與完整性

- Day27 source manifest SHA-256：`7604c9f1bbdf48185d48c6ea79bbcd5a5cdaf1dbecde780ea3a9f2eb662e005f`
- Day27 summary SHA-256：`3c9e77687e3bc00643d63232dec7bc227e508cac6845cd70b8f3d56b81915bfa`
- Day27 run manifest SHA-256：`a2d54c70e189dc7dd50691d713c8b6490f9ae72602e094fc061a50bda72e6ff7`
- 本地輸出完整性：**23 個 artifact、MATCH**；未使用數位簽章。

**SHA-256 用來檢查內容有沒有變，不代表來源真實，也不能證明科學結論。**

## 二、QC 與缺值

人工候選 8 筆；通過 QC 3、QC 待審 2、QC 排除 3。

| 人工變異 | 人工分數 | 分數狀態 | 最終審閱 |
|---|---:|---|---|
| `chr1:100:A>C` | 0.42 | `SYNTHETIC_SCORE` | `REVIEW_REQUIRED` |
| `chr1:110:A>G` | -0.17 | `SYNTHETIC_SCORE` | `REVIEW_REQUIRED` |
| `chr1:110:A>C` | `null` | `MISSING_SCORE` | `REVIEW_REQUIRED` |

**`null` 表示此處沒有分數，不是 0；即使人工分數有值也需要人工審閱。**

## 三、模型限制與不確定性

- **SYNTHETIC_ONLY：**所有位置、REF/ALT 軌跡、分數與 Motif 均是人工教學資料；不是 AlphaGenome API 或實驗結果。
- **ASSEMBLY_UNVERIFIED：**TOY_NOT_GRCH38；未獨立以 GRCh38 FASTA 核對 REF，也未進行臨床變異正規化。
- **MODEL_TISSUE_NOT_RUN：**未執行真實模型或組織情境預測；不能把差分數值當成剪接或基因調控機制。
- **POPULATION_TRANS_NOT_VALIDATED：**沒有族群泛化、跨區域 trans 調控、臨床疾病風險或效能校準的驗證。
- **PROVENANCE_NOT_SIGNED：**SHA-256 只能檢查相對於本地 manifest 的內容一致性，不是數位簽章或來源真實性的證明。
- **MANUAL_REVIEW_REQUIRED：**SELECTED 只代表通過人工教學 QC；不能直接生成 ACMG/AMP 分類、診斷或治療建議。

## 四、可以與不可以解讀什麼？

- Day24 的 −0.30 只是人工軌跡差分。
- Day26 的 −0.12／+0.12 是另一套人工 PPM-like 序列偏好。兩種單位**不能相加**。
- 沒有執行真實模型、族群泛化測試或任何臨床分類。
- `manual_review=REQUIRED` 代表**仍待人工檢查**，不是已獲核准。

## 五、總結

一份負責任的報告，不只是呈現知道的答案，還要說清楚資料從哪裡來、哪裡是未知，以及哪些推論尚無證據。
