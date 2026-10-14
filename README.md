# AlphaGenome Atlas — 30 Days of Offline Python Examples

**RESEARCH USE ONLY · NO CLINICAL USE**

**作者：ARiN · 2026 iThome 鐵人賽 · Python 3.12+**

這份 repository 收錄以最後修訂教材為基準的 **原創教學 Python Lab、合成 Fixture、離線測試與證據追蹤範例**。不收錄 30 天文章全文，也不包含 Google AI Studio 的 React/Node 原始專案。所有範例僅供學習、研究方法示範，**不得用於臨床診斷**。

> **嚴格區分資料來源：** `SYNTHETIC / FIXTURE` 不是真實 AlphaGenome 推論；AVI 排序分數不等於致病機率；ACMG 分類必須經合格專家審核。請勿輸入真實病患資料或 API 金鑰。

**30 天文章：**[iThome 系列總入口](https://ithelp.ithome.com.tw/users/20183596/ironman/9657) · **Research Agent 線上展示：**[NeuroGene Atlas Demo](https://neurogene-atlas-demo-jqvfylqbia-de.a.run.app/)（獨立的應用程式，並非本 Python repo 的同一套來源碼；畫面上的模擬數值不可視為正式模型證據）。

## 1. 快速開始

主要範例使用 Python 標準函式庫，不需要設定 AlphaGenome 或 Gemini 金鑰。請使用 Python 3.12 或更新版本。

```bash
python3.12 scripts/check_repo_safety.py
python3.12 scripts/run_tests.py

# Day 04：離線 Research Agent 資料流範例
cd examples/day04
python3.12 day04_offline_mvp_bridge_lab.py run --mode offline \
  --fixture day04_synthetic_mvp.json --output-dir day04_output
```

Day 17–18 **Synthetic-only 離線解讀**（從 repository 根目錄執行，無須 API Key）：

```bash
python3.12 examples/historical_readers/day17_historical_avi_reader.py --verify-only
python3.12 examples/historical_readers/day18_historical_splicing_reader.py --verify-only
python3.12 examples/historical_readers/day17_historical_avi_reader.py --output-dir day17_output
python3.12 examples/historical_readers/day18_historical_splicing_reader.py --output-dir day18_output
```

`examples/historical_readers/README.md` 說明虛構數據來源、SHA-256 fixture manifest 與詮釋界線；CLI 檔名保留 `historical_` 是遷移相容策略，**不能**藉此認為結果來自真實 API。請勿將 `day17_output` 或 `day18_output` 提交為模型研究結果。

以下是另一個獨立示例；先回到 repository 根目錄（`cd ../..`）。

Day 21–28 採用**同一套共用模組**，避免複製 8 份版本不一致的 CLI：

```bash
cd examples/variant_pipeline
python3.12 day27_alpha_variant_cli.py run --mode demo --output-dir day27_output
python3.12 day27_alpha_variant_cli.py verify --run-dir day27_output
python3.12 day28_offline_evidence_report_lab.py run --mode offline \
  --run-dir day27_output --output-dir day28_output
python3.12 day28_offline_evidence_report_lab.py verify \
  --run-dir day27_output --report-dir day28_output
```

## 2. 程式結構

```text
alphagenome-atlas-30day-examples/
├── README.md
├── article_links.csv                # 30 篇文章 URL 與範例位置
├── LICENSE / NOTICE.md              # MIT 原創程式碼與權利界線
├── examples/
│   ├── day02/ ... day16/            # 單日離線 Lab
│   ├── historical_readers/          # Day 17–18，SHA-256 合成 Fixture 的離線解讀
│   ├── day19/ day20/                # 重試／區間註釋
│   ├── variant_pipeline/           # Day 21–28，共用 VCF → Evidence workflow
│   └── day29/ day30/               # 隱私與 Capstone
├── scripts/
│   ├── check_repo_safety.py
│   └── run_tests.py
├── docs/                             # 科研限制與發布紀錄
└── .github/workflows/offline-tests.yml
```

## 3. 30 天文章與教學範例對照

> **文章連結**是原 iThome 發表頁；此表的標題以最後重整的教學版為主。若 iThome 尚未同步修訂，兩者用字和 Lab 內容可能不同。**不可把「範例測試通過」解讀為 AlphaGenome/Gemini 真實 API 成功。**

| Day | 最終整理主題 | iThome 文章 | 範例所在位置 |
|:---:|---|---|---|
| 01 | 破解 98% 基因暗物質 —— 從 DNA 變異找到可研究的線索 | [Day 01](https://ithelp.ithome.com.tw/articles/10410043) | 概念篇（無程式） |
| 02 | 基因調控與組織情境 —— 同一變異為何要看不同細胞 | [Day 02](https://ithelp.ithome.com.tw/articles/10410050) | [`examples/day02/`](examples/day02/) |
| 03 | 保守性與分子預測 —— 兩種互補的變異線索 | [Day 03](https://ithelp.ithome.com.tw/articles/10410841) | [`examples/day03/`](examples/day03/) |
| 04 | In Silico 預先計算 —— 看懂 Variant View 與 Locus View | [Day 04](https://ithelp.ithome.com.tw/articles/10412124) | [`examples/day04/`](examples/day04/) |
| 05 | DeepMind 生命科學版圖 —— 從 AlphaFold 結構到 AlphaGenome 調控 | [Day 05](https://ithelp.ithome.com.tw/articles/10412212) | [`examples/day05/`](examples/day05/) |
| 06 | 百萬鹼基視野 —— 使用上下文捕捉三維空間中的遠距調控 | [Day 06](https://ithelp.ithome.com.tw/articles/10412710) | [`examples/day06/`](examples/day06/) |
| 07 | 演算法核心拆解 —— CNN 與 Transformer 的多尺度融合 | [Day 07](https://ithelp.ithome.com.tw/articles/10413647) | [`examples/day07/`](examples/day07/) |
| 08 | 單鹼基解析度圖譜 —— 一次讀懂剪接與染色質訊號 | [Day 08](https://ithelp.ithome.com.tw/articles/10414152) | [`examples/day08/`](examples/day08/) |
| 09 | 突變效應預測本質 —— REF/ALT 差分與視窗位移 | [Day 09](https://ithelp.ithome.com.tw/articles/10414660) | [`examples/day09/`](examples/day09/) |
| 10 | 1 PB 巨量資料庫 —— 用座標索引理解 Atlas 查詢 | [Day 10](https://ithelp.ithome.com.tw/articles/10415326) | [`examples/day10/`](examples/day10/) |
| 11 | 打破維度災難 —— AVI 的 18 項輸入特徵怎麼讀 | [Day 11](https://ithelp.ithome.com.tw/articles/10415888) | [`examples/day11/`](examples/day11/) |
| 12 | 讀懂 AVI 評分指標 —— 研究優先級不是發病機率 | [Day 12](https://ithelp.ithome.com.tw/articles/10416613) | [`examples/day12/`](examples/day12/) |
| 13 | DNM1 個案導讀 —— 從模型線索走到可追溯證據 | [Day 13](https://ithelp.ithome.com.tw/articles/10416857) | [`examples/day13/`](examples/day13/) |
| 14 | DNM1 剪接假說 —— 如何設計獨立驗證 | [Day 14](https://ithelp.ithome.com.tw/articles/10417350) | [`examples/day14/`](examples/day14/) |
| 15 | 高通量飽和突變驗證 —— 從計算候選到實驗設計 | [Day 15](https://ithelp.ithome.com.tw/articles/10417946) | [`examples/day15/`](examples/day15/) |
| 16 | Python SDK 使用前檢查 —— 版本、金鑰與授權 | [Day 16](https://ithelp.ithome.com.tw/articles/10418468) | [`examples/day16/`](examples/day16/) |
| 17 | 單點變異查詢入門 —— 從 SNV 輸出看懂模型回應 | [Day 17](https://ithelp.ithome.com.tw/articles/10418899) | [`examples/historical_readers/`](examples/historical_readers/) |
| 18 | 多模態軌跡練習 —— 比較不同組織的訊號變化 | [Day 18](https://ithelp.ithome.com.tw/articles/10419245) | [`examples/historical_readers/`](examples/historical_readers/) |
| 19 | 批次查詢與限流概念 —— asyncio 與重試策略實作 | [Day 19](https://ithelp.ithome.com.tw/articles/10419813) | [`examples/day19/`](examples/day19/) |
| 20 | 基因座區間探索 —— 從座標系統理解 Motif 標註 | [Day 20](https://ithelp.ithome.com.tw/articles/10420199) | [`examples/day20/`](examples/day20/) |
| 21 | 從 VCF 到研究審閱封包 —— QC、註釋與證據界線 | [Day 21](https://ithelp.ithome.com.tw/articles/10420536) | [`examples/variant_pipeline/`](examples/variant_pipeline/) |
| 22 | VCF 資料前處理 —— FILTER、QUAL、REF 與 BED QC | [Day 22](https://ithelp.ithome.com.tw/articles/10420780) | [`examples/variant_pipeline/`](examples/variant_pipeline/) |
| 23 | 自動化批次標註 —— 從缺失值到結構化欄位設計 | [Day 23](https://ithelp.ithome.com.tw/articles/10421282) | [`examples/variant_pipeline/`](examples/variant_pipeline/) |
| 24 | REF/ALT 軌跡差分分析 —— 從訊號變化提出機制問題 | [Day 24](https://ithelp.ithome.com.tw/articles/10421741) | [`examples/variant_pipeline/`](examples/variant_pipeline/) |
| 25 | 合成基因座軌跡圖 —— 用 Python 畫出 REF/ALT 對比 | [Day 25](https://ithelp.ithome.com.tw/articles/10422140) | [`examples/variant_pipeline/`](examples/variant_pipeline/) |
| 26 | Motif 文法與變異 —— 從序列偏好理解 REF/ALT 差異 | [Day 26](https://ithelp.ithome.com.tw/articles/10422542) | [`examples/variant_pipeline/`](examples/variant_pipeline/) |
| 27 | 封裝分析流程 —— 打造自己的 AlphaVariant CLI | [Day 27](https://ithelp.ithome.com.tw/articles/10422802) | [`examples/variant_pipeline/`](examples/variant_pipeline/) |
| 28 | AlphaGenome 的限制與證據追蹤 —— 學會說清楚我們不知道什麼 | [Day 28](https://ithelp.ithome.com.tw/articles/10423103) | [`examples/variant_pipeline/`](examples/variant_pipeline/) |
| 29 | AI 基因解密時代的倫理防線、數據隱私與未來十年精準醫療藍圖 | [Day 29](https://ithelp.ithome.com.tw/articles/10423247) | [`examples/day29/`](examples/day29/) |
| 30 | 30 天實戰總結 —— 從 AI Studio MVP 到可追溯的 Python 研究工具 | [Day 30](https://ithelp.ithome.com.tw/articles/10423663) | [`examples/day30/`](examples/day30/) |

## 4. 科學證據界線

- **Day 01：**純概念篇，沒有冒充可執行程式的占位檔。
- **Day 02–16、19–30：**合成資料／離線教學；程式測試僅驗證軟體邏輯，不是臨床、濕實驗或官方服務效能評估。
- **Day 17–18：**使用 `chrToy` 合成教學 Fixture 的完全離線讀取工具，解讀虛構 AVI／分位數／PHRED-like 尺度、RNA splicing 與 REF/ALT 差分。SHA-256 只驗證**教學 Fixture 的檔案完整性**，不代表 Google DeepMind 模型預測或官方來源驗證；兩支 `.py` 檔沿用 `historical_` 舊名稱僅為相容性，內容不含歷史模型數值。
- **Day 11–12：**現在全部使用 `SYNTHETIC_EDUCATIONAL_DEMO`，包含獨立的虛構 AVI-like 對照組；沒有收錄先前 Reboot 的歷史模型分數。
- **Day 21–28：**資料預處理、差分圖表、Motif 類比與 CLI 是合成示範；不得視為「官方 AlphaGenome API 已連線」。
- **個資：**不得上傳實際患者 VCF、身分識別、病歷表型、院內資料或憑證。模擬 Token 與假名不等於去識別保證。

## 5. 發布與授權

本 repository 的**原創教學程式碼**以 MIT 授權；獲授權的原創教學圖文可另以 CC BY 4.0 標示。第三方圖片、Logo、API/SDK 原始碼、模型輸出與託管服務不因本 repo 授權而改變權利歸屬。詳見 `LICENSE` 與 `NOTICE.md`。

**使用提醒：** 這是從經 SHA-256 鎖定的離線教學候選版整理的公開程式碼範例。讀者可自行執行 `scripts/check_repo_safety.py` 與 `scripts/run_tests.py`；測試通過僅代表示範程式的離線行為符合預期，並非官方 AlphaGenome API、臨床或濕實驗驗證。
