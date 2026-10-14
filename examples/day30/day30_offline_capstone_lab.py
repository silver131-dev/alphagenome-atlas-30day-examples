"""Day 30 | 30-day evidence capstone (OFFLINE, synthetic data only).

Read existing Day27/28/29 educational outputs, verify local manifests and
source links, then produce a minimal, non-clinical teaching dashboard.
No import of alphagenome, no API, no patient data. Hashes are NOT signatures.
Python 3.12+ standard library only.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import re
import sys
from pathlib import Path, PurePosixPath

SCHEMA = "day30.synthetic.capstone.v1"
EVIDENCE = "SYNTHETIC_EDUCATIONAL_DEMO"
STAGES = (
    ("01–05", "從生命暗物質出發", "調控、組織情境、研究問題"),
    ("06–10", "看懂模型與資料", "1 Mb、CNN/Transformer、REF/ALT"),
    ("11–15", "評分與實驗假說", "AVI、DNM1、驗證設計"),
    ("16–19", "從 MVP 走到 Python", "AI Studio、SDK、API、重試"),
    ("20–27", "離線管線與視覺化", "VCF、QC、Motif、CLI"),
    ("28–30", "證據與研究治理", "Provenance、隱私、完整性"),
)
MANIFESTS = {
    "day27": ("day27_run_manifest.json", 23),
    "day28": ("day28_report_manifest.json", 5),
    "day29": ("day29_report_manifest.json", 6),
}


class EvidenceError(ValueError):
    """Raised when an educational artifact is missing or inconsistent."""


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def read_json(path: Path) -> dict:
    if not path.is_file() or path.is_symlink() or path.stat().st_size > 10_000_000:
        raise EvidenceError(f"Missing, linked or oversized JSON file: {path.name}")
    try:
        with path.open(encoding="utf-8") as f:
            result = json.load(f, parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))
    except (UnicodeError, ValueError, json.JSONDecodeError) as exc:
        raise EvidenceError(f"Invalid JSON: {path.name}") from exc
    if not isinstance(result, dict):
        raise EvidenceError(f"Expected JSON object: {path.name}")
    return result


def require(value: bool, reason: str) -> None:
    if not value:
        raise EvidenceError(reason)


def check_flag(obj: dict, name: str, expected) -> None:
    require(type(obj.get(name)) is type(expected) and obj.get(name) == expected,
            f"Invalid {name}: expected {expected!r}")


def verify_manifest(folder: Path, manifest_name: str, count: int) -> tuple[dict, str]:
    folder = Path(folder)
    manifest_path = folder / manifest_name
    manifest = read_json(manifest_path)
    artifacts = manifest.get("artifacts")
    require(type(artifacts) is dict and len(artifacts) == count, "Manifest artifact set invalid")
    check_flag(manifest, "artifact_count", count)
    check_flag(manifest, "mode", "OFFLINE_ONLY")
    check_flag(manifest, "real_api_calls_this_run", 0)
    root = folder.resolve(strict=True)
    for raw_name, info in artifacts.items():
        require(isinstance(raw_name, str), "Manifest path must be text")
        p = PurePosixPath(raw_name)
        require(not p.is_absolute() and len(p.parts)>0 and
                all(x not in (".", "..", "") for x in p.parts), "Unsafe relative artifact path")
        require(type(info) is dict and type(info.get("bytes")) is int and
                type(info.get("sha256")) is str and
                re.fullmatch(r"[0-9a-f]{64}", info["sha256"]) is not None,
                "Malformed artifact checksum record")
        candidate = folder.joinpath(*p.parts)
        require(candidate.is_file() and not candidate.is_symlink(), "Artifact file missing or linked")
        require(candidate.resolve().is_relative_to(root), "Artifact escapes source root")
        require(candidate.stat().st_size == info["bytes"], "Artifact byte count mismatch")
        require(digest(candidate) == info["sha256"], "Artifact SHA-256 mismatch")
    return manifest, digest(manifest_path)


def verify_source_chain(day27: Path, day28: Path, day29: Path) -> dict:
    day27, day28, day29 = map(Path, (day27, day28, day29))
    manifest27, sha27 = verify_manifest(day27, *MANIFESTS["day27"])
    manifest28, sha28 = verify_manifest(day28, *MANIFESTS["day28"])
    manifest29, sha29 = verify_manifest(day29, *MANIFESTS["day29"])
    require(manifest28.get("source_day27_manifest_sha256") == sha27,
            "Day28 does not refer to the verified Day27 manifest")
    require(manifest29.get("source_manifest_sha256") == sha28,
            "Day29 manifest does not refer to Day28 manifest")

    s27 = read_json(day27 / "day27_pipeline_summary.json")
    s28 = read_json(day28 / "day28_evidence_report.json")
    s29 = read_json(day29 / "day29_privacy_audit.json")
    preview29 = read_json(day29 / "day29_local_preview.json")
    for s in (s27, s28, s29):
        check_flag(s, "evidence_class", EVIDENCE)
        check_flag(s, "mode", "OFFLINE_ONLY")
        check_flag(s, "real_api_calls_this_run", 0)
    check_flag(s27, "schema_version", "day27.synthetic.pipeline.v1")
    check_flag(s28, "schema_version", "day28.synthetic.evidence_report.v1")
    check_flag(s29, "schema_version", "day29.synthetic.privacy_gate.v1")
    check_flag(s27, "status", "PASS")
    check_flag(s28, "status", "PASS")
    check_flag(s29, "status", "PASS_FOR_LOCAL_SYNTHETIC_EDUCATION_ONLY")
    check_flag(s27, "assembly", "TOY_NOT_GRCH38")
    check_flag(s28, "assembly", "TOY_NOT_GRCH38")
    for name in ("model_prediction", "official_atlas_data", "atlas_query_performed", "clinical_interpretation"):
        if name in s27: check_flag(s27, name, False)
    for name in ("model_prediction", "atlas_query_performed", "clinical_interpretation"):
        check_flag(s28, name, False)
    check_flag(s27, "pipeline_stage_count", 6)
    stages = s27.get("stages")
    require(isinstance(stages,list) and [x.get("day") for x in stages]==list(range(21,27)),
            "Unexpected pipeline stage sequence")
    for i, s in enumerate(stages):
        require(s.get("status")=="PASS", f"Day {21+i} did not pass")
    check_flag(s27, "score_combination_performed", False)
    check_flag(s27, "source_manifest_sha256", s28["source_provenance"]["source_bundle_sha256"])
    src28 = s28.get("source_provenance", {})
    require(isinstance(src28, dict) and src28.get("day27_run_manifest_sha256")==sha27,
            "Day28 report source Day27 manifest mismatch")
    require(src28.get("day27_summary_sha256")==digest(day27 / "day27_pipeline_summary.json"),
            "Day28 report does not refer to this Day27 summary")
    check_flag(src28, "artifact_count_verified", 23)
    check_flag(src28, "hash_check", "MATCH")
    check_flag(src28, "signature_verified", False)
    check_flag(src28, "native_alpha_genome_response_retained", False)
    check_flag(s28, "manual_review", "REQUIRED")
    check_flag(s28, "confidence_calibrated", False)
    require(s28.get("clinical_classification") is None and s28.get("disease_probability") is None,
            "Clinical result must remain null")
    counts=s28.get("pipeline_counts")
    require(type(counts) is dict, "Missing Day28 counts")
    numeric = dict(candidate_count=8, qc_selected=3, qc_review_required=2,
                   qc_rejected=3, annotated_count=2, missing_score_count=1, stage_count=6)
    for key, val in numeric.items():check_flag(counts,key,val)
    require([s27.get("candidate_count"),s27.get("selected_count"),s27.get("missing_score_count")]
            ==[8,3,1],"Day27 counts mismatch")
    ann=s28.get("annotations")
    require(type(ann) is list and len(ann)==3,"Wrong selected candidate count")
    require(sum(x.get("score_status")=="MISSING_SCORE" and x.get("synthetic_score") is None for x in ann)==1,
            "Missing score must stay JSON null")
    require(all(x.get("review_status")=="REVIEW_REQUIRED" for x in ann),
            "Candidates must remain pending review")

    source29=s29.get("source", {})
    require(type(source29) is dict,"Missing Day29 source provenance")
    require(source29.get("day28_manifest_sha256")==sha28 and
            source29.get("day28_report_sha256")==digest(day28/"day28_evidence_report.json"),
            "Day29 evidence does not refer to verified Day28 data")
    require(source29.get("day28_manifest_report_entry_matches") is True and
            source29.get("day28_manifest_signed") is False,
            "Day29 report integrity/signature status invalid")
    require(manifest29.get("source_report_sha256")==digest(day28/"day28_evidence_report.json"),
            "Day29 source report digest mismatch")
    for key, val in dict(release_gate="BLOCKED", governance_review_status="MANUAL_REVIEW_REQUIRED",
                         anonymization_claim="NOT_ESTABLISHED",pseudonymization_claim="NOT_ESTABLISHED",
                         token_is_anonymization=False, reidentification_risk_eliminated=False).items():
        check_flag(s29,key,val)
    policy=s29.get("demo_policy",{})
    require(type(policy) is dict,"Missing governance policy")
    for key in ("publish_allowed","external_transfer_allowed","human_approval_completed",
                "consent_verified","irb_or_ethics_review_verified","real_patient_processing_allowed"):
        check_flag(policy,key,False)
    ag=s29.get("aggregate_counts",{})
    require(type(ag) is dict,"Missing Day29 aggregate counts")
    for key,val in dict(candidate_count=8,selected_count=3,missing_score_count=1,
                        day27_artifact_count_claimed_by_day28=23).items():
        check_flag(ag,key,val)
    for key in ("random_token_included","per_variant_data_included","publish_allowed","external_transfer_allowed"):
        check_flag(preview29,key,False)
    for key in ("local_demo_token","variant_id","variant_query","patient_id","full_name"):
        require(key not in preview29,"Private field appears in local preview")
    check_flag(preview29,"real_api_calls_this_run",0)
    require(preview29.get("governance_review_status")=="MANUAL_REVIEW_REQUIRED",
            "Day29 local preview must require review")
    require(digest(day28/"day28_evidence_report.json")==manifest28["artifacts"]["day28_evidence_report.json"]["sha256"],
            "Day28 manifest mismatch")
    return {
      "schema_version":SCHEMA,"status":"PASS_FOR_SYNTHETIC_OFFLINE_TEACHING_ONLY",
      "evidence_class":EVIDENCE,"mode":"OFFLINE_ONLY",
      "source_manifest_digests":{"day27":sha27,"day28":sha28,"day29":sha29},
      "source_artifact_counts":{"day27":23,"day28":5,"day29":6,"total":34},
      "day27_pipeline_stages":6,
      "candidate_count":8,"qc_selected":3,"qc_review_required":2,"qc_rejected":3,
      "annotated_count":2,"missing_score_count":1,
      "source_chains_verified":True,"source_manifest_digital_signatures_verified":False,
      "models_called_in_this_capstone":False,"historical_live_evidence_imported":False,
      "real_api_calls_this_run":0,"genome_assembly":"TOY_NOT_GRCH38",
      "privacy_release_gate":"BLOCKED","external_transfer_allowed":False,
      "human_approval_completed":False,"clinical_interpretation":False,
      "anonymization_claim":"NOT_ESTABLISHED",
      "reidentification_risk_eliminated":False,
      "manual_review":"REQUIRED",
      "github_live_visibility_independently_checked":False,
      "science_limitations":[
        "本驗收僅核對 Day27–29 人工合成資料與本地雜湊鏈，不驗證真實生物效應。",
        "SHA-256 不等於數位簽章，也無法單獨證明資料來源可靠。",
        "歷史 Reboot SDK 0.9.0 曾有 LIVE 成功紀錄，這次未匯入或重新執行該 API。",
        "新 Bridge 0.10.0 的 LIVE 相容性與當前 API 權限仍待驗證。",
        "隨機 Token、假名化及合成標記不能證明真實基因資料已匿名化。",
        "Day29 禁止發佈與外傳，未完成人工審閱、IRB 或臨床分類。",
      ],
    }


def read_approved_titles(path: Path) -> list[dict]:
    if not path.is_file() or path.is_symlink():raise EvidenceError("Missing approved title CSV")
    with path.open(encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        rows=list(reader)
    require(len(rows)==28,"Expected precisely 28 author-approved titles")
    for n,row in enumerate(rows,1):
        require(row.get("day")==f"{n:02d}" and row.get("title_status")=="AUTHOR_APPROVED",
                "Approval matrix out of sequence")
        title=row.get("approved_title")
        require(isinstance(title,str) and title.strip() and "\n" not in title and "|" not in title,
                "Invalid title text")
    return rows


def make_map(report:dict,titles:list[dict])->str:
    lines=["# AlphaGenome Atlas｜30 天課程導覽（作者審稿候選）", "",
        "Day01–28 為作者已確認的**標題**；其餘各日正文仍按作者審稿狀態處理。",
        "Day29 已發表，但此處的 Day29 隱私程式為附加式離線教學，不是正式文章已更新。",
        "Day30 為本次新增的候選稿，尚未發表。GitHub 維持私人狀態。", "",
        "## 六階段學習地圖", "",
        "| 範圍 | 學習主題 | 關鍵技能 |", "|---|---|---|"]
    for interval,subject,skill in STAGES:
        lines.append(f"| Day {interval} | {subject} | {skill} |")
    lines.extend(["","## Day 01–30 每日標題（審稿索引）", "",
                  "| Day | 標題 | 題名來源／內容狀態 |", "|---|---|---|"])
    for row in titles:
        day=int(row["day"])
        content="Day01 正文作者已核准" if day==1 else "正文待作者確認"
        lines.append(f"| {day:02d} | {row['approved_title']} | 作者核准標題；{content} |")
    lines.append("| 29 | 前瞻展望 —— AI 基因解密時代的倫理防線、數據隱私與未來十年精準醫療藍圖 | 已發表；附加實戰包未上線 |")
    lines.append("| 30 | 30 天實戰總結 —— 從 AI Studio MVP 到可追溯的 Python 研究工具 | 新稿；待作者審查 |")
    lines.extend(["","## 程式與真實研究證據界線", "",
      "- OFFLINE：Day21–29 人工合成資料之 Python 示範可離線重現。",
      "- HISTORICAL LIVE：Reboot SDK 0.9.0 曾有歷史 Atlas AVI 與 Model 預測衍生摘要；**不在本次 Day30 可執行資料鏈內**。",
      "- CURRENT LIVE：沒有授權也沒有執行。SDK 0.10.0 實際服務相容性仍待驗證。",
      "- Release：本次檢查結果 **BLOCKED**。PRIVATE GitHub 不能宣稱為公開資源。",
      "- 以上 SHA-256 僅驗證相對本地 manifest 的位元一致性，不構成來源認證。",
      "", f"**Day27／28／29 完整性核對：{report['source_artifact_counts']['total']} 個本地來源檔。**", ""])
    return "\n".join(lines)


def make_checklist(report:dict)->str:
    return "\n".join([
      "# Day30｜研究成果交付前 Checklist", "",
      "- [x] 已檢查 Day27 23 個人工合成輸出及 SHA-256。",
      "- [x] 已檢查 Day28 5 個報告產物及來源關聯。",
      "- [x] 已檢查 Day29 6 個隱私審查產物及來源關聯。",
      "- [x] 保留 JSON null 與 MISSING_SCORE；沒有把缺值補成零。",
      "- [x] 只做 OFFLINE 合成教學，real_api_calls_this_run=0。",
      "- [ ] 正式研究資料來源／參考基因組 REF 獨立核實（未進行）。",
      "- [ ] SDK 0.10.0 的安裝版簽名與經授權 LIVE 查詢（未進行）。",
      "- [ ] 真實病人／倫理／隱私與資料使用許可（未確認）。",
      "- [ ] 人工批准新教學文章及任何公開發布。",
      "- [ ] GitHub 公開／授權／發行版（PUBLICATION_ON_HOLD）。",
      "", "**結論：PASS_FOR_SYNTHETIC_OFFLINE_TEACHING_ONLY，不代表臨床、法遵或公開發佈通過。**", ""])


def make_svg(report:dict)->str:
    # Every label comes from fixed, trusted text. Counts from validated synthetic fixtures only.
    colors={"bg":"#E5DEDA","ink":"#293044","muted":"#65677A","edge":"#CFC3BE", "white":"#FFFDFC","rose":"#CB8F94","gold":"#B99B67","blue":"#85A8BE","violet":"#AAA2C4"}
    pieces=[f'<svg xmlns="http://www.w3.org/2000/svg" width="1920" height="1080" viewBox="0 0 1920 1080">',
            f'<rect width="1920" height="1080" fill="{colors["bg"]}"/>',
            '<text x="90" y="100" fill="#293044" font-family="Noto Sans CJK TC, Arial, sans-serif" font-weight="700" font-size="53">AlphaGenome Atlas｜30 天實戰總結</text>',
            '<text x="94" y="150" fill="#65677A" font-family="Noto Sans CJK TC, Arial, sans-serif" font-size="26">Google AI Studio MVP → Python → 研究證據 → 倫理治理</text>']
    for i,(span,subject,detail) in enumerate(STAGES):
        col,row=i%3,i//3
        x=92+col*581; y=218+row*267
        color=[colors['blue'],colors['violet'],colors['gold'],colors['rose'],colors['blue'],colors['violet']][i]
        pieces.extend([
          f'<rect x="{x}" y="{y}" width="543" height="220" rx="28" fill="{colors["white"]}" stroke="{colors["edge"]}" stroke-width="2"/>',
          f'<rect x="{x+24}" y="{y+24}" width="115" height="43" rx="20" fill="{color}" opacity="0.22"/>',
          f'<text x="{x+37}" y="{y+53}" fill="#293044" font-size="23" font-family="Arial, sans-serif" font-weight="700">Day {span}</text>',
          f'<text x="{x+27}" y="{y+116}" fill="#293044" font-size="33" font-family="Noto Sans CJK TC, Arial, sans-serif" font-weight="700">{html.escape(subject)}</text>',
          f'<text x="{x+28}" y="{y+163}" fill="#65677A" font-size="22" font-family="Noto Sans CJK TC, Arial, sans-serif">{html.escape(detail)}</text>',
        ])
    pieces.extend([
      '<rect x="91" y="777" width="1740" height="114" rx="26" fill="#FFFDFC" fill-opacity="0.93"/>',
      '<text x="121" y="821" fill="#293044" font-size="29" font-weight="700" font-family="Noto Sans CJK TC, Arial, sans-serif">可驗證的離線證據鏈</text>',
      '<text x="121" y="862" fill="#65677A" font-size="23" font-family="Noto Sans CJK TC, Arial, sans-serif">Day27：23 檔  →  Day28：5 檔  →  Day29：6 檔  ｜  共 34 個來源檔 SHA 相符</text>',
      '<text x="93" y="958" fill="#293044" font-size="25" font-family="Noto Sans CJK TC, Arial, sans-serif" font-weight="700">SYNTHETIC • OFFLINE ONLY • MANUAL REVIEW REQUIRED</text>',
      '<text x="93" y="998" fill="#65677A" font-size="21" font-family="Noto Sans CJK TC, Arial, sans-serif">未執行真實 API｜不可視為臨床效能或資料匿名化驗證｜GitHub PRIVATE • PUBLICATION ON HOLD</text>',
      '</svg>'])
    return "\n".join(pieces)


def make_html(report:dict, svg:str)->str:
    # Minimal dashboard; no external links, JS, remote fonts, or patient-level fields.
    summary=f"{report['source_artifact_counts']['total']} 個來源檔完整性檢查通過"
    return "\n".join([
      "<!DOCTYPE html>","<html lang=\"zh-Hant\"><head><meta charset=\"UTF-8\"><meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
      "<title>Day30 | Offline Evidence Capstone</title>",
      "<style>body{margin:0;background:#E5DEDA;color:#293044;font-family:system-ui,'Microsoft JhengHei',sans-serif}main{max-width:1200px;margin:auto;padding:24px}svg{max-width:100%;height:auto;border-radius:16px}section{background:#fffdfa;padding:22px;border-radius:16px;margin:18px 0}li{margin:8px 0}</style></head><body><main>",
      "<h1>Day30｜30 天實戰總結</h1>",svg,
      f"<section><h2>{html.escape(summary)}</h2><p>這是本地人工合成資料的雜湊一致性核對，非伺服器簽章，也不是生醫模型效度證明。</p></section>",
      "<section><h2>真正完成什麼？</h2><ul><li>離線六階段分析與研究審閱封包</li><li>保留缺值、逐階段來源與檔案完整性</li><li>隱私審查狀態仍 BLOCKED，不允許外傳</li></ul></section>",
      "<section><h2>仍待驗證什麼？</h2><p>LIVE API、真實模型與組織效能、病人資料治理、臨床／法遵及公開授權。</p></section>",
      "</main></body></html>"
    ])


def write_new(folder: Path, filename: str, content: str) -> dict:
    path=folder / filename
    path.write_text(content,encoding="utf-8")
    return {"bytes":path.stat().st_size,"sha256":digest(path)}


def create_report(output_dir:Path, day27:Path, day28:Path, day29:Path,titles_file:Path)->dict:
    report=verify_source_chain(day27,day28,day29)
    titles=read_approved_titles(titles_file)
    svg=make_svg(report)
    contents={
       "day30_capstone_report.json":json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n",
       "day30_learning_map.md":make_map(report,titles),
       "day30_release_checklist.md":make_checklist(report),
       "day30_capstone_map.svg":svg,
       "day30_capstone_dashboard.html":make_html(report,svg),
    }
    output_dir=Path(output_dir)
    if output_dir.exists():raise EvidenceError("Output directory exists; choose a NEW output location")
    output_dir.mkdir(parents=True)
    files={k:write_new(output_dir,k,v) for k,v in contents.items()}
    manifest={"schema_version":"day30.synthetic.capstone_manifest.v1", "mode":"OFFLINE_ONLY",
              "artifact_count":len(files),"artifacts":files,"signed":False,
              "source_manifest_digests":report["source_manifest_digests"],
              "titles_source_sha256":digest(Path(titles_file)), "real_api_calls_this_run":0}
    (output_dir/"day30_report_manifest.json").write_text(
       json.dumps(manifest,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return {"status":"PASS","mode":"OFFLINE_ONLY","source_files_verified":34,
            "output_artifact_count":len(files),"release_gate":"BLOCKED",
            "real_api_calls_this_run":0,
            "report_file":str(output_dir/"day30_capstone_report.json"),
            "html_file":str(output_dir/"day30_capstone_dashboard.html")}


def verify_report(report_dir:Path,day27:Path,day28:Path,day29:Path,titles_file:Path)->dict:
    source=verify_source_chain(day27,day28,day29)
    m,_=verify_manifest(report_dir,"day30_report_manifest.json",5)
    require(m.get("source_manifest_digests")==source["source_manifest_digests"],
            "Output source-chain manifest digest mismatch")
    require(m.get("titles_source_sha256")==digest(Path(titles_file)),
            "Approved-title source has changed")
    stored=read_json(Path(report_dir)/"day30_capstone_report.json")
    require(stored==source,"Report content does not match current sources")
    return {"status":"PASS","source_files_verified":34,
            "output_artifact_count":5,"release_gate":"BLOCKED",
            "real_api_calls_this_run":0,"checksums":"MATCH"}


def main(argv:list[str]|None=None)->int:
    ap=argparse.ArgumentParser(description="30 Day AlphaGenome offline capstone")
    sub=ap.add_subparsers(dest="command",required=True)
    for command in ("run","verify"):
        p=sub.add_parser(command)
        p.add_argument("--day27-dir",default="example_output/day27_output")
        p.add_argument("--day28-dir",default="example_output/day28_output")
        p.add_argument("--day29-dir",default="example_output/day29_output")
        p.add_argument("--titles-file",default="author_approved_day01_day28_titles.csv")
        if command=="run":
            p.add_argument("--mode",required=True,choices=["offline"])
            p.add_argument("--output-dir",default="day30_output")
        else:
            p.add_argument("--report-dir",required=True)
    args=ap.parse_args(argv)
    try:
        if args.command=="run":
            result=create_report(Path(args.output_dir),Path(args.day27_dir),Path(args.day28_dir),
                                 Path(args.day29_dir),Path(args.titles_file))
        else:
            result=verify_report(Path(args.report_dir),Path(args.day27_dir),Path(args.day28_dir),
                                 Path(args.day29_dir),Path(args.titles_file))
    except (EvidenceError,OSError,KeyError,TypeError) as exc:
        print(f"ERROR: {exc}",file=sys.stderr)
        return 2
    print(json.dumps(result,ensure_ascii=False,sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
