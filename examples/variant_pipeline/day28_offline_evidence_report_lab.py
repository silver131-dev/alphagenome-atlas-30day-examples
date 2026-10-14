"""Day 28 | Explainable, auditable OFFLINE evidence-report teaching CLI.

Reads a verified Day27 synthetic pipeline. Produces a cautious review packet.
No LIVE mode, network calls, API keys, AlphaGenome SDK, clinical inference.
Python 3.12+ standard library only. NOT the separate Day28 Reporter MVP or GitHub v0.1.0.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from html import escape
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import sys
import tempfile
import xml.etree.ElementTree as ET

import day27_alpha_variant_cli as d27

SCHEMA = "day28.synthetic.evidence_report.v1"
MODE = "OFFLINE_ONLY"
EVIDENCE = "SYNTHETIC_EDUCATIONAL_DEMO"
ASSEMBLY = "TOY_NOT_GRCH38"
ARTIFACTS = ("day28_evidence_report.json", "day28_evidence_report.md",
             "day28_evidence_dashboard.html", "day28_evidence_map.svg",
             "day28_evidence_table.csv")
SHA_RE = re.compile(r"[0-9a-f]{64}\Z")
LIMITATIONS = (
    {"code":"SYNTHETIC_ONLY", "text":"所有位置、REF/ALT 軌跡、分數與 Motif 均是人工教學資料；不是 AlphaGenome API 或實驗結果。"},
    {"code":"ASSEMBLY_UNVERIFIED", "text":"TOY_NOT_GRCH38；未獨立以 GRCh38 FASTA 核對 REF，也未進行臨床變異正規化。"},
    {"code":"MODEL_TISSUE_NOT_RUN", "text":"未執行真實模型或組織情境預測；不能把差分數值當成剪接或基因調控機制。"},
    {"code":"POPULATION_TRANS_NOT_VALIDATED", "text":"沒有族群泛化、跨區域 trans 調控、臨床疾病風險或效能校準的驗證。"},
    {"code":"PROVENANCE_NOT_SIGNED", "text":"SHA-256 只能檢查相對於本地 manifest 的內容一致性，不是數位簽章或來源真實性的證明。"},
    {"code":"MANUAL_REVIEW_REQUIRED", "text":"SELECTED 只代表通過人工教學 QC；不能直接生成 ACMG/AMP 分類、診斷或治療建議。"},
)
SVG_NS="http://www.w3.org/2000/svg"
ET.register_namespace("", SVG_NS)


def sha256(path:Path)->str:
    return d27.file_hash(Path(path))


def annotations_with_missing(rows:list)->tuple[list[dict],list[str]]:
    """Small schema guard, distinguishing a real numeric 0.0 from missing None."""
    if not isinstance(rows,list) or len(rows)!=3:
        raise ValueError("Expected 3 selected synthetic annotations")
    seen=set(); output=[]; missing=[]
    for row in rows:
        if not isinstance(row,dict): raise ValueError("Annotation is not a JSON object")
        variant=row.get("variant_id")
        if not isinstance(variant,str) or not re.fullmatch(r"chr1:\d+:[ACGT]>[ACGT]",variant) or variant in seen:
            raise ValueError("Duplicate or non-synthetic annotation identifier")
        seen.add(variant)
        value=row.get("synthetic_score")
        if value is None:
            if row.get("score_status")!="MISSING_SCORE":
                raise ValueError("Null score must be MISSING_SCORE")
            missing.append(variant)
        elif (not isinstance(value,(float,int)) or isinstance(value,bool)
              or not math.isfinite(value) or row.get("score_status")!="SYNTHETIC_SCORE"):
            raise ValueError("Non-missing score must be finite numeric synthetic")
        if row.get("review_status")!="REVIEW_REQUIRED" or row.get("score_source")!="SYNTHETIC_FIXTURE_ONLY":
            raise ValueError("Annotation source/review contract differs")
        output.append({"variant_id":variant,"synthetic_score":value,
                       "score_status":row["score_status"],"review_status":"REVIEW_REQUIRED"})
    if len(missing)!=1:
        raise ValueError("Synthetic teaching packet must preserve one missing score")
    return output,missing


def report_from_day27(run_dir:Path)->dict:
    """Read and verify sealed Day27 offline output before deriving any report."""
    run_dir=Path(run_dir)
    if run_dir.is_symlink() or not run_dir.is_dir():
        raise ValueError("Day27 run directory missing or symlinked")
    verify=d27.verify_run(run_dir)
    if verify.get("checksums")!="MATCH" or verify.get("artifact_count")!=23:
        raise ValueError("Day27 output is not the expected 23-artifact run")
    summary_path=run_dir/"day27_pipeline_summary.json"
    run_manifest_path=run_dir/"day27_run_manifest.json"
    summary=d27.load_json(summary_path)
    d23=d27.load_json(run_dir/"day23"/"day23_annotation_report.json")
    for key,value in {"evidence_class":EVIDENCE,"mode":MODE,"assembly":ASSEMBLY,
                      "real_api_calls_this_run":0,"model_prediction":False,
                      "official_atlas_data":False,"atlas_query_performed":False,
                      "clinical_interpretation":False,"score_source":"SYNTHETIC_FIXTURE_ONLY",
                      "manual_review":"REQUIRED","biological_mechanism_status":"NOT_VALIDATED"}.items():
        if summary.get(key)!=value: raise ValueError("Unexpected Day27 summary evidence: "+key)
    expected={"candidate_count":8,"selected_count":3,"qc_review_count":2,
              "qc_rejected_count":3,"annotated_count":2,"missing_score_count":1,
              "largest_absolute_delta_position_1based":110,"largest_absolute_delta":-0.3,
              "ppm_anchor_deltas":{"MOTIF_A":-0.12,"MOTIF_B":0.12},
              "score_combination_performed":False}
    for key,value in expected.items():
        if summary.get(key)!=value: raise ValueError("Unexpected Day27 tutorial count: "+key)
    if any(d23.get(k)!=v for k,v in {"evidence_class":EVIDENCE,"mode":MODE,"assembly":ASSEMBLY,
                                    "real_api_calls_this_run":0,"clinical_interpretation":False}.items()):
        raise ValueError("Day23 packet provenance or mode mismatch")
    annotations,missing=annotations_with_missing(d23.get("annotations"))
    if (d23.get("missing_score_count")!=len(missing)
        or d23.get("candidate_count")!=8 or d23.get("selected_count")!=3):
        raise ValueError("Day23 annotations and summary counts differ")
    # These scores are intentionally NOT added across modalities.
    return {
        "schema_version":SCHEMA,"status":"PASS","evidence_class":EVIDENCE,
        "mode":MODE,"assembly":ASSEMBLY,"variant_id":summary["variant_id"],
        "source_provenance":{
            "source_kind":"DAY27_OFFLINE_SYNTHETIC_PIPELINE",
            "source_schema_version":summary["schema_version"],
            "source_bundle_sha256":summary["source_manifest_sha256"],
            "day27_summary_sha256":sha256(summary_path),
            "day27_run_manifest_sha256":sha256(run_manifest_path),
            "artifact_count_verified":verify["artifact_count"],
            "hash_check":"MATCH","signature_verified":False,
            "native_alpha_genome_response_retained":False,
        },
        "pipeline_counts":{
            "candidate_count":8,"qc_selected":3,"qc_review_required":2,
            "qc_rejected":3,"annotated_count":2,"missing_score_count":len(missing),
            "stage_count":6,
        },
        "annotations":annotations,
        "missing_values":{"missing_score_variant_ids":missing,
                          "missing_is_zero":False,"zero_is_a_valid_numeric_value":True},
        "teaching_observations":{
            "ref_alt_delta":{"position_1based":110,"value":-0.3,"unit":"ARTIFICIAL_PROXY"},
            "motif_ppm_like":{"MOTIF_A":-0.12,"MOTIF_B":0.12,
                              "unit":"ARTIFICIAL_PPM_MEAN_PROBABILITY"},
            "combined_biological_score":None,
            "mechanism_status":"NOT_VALIDATED",
        },
        "model_limitations":list(LIMITATIONS),
        "manual_review":"REQUIRED","confidence_calibrated":False,
        "clinical_classification":None,"disease_probability":None,
        "model_prediction":False,"atlas_query_performed":False,
        "clinical_interpretation":False,"real_api_calls_this_run":0,
    }


def markdown_report(report:dict)->str:
    s=report["source_provenance"]; c=report["pipeline_counts"]
    rows="\n".join(
        f"| `{r['variant_id']}` | {('`null`' if r['synthetic_score'] is None else r['synthetic_score'])} | "
        f"`{r['score_status']}` | `{r['review_status']}` |" for r in report["annotations"])
    limits="\n".join(f"- **{x['code']}：**{x['text']}" for x in report["model_limitations"])
    return ("# Day 28｜合成資料的限制與證據追蹤報告\n\n"
            "> **研究教學示例（OFFLINE / TOY_NOT_GRCH38）**。不是 AlphaGenome API 回應、病人資料或臨床結論。\n\n"
            "## 一、來源與完整性\n\n"
            f"- Day27 source manifest SHA-256：`{s['source_bundle_sha256']}`\n"
            f"- Day27 summary SHA-256：`{s['day27_summary_sha256']}`\n"
            f"- Day27 run manifest SHA-256：`{s['day27_run_manifest_sha256']}`\n"
            f"- 本地輸出完整性：**{s['artifact_count_verified']} 個 artifact、{s['hash_check']}**；未使用數位簽章。\n\n"
            "**SHA-256 用來檢查內容有沒有變，不代表來源真實，也不能證明科學結論。**\n\n"
            "## 二、QC 與缺值\n\n"
            f"人工候選 {c['candidate_count']} 筆；通過 QC {c['qc_selected']}、QC 待審 {c['qc_review_required']}、QC 排除 {c['qc_rejected']}。\n\n"
            "| 人工變異 | 人工分數 | 分數狀態 | 最終審閱 |\n|---|---:|---|---|\n"+rows+"\n\n"
            "**`null` 表示此處沒有分數，不是 0；即使人工分數有值也需要人工審閱。**\n\n"
            "## 三、模型限制與不確定性\n\n"+limits+"\n\n"
            "## 四、可以與不可以解讀什麼？\n\n"
            "- Day24 的 −0.30 只是人工軌跡差分。\n"
            "- Day26 的 −0.12／+0.12 是另一套人工 PPM-like 序列偏好。兩種單位**不能相加**。\n"
            "- 沒有執行真實模型、族群泛化測試或任何臨床分類。\n"
            "- `manual_review=REQUIRED` 代表**仍待人工檢查**，不是已獲核准。\n\n"
            "## 五、總結\n\n"
            "一份負責任的報告，不只是呈現知道的答案，還要說清楚資料從哪裡來、哪裡是未知，以及哪些推論尚無證據。\n")


def draw_evidence_svg(report:dict)->str:
    r=ET.Element(f"{{{SVG_NS}}}svg",{"width":"1920","height":"1080","viewBox":"0 0 1920 1080",
                                       "role":"img","aria-label":"Four pillars of synthetic evidence reporting"})
    def elem(tag,attributes,text=None):
        e=ET.SubElement(r,f"{{{SVG_NS}}}{tag}",attributes)
        if text is not None:e.text=text
    elem("rect",{"width":"1920","height":"1080","fill":"#E5DEDA"})
    elem("text",{"x":"110","y":"130","font-size":"59","font-weight":"700","fill":"#344250"},"Day 28  |  Evidence & Limitations Reporter")
    elem("text",{"x":"112","y":"191","font-size":"31","fill":"#826670"},"From Day 27 offline pipeline  |  Research review packet")
    cards=[("01","PROVENANCE","23 output artifacts","SHA-256: MATCH"),
           ("02","MISSING VALUES","2 toy scores, 1 null","Unknown is NOT zero"),
           ("03","LIMITATIONS","Synthetic / No Model","No clinical calibration"),
           ("04","MANUAL REVIEW","REQUIRED","Not an approval")]
    for i,(num,title,sub,more) in enumerate(cards):
        x=108+(i%2)*885; y=266+(i//2)*308
        elem("rect",{"x":str(x),"y":str(y),"width":"819","height":"258","rx":"24","fill":"#FFFDFC","stroke":"#D7BEBF","stroke-width":"3"})
        elem("circle",{"cx":str(x+72),"cy":str(y+70),"r":"38","fill":"#D7B5B4"})
        elem("text",{"x":str(x+46),"y":str(y+84),"font-size":"35","font-weight":"700","fill":"#4E4056"},num)
        elem("text",{"x":str(x+135),"y":str(y+86),"font-size":"46","font-weight":"700","fill":"#323C50"},title)
        elem("text",{"x":str(x+49),"y":str(y+157),"font-size":"32","fill":"#5C6370"},sub)
        elem("text",{"x":str(x+49),"y":str(y+204),"font-size":"28","fill":"#36796C"},more)
    elem("rect",{"x":"107","y":"897","width":"1705","height":"104","rx":"15","fill":"#F9EBE7"})
    elem("text",{"x":"143","y":"952","font-size":"33","font-weight":"700","fill":"#884A60"},
         "TOY_NOT_GRCH38  |  No Live API  |  No Patient Data  |  No Diagnosis")
    return ET.tostring(r,encoding="unicode")+"\n"


def html_report(report:dict)->str:
    h=lambda v:escape(str(v),quote=True)
    rows="".join(f"<tr><td><code>{h(r['variant_id'])}</code></td><td>{('null' if r['synthetic_score'] is None else h(r['synthetic_score']))}</td>"
                 f"<td>{h(r['score_status'])}</td><td>REVIEW_REQUIRED</td></tr>" for r in report["annotations"])
    limits="".join(f"<li><b>{h(x['code'])}</b> — {h(x['text'])}</li>" for x in report["model_limitations"])
    return f'''<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Day28 | Evidence Reporter</title>
<style>body{{margin:0;background:#E5DEDA;color:#344256;font:17px/1.65 system-ui,"Noto Sans TC",sans-serif}}
main{{max-width:1120px;margin:auto;padding:34px 20px 60px}}section,header{{background:#fffaf8;border-radius:18px;padding:23px;margin:14px 0}}
h1{{margin:0}}.warn{{border-left:6px solid #b06b72;background:#fff2ed;padding:17px;border-radius:10px}}
.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(165px,1fr));gap:13px}}.tile{{border:1px solid #eee0db;border-radius:12px;padding:16px}}
.big{{font-size:2em;font-weight:800;color:#316d6a}}table{{width:100%;border-collapse:collapse}}td,th{{padding:10px;border-bottom:1px solid #e8dbd7;text-align:left}}
th{{background:#f2e7e2}}img{{width:100%;height:auto}}a{{color:#365b86}}small{{color:#666}}</style></head><body><main>
<header><h1>Day 28｜限制與證據追蹤</h1><p>學會說清楚我們不知道什麼 — OFFLINE 教學報告</p></header>
<div class="warn"><b>SYNTHETIC EDUCATIONAL DEMO · TOY_NOT_GRCH38</b><br>這些不是 AlphaGenome 真實輸出、病人資料或診斷。SHA-256 僅核對本地檔案完整性，並非權威來源認證。</div>
<section><div class="grid"><div class="tile"><div class="big">{h(report['source_provenance']['artifact_count_verified'])}</div>已核對檔案</div>
<div class="tile"><div class="big">{h(report['pipeline_counts']['candidate_count'])}</div>候選</div>
<div class="tile"><div class="big">{h(report['pipeline_counts']['missing_score_count'])}</div>缺人工分數</div>
<div class="tile"><div class="big">REQUIRED</div>人工審閱</div></div></section>
<section><h2>從 Day 27 到證據報告</h2><img src="day28_evidence_map.svg" alt="Provenance, missing values, limitations and manual review"><p>Day 27 來源完整性：{h(report['source_provenance']['hash_check'])}；但未驗證原生模型或資料來源真實性。</p></section>
<section><h2>QC 後的人工候選</h2><table><thead><tr><th>Variant</th><th>人工分數</th><th>分數狀態</th><th>審閱</th></tr></thead><tbody>{rows}</tbody></table></section>
<section><h2>模型／證據限制</h2><ul>{limits}</ul></section>
<section><a href="day28_evidence_report.json">完整 JSON</a> · <a href="day28_evidence_report.md">繁中 Markdown</a> · <a href="day28_evidence_table.csv">CSV</a></section>
<small>Research Use Only · OFFLINE_ONLY · No clinical classification · Real API calls: 0</small>
</main></body></html>'''


def report_csv(report:dict)->str:
    import io
    buf=io.StringIO()
    writer=csv.DictWriter(buf,fieldnames=["variant_id","synthetic_score","score_status","review_status","evidence_class"],lineterminator="\n")
    writer.writeheader()
    for r in report["annotations"]:
        writer.writerow({**r,"synthetic_score":"" if r["synthetic_score"] is None else r["synthetic_score"],
                         "evidence_class":"SYNTHETIC_FIXTURE"})
    return buf.getvalue()


def sha_manifest(folder:Path)->dict:
    return {name:{"sha256":sha256(folder/name),"bytes":(folder/name).stat().st_size} for name in ARTIFACTS}


def write_outputs(run_dir:Path,output_dir:Path)->dict:
    report=report_from_day27(Path(run_dir))
    target=Path(output_dir).expanduser().absolute()
    if target.exists() or target.is_symlink():
        raise ValueError("Output directory already exists; choose a new folder, never overwrite")
    target.parent.mkdir(parents=True,exist_ok=True)
    temp=Path(tempfile.mkdtemp(prefix=".day28-working-",dir=str(target.parent)))
    try:
        (temp/ARTIFACTS[0]).write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
        (temp/ARTIFACTS[1]).write_text(markdown_report(report),encoding="utf-8")
        (temp/ARTIFACTS[2]).write_text(html_report(report),encoding="utf-8")
        (temp/ARTIFACTS[3]).write_text(draw_evidence_svg(report),encoding="utf-8")
        (temp/ARTIFACTS[4]).write_text(report_csv(report),encoding="utf-8")
        manifest={"schema_version":"day28.synthetic.report_manifest.v1","mode":MODE,
                  "artifact_count":len(ARTIFACTS),"artifacts":sha_manifest(temp),
                  "source_day27_manifest_sha256":report["source_provenance"]["day27_run_manifest_sha256"],
                  "signed":False,"real_api_calls_this_run":0}
        (temp/"day28_report_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
        if target.exists():raise ValueError("Target appeared during report creation; refusing overwrite")
        temp.rename(target)
    finally:
        if temp.exists():shutil.rmtree(temp)
    return {"status":"PASS","mode":MODE,"variant_id":report["variant_id"],
            "source_artifacts_verified":report["source_provenance"]["artifact_count_verified"],
            "candidate_count":report["pipeline_counts"]["candidate_count"],
            "missing_score_count":report["pipeline_counts"]["missing_score_count"],
            "output_artifact_count":len(ARTIFACTS),"manual_review":"REQUIRED",
            "report_file":str(target/ARTIFACTS[0]),"html_file":str(target/ARTIFACTS[2]),
            "real_api_calls_this_run":0}


def verify_report(run_dir:Path,report_dir:Path)->dict:
    report_dir=Path(report_dir)
    if not report_dir.is_dir() or report_dir.is_symlink():raise ValueError("Report dir missing/symlinked")
    manifest=d27.load_json(report_dir/"day28_report_manifest.json")
    if (manifest.get("schema_version")!="day28.synthetic.report_manifest.v1" or manifest.get("mode")!=MODE
        or manifest.get("signed") is not False or manifest.get("real_api_calls_this_run")!=0
        or manifest.get("artifact_count")!=len(ARTIFACTS)
        or not isinstance(manifest.get("artifacts"),dict)
        or set(manifest["artifacts"])!=set(ARTIFACTS)):
        raise ValueError("Invalid Day28 report manifest")
    for name,record in manifest["artifacts"].items():
        p=PurePosixPath(name); path=report_dir/name
        if (p.is_absolute() or ".." in p.parts or path.is_symlink()
            or not path.is_file() or record.get("bytes")!=path.stat().st_size
            or not SHA_RE.fullmatch(str(record.get("sha256"))) or sha256(path)!=record["sha256"]):
            raise ValueError("Day28 report checksum mismatch: "+name)
    original=report_from_day27(Path(run_dir))
    saved=d27.load_json(report_dir/"day28_evidence_report.json")
    if saved!=original or manifest.get("source_day27_manifest_sha256")!=original["source_provenance"]["day27_run_manifest_sha256"]:
        raise ValueError("Report no longer matches verified Day27 evidence")
    return {"status":"PASS","mode":MODE,"report_artifact_count":len(ARTIFACTS),
            "source_artifact_count":original["source_provenance"]["artifact_count_verified"],
            "checksums":"MATCH","source_alignment":"MATCH","real_api_calls_this_run":0}


def main(argv=None)->int:
    parser=argparse.ArgumentParser(description="Day28 OFFLINE Evidence Reporter (synthetic only; no clinical use)")
    subs=parser.add_subparsers(dest="command",required=True)
    run=subs.add_parser("run",help="Create JSON, Markdown, CSV, SVG, HTML evidence reports")
    run.add_argument("--mode",choices=["offline"],required=True)
    run.add_argument("--run-dir",default="example_output/day27_output")
    run.add_argument("--output-dir",default="day28_output")
    check=subs.add_parser("verify",help="Verify Day28 artifacts against the Day27 source")
    check.add_argument("--run-dir",default="example_output/day27_output")
    check.add_argument("--report-dir",default="day28_output")
    args=parser.parse_args(argv)
    try:
        result=(write_outputs(Path(args.run_dir),Path(args.output_dir)) if args.command=="run"
                else verify_report(Path(args.run_dir),Path(args.report_dir)))
    except (ValueError,KeyError,TypeError,OSError,RuntimeError,IndexError) as exc:
        parser.exit(2,f"DAY28_OFFLINE_ERROR: {exc}\n")
    print(json.dumps(result,ensure_ascii=False,sort_keys=True,allow_nan=False))
    return 0


if __name__=="__main__":
    sys.exit(main())
