"""Day 29: local-only, synthetic-only privacy review learning exercise.

Consumes the validated Day 28 synthetic report + UNSIGNED hash manifest.
Does NOT anonymize real genomic data, certify consent, or authorize disclosure.
No third-party packages, AlphaGenome SDK, credentials or network calls.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import secrets
import sys
from pathlib import Path
from xml.etree import ElementTree as ET

DAY28_SCHEMA = "day28.synthetic.evidence_report.v1"
DAY28_MANIFEST_SCHEMA = "day28.synthetic.report_manifest.v1"
REPORT_SCHEMA = "day29.synthetic.privacy_gate.v1"
VARIANT_RE = re.compile(r"chr(?:[1-9]|1[0-9]|2[0-2]|X|Y|M):[1-9][0-9]*:[ACGT]+>[ACGT]+\Z")
DEMO_TOKEN_RE = re.compile(r"DEMO_[0-9a-f]{32}\Z")
SOURCE_KEYS = frozenset({
    "schema_version", "status", "evidence_class", "mode", "assembly", "variant_id",
    "source_provenance", "pipeline_counts", "annotations", "missing_values",
    "teaching_observations", "model_limitations", "manual_review",
    "confidence_calibrated", "clinical_classification", "disease_probability",
    "model_prediction", "atlas_query_performed", "clinical_interpretation",
    "real_api_calls_this_run",
})
EXPECTED_PROVENANCE = frozenset({
    "source_kind", "source_schema_version", "source_bundle_sha256",
    "day27_summary_sha256", "day27_run_manifest_sha256", "artifact_count_verified",
    "hash_check", "signature_verified", "native_alpha_genome_response_retained",
})
EXPECTED_COUNTS = frozenset({
    "candidate_count", "qc_selected", "qc_review_required", "qc_rejected",
    "annotated_count", "missing_score_count", "stage_count",
})
EXPECTED_ANNOTATION = frozenset({"variant_id", "synthetic_score", "score_status", "review_status"})
BANNED_KEYS = frozenset({
    "name", "full_name", "patient_name", "patient_id", "mrn", "medical_record_number",
    "date_of_birth", "dob", "birth_date", "age", "phone", "phone_number", "email",
    "address", "hospital_id", "identity_number", "national_id", "family_history",
    "clinical_notes", "sex", "gender", "consent", "api_key", "token", "password",
})
ARTIFACT_FILES = (
    "day29_privacy_audit.json", "day29_local_preview.json",
    "day29_privacy_report.md", "day29_privacy_dashboard.html",
    "day29_privacy_flow.svg", "day29_release_decision.csv",
)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _reject_direct_identifier_keys(obj: object) -> None:
    """Demonstration: reject some known column names, not a PHI detector."""
    if isinstance(obj, dict):
        for key, value in obj.items():
            if isinstance(key, str) and key.casefold() in BANNED_KEYS:
                raise ValueError("DIRECT_IDENTIFIER_FIELD_NOT_ALLOWED")
            _reject_direct_identifier_keys(value)
    elif isinstance(obj, list):
        for item in obj:
            _reject_direct_identifier_keys(item)


def _parse_json(path: Path) -> tuple[dict, bytes]:
    if not path.is_file() or path.is_symlink() or path.stat().st_size > 1024 * 1024:
        raise ValueError("MISSING_SYMLINK_OR_OVERSIZED_INPUT")
    data = path.read_bytes()
    try:
        parsed = json.loads(data)
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError("INVALID_JSON") from exc
    if not isinstance(parsed, dict):
        raise ValueError("EXPECTED_JSON_OBJECT")
    return parsed, data


def _number(value: object) -> bool:
    return type(value) in (int, float) and math.isfinite(value)


def validate_day28(report: dict) -> dict:
    """Fail closed for the *specified synthetic teaching schema*.

    A self-described synthetic flag cannot establish whether source data are
    genuinely fictitious. Do NOT use this function on real patient data.
    """
    _reject_direct_identifier_keys(report)
    if set(report) != SOURCE_KEYS:
        raise ValueError("UNEXPECTED_DAY28_FIELDS")
    expected = {
        "schema_version": DAY28_SCHEMA,
        "status": "PASS", "evidence_class": "SYNTHETIC_EDUCATIONAL_DEMO",
        "mode": "OFFLINE_ONLY", "assembly": "TOY_NOT_GRCH38",
        "manual_review": "REQUIRED", "confidence_calibrated": False,
        "clinical_classification": None, "disease_probability": None,
        "model_prediction": False, "atlas_query_performed": False,
        "clinical_interpretation": False, "real_api_calls_this_run": 0,
    }
    for key, expected_value in expected.items():
        if type(report[key]) is not type(expected_value) or report[key] != expected_value:
            raise ValueError("NOT_STRICTLY_SYNTHETIC_OFFLINE_DAY28")
    if not isinstance(report["variant_id"], str) or VARIANT_RE.fullmatch(report["variant_id"]) is None:
        raise ValueError("BAD_TEACHING_VARIANT_ID")

    p = report["source_provenance"]
    if not isinstance(p, dict) or set(p) != EXPECTED_PROVENANCE:
        raise ValueError("BAD_SOURCE_PROVENANCE_SCHEMA")
    for key, val in {
        "source_kind": "DAY27_OFFLINE_SYNTHETIC_PIPELINE",
        "source_schema_version": "day27.synthetic.pipeline.v1",
        "artifact_count_verified": 23, "hash_check": "MATCH",
        "signature_verified": False, "native_alpha_genome_response_retained": False,
    }.items():
        if type(p[key]) is not type(val) or p[key] != val:
            raise ValueError("BAD_SYNTHETIC_SOURCE_PROVENANCE")
    for key in ("source_bundle_sha256", "day27_summary_sha256", "day27_run_manifest_sha256"):
        if not isinstance(p[key], str) or re.fullmatch(r"[0-9a-f]{64}", p[key]) is None:
            raise ValueError("INVALID_SHA256_LABEL")

    c = report["pipeline_counts"]
    if not isinstance(c, dict) or set(c) != EXPECTED_COUNTS:
        raise ValueError("BAD_COUNTS_SCHEMA")
    if any(type(val) is not int or val < 0 for val in c.values()):
        raise ValueError("BAD_COUNT_VALUE")
    if (c["candidate_count"] != c["qc_selected"] + c["qc_review_required"] + c["qc_rejected"]
            or c["qc_selected"] != c["annotated_count"] + c["missing_score_count"]
            or c["stage_count"] != 6):
        raise ValueError("INCONSISTENT_SYNTHETIC_COUNTS")
    a = report["annotations"]
    if not isinstance(a, list) or len(a) != c["qc_selected"]:
        raise ValueError("INCORRECT_ANNOTATION_ROWS")
    identifiers: list[str] = []
    missing = 0
    annotated = 0
    for entry in a:
        if not isinstance(entry, dict) or set(entry) != EXPECTED_ANNOTATION:
            raise ValueError("BAD_ANNOTATION_SCHEMA")
        variant = entry["variant_id"]
        if not isinstance(variant, str) or VARIANT_RE.fullmatch(variant) is None:
            raise ValueError("BAD_ANNOTATION_VARIANT")
        identifiers.append(variant)
        if entry["review_status"] != "REVIEW_REQUIRED":
            raise ValueError("REVIEW_FLAG_MISSING")
        score = entry["synthetic_score"]
        if score is None:
            missing += 1
            if entry["score_status"] != "MISSING_SCORE":
                raise ValueError("MISSING_SCORE_MISMATCH")
        elif _number(score):
            annotated += 1
            if entry["score_status"] != "SYNTHETIC_SCORE":
                raise ValueError("SCORE_STATUS_MISMATCH")
        else:
            raise ValueError("INVALID_SYNTHETIC_SCORE")
    if len(set(identifiers)) != len(identifiers) or report["variant_id"] not in identifiers:
        raise ValueError("INVALID_ANNOTATION_IDENTIFIERS")
    if missing != c["missing_score_count"] or annotated != c["annotated_count"]:
        raise ValueError("SCORE_COUNTS_MISMATCH")
    m = report["missing_values"]
    if not isinstance(m, dict) or set(m) != {
        "missing_score_variant_ids", "missing_is_zero", "zero_is_a_valid_numeric_value"
    }:
        raise ValueError("MISSING_VALUES_SCHEMA")
    if (type(m["missing_is_zero"]) is not bool or m["missing_is_zero"] is not False
            or m["zero_is_a_valid_numeric_value"] is not True
            or not isinstance(m["missing_score_variant_ids"], list)
            or set(m["missing_score_variant_ids"]) != {
                x["variant_id"] for x in a if x["synthetic_score"] is None
            }):
        raise ValueError("MISSING_VALUES_CONTRACT")

    t = report["teaching_observations"]
    if (not isinstance(t, dict) or t.get("combined_biological_score", "BAD") is not None
            or t.get("mechanism_status") != "NOT_VALIDATED"):
        raise ValueError("UNVERIFIED_MECHANISM_NOT_LABELLED")
    limitations = report["model_limitations"]
    if (not isinstance(limitations, list) or len(limitations) < 4
            or any(not isinstance(item, dict)
                   or set(item) != {"code", "text"}
                   or not isinstance(item["code"], str)
                   or not isinstance(item["text"], str) for item in limitations)):
        raise ValueError("INVALID_LIMITATIONS")
    required = {"SYNTHETIC_ONLY", "PROVENANCE_NOT_SIGNED", "MANUAL_REVIEW_REQUIRED"}
    if not required.issubset({item["code"] for item in limitations}):
        raise ValueError("REQUIRED_LIMITATIONS_MISSING")
    return {"candidate_count": c["candidate_count"],
            "selected_count": c["qc_selected"],
            "missing_score_count": c["missing_score_count"],
            "day27_artifact_count_claimed_by_day28": p["artifact_count_verified"]}


def load_and_validate(report_path: Path, manifest_path: Path) -> tuple[dict, dict]:
    manifest, _ = _parse_json(manifest_path)
    if (manifest.get("schema_version") != DAY28_MANIFEST_SCHEMA
            or manifest.get("mode") != "OFFLINE_ONLY"
            or manifest.get("signed") is not False
            or manifest.get("real_api_calls_this_run") != 0
            or manifest.get("artifact_count") != 5):
        raise ValueError("BAD_DAY28_UNSIGNED_MANIFEST")
    entries = manifest.get("artifacts")
    if not isinstance(entries, dict) or len(entries) != 5:
        raise ValueError("BAD_DAY28_ARTIFACT_MANIFEST")
    entry = entries.get("day28_evidence_report.json")
    if (not isinstance(entry, dict) or set(entry) != {"sha256", "bytes"}
            or type(entry["bytes"]) is not int):
        raise ValueError("MISSING_DAY28_REPORT_DIGEST")
    report, data = _parse_json(report_path)
    if sha256_bytes(data) != entry["sha256"] or len(data) != entry["bytes"]:
        raise ValueError("DAY28_INPUT_HASH_MISMATCH")
    counts = validate_day28(report)
    return counts, {
        "day28_report_sha256": sha256_bytes(data),
        "day28_manifest_sha256": sha256_bytes(manifest_path.read_bytes()),
        "day28_manifest_signed": False,
        "day27_source_artifacts_independently_checked_this_run": False,
        "day28_manifest_report_entry_matches": True,
    }


def build_privacy_audit(counts: dict, provenance: dict) -> dict:
    return {
        "schema_version": REPORT_SCHEMA,
        "status": "PASS_FOR_LOCAL_SYNTHETIC_EDUCATION_ONLY",
        "evidence_class": "SYNTHETIC_EDUCATIONAL_DEMO",
        "mode": "OFFLINE_ONLY",
        "local_demo_token": "DEMO_" + secrets.token_hex(16),
        "token_is_anonymization": False,
        "anonymization_claim": "NOT_ESTABLISHED",
        "pseudonymization_claim": "NOT_ESTABLISHED",
        "reidentification_risk_eliminated": False,
        "source": provenance,
        "aggregate_counts": dict(counts),
        "suppressed_from_preview": [
            "variant_ids", "per_variant_scores", "phenotype_details", "local_demo_token"
        ],
        "demo_policy": {
            "purpose": "LOCAL_EDUCATION",
            "external_transfer_allowed": False,
            "publish_allowed": False,
            "real_patient_processing_allowed": False,
            "consent_verified": False,
            "irb_or_ethics_review_verified": False,
            "human_approval_completed": False,
        },
        "governance_review_status": "MANUAL_REVIEW_REQUIRED",
        "release_gate": "BLOCKED",
        "risk_notes": [
            "Self-declared synthetic labels do not prove that a dataset is fictitious.",
            "Random tokens remove neither genetic uniqueness nor potential family linkability.",
            "SHA-256 can check local file consistency, not author identity or source authenticity.",
            "No consent, ethics approval, clinical use, external transfer or legal compliance is established.",
        ],
        "real_api_calls_this_run": 0,
    }


def build_minimal_preview(audit: dict) -> dict:
    counts = audit["aggregate_counts"]
    return {
        "data_kind": "SYNTHETIC_EDUCATIONAL_DEMO",
        "mode": "OFFLINE_ONLY",
        "candidate_count": counts["candidate_count"],
        "selected_count": counts["selected_count"],
        "missing_score_count": counts["missing_score_count"],
        "per_variant_data_included": False,
        "random_token_included": False,
        "anonymization_claim": "NOT_ESTABLISHED",
        "external_transfer_allowed": False,
        "publish_allowed": False,
        "governance_review_status": "MANUAL_REVIEW_REQUIRED",
        "real_api_calls_this_run": 0,
    }


def format_markdown(audit: dict, preview: dict) -> str:
    c = audit["aggregate_counts"]
    return (
        "# Day 29｜本機合成資料隱私與發布審查紀錄\n\n"
        "> **OFFLINE SYNTHETIC DEMO — 外傳與公開發布仍被阻擋。**\n\n"
        "## 1. 這份報告能證明什麼？\n\n"
        f"- 來源 Day 28 JSON SHA-256：`{audit['source']['day28_report_sha256']}`。\n"
        "- 只核對本機、未簽章的 Day 28 manifest 裡該 JSON 的雜湊；**不能證明資料來源真實**。\n"
        "- Day 27 的 23 個輸出來自 Day 28 記錄，**本次並未逐一重驗**。\n\n"
        "## 2. 資料最小化\n\n"
        f"- 合成候選：{c['candidate_count']}；QC 通過：{c['selected_count']}；缺分數：{c['missing_score_count']}。\n"
        "- 預覽省略逐筆變異、個別分數與本機隨機 Token。\n"
        "- `null` 是缺值，不應寫成 0 分。\n\n"
        "## 3. 隱私與研究治理關卡\n\n"
        "- `external_transfer_allowed=false`；`publish_allowed=false`。\n"
        "- `MANUAL_REVIEW_REQUIRED`；`release_gate=BLOCKED`。\n"
        "- 隨機 Token 不衍生自姓名或病歷號，但**不保證基因資料無法被重新識別**。\n"
        "- `SYNTHETIC` 標籤、SHA、布林旗標皆不能代替同意、倫理審查與法遵判定。\n\n"
        "## 4. 對照 Day 28 的報告\n\n"
        "- Day 28 談模型證據、缺值與解讀限制。\n"
        "- Day 29 加入資料是否應進入外部系統的治理問題，但**本教學不執行外傳**。\n\n"
        "## 5. 未來十年：一張有待驗證的研究路線圖\n\n"
        "- 近期：資料最小化、來源紀錄、可審計的人工判斷。\n"
        "- 中期：探索跨院可信任運算、聯邦學習與其隱私攻擊防護。\n"
        "- 長期：在功能實驗、族群驗證與臨床研究支持下，讓模型結果成為可審查的研究線索。\n\n"
        "這些都是**願景**，不代表產品已部署、符合任何法規或可用於醫療決策。\n\n"
        "---\n\n*僅本機合成教學；非匿名化服務、臨床工具或資安認證。*\n"
    )


def format_html(preview: dict) -> str:
    # Only static English/Chinese literals and validated nonnegative integers enter HTML.
    count, selected, missing = (preview[k] for k in
        ("candidate_count", "selected_count", "missing_score_count"))
    return f'''<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8">
<title>Day 29｜Synthetic Privacy Gate</title>
<style>
body{{background:#E5DEDA;color:#2f3036;font-family:system-ui,"Noto Sans TC",sans-serif;margin:0;padding:45px}}
main{{max-width:1150px;margin:auto;background:#fffaf8;border-radius:25px;padding:40px;box-shadow:0 10px 32px #593c3120}}
h1{{font-size:32px;margin:0}} p{{line-height:1.75}} .pill{{display:inline-block;border-radius:100px;background:#f5e8db;padding:7px 15px;margin:14px 0}}
.cards{{display:grid;grid-template-columns:repeat(3,1fr);gap:16px;margin:26px 0}}
.card{{background:#f7f3f0;border-radius:16px;padding:22px}} .big{{font-size:38px;font-weight:700}}
.block{{background:#fff0ee;padding:24px;border-left:5px solid #b76764;border-radius:12px}}
small{{color:#5f6067}}
</style></head><body><main><h1>Day 29｜基因隱私與證據發布關卡</h1>
<span class="pill">SYNTHETIC EDUCATIONAL DEMO · OFFLINE ONLY</span>
<p>接續 Day 28 的證據報告：問的不只是「分析得到了什麼」，還有「資料應否被傳輸、分享或拿來做結論」。</p>
<section class="cards"><div class="card"><div class="big">{count}</div>人工候選</div>
<div class="card"><div class="big">{selected}</div>QC 通過（非診斷）</div>
<div class="card"><div class="big">{missing}</div>保留 null 的缺值</div></section>
<div class="block"><strong>外傳：禁止｜公開：禁止｜人工審閱：尚未完成</strong>
<p>這是程式設計示範，<strong>不是匿名化或法遵認證</strong>。隨機 Token 不能消除真實基因資料的再識別風險。</p></div>
<p><strong>目前顯示的只是合成彙總。</strong>不顯示變異 ID、逐筆分數、表型或隨機 Token。</p>
<small>SHA-256 只協助核對本機檔案與未簽章 manifest 的內容一致性。Day 27 的 23 個檔案並未由本篇再驗證。</small>
</main></body></html>'''


def format_svg(preview: dict) -> str:
    ns = "http://www.w3.org/2000/svg"
    svg = ET.Element("svg", xmlns=ns, width="1920", height="1080", viewBox="0 0 1920 1080")
    ET.SubElement(svg, "rect", x="0", y="0", width="1920", height="1080", fill="#E5DEDA")
    ET.SubElement(svg, "rect", x="90", y="75", width="1740", height="900", rx="50", fill="#fbf9f8")
    def label(x: int, y: int, value: str, size: int = 35, fill: str = "#363941", weight: str = "500") -> None:
        node = ET.SubElement(svg,"text",x=str(x),y=str(y),fill=fill,
            attrib={"font-family":"Noto Sans CJK TC", "font-size":str(size),"font-weight":weight})
        node.text=value
    label(145,160,"DAY 29  |  Privacy & Evidence Gate",54,"#344454","700")
    label(145,217,"AlphaGenome Atlas · SYNTHETIC EDUCATIONAL DEMO",29,"#787078")
    stages=[("01","合成資料標示","只接受教學標記"),("02","完整性核對","本地 SHA 相符"),
      ("03","最小化預覽","不輸出個別變異"),("04","人工審閱","仍未獲得批准"),("05","外部傳輸","預設拒絕")]
    for i,(number,title,description) in enumerate(stages):
        x=146+i*334
        ET.SubElement(svg,"rect",x=str(x),y="326",width="309",height="286",rx="25",fill=["#e9eff3","#e9f0e9","#f4ede7","#eeeafb","#f8e5e2"][i])
        label(x+26,379,number,27,"#62707d","700")
        label(x+26,475,title,37,"#2f3944","700")
        label(x+26,535,description,23,"#61636c")
        if i<4: label(x+315,470,"›",36,"#ba8d80","700")
    ET.SubElement(svg,"rect",x="145",y="672",width="1630",height="170",rx="25",fill="#f7e4df")
    label(190,743,"NOT ANONYMIZED · NO REAL API · NO CLINICAL DECISION",38,"#9e4c4d","700")
    label(190,805,"Random token ≠ zero re-identification risk; human review is still required.",27,"#6b5a5b")
    label(145,922,"All values are synthetic. No patient data are accepted as valid evidence of safety.",25,"#777177")
    return ET.tostring(svg, encoding="unicode")


def run(report_path: Path, day28_manifest: Path, output_dir: Path) -> dict:
    counts, prov = load_and_validate(report_path, day28_manifest)
    audit = build_privacy_audit(counts, prov)
    preview = build_minimal_preview(audit)
    # Refuse to overwrite existing output, even if it contains previous synthetic results.
    if output_dir.exists():
        raise ValueError("OUTPUT_DIR_EXISTS_NO_OVERWRITE")
    output_dir.mkdir(parents=True, exist_ok=False)
    texts = {
        "day29_privacy_audit.json": json.dumps(audit, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        "day29_local_preview.json": json.dumps(preview, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        "day29_privacy_report.md": format_markdown(audit, preview),
        "day29_privacy_dashboard.html": format_html(preview),
        "day29_privacy_flow.svg": format_svg(preview),
        "day29_release_decision.csv": (
            "data_kind,mode,release_gate,external_transfer_allowed,human_approval_completed,anonymization_claim\n"
            "SYNTHETIC_EDUCATIONAL_DEMO,OFFLINE_ONLY,BLOCKED,false,false,NOT_ESTABLISHED\n"
        ),
    }
    for filename, content in texts.items():
        (output_dir / filename).write_text(content, encoding="utf-8")
    manifest = {
        "schema_version": "day29.synthetic.report_manifest.v1",
        "signed": False,
        "mode": "OFFLINE_ONLY",
        "source_report_sha256": prov["day28_report_sha256"],
        "source_manifest_sha256": prov["day28_manifest_sha256"],
        "artifact_count": len(texts),
        "artifacts": {
            fn: {"sha256": sha256_bytes((output_dir/fn).read_bytes()),
                 "bytes": (output_dir/fn).stat().st_size}
            for fn in ARTIFACT_FILES
        },
        "real_api_calls_this_run": 0,
    }
    (output_dir / "day29_report_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    return {
        "status": "PASS", "mode": "OFFLINE_ONLY", "source_day28_sha_match": True,
        "candidate_count": counts["candidate_count"], "selected_count": counts["selected_count"],
        "missing_score_count": counts["missing_score_count"],
        "output_artifact_count": len(ARTIFACT_FILES),
        "external_transfer_allowed": False, "publish_allowed": False,
        "governance_review_status": "MANUAL_REVIEW_REQUIRED", "release_gate": "BLOCKED",
        "real_api_calls_this_run": 0,
    }


def verify(report_path: Path, day28_manifest: Path, output_dir: Path) -> dict:
    counts, prov = load_and_validate(report_path, day28_manifest)
    manifest, _ = _parse_json(output_dir / "day29_report_manifest.json")
    if (manifest.get("schema_version") != "day29.synthetic.report_manifest.v1"
            or manifest.get("signed") is not False
            or manifest.get("mode") != "OFFLINE_ONLY"
            or manifest.get("real_api_calls_this_run") != 0
            or manifest.get("source_report_sha256") != prov["day28_report_sha256"]
            or manifest.get("source_manifest_sha256") != prov["day28_manifest_sha256"]
            or manifest.get("artifact_count") != len(ARTIFACT_FILES)
            or set(manifest.get("artifacts", {})) != set(ARTIFACT_FILES)):
        raise ValueError("DAY29_MANIFEST_MISMATCH")
    for filename in ARTIFACT_FILES:
        path = output_dir / filename
        if not path.is_file() or path.is_symlink():
            raise ValueError("DAY29_ARTIFACT_MISSING")
        digest = sha256_bytes(path.read_bytes())
        entry = manifest["artifacts"][filename]
        if digest != entry["sha256"] or path.stat().st_size != entry["bytes"]:
            raise ValueError("DAY29_ARTIFACT_HASH_MISMATCH")
    audit, _ = _parse_json(output_dir / "day29_privacy_audit.json")
    preview, _ = _parse_json(output_dir / "day29_local_preview.json")
    if (audit.get("source") != prov or audit.get("aggregate_counts") != counts
            or not isinstance(audit.get("local_demo_token"), str)
            or DEMO_TOKEN_RE.fullmatch(audit["local_demo_token"]) is None
            or audit.get("release_gate") != "BLOCKED"
            or audit.get("demo_policy", {}).get("external_transfer_allowed") is not False
            or audit.get("demo_policy", {}).get("publish_allowed") is not False
            or preview != build_minimal_preview(audit)
            or preview.get("per_variant_data_included") is not False
            or preview.get("random_token_included") is not False):
        raise ValueError("DAY29_PRIVACY_GATE_OUTPUT_MISMATCH")
    for forbidden in ("chr1:", "DEMO_", "0.42", "-0.17"):
        if forbidden in (output_dir / "day29_local_preview.json").read_text(encoding="utf-8"):
            raise ValueError("PREVIEW_CONTAINS_PER_VARIANT_CONTENT")
    return {"status":"PASS","verified_artifacts":len(ARTIFACT_FILES),
            "source_day28_sha_match":True,"governance_review_status":"MANUAL_REVIEW_REQUIRED",
            "release_gate":"BLOCKED","real_api_calls_this_run":0}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Local, synthetic-only Day 29 privacy and disclosure lab")
    sub = parser.add_subparsers(dest="command", required=True)
    for command in ("run", "verify"):
        sp = sub.add_parser(command)
        sp.add_argument("--mode",choices=["offline"], default="offline")
        sp.add_argument("--day28-report",default="fixtures/day28_evidence_report.json")
        sp.add_argument("--day28-manifest",default="fixtures/day28_report_manifest.json")
        sp.add_argument("--output-dir",default="day29_output")
    args = parser.parse_args(argv)
    try:
        if args.command == "run":
            result = run(Path(args.day28_report), Path(args.day28_manifest), Path(args.output_dir))
        else:
            result = verify(Path(args.day28_report), Path(args.day28_manifest), Path(args.output_dir))
    except (ValueError, KeyError, TypeError, OSError) as exc:
        # Only an internal controlled error code, never source patient contents.
        print(f"Day29 privacy lab failed: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result,ensure_ascii=False,sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
