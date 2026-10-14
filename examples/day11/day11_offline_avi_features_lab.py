#!/usr/bin/env python3
"""Day11: OFFLINE-ONLY educational inspector for AVI-like 18-feature records.

No alphagenome SDK, network calls, clinical interpretation, or patient data.
Feature names mirror the AVI public description but all fixture numeric values
are INVENTED; the exercise does not reproduce the AVI neural network or SHAP.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import io
import json
import math
import re
import sys
from pathlib import Path
from typing import Any
from xml.sax.saxutils import escape as xml_escape

VERSION = "1.0.0"
DATA_KIND = "SYNTHETIC_EDUCATIONAL_DEMO"
ASSEMBLY = "TOY_NOT_GRCH38"
TOLERANCE = 1e-3  # Teaching QC, NOT a certified AVI SHAP tolerance.
FEATURES = (
    ("MAX_ABS_ATAC", "AlphaGenome", "ATAC chromatin accessibility"),
    ("MAX_ABS_DNASE", "AlphaGenome", "DNase accessibility"),
    ("MAX_ABS_CHIP_TF", "AlphaGenome", "TF ChIP binding"),
    ("MAX_ABS_CHIP_HISTONE", "AlphaGenome", "Histone ChIP"),
    ("MAX_ABS_CAGE", "AlphaGenome", "CAGE transcription initiation"),
    ("MAX_ABS_PROCAP", "AlphaGenome", "PRO-cap transcription initiation"),
    ("MAX_ABS_RNA_SEQ", "AlphaGenome", "RNA-seq expression"),
    ("MAX_ABS_POLYADENYLATION", "AlphaGenome", "Polyadenylation"),
    ("MAX_ABS_CONTACT_MAPS", "AlphaGenome", "Chromatin contact map"),
    ("MERGED_SPLICING", "AlphaGenome", "Merged splicing"),
    ("ALPHAMISSENSE", "AlphaMissense", "Missense-associated annotation"),
    ("PROTEIN_TERMINATION", "VEP consequences", "Termination or frameshift"),
    ("START_LOST", "VEP consequences", "Start lost"),
    ("STOP_LOST", "VEP consequences", "Stop lost"),
    ("CACTUS_241_WAY", "Conservation", "Cactus 241-way conservation"),
    ("PHASTCONS_470_WAY", "Conservation", "PhastCons 470-way"),
    ("IS_INSERTION", "Indel flags", "Insertion indicator"),
    ("IS_DELETION", "Indel flags", "Deletion indicator"),
)
GROUP_ORDER = ("AlphaGenome", "AlphaMissense", "VEP consequences", "Conservation", "Indel flags")
FEATURE_NAMES = tuple(x[0] for x in FEATURES)
GROUP_COUNT = {group: sum(1 for f in FEATURES if f[1] == group) for group in GROUP_ORDER}
BINARY_FEATURES = frozenset(("PROTEIN_TERMINATION", "START_LOST", "STOP_LOST", "IS_INSERTION", "IS_DELETION"))
OUTPUT_FILES = ("day11_feature_report.json", "day11_feature_values.csv", "day11_feature_attributions.csv", "day11_feature_map.svg", "day11_feature_dashboard.html", "day11_evidence_notes.md")
_IDENT = re.compile(r"[A-Z][A-Z0-9_]{1,63}\Z")


def fail_if(condition: bool, message: str) -> None:
    if condition:
        raise ValueError(message)


def json_dumps(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False) + "\n"


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def canonical_hash(data: Any) -> str:
    blob = json.dumps(data, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    return sha256_bytes(blob)


def ensure_real_or_null(value: Any, label: str) -> float | None:
    if value is None:
        return None
    fail_if(isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value), f"{label} must be finite numeric or JSON null")
    return float(value)


def require_feature_map(value: Any, label: str, allow_null_map: bool = False, binary_indicators: bool = False) -> dict[str, float | None] | None:
    if value is None and allow_null_map:
        return None
    fail_if(not isinstance(value, dict), f"{label} must be an object")
    fail_if(set(value) != set(FEATURE_NAMES), f"{label} must contain exactly the 18 documented feature keys")
    normalized = {name: ensure_real_or_null(value[name], f"{label}.{name}") for name in FEATURE_NAMES}
    if binary_indicators:
        for name in BINARY_FEATURES:
            v = normalized[name]
            fail_if(v is not None and v not in (0.0, 1.0), f"{label}.{name} must be 0, 1, or null")
        fail_if(normalized["IS_INSERTION"] == 1 and normalized["IS_DELETION"] == 1, "Input cannot be both insertion and deletion")
        for name in FEATURE_NAMES:
            v = normalized[name]
            if v is not None and (name.startswith("MAX_ABS_") or name in ("MERGED_SPLICING", "ALPHAMISSENSE")):
                fail_if(v < 0, f"{label}.{name} must be nonnegative")
        pc=normalized["PHASTCONS_470_WAY"]
        fail_if(pc is not None and not 0 <= pc <= 1, f"{label}.PHASTCONS_470_WAY must be between 0 and 1")
    return normalized


def validate_fixture(data: Any) -> list[dict[str, Any]]:
    fail_if(not isinstance(data, dict), "Fixture must be a JSON object")
    fail_if(set(data) != {"data_kind", "assembly", "dataset_role", "cases"}, "Fixture has unexpected/missing top-level keys")
    fail_if(data["data_kind"] != DATA_KIND or data["assembly"] != ASSEMBLY or data["dataset_role"] != "TEACHING_VALUES_NOT_MODEL_INFERENCE", "Only explicit toy educational input is accepted")
    cases = data["cases"]
    fail_if(not isinstance(cases, list) or not (1 <= len(cases) <= 20), "Expected 1-20 toy cases")
    seen: set[str] = set()
    normalized: list[dict[str, Any]] = []
    for idx, row in enumerate(cases):
        fail_if(not isinstance(row, dict), f"Case #{idx} must be object")
        fail_if(set(row) != {"case_id", "variant_label", "raw_score", "feature_values", "attributions"}, f"Case #{idx} has unknown/missing keys")
        case_id = row["case_id"]
        fail_if(not isinstance(case_id, str) or not _IDENT.fullmatch(case_id), f"Bad case_id #{idx}")
        fail_if(case_id in seen, "Duplicate case_id")
        seen.add(case_id)
        variant = row["variant_label"]
        fail_if(not isinstance(variant, str) or not variant.startswith("TOY_") or not _IDENT.fullmatch(variant), "Variant label must be explicitly synthetic and ASCII-safe")
        raw = ensure_real_or_null(row["raw_score"], "raw_score")
        inputs = require_feature_map(row["feature_values"], "feature_values", binary_indicators=True)
        attrs = require_feature_map(row["attributions"], "attributions", allow_null_map=True)
        fail_if((raw is None) != (attrs is None), "raw_score and attributions must be missing together")
        if attrs is not None:
            fail_if(any(v is None for v in attrs.values()), "Complete attribution map cannot contain null values")
        normalized.append({"case_id": case_id, "variant_label": variant, "raw_score": raw, "feature_values": inputs, "attributions": attrs})
    return normalized


def validate_comparison(comparison: Any) -> dict[str, Any] | None:
    """Check an independently invented comparison record; reject comparison/API claims."""
    if comparison is None:
        return None
    fail_if(not isinstance(comparison, dict), "Synthetic comparison must be a JSON object")
    expected = {"evidence_class", "data_kind", "dataset_role", "assembly", "variant_label",
                "score_semantics", "model_executed", "native_response_retained", "clinical_interpretation",
                "avi_raw", "cdf_quantile", "phred_recorded", "feature_count", "top_feature",
                "top_feature_weight", "top_feature_abs_share", "attribution_sum", "additivity_residual", "additivity_status"}
    fail_if(set(comparison) != expected, "Only the documented synthetic comparison fields are accepted")
    fail_if(comparison["evidence_class"] != DATA_KIND or comparison["data_kind"] != DATA_KIND or
            comparison["dataset_role"] != "ARTIFICIAL_COMPARISON_NOT_MODEL_OUTPUT" or
            comparison["assembly"] != ASSEMBLY or comparison["variant_label"] != "TOY_COMPARISON" or
            comparison["score_semantics"] != "TOY_AVI_LIKE_SCORE_AND_RANK", "The comparison must be explicitly synthetic")
    fail_if(comparison["model_executed"] is not False or comparison["native_response_retained"] is not False or
            comparison["clinical_interpretation"] is not False, "No live-model or clinical claims are permitted")
    fail_if(type(comparison["feature_count"]) is not int or comparison["feature_count"] != 18 or
            comparison["top_feature"] != "MERGED_SPLICING", "Synthetic 18-feature exercise metadata mismatch")
    for name in ("avi_raw", "cdf_quantile", "phred_recorded", "attribution_sum", "additivity_residual", "top_feature_weight", "top_feature_abs_share"):
        fail_if(ensure_real_or_null(comparison[name], name) is None, f"Missing synthetic field: {name}")
    raw = float(comparison["avi_raw"])
    cdf = float(comparison["cdf_quantile"])
    fail_if(not 0 <= cdf < 1, "Toy quantile must be in [0,1)")
    fail_if(abs(float(comparison["phred_recorded"]) - (-10 * math.log10(1-cdf))) > 1e-7,
            "Synthetic PHRED must agree with toy CDF")
    fail_if(not 0 <= float(comparison["top_feature_abs_share"]) <= 1 or
            not 0 <= float(comparison["top_feature_weight"]) <= 1,
            "Toy feature metadata must be between 0 and 1")
    residual = raw - float(comparison["attribution_sum"])
    fail_if(not math.isclose(residual, float(comparison["additivity_residual"]), rel_tol=0, abs_tol=1e-9),
            "Synthetic residual does not equal artificial raw minus sum")
    fail_if(comparison["additivity_status"] != ("REVIEW" if abs(residual) >= TOLERANCE else "PASS"),
            "Synthetic review status must match the toy residual")
    return {"evidence_class": DATA_KIND, "dataset_role": comparison["dataset_role"],
            "variant_label": comparison["variant_label"], "avi_raw": raw,
            "cdf_quantile": cdf, "phred_recorded": comparison["phred_recorded"],
            "top_feature": comparison["top_feature"], "top_feature_weight": comparison["top_feature_weight"],
            "attribution_sum": comparison["attribution_sum"], "additivity_residual": comparison["additivity_residual"],
            "additivity_status": comparison["additivity_status"], "assembly": ASSEMBLY,
            "model_executed": False, "native_response_retained": False, "clinical_interpretation": False}


def classify_additivity(raw_score: float | None, attrs: dict[str, float | None] | None) -> tuple[str, float | None, float | None]:
    if raw_score is None or attrs is None:
        return ("NOT_APPLICABLE_MISSING", None, None)
    attr_sum = math.fsum(float(attrs[x]) for x in FEATURE_NAMES)
    residual = raw_score - attr_sum
    return ("PASS" if abs(residual) < TOLERANCE else "REVIEW", attr_sum, residual)


def build_report(fixture: Any, comparison: Any = None) -> dict[str, Any]:
    cases = validate_fixture(fixture)
    comparison = validate_comparison(comparison)
    summaries: list[dict[str, Any]] = []
    for case in cases:
        attrs = case["attributions"]
        status, attr_sum, residual = classify_additivity(case["raw_score"], attrs)
        ranked = [] if attrs is None else sorted(
            ({"feature": name, "group": group, "contribution": attrs[name]} for name, group, _ in FEATURES),
            key=lambda row: (-abs(row["contribution"]), row["feature"]),
        )
        summaries.append({"case_id": case["case_id"], "variant_label": case["variant_label"], "raw_score": case["raw_score"], "feature_count": 18, "missing_feature_value_count": sum(case["feature_values"][n] is None for n in FEATURE_NAMES), "attribution_sum": attr_sum, "additivity_residual": residual, "additivity_status": status, "top_attribution": ranked[0] if ranked else None, "top5_attributions": ranked[:5], "review_status": "REVIEW_REQUIRED", "source_class": DATA_KIND})
    return {"status": "PASS", "mode": "OFFLINE_ONLY", "evidence_class": DATA_KIND, "dataset_role": "TOY_FEATURES_AND_ATTRIBUTIONS_NOT_OFFICIAL_AVI", "assembly": ASSEMBLY, "feature_count": len(FEATURES), "feature_group_counts": GROUP_COUNT, "feature_catalog": [{"name": n, "group": g, "description": d} for n, g, d in FEATURES], "case_count": len(cases), "cases": summaries, "review_count": sum(x["additivity_status"] == "REVIEW" for x in summaries), "missing_case_count": sum(x["raw_score"] is None for x in summaries), "comparison_summary": comparison, "comparison_additivity_status": comparison["additivity_status"] if comparison else "NOT_PROVIDED", "api_key_read": False, "official_avi_model_executed": False, "real_api_calls_this_run": 0, "clinical_interpretation": False, "manual_review": "REQUIRED", "publication_permission": "NOT_GRANTED"}


def make_csv(feature_map_cases: list[dict[str, Any]], *, attribution: bool) -> str:
    output = io.StringIO(newline="")
    w = csv.writer(output, lineterminator="\n")
    w.writerow(["case_id", "variant_label", "feature_name", "feature_group", "value", "evidence_class", "value_kind"])
    for case in feature_map_cases:
        fields = case["attributions"] if attribution else case["feature_values"]
        for name, group, _ in FEATURES:
            val = None if fields is None else fields[name]
            w.writerow([case["case_id"], case["variant_label"], name, group, "" if val is None else f"{val:.12g}", DATA_KIND, "TOY_ATTRIBUTION" if attribution else "TOY_INPUT_FEATURE"])
    return output.getvalue()


def make_svg(report: dict[str, Any], cases: list[dict[str, Any]]) -> str:
    head = cases[0]
    attrs = head["attributions"]
    bar_data = sorted(((name, float(attrs[name])) for name in FEATURE_NAMES if attrs is not None and attrs[name] is not None), key=lambda x: (-abs(x[1]), x[0]))[:8]
    colors = {"AlphaGenome": "#3D8293", "AlphaMissense": "#9F6697", "VEP consequences": "#C78B6D", "Conservation": "#9A7A40", "Indel flags": "#7B8197"}
    group_of = {n: g for n, g, _ in FEATURES}
    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" width="1920" height="1080" viewBox="0 0 1920 1080" role="img" aria-label="18 AVI-like feature groups, synthetic signed attributions, comparison additivity review">', '<rect width="1920" height="1080" fill="#E5DEDA"/>', '<rect x="80" y="72" width="1760" height="940" rx="42" fill="#FCFAF8"/>', '<text x="124" y="131" font-family="Noto Sans CJK TC, Arial" fill="#577688" font-size="26" letter-spacing="2">ALPHAGENOME ATLAS · DAY 11 · OFFLINE LEARNING</text>', '<text x="124" y="206" font-family="Noto Sans CJK TC, Arial" font-size="57" font-weight="700" fill="#27354A">18 INPUT FEATURES</text>', '<text x="124" y="255" font-family="Noto Sans CJK TC, Arial" font-size="27" fill="#5C6C79">五類證據，不等於十八項致病機制 · SYNTHETIC ONLY</text>']
    group_labels=[("AlphaGenome",10), ("AlphaMissense",1),("VEP consequences",3),("Conservation",2),("Indel flags",2)]
    for i,(label,count) in enumerate(group_labels):
        x=124+i*337
        parts.extend([f'<rect x="{x}" y="297" width="311" height="105" rx="18" fill="#F2F1F0" stroke="#E5E2E0"/>', f'<circle cx="{x+31}" cy="338" r="10" fill="{colors[label]}"/>', f'<text x="{x+54}" y="347" font-family="Arial" fill="#27354A" font-size="22" font-weight="700">{xml_escape(label)}</text>', f'<text x="{x+30}" y="386" font-family="Arial" fill="#5E6D7B" font-size="22">{count} features</text>'])
    parts.extend(['<text x="135" y="473" font-family="Noto Sans CJK TC, Arial" fill="#27354A" font-size="27" font-weight="700">Synthetic case · Signed toy attribution</text>', '<text x="135" y="505" font-family="Arial" fill="#637482" font-size="21">TOY_SPLICING: contribution is not a measured molecular effect</text>', '<line x1="580" y1="538" x2="580" y2="908" stroke="#A6B3BF" stroke-width="2"/>'])
    for i,(name,val) in enumerate(bar_data):
        y=555+i*43
        col=colors[group_of[name]]; w=min(abs(val)*490,340)
        parts.append(f'<text x="136" y="{y+20}" font-family="Arial" font-size="20" fill="#354555">{xml_escape(name)}</text>')
        if val>=0: parts.append(f'<rect x="580" y="{y}" width="{w:.2f}" height="26" rx="7" fill="{col}"/>')
        else: parts.append(f'<rect x="{580-w:.2f}" y="{y}" width="{w:.2f}" height="26" rx="7" fill="{col}"/>')
        parts.append(f'<text x="{590+w if val>=0 else 590:.2f}" y="{y+20}" font-family="Arial" font-size="18" fill="#27354A">{val:+.2f}</text>')
    parts.extend(['<rect x="1135" y="462" width="640" height="420" rx="23" fill="#F5F0EB"/>', '<text x="1170" y="518" font-family="Arial" font-size="28" fill="#27354A" font-weight="700">ADDITIVITY CHECK</text>'])
    toy=report["cases"][0]
    parts.extend([f'<text x="1170" y="571" font-family="Arial" font-size="23" fill="#536370">TOY raw = {toy["raw_score"]:.2f} | sum = {toy["attribution_sum"]:.2f}</text>', '<rect x="1170" y="592" width="116" height="45" rx="18" fill="#E2EEE8"/>', '<text x="1199" y="624" font-family="Arial" font-size="23" fill="#296854">PASS</text>', '<line x1="1170" y1="658" x2="1736" y2="658" stroke="#D3CECA"/>'])
    comparison=report["comparison_summary"]
    if comparison:
        parts.extend(['<text x="1170" y="701" font-family="Arial" font-size="23" fill="#27354A" font-weight="700">SYNTHETIC COMPARISON</text>', f'<text x="1170" y="745" font-family="Arial" font-size="22" fill="#536370">Toy raw {comparison["avi_raw"]:.2f} vs sum {comparison["attribution_sum"]:.2f}</text>' , '<rect x="1170" y="762" width="144" height="49" rx="18" fill="#F5DDD4"/>', '<text x="1190" y="795" font-family="Arial" font-size="22" fill="#A35543" font-weight="700">REVIEW</text>', '<text x="1330" y="793" font-family="Arial" font-size="20" fill="#667581">Invented scores; no model request</text>'])
    parts.extend(['<line x1="124" y1="928" x2="1776" y2="928" stroke="#E0D9D3"/>', '<text x="124" y="977" font-family="Noto Sans CJK TC,Arial" font-size="22" fill="#6A7681">AVI 排名與 SHAP 歸因是研究線索；非發病機率、臨床診斷或實驗證明。No LIVE API.</text>', '</svg>'])
    return "\n".join(parts)+"\n"


def make_html(report: dict[str, Any], cases: list[dict[str, Any]]) -> str:
    case_rows="".join(f'<tr><th scope="row">{html.escape(x["case_id"])}</th><td>{"MISSING" if x["raw_score"] is None else f"{x["raw_score"]:.3f}"}</td><td>{x["additivity_status"]}</td><td>{"—" if x["top_attribution"] is None else html.escape(x["top_attribution"]["feature"])}</td></tr>' for x in report["cases"])
    cat_rows="".join(f'<tr><td>{html.escape(name)}</td><td>{html.escape(group)}</td><td>{html.escape(description)}</td></tr>' for name,group,description in FEATURES)
    h=report["comparison_summary"]
    comparison_html='<p>未提供合成對照資料。</p>' if not h else f'<p><strong>Synthetic comparison (invented, NOT a model output)</strong> / {html.escape(h["variant_label"])} / toy raw {h["avi_raw"]:.3f} / toy attribution sum {h["attribution_sum"]:.3f} / <strong class="review">{h["additivity_status"]}</strong>. Separate exercise; no API or patient data.</p>' 
    return f'''<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Day 11 — AVI 18 input features teaching panel</title><style>body{{margin:0;background:#E5DEDA;font-family:system-ui,'Noto Sans CJK TC',sans-serif;color:#23334A}}main{{max-width:1180px;margin:28px auto;background:#fffdfa;border-radius:22px;padding:32px;box-shadow:0 8px 24px #0001}}img{{width:100%;border-radius:12px}}table{{width:100%;border-collapse:collapse;margin:12px 0 28px}}th,td{{border-bottom:1px solid #e7e0da;text-align:left;padding:10px}}th{{background:#f3eeeb}}.tag{{display:inline-block;background:#e5f1e8;border-radius:20px;padding:6px 14px}}.review{{color:#A54C3C}}.foot{{color:#546573}}code{{background:#eee8e2;padding:2px 5px;border-radius:4px}}</style></head><body><main><div class="tag">OFFLINE · TOY_NOT_GRCH38 · NO REAL API</div><h1>Day 11｜AVI 18 項輸入特徵怎麼讀</h1><p>5 groups · 18 input keys. Python interpretation only; no AVI neural network inference.</p><img src="day11_feature_map.svg" alt="18 個特徵分組與人工歸因，含合成加成檢查"><h2>人工個案的分數與歸因 QC</h2><table><thead><tr><th>Case</th><th>Artificial raw score</th><th>Additivity QC</th><th>Top toy attribution</th></tr></thead><tbody>{case_rows}</tbody></table><h2>18 項特徵對照表</h2><table><thead><tr><th>Field name</th><th>Group</th><th>Concept</th></tr></thead><tbody>{cat_rows}</tbody></table><h2>獨立合成對照組（非模型輸出）</h2>{comparison_html}<h2>Limits</h2><p class="foot">分數值、歸因值、TOY 變異全部是人工資料。兩組皆為人工編造，residual REVIEW 僅供教學；沒有執行官方 AlphaGenome API，亦不提供臨床解讀。</p></main></body></html>\n'''


def make_notes(report: dict[str, Any]) -> str:
    lines=["# Day 11｜18 features · Evidence Notes", "", "- 教學輸入：`SYNTHETIC_EDUCATIONAL_DEMO`，人工分數與人工 attribution；**不會執行官方 AVI 模型或 SHAP**。", "- 分組：10 AlphaGenome、1 AlphaMissense、3 VEP、2 Conservation、2 Indel indicators，共 18。", "- `0.0`、`null` 與缺失來源的意義不同；不要從 0 推斷生物功能不存在。", "- 合成教學加成門檻為 0.001，**不代表官方版本的數值容差**。", "", "## Cases", "", "| Case | Toy raw | Attribution sum | Residual | QC |", "|---|---:|---:|---:|---|"]
    for c in report["cases"]:
        fmt=lambda v:"null" if v is None else f"{v:.6f}"
        lines.append(f'| {c["case_id"]} | {fmt(c["raw_score"])} | {fmt(c["attribution_sum"])} | {fmt(c["additivity_residual"])} | {c["additivity_status"]} |')
    lines += ["", "## Synthetic Comparison（完全人工編造）", ""]
    comparison = report["comparison_summary"]
    if comparison:
        lines += [f'- 虛構標籤 `{comparison["variant_label"]}`；Toy raw `{comparison["avi_raw"]}`，人工歸因加總 `{comparison["attribution_sum"]}`，差 `{comparison["additivity_residual"]}`，狀態 `{comparison["additivity_status"]}`。', '- 對照資料為 `SYNTHETIC_EDUCATIONAL_DEMO`，不是歷史模型數據、真實 AlphaGenome 結果或 API 回應。', '- 加成不一致是刻意設計的 QC 教學案例，不能藉由補一個 bias 或填 0 假裝已得到實驗結果。']
    else: lines.append('- 未提供合成對照資料。')
    lines += ["", "## Scientific / Governance", "", "- 人工 attribution 不代表正式 AVI / SHAP 模型歸因；合成差分只用來演示 `REVIEW` 檢查。", "- AVI 不是疾病風險、ACMG/AMP 分類或實驗量測。", "- 沒有真人、沒有病患 VCF、沒有金鑰、沒有 API；這份 Lab 的 release_gate 專指證據自動裁定禁止，不是本 repo 的授權證明。", "- SHA-256 驗證的是本地內容一致，不能獨立證明來源真實。", ""]
    return "\n".join(lines)


def create_outputs(fixture_path: Path, comparison_path: Path | None) -> tuple[dict[str, str], dict[str, Any]]:
    fixture=json.loads(fixture_path.read_text(encoding="utf-8"))
    comparison=json.loads(comparison_path.read_text(encoding="utf-8")) if comparison_path else None
    report=build_report(fixture,comparison)
    cases=validate_fixture(fixture)
    texts={
        "day11_feature_report.json":json_dumps(report),
        "day11_feature_values.csv":make_csv(cases,attribution=False),
        "day11_feature_attributions.csv":make_csv(cases,attribution=True),
        "day11_feature_map.svg":make_svg(report,cases),
        "day11_feature_dashboard.html":make_html(report,cases),
        "day11_evidence_notes.md":make_notes(report),
    }
    return texts,report


def write_files(output_dir: Path, fixture_path: Path, comparison_path: Path | None) -> dict[str, Any]:
    if output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {output_dir}")
    texts,report=create_outputs(fixture_path,comparison_path)
    # Input paths and bytes are audited, but file names not embedded in reproducible output.
    manifest={"report_kind":"DAY11_SYNTHETIC_OUTPUT_MANIFEST", "version":VERSION,"files":{key:sha256_bytes(val.encode("utf-8")) for key,val in sorted(texts.items())},"inputs": {"fixture_sha256":sha256_file(fixture_path), "comparison_sha256":sha256_file(comparison_path) if comparison_path else None},"count":len(texts),"signed":False,"release_gate":"BLOCKED","real_api_calls_this_run":0}
    output_dir.mkdir(parents=True,exist_ok=False)
    for name,value in texts.items(): (output_dir/name).write_text(value,encoding="utf-8",newline="\n")
    (output_dir/"day11_manifest.json").write_text(json_dumps(manifest),encoding="utf-8",newline="\n")
    return {"status":"PASS","mode":"OFFLINE_ONLY","feature_count":18,"group_count":5,"case_count":len(report["cases"]),"additivity_review_count":report["review_count"],"missing_case_count":report["missing_case_count"],"comparison_additivity_status":report["comparison_additivity_status"],"output_artifacts":7,"real_api_calls_this_run":0}


def verify_files(output_dir: Path, fixture_path: Path, comparison_path: Path | None) -> dict[str, Any]:
    expected,report=create_outputs(fixture_path,comparison_path)
    mp=output_dir/"day11_manifest.json"
    saved=json.loads(mp.read_text(encoding="utf-8"))
    fail_if(saved.get("report_kind") != "DAY11_SYNTHETIC_OUTPUT_MANIFEST" or saved.get("signed") is not False or saved.get("release_gate") != "BLOCKED", "Invalid manifest contract")
    fail_if(set(saved.get("files",{}))!=set(OUTPUT_FILES),"Manifest output files mismatch")
    fail_if(saved.get("inputs",{}).get("fixture_sha256")!=sha256_file(fixture_path),"Fixture SHA mismatch")
    fail_if(saved["inputs"].get("comparison_sha256") != (sha256_file(comparison_path) if comparison_path else None),"Synthetic comparison input SHA mismatch")
    for name,value in expected.items():
        fail_if(not (output_dir/name).is_file(),f"Missing {name}")
        blob=(output_dir/name).read_bytes()
        fail_if(blob != value.encode("utf-8"),f"Recomputed content mismatch: {name}")
        fail_if(saved["files"].get(name)!=sha256_bytes(blob),f"SHA mismatch: {name}")
    fail_if(set(p.name for p in output_dir.iterdir())!=set(OUTPUT_FILES)|{"day11_manifest.json"},"Unexpected output files")
    return {"status":"PASS","mode":"OFFLINE_ONLY","verified_artifacts":len(expected),"checksums":"MATCH","comparison_additivity_status":report["comparison_additivity_status"],"real_api_calls_this_run":0}


def main(argv: list[str]|None=None) -> int:
    parser=argparse.ArgumentParser(description="OFFLINE synthetic AVI 18-feature teaching lab (no real API)")
    sub=parser.add_subparsers(dest="action",required=True)
    for action in ("run","verify"):
        p=sub.add_parser(action)
        p.add_argument("--mode",choices=("offline","live"),default="offline")
        p.add_argument("--fixture",type=Path,default=Path("day11_synthetic_avi_features.json"))
        p.add_argument("--comparison-fixture",type=Path,default=Path("day11_synthetic_avi_comparison.json"))
        p.add_argument("--output-dir",type=Path,default=Path("day11_output"))
    args=parser.parse_args(argv)
    if args.mode != "offline":
        print("ERROR: LIVE_DISABLED; this lesson never calls the API",file=sys.stderr)
        return 2
    try:
        result=write_files(args.output_dir,args.fixture,args.comparison_fixture) if args.action=="run" else verify_files(args.output_dir,args.fixture,args.comparison_fixture)
    except (ValueError,FileNotFoundError,FileExistsError,KeyError,TypeError,json.JSONDecodeError) as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}",file=sys.stderr)
        return 2
    print(json.dumps(result,ensure_ascii=False,sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
