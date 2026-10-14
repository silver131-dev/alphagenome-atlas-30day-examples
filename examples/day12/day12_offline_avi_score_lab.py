"""Day 12 | Local-only AVI Phred ranking and attribution interpretation.

Requires Python 3.12+ standard library. The synthetic toy scores are NOT
AlphaGenome predictions; the separate comparison record is also entirely
invented (NOT historical data or a native API response). Never diagnose.
No API calls, no network clients, no patient inputs, no live mode.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import re
import sys
from pathlib import Path
from xml.sax.saxutils import escape as xml_escape

TOY_ID = re.compile(r"TOY_[A-Z0-9]{1,8}\Z", re.ASCII)
TOY_FEATURE = re.compile(r"TEACHING_[A-Z_]{2,40}\Z", re.ASCII)
ARTIFACTS = (
    "day12_avi_interpretation_report.json",
    "day12_avi_candidates.csv",
    "day12_avi_education.svg",
    "day12_avi_dashboard.html",
    "day12_evidence_notes.md",
)


def sha256_file(path: Path) -> str:
    """Hash bytes without modifying their file."""
    h = hashlib.sha256()
    with path.open("rb") as file:
        for part in iter(lambda: file.read(65536), b""):
            h.update(part)
    return h.hexdigest()


def _unique_pairs(pairs):
    values = {}
    for key, value in pairs:
        if key in values:
            raise ValueError(f"Duplicate JSON key: {key}")
        values[key] = value
    return values


def load_json(path: Path):
    with path.open("r", encoding="utf-8") as file:
        return json.load(
            file, object_pairs_hook=_unique_pairs,
            parse_constant=lambda item: (_ for _ in ()).throw(ValueError("Nonfinite JSON value")),
        )


def save_json(path: Path, data) -> None:
    path.write_text(json.dumps(data, sort_keys=True, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _finite_number(value, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (float, int)):
        raise ValueError(f"{name} must be finite numeric")
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return float(value)


def phred_from_cdf(cdf: float) -> float:
    """Convert a reference-distribution CDF to a rank scale, not disease risk.

    P = -10 log10(1-CDF). CDF=1 is intentionally rejected; an infinite
    Phred is not represented as a finite result in this teaching lab.
    """
    cdf = _finite_number(cdf, "cdf_quantile")
    if not (0.0 <= cdf < 1.0):
        raise ValueError("cdf_quantile requires 0 <= cdf < 1")
    return -10.0 * math.log10(1.0 - cdf)


def top_percent_from_cdf(cdf: float) -> float:
    cdf = _finite_number(cdf, "cdf_quantile")
    if not (0.0 <= cdf < 1.0):
        raise ValueError("cdf_quantile requires 0 <= cdf < 1")
    return 100.0 * (1.0 - cdf)


def validate_synthetic(data: dict) -> list[dict]:
    if not isinstance(data, dict):
        raise ValueError("Synthetic fixture must be object")
    if data.get("evidence_class") != "SYNTHETIC_EDUCATIONAL_DEMO":
        raise ValueError("Only explicit synthetic teaching data accepted")
    if data.get("assembly") != "TOY_NOT_GRCH38":
        raise ValueError("Synthetic assembly must be TOY_NOT_GRCH38")
    if data.get("score_semantics") != "PHRED_RANK_TEACHING_ONLY":
        raise ValueError("Invalid teaching score semantics")
    if data.get("clinical_interpretation") is not False:
        raise ValueError("Clinical interpretation must be false")
    candidates = data.get("candidates")
    if not isinstance(candidates, list) or not (1 <= len(candidates) <= 32):
        raise ValueError("1-32 synthetic candidates are required")
    seen = set()
    for item in candidates:
        if not isinstance(item, dict):
            raise ValueError("Candidate must be object")
        name = item.get("variant_id")
        if not isinstance(name, str) or not TOY_ID.fullmatch(name):
            raise ValueError("Only TOY_* artificial identifiers allowed")
        if name in seen:
            raise ValueError("Duplicate candidate")
        seen.add(name)
        cdf = item.get("cdf_quantile")
        attrs = item.get("toy_attributions")
        if not isinstance(attrs, list) or len(attrs) > 12:
            raise ValueError("Invalid teaching attributions")
        if cdf is None and attrs:
            raise ValueError("Missing score must not have pretend attribution")
        if cdf is not None:
            phred_from_cdf(cdf)
            if not attrs:
                raise ValueError("Scored example requires attributions")
        attr_names = set()
        for attr in attrs:
            if not isinstance(attr, dict):
                raise ValueError("Attribution must be an object")
            label = attr.get("feature")
            if not isinstance(label, str) or not TOY_FEATURE.fullmatch(label):
                raise ValueError("Only TEACHING_* feature names permitted")
            if label in attr_names:
                raise ValueError("Duplicate teaching feature")
            attr_names.add(label)
            magnitude = _finite_number(attr.get("value"), "toy_attribution")
            if not (-1 <= magnitude <= 1):
                raise ValueError("Teaching attribution must be between -1 and 1")
    return candidates


def validate_comparison(data: dict) -> dict:
    """Validate a fully invented comparison record, never an API/historical response."""
    required = {"evidence_class", "data_kind", "dataset_role", "assembly", "variant_label",
                "score_semantics", "model_executed", "native_response_retained", "clinical_interpretation",
                "avi_raw", "cdf_quantile", "phred_recorded", "feature_count", "top_feature",
                "top_feature_weight", "top_feature_abs_share", "attribution_sum", "additivity_residual", "additivity_status"}
    if not isinstance(data, dict) or set(data) != required:
        raise ValueError("Only the documented synthetic comparison fields are accepted")
    if any(data[k] != 'SYNTHETIC_EDUCATIONAL_DEMO' for k in ('evidence_class','data_kind')):
        raise ValueError("This comparison must be explicitly synthetic")
    if data['dataset_role'] != 'ARTIFICIAL_COMPARISON_NOT_MODEL_OUTPUT' or data['assembly'] != 'TOY_NOT_GRCH38' or data['variant_label'] != 'TOY_COMPARISON' or data['score_semantics'] != 'TOY_AVI_LIKE_SCORE_AND_RANK':
        raise ValueError("No real variants, assemblies or model scores are accepted")
    if any(data[k] is not False for k in ('model_executed','native_response_retained','clinical_interpretation')):
        raise ValueError("Historical, live, native and clinical claims are forbidden")
    if type(data['feature_count']) is not int or data['feature_count'] != 18 or data['top_feature'] != 'MERGED_SPLICING':
        raise ValueError("Toy feature metadata mismatch")
    raw = _finite_number(data['avi_raw'], 'toy_raw')
    cdf = _finite_number(data['cdf_quantile'], 'toy_cdf')
    phred = _finite_number(data['phred_recorded'], 'toy_phred')
    total = _finite_number(data['attribution_sum'], 'toy_attribution_sum')
    residual = _finite_number(data['additivity_residual'], 'toy_residual')
    weight = _finite_number(data['top_feature_weight'], 'toy_top_weight')
    share = _finite_number(data['top_feature_abs_share'], 'toy_top_share')
    if not (0 <= weight <= 1 and 0 <= share <= 1):
        raise ValueError("Toy feature weights must be between 0 and 1")
    if abs(phred_from_cdf(cdf) - phred) > 1e-7:
        raise ValueError("Toy PHRED/CDF mismatch")
    if abs(raw - total - residual) > 1e-9:
        raise ValueError("Toy additivity residual mismatch")
    if data['additivity_status'] != ('REVIEW' if abs(residual) >= 1e-3 else 'PASS'):
        raise ValueError("Toy additivity status must reflect the difference")
    return data


def build_report(toy_data: dict, comparison_data: dict, toy_sha: str, comparison_sha: str) -> dict:
    candidates = validate_synthetic(toy_data)
    comparison = validate_comparison(comparison_data)
    rows = []
    for candidate in candidates:
        cdf = candidate["cdf_quantile"]
        attrs = candidate["toy_attributions"]
        top = max(attrs, key=lambda x: abs(x["value"])) if attrs else None
        rows.append({
            "variant_id": candidate["variant_id"],
            "cdf_quantile": cdf,
            "phred": round(phred_from_cdf(cdf), 6) if cdf is not None else None,
            "top_percent": round(top_percent_from_cdf(cdf), 6) if cdf is not None else None,
            "score_status": "SYNTHETIC_RANK" if cdf is not None else "MISSING_SCORE",
            "top_teaching_feature": top["feature"] if top else None,
            "top_teaching_weight": top["value"] if top else None,
            "toy_attributions": attrs,
            "evidence_class": "SYNTHETIC_FIXTURE",
            "clinical_interpretation": False,
        })
    return {
        "report_kind": "DAY12_AVI_RANK_INTERPRETATION_OFFLINE",
        "evidence_class": "SYNTHETIC_EDUCATIONAL_DEMO",
        "mode": "OFFLINE_ONLY",
        "assembly": "TOY_NOT_GRCH38",
        "candidate_count": len(rows),
        "scored_count": sum(x["phred"] is not None for x in rows),
        "missing_score_count": sum(x["phred"] is None for x in rows),
        "reference_phred_examples": [
            {"phred": score, "top_percent": round(100 * 10 ** (-score / 10), 6)}
            for score in (10, 20, 30)
        ],
        "toy_candidates": rows,
        "synthetic_comparison_reference": {
            "evidence_class": "SYNTHETIC_EDUCATIONAL_DEMO",
            "dataset_role": comparison["dataset_role"],
            "variant_label": comparison["variant_label"],
            "assembly": "TOY_NOT_GRCH38",
            "toy_raw": comparison["avi_raw"],
            "cdf_quantile": comparison["cdf_quantile"],
            "phred": round(phred_from_cdf(comparison["cdf_quantile"]), 6),
            "top_percent": round(top_percent_from_cdf(comparison["cdf_quantile"]), 6),
            "top_feature": comparison["top_feature"],
            "top_feature_weight": comparison["top_feature_weight"],
            "top_feature_abs_share": comparison["top_feature_abs_share"],
            "feature_count": comparison["feature_count"],
            "attribution_sum": comparison["attribution_sum"],
            "toy_raw_minus_attribution_sum": comparison["additivity_residual"],
            "additivity_status": comparison["additivity_status"],
            "model_executed": False,
            "clinical_interpretation": False,
        },
        "input_provenance": {
            "toy_fixture_sha256": toy_sha,
            "synthetic_comparison_sha256": comparison_sha,
            "comparison_is_synthetic": True,
            "real_api_calls_this_run": 0,
        },
        "missing_is_not_zero": True,
        "toy_attributions_are_not_official_shap": True,
        "clinical_interpretation": False,
        "model_inference_performed": False,
        "real_api_calls_this_run": 0,
        "manual_review": "REQUIRED",
        "release_gate": "BLOCKED",
    }


def _to_svg(report: dict) -> str:
    """Accessible, deterministic 16:9 teaching visualization (no external fonts)."""
    bg, ink, sub = "#E5DEDA", "#292D35", "#62616A"
    colors = ["#B27F78", "#739EAF", "#9E90B6", "#C9C0BE"]
    rows = report["toy_candidates"]
    def label(x):
        return xml_escape(str(x))
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="1920" height="1080" viewBox="0 0 1920 1080" role="img" aria-labelledby="st sd">',
           '<title id="st">Day 12 | AVI Phred is a ranking, not disease probability</title>',
           '<desc id="sd">All candidate scores and the separate comparison are invented teaching values; no model inference or clinical prediction.</desc>',
           f'<rect width="1920" height="1080" fill="{bg}"/>',
           '<defs><linearGradient id="gl" x1="0" x2="1"><stop stop-color="#FFF9F6"/><stop offset="1" stop-color="#F1E6E1"/></linearGradient></defs>',
           '<rect x="68" y="62" width="1784" height="956" rx="40" fill="url(#gl)" stroke="#D2C2BC"/>',
           f'<text x="132" y="151" font-size="58" font-weight="700" fill="{ink}" font-family="Arial,sans-serif">Day 12  |  AVI PHRED</text>',
           f'<text x="134" y="209" font-size="29" fill="{sub}" font-family="Arial,sans-serif">Predicted-impact ranking  |  NOT a disease probability</text>']
    # Three PHRED reference examples
    for i, (score, pct) in enumerate(((10, "Top 10%"), (20, "Top 1%"), (30, "Top 0.1%"))):
        x = 135 + i * 566
        out.append(f'<rect x="{x}" y="252" width="522" height="148" rx="22" fill="#FFFFFF" stroke="#DDCECA"/>')
        out.append(f'<text x="{x+27}" y="319" font-size="57" font-weight="700" fill="{colors[i]}" font-family="Arial">{score}</text>')
        out.append(f'<text x="{x+163}" y="310" font-size="27" fill="{ink}" font-family="Arial">{pct}</text>')
        out.append(f'<text x="{x+163}" y="344" font-size="19" fill="{sub}" font-family="Arial">Reference rank tail</text>')
    out.append(f'<text x="139" y="483" font-size="31" font-weight="700" fill="{ink}" font-family="Arial">SYNTHETIC: same rank, different proxy features</text>')
    start_x = 135
    bar_max = 820
    for i, row in enumerate(rows):
        y = 535 + i * 99
        out.append(f'<text x="{start_x}" y="{y+17}" font-size="27" font-weight="700" fill="{ink}" font-family="Arial">{label(row["variant_id"])}</text>')
        if row["phred"] is not None:
            width = row["phred"] / 30 * bar_max
            out.extend([f'<rect x="330" y="{y-15}" width="{bar_max}" height="37" rx="12" fill="#E9E8EA"/>',
                        f'<rect x="330" y="{y-15}" width="{width:.1f}" height="37" rx="12" fill="{colors[i]}"/>',
                        f'<text x="1180" y="{y+16}" font-size="26" font-weight="700" fill="{ink}" font-family="Arial">{row["phred"]:.1f}</text>',
                        f'<text x="1310" y="{y+16}" font-size="20" fill="{sub}" font-family="Arial">Top {row["top_percent"]:.3f}%</text>'])
        else:
            out.append(f'<text x="330" y="{y+15}" font-size="27" fill="{sub}" font-family="Arial">MISSING_SCORE (not zero)</text>')
    h = report["synthetic_comparison_reference"]
    out.extend(['<rect x="132" y="933" width="1652" height="62" rx="16" fill="#EAE3DE"/>',
                f'<text x="158" y="972" font-size="21" font-weight="700" fill="{ink}" font-family="Arial">SYNTHETIC COMPARISON  |  Toy PHRED {h["phred"]:.3f}  |  Top {h["top_percent"]:.3f}%  |  Toy additivity: REVIEW</text>' ,
                '<text x="134" y="899" font-size="21" fill="#754D49" font-family="Arial">All values are invented; toy attribution differences do not establish mechanisms or model behavior.</text>',
                '</svg>'])
    return "\n".join(out) + "\n"


def _to_csv(report: dict, path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as file:
        w = csv.writer(file, lineterminator="\n")
        w.writerow(["toy_variant", "phred", "top_percent", "score_status", "top_teaching_feature", "top_teaching_weight", "evidence_class", "clinical_interpretation"])
        for row in report["toy_candidates"]:
            w.writerow([row["variant_id"], "" if row["phred"] is None else row["phred"],
                        "" if row["top_percent"] is None else row["top_percent"],
                        row["score_status"], row["top_teaching_feature"] or "",
                        "" if row["top_teaching_weight"] is None else row["top_teaching_weight"],
                        row["evidence_class"], "false"])


def _to_html(report: dict, svg: str) -> str:
    def h(value):
        return html.escape(str(value), quote=True)
    body = []
    for row in report["toy_candidates"]:
        score = "MISSING (JSON null)" if row["phred"] is None else f'{row["phred"]:.3f}'
        pct = "not available" if row["top_percent"] is None else f'{row["top_percent"]:.3f}%'
        body.append("<tr>" + "".join(f"<td>{h(item)}</td>" for item in [row["variant_id"],score,pct,row["score_status"],row["top_teaching_feature"] or "none"]) + "</tr>")
    comparison = report["synthetic_comparison_reference"]
    return ("<!doctype html>\n<html lang=\"zh-Hant\"><head><meta charset=\"utf-8\">"
            "<title>Day 12 | AVI PHRED Ranking — Offline</title>"
            "<style>body{font:19px/1.6 system-ui,sans-serif;background:#E5DEDA;color:#292D35;max-width:1200px;margin:auto;padding:2rem}"
            "svg{width:100%;height:auto;border-radius:18px}table{border-collapse:collapse;width:100%;background:#fff}td,th{padding:.7rem;border-bottom:1px solid #ddd;text-align:left}"
            ".notice{padding:1rem;background:#fff4eb;border-left:5px solid #B27F78;margin:1.5rem 0}code{background:#eee;padding:.1rem .3rem}</style></head><body>"
            "<h1>Day 12：AVI 是研究排序，不是發病機率</h1>"
            "<div class=\"notice\">SYNTHETIC EDUCATIONAL DEMO / OFFLINE ONLY. No patient data, no model inference, no real API calls.</div>"
            + svg + "<h2>人工候選：同分不同來源</h2><table><thead><tr><th>候選</th><th>PHRED</th><th>Top percentile</th><th>狀態</th><th>人工特徵</th></tr></thead><tbody>"
            + "".join(body) + "</tbody></table>"
            f'<h2>分開呈現：合成對照範例（非 API）</h2><p>{h(comparison["variant_label"])}；Toy PHRED {comparison["phred"]:.4f}；Top {comparison["top_percent"]:.4f}%；人工主要歸因欄位 {h(comparison["top_feature"])}。</p>' 
            '<div class="notice">additivity_status=REVIEW；刻意設定的人工歸因加總不一致，僅用於驗證 QC；不是任何 AlphaGenome 模型或臨床結果。</div>'
            "<p>Missing is not zero. Higher Phred indicates more extreme predicted-impact ranking relative to its reference distribution; it does not measure a person's chance of disease.</p>"
            "<p>All toy variant identifiers and scores are fictional. This exercise never calls a model.</p></body></html>\n")


def _to_notes(report: dict) -> str:
    comparison = report["synthetic_comparison_reference"]
    return f"""# Day12｜AVI 教學證據與限制

- `SYNTHETIC_EDUCATIONAL_DEMO`：所有候選與獨立對照組完全人工編造，不是 AlphaGenome 輸出。
- 對照標籤 `{comparison['variant_label']}`：Toy raw {comparison['toy_raw']:.3f}；Toy CDF {comparison['cdf_quantile']:.3f}；Toy PHRED {comparison['phred']:.6f}；Top {comparison['top_percent']:.3f}% 。
- 人工歸因欄位 `{comparison['top_feature']}`：其相對大小不是疾病／剪接發生機率。
- `additivity_status={comparison['additivity_status']}`；`toy_raw - attribution_sum = {comparison['toy_raw_minus_attribution_sum']:.3f}`；故意的不一致提供 QC 練習，並非模型誤差測量。
- `TOY_NOT_GRCH38`：沒有使用真實基因座或人類基因組組裝參照。
- PHRED 教學轉換非風險、致病性或基因功能驗證；分數缺值保留 JSON null。
- SHA-256 僅驗證本地合成檔案內容一致，不能證明模型或臨床有效性。
- `OFFLINE_ONLY`；`real_api_calls_this_run=0`；教學 `release_gate=BLOCKED` 禁止自動科研／臨床確證。
"""


def run(toy_path: Path, comparison_path: Path, output_dir: Path, mode: str = "offline") -> dict:
    if mode != "offline":
        raise ValueError("Only offline mode is implemented; live API requires separate authorization")
    if output_dir.exists():
        raise FileExistsError("Output directory exists; choose a new --output-dir")
    report = build_report(load_json(toy_path), load_json(comparison_path), sha256_file(toy_path), sha256_file(comparison_path))
    output_dir.mkdir(parents=True)
    save_json(output_dir / ARTIFACTS[0], report)
    _to_csv(report, output_dir / ARTIFACTS[1])
    svg = _to_svg(report)
    (output_dir / ARTIFACTS[2]).write_text(svg, encoding="utf-8")
    (output_dir / ARTIFACTS[3]).write_text(_to_html(report, svg), encoding="utf-8")
    (output_dir / ARTIFACTS[4]).write_text(_to_notes(report), encoding="utf-8")
    manifest = {
        "mode": "OFFLINE_ONLY",
        "evidence_class": "SYNTHETIC_EDUCATIONAL_DEMO",
        "source_sha256": {"toy": sha256_file(toy_path), "synthetic_comparison": sha256_file(comparison_path)},
        "artifact_sha256": {name: sha256_file(output_dir / name) for name in ARTIFACTS},
        "manifest_is_not_a_signature": True,
        "release_gate": "BLOCKED",
        "real_api_calls_this_run": 0,
    }
    save_json(output_dir / "day12_manifest.json", manifest)
    return {
        "status": "PASS", "mode": "OFFLINE_ONLY",
        "candidate_count": report["candidate_count"],
        "scored_count": report["scored_count"],
        "missing_score_count": report["missing_score_count"],
        "synthetic_comparison_phred": report["synthetic_comparison_reference"]["phred"],
        "synthetic_comparison_additivity_status": report["synthetic_comparison_reference"]["additivity_status"],
        "output_artifact_count": len(ARTIFACTS),
        "output_dir": str(output_dir), "real_api_calls_this_run": 0,
    }


def verify(toy_path: Path, comparison_path: Path, report_dir: Path) -> dict:
    manifest = load_json(report_dir / "day12_manifest.json")
    if manifest.get("mode") != "OFFLINE_ONLY" or manifest.get("release_gate") != "BLOCKED" or manifest.get("real_api_calls_this_run") != 0:
        raise ValueError("Manifest status mismatch")
    expected_sources = {"toy": sha256_file(toy_path), "synthetic_comparison": sha256_file(comparison_path)}
    if manifest.get("source_sha256") != expected_sources:
        raise ValueError("Source SHA mismatch")
    d = manifest.get("artifact_sha256")
    if not isinstance(d, dict) or set(d) != set(ARTIFACTS):
        raise ValueError("Artifact manifest missing or unexpected names")
    for name in ARTIFACTS:
        if not re.fullmatch(r"[0-9a-f]{64}", str(d[name])) or sha256_file(report_dir / name) != d[name]:
            raise ValueError(f"Output SHA mismatch: {name}")
    expected = build_report(load_json(toy_path), load_json(comparison_path), **{
        "toy_sha": expected_sources["toy"], "comparison_sha": expected_sources["synthetic_comparison"]
    })
    if load_json(report_dir / ARTIFACTS[0]) != expected:
        raise ValueError("Report content does not match source calculations")
    with (report_dir / ARTIFACTS[1]).open(newline="", encoding="utf-8") as file:
        records = list(csv.DictReader(file))
    if len(records) != expected["candidate_count"] or any(records[i]["toy_variant"] != expected["toy_candidates"][i]["variant_id"] for i in range(len(records))):
        raise ValueError("CSV row contract mismatch")
    return {"status": "PASS", "checksums": "MATCH", "source_files_verified": 2,
            "output_artifacts_verified": len(ARTIFACTS), "manual_review": "REQUIRED",
            "release_gate": "BLOCKED", "real_api_calls_this_run": 0}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Offline AVI PHRED ranking lesson; not clinical software")
    sub = parser.add_subparsers(dest="command", required=True)
    for cmd in ("run", "verify"):
        part = sub.add_parser(cmd)
        part.add_argument("--mode", default="offline", choices=["offline", "live"])
        part.add_argument("--fixture", default="day12_synthetic_avi_fixture.json")
        part.add_argument("--comparison-fixture", default="day12_synthetic_avi_comparison.json")
        part.add_argument("--output-dir" if cmd == "run" else "--report-dir", default="day12_output")
    args = parser.parse_args(argv)
    try:
        if args.mode != "offline":
            raise ValueError("LIVE mode intentionally blocked")
        result = (run(Path(args.fixture), Path(args.comparison_fixture), Path(args.output_dir)) if args.command == "run"
                  else verify(Path(args.fixture), Path(args.comparison_fixture), Path(args.report_dir)))
        print(json.dumps(result, ensure_ascii=False, sort_keys=True))
        return 0
    except (ValueError, OSError, KeyError, TypeError) as ex:
        print(json.dumps({"status": "ERROR", "error_type": type(ex).__name__, "message": str(ex), "real_api_calls_this_run": 0}, ensure_ascii=False))
        return 2


if __name__ == "__main__":
    sys.exit(main())
