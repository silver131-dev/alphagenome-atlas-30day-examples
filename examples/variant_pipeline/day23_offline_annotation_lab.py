"""Day 23 | OFFLINE synthetic VCF annotation schema and missing-value lesson.

Reuses Day21's ORIGINAL synthetic parser/QC and four byte-identical fixtures.
Does not load an AlphaGenome SDK, use a key, make a request, or interpret disease.
Python 3.12+ standard library only. NOT a clinical or production VCF engine.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from html import escape
import json
import math
from pathlib import Path
import re

from day21_offline_vcf_review_lab import (
    ASSEMBLY, QC_ORDER, VARIANT_RE, build_review_packet,
    load_toy_reference, load_toy_scores, parse_toy_bed,
)

EVIDENCE_CLASS = "SYNTHETIC_EDUCATIONAL_DEMO"
SCHEMA_VERSION = "day23.synthetic.annotation.v1"
ANNOTATION_FIELDS = (
    "variant_id", "chrom", "pos_1based", "ref", "alt", "qc_status",
    "qc_reason", "synthetic_score", "score_status", "score_source",
    "review_status", "source_line",
)
SHA256_RE = re.compile(r"[a-f0-9]{64}\Z")
EXPECTED_QC_ORDER = ("FILTER", "BED", "REF", "INDEL", "MISSING_QUAL", "LOW_QUAL")
ALLOWED_QC = {"SELECTED", "REVIEW_REQUIRED", "REJECTED"}
REASONS_BY_STATUS = {
    "SELECTED": {"TEACHING_QC_PASS"},
    "REVIEW_REQUIRED": {"INDEL_NOT_NORMALIZED", "MISSING_QUAL"},
    "REJECTED": {"FILTER_NOT_PASS", "OUTSIDE_TEACHING_BED", "REF_MISMATCH", "LOW_QUAL"},
}


def score_state(value):
    """None means absent; 0.0 is a PRESENT numeric teaching value."""
    if value is None:
        return "MISSING_SCORE"
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError("Synthetic score must be finite numeric or null")
    return "SYNTHETIC_SCORE"


def split_variant(variant_id):
    if not isinstance(variant_id, str) or VARIANT_RE.fullmatch(variant_id) is None:
        raise ValueError("Invalid artificial variant ID")
    chrom, pos, alleles = variant_id.split(":")
    ref, alt = alleles.split(">")
    if ref == alt:
        raise ValueError("REF and ALT must differ")
    return chrom, int(pos), ref, alt


def validate_day21_packet(packet):
    """Accept only the narrow, explicit Day21 synthetic contract, never live/patient data."""
    if not isinstance(packet, dict):
        raise ValueError("Day21 packet must be a JSON object")
    if (packet.get("evidence_class") != EVIDENCE_CLASS
            or packet.get("mode") != "OFFLINE_ONLY"
            or packet.get("assembly") != ASSEMBLY
            or packet.get("score_source") != "SYNTHETIC_FIXTURE_ONLY"
            or packet.get("atlas_query_performed") is not False
            or packet.get("clinical_interpretation") is not False
            or packet.get("real_api_calls_this_run") != 0):
        raise ValueError("Refusing input not explicitly marked as offline synthetic teaching data")
    if tuple(QC_ORDER) != EXPECTED_QC_ORDER or packet.get("qc_order") != list(EXPECTED_QC_ORDER):
        raise ValueError("QC precedence differs from validated Day21")
    qc_rows, review_rows = packet.get("qc_rows"), packet.get("review_packet")
    if not isinstance(qc_rows, list) or not qc_rows or not isinstance(review_rows, list):
        raise ValueError("Day21 packet requires QC and annotation lists")
    qc_by_id = {}
    for row in qc_rows:
        if not isinstance(row, dict):
            raise ValueError("Malformed QC row")
        variant_id = row.get("variant_id")
        split_variant(variant_id)
        if variant_id in qc_by_id or row.get("qc_status") not in ALLOWED_QC:
            raise ValueError("Duplicate QC variant or invalid status")
        if row.get("reason") not in REASONS_BY_STATUS[row["qc_status"]]:
            raise ValueError("QC reason is not valid for its status")
        qual = row.get("qual")
        if qual is not None and (isinstance(qual, bool) or not isinstance(qual, (int, float)) or not math.isfinite(qual) or qual < 0):
            raise ValueError("Invalid QUAL in QC row")
        if (not isinstance(row.get("filter"), str) or not row["filter"]
                or not isinstance(row.get("source_line"), int)
                or isinstance(row.get("source_line"), bool)
                or row["source_line"] < 1):
            raise ValueError("QC filter and original source line required")
        qc_by_id[variant_id] = row
    counts = {
        "SELECTED": sum(r["qc_status"] == "SELECTED" for r in qc_rows),
        "REVIEW_REQUIRED": sum(r["qc_status"] == "REVIEW_REQUIRED" for r in qc_rows),
        "REJECTED": sum(r["qc_status"] == "REJECTED" for r in qc_rows),
    }
    if any(packet.get(field) != expected for field, expected in (
        ("candidate_count", len(qc_rows)), ("selected_count", counts["SELECTED"]),
        ("review_required_qc_count", counts["REVIEW_REQUIRED"]),
        ("rejected_count", counts["REJECTED"]),
    )):
        raise ValueError("QC count contract mismatch")
    reviewed = set()
    missing = 0
    for row in review_rows:
        if not isinstance(row, dict):
            raise ValueError("Malformed annotation row")
        variant_id = row.get("variant_id")
        split_variant(variant_id)
        if variant_id in reviewed or variant_id not in qc_by_id:
            raise ValueError("Duplicate or unmatched annotation")
        if qc_by_id[variant_id]["qc_status"] != "SELECTED":
            raise ValueError("Non-selected variants must not be annotated")
        value = row.get("synthetic_score")
        status = score_state(value)
        if (row.get("score_status") != status
                or row.get("review_status") != "REVIEW_REQUIRED"
                or row.get("evidence_class") != "SYNTHETIC_FIXTURE"):
            raise ValueError("Annotation status or source is inconsistent")
        missing += value is None
        reviewed.add(variant_id)
    if (reviewed != {vid for vid, row in qc_by_id.items() if row["qc_status"] == "SELECTED"}
            or packet.get("missing_score_count") != missing):
        raise ValueError("Missing, extra, or wrongly counted annotations")
    return qc_by_id


def annotation_schema(packet, input_provenance):
    """Transform verified Day21 rows into Day23's explicit analysis-ready schema."""
    qc_by_id = validate_day21_packet(packet)
    if not isinstance(input_provenance, dict) or not input_provenance:
        raise ValueError("Input provenance is required")
    for name, record in input_provenance.items():
        if not isinstance(name, str) or not isinstance(record, dict):
            raise ValueError("Malformed provenance")
        if (not isinstance(record.get("filename"), str)
                or not record["filename"] or "/" in record["filename"] or "\\" in record["filename"]
                or not isinstance(record.get("sha256"), str)
                or SHA256_RE.fullmatch(record["sha256"]) is None):
            raise ValueError("Each provenance record needs basename and SHA-256")
    annotations = []
    for row in packet["review_packet"]:
        variant_id = row["variant_id"]
        chrom, pos, ref, alt = split_variant(variant_id)
        qc = qc_by_id[variant_id]
        annotations.append({
            "variant_id": variant_id,
            "chrom": chrom,
            "pos_1based": pos,
            "ref": ref,
            "alt": alt,
            "qc_status": "SELECTED",
            "qc_reason": qc["reason"],
            "synthetic_score": row["synthetic_score"],
            "score_status": row["score_status"],
            "score_source": "SYNTHETIC_FIXTURE_ONLY",
            "review_status": "REVIEW_REQUIRED",
            "source_line": qc["source_line"],
        })
    if len(annotations) != packet["selected_count"]:
        raise ValueError("Day23 annotation count mismatch")
    return {
        "schema_version": SCHEMA_VERSION,
        "evidence_class": EVIDENCE_CLASS,
        "mode": "OFFLINE_ONLY",
        "assembly": ASSEMBLY,
        "annotation_columns": list(ANNOTATION_FIELDS),
        "candidate_count": packet["candidate_count"],
        "selected_count": packet["selected_count"],
        "review_required_qc_count": packet["review_required_qc_count"],
        "rejected_count": packet["rejected_count"],
        "annotated_count": sum(r["score_status"] == "SYNTHETIC_SCORE" for r in annotations),
        "missing_score_count": sum(r["score_status"] == "MISSING_SCORE" for r in annotations),
        "qc_order": list(QC_ORDER),
        "qc_rows": packet["qc_rows"],
        "annotations": annotations,
        "missing_data_policy": "JSON null + MISSING_SCORE; blank CSV numeric cell + MISSING_SCORE; never substitute zero",
        "input_provenance": input_provenance,
        "pipeline_source": "REUSED_UNMODIFIED_DAY21_SYNTHETIC_QC",
        "official_avi": False,
        "atlas_query_performed": False,
        "clinical_interpretation": False,
        "real_api_calls_this_run": 0,
    }


def render_html(report):
    """Create a standalone HTML table; escape all data-derived labels."""
    def h(value):
        return escape("null" if value is None else str(value), quote=True)
    annotation_rows = "\n".join(
        "<tr>" + "".join(f"<td>{h(row[column])}</td>" for column in
                           ("variant_id", "qc_status", "synthetic_score", "score_status", "score_source", "review_status")) + "</tr>"
        for row in report["annotations"]
    )
    qc_rows = "\n".join(
        "<tr>" + "".join(f"<td>{h(row.get(column))}</td>" for column in
                           ("variant_id", "qc_status", "reason")) + "</tr>"
        for row in report["qc_rows"]
    )
    return f'''<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Day 23 | Synthetic Annotation Schema</title>
<style>:root{{font:16px/1.55 system-ui,"Noto Sans TC",sans-serif;color:#2e3649;background:#f1f1f2}}
body{{margin:auto;max-width:1140px;padding:30px 18px 65px}}h1{{font-size:1.85rem}}
.warning{{background:#fff1de;border-left:5px solid #b8743e;padding:15px;border-radius:8px}}
.stats{{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px}}
.stats div,section{{background:white;border:1px solid #dfdfdf;border-radius:12px;padding:16px;margin-top:18px}}
.num{{font-size:2rem;font-weight:800;color:#286683}}.label{{font-size:.85rem;color:#697386}}
table{{width:100%;border-collapse:collapse}}th,td{{text-align:left;padding:10px;border-bottom:1px solid #eaeaea}}
th{{background:#e9f0f5}}section{{overflow-x:auto}}.missing{{color:#a83b43}}
footer{{font-size:.85rem;color:#687386;margin-top:25px}}</style></head><body>
<h1>Day 23｜從缺失值到結構化標註欄位</h1>
<div class="warning"><strong>SYNTHETIC · OFFLINE_ONLY · TOY_NOT_GRCH38</strong><br>
所有分數均為人工教學值，不是 AlphaGenome AVI，沒有患者資料或真實 API 呼叫。所有 SELECTED 仍待研究人工審閱。</div>
<div class="stats"><div><span class="num">{report['candidate_count']}</span><br><span class="label">候選</span></div>
<div><span class="num">{report['selected_count']}</span><br><span class="label">通過教學 QC</span></div>
<div><span class="num">{report['annotated_count']}</span><br><span class="label">有人工分數</span></div>
<div><span class="num">{report['missing_score_count']}</span><br><span class="label">分數缺值</span></div></div>
<section><h2>標註表｜一筆候選＝一列資料</h2><p>請注意 <code>synthetic_score=null</code> 表示缺分數，並非 0；分數正負號沒有臨床風險意義。</p>
<table><thead><tr><th>Variant</th><th>QC</th><th>人工分數</th><th>分數狀態</th><th>來源</th><th>審閱</th></tr></thead><tbody>{annotation_rows}</tbody></table></section>
<section><h2>全部 QC 紀錄（保留失敗原因）</h2><table><thead><tr><th>Variant</th><th>QC Status</th><th>Reason</th></tr></thead><tbody>{qc_rows}</tbody></table></section>
<footer>schema_version: {h(report['schema_version'])} · input provenance SHA-256 另見 JSON · official_avi=false · real_api_calls_this_run=0 · 研究教育用途。</footer></body></html>'''


def emit_artifacts(report, output_dir):
    """Write deterministic JSON/CSV/HTML without external services."""
    folder = Path(output_dir)
    folder.mkdir(parents=True, exist_ok=True)
    json_file = folder / "day23_annotation_report.json"
    csv_file = folder / "day23_annotation_table.csv"
    html_file = folder / "day23_annotation_panel.html"
    json_file.write_text(json.dumps(report, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    with csv_file.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=ANNOTATION_FIELDS, extrasaction="raise")
        writer.writeheader()
        writer.writerows(report["annotations"])
    html_file.write_text(render_html(report), encoding="utf-8")
    return {"status": "PASS", "mode": "OFFLINE_ONLY",
            "candidate_count": report["candidate_count"],
            "selected_count": report["selected_count"],
            "annotated_count": report["annotated_count"],
            "missing_score_count": report["missing_score_count"],
            "json_file": str(json_file), "csv_file": str(csv_file), "html_file": str(html_file),
            "real_api_calls_this_run": 0}


def execute_offline(vcf, reference, regions, scores, output_dir, input_day21=None):
    if input_day21 is not None:
        path = Path(input_day21)
        packet = json.loads(path.read_text(encoding="utf-8"))
        provenance = {"day21_review_packet": {
            "filename": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }}
    else:
        paths = {"vcf": Path(vcf), "reference": Path(reference),
                 "regions": Path(regions), "scores": Path(scores)}
        data = {key: path.read_text(encoding="utf-8") for key, path in paths.items()}
        packet = build_review_packet(
            data["vcf"], load_toy_reference(data["reference"]),
            parse_toy_bed(data["regions"]), load_toy_scores(data["scores"])
        )
        provenance = {
            key: {"filename": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            for key, path in paths.items()
        }
    return emit_artifacts(annotation_schema(packet, provenance), output_dir)


def main(argv=None):
    p = argparse.ArgumentParser(description="Day23 offline synthetic annotation and schema lab")
    p.add_argument("--mode", choices=["offline"], required=True, help="LIVE is intentionally disabled")
    p.add_argument("--input-day21", help="Optional prior Day21 synthetic JSON; excludes non-synthetic data")
    p.add_argument("--vcf", default="day21_synthetic.vcf")
    p.add_argument("--reference", default="day21_toy_reference.json")
    p.add_argument("--regions", default="day21_synthetic_regions.bed")
    p.add_argument("--scores", default="day21_synthetic_scores.json")
    p.add_argument("--output-dir", default="day23_output")
    args = p.parse_args(argv)
    result = execute_offline(args.vcf, args.reference, args.regions, args.scores,
                             args.output_dir, input_day21=args.input_day21)
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))
    return result


if __name__ == "__main__":
    main()
