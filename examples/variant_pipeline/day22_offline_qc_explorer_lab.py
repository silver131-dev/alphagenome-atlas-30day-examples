"""Day 22 — inspect the first-decision FILTER/BED/REF/indel/QUAL QC chain.

This offline educational program REUSES the verified Day 21 parsing and QC
functions, on the byte-identical synthetic fixtures shipped alongside it.
It does not implement the full VCF standard, variant normalization, genome
alignment, clinical interpretation, or any AlphaGenome network call.
Python 3.12+ standard library only.
"""
from __future__ import annotations

import argparse
from collections import Counter
from html import escape
import hashlib
import json
import math
from pathlib import Path

from day21_offline_vcf_review_lab import (
    ASSEMBLY, QC_ORDER, check_qc, in_bed,
    load_toy_reference, parse_synthetic_vcf, parse_toy_bed,
)

EVIDENCE_CLASS = "SYNTHETIC_EDUCATIONAL_DEMO"
FAIL_REASON = {
    "FILTER": ("REJECTED", "FILTER_NOT_PASS"),
    "BED": ("REJECTED", "OUTSIDE_TEACHING_BED"),
    "REF": ("REJECTED", "REF_MISMATCH"),
    "INDEL": ("REVIEW_REQUIRED", "INDEL_NOT_NORMALIZED"),
    "MISSING_QUAL": ("REVIEW_REQUIRED", "MISSING_QUAL"),
    "LOW_QUAL": ("REJECTED", "LOW_QUAL"),
}


def validate_min_qual(value):
    """Validate a teaching threshold; it is NOT a clinical decision cutoff."""
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
        raise ValueError("min_qual must be a finite nonnegative number")
    return float(value)


def trace_qc(row, reference, intervals, min_qual=20.0):
    """Explain *only* the checks reached before Day21's deciding QC result.

    Checks later in the chain are deliberately NOT evaluated once a prior check
    makes the decision. This is central to the Day22 lesson.
    """
    min_qual = validate_min_qual(min_qual)
    steps = []
    status = "SELECTED"
    reason = "TEACHING_QC_PASS"

    for step in QC_ORDER:
        if step == "FILTER":
            passed = row["filter"] == "PASS"
            detail = "FILTER=PASS" if passed else "FILTER is not PASS"
        elif step == "BED":
            passed = in_bed(row["chrom"], row["pos_1based"], intervals)
            detail = f"VCF POS {row['pos_1based']} -> BED anchor {row['pos_1based'] - 1}"
        elif step == "REF":
            seq = reference.get(row["chrom"], "")
            start = row["pos_1based"] - 1
            expected = seq[start:start + len(row["ref"])]
            passed = expected == row["ref"]
            detail = "REF matches artificial reference" if passed else "REF mismatch / absent toy reference"
        elif step == "INDEL":
            passed = len(row["ref"]) == 1 and len(row["alt"]) == 1
            detail = "SNV teaching subset" if passed else "indel requires independent normalization/review"
        elif step == "MISSING_QUAL":
            passed = row["qual"] is not None
            detail = "QUAL present" if passed else "QUAL missing (dot in synthetic VCF)"
        elif step == "LOW_QUAL":
            passed = row["qual"] >= min_qual
            detail = f"QUAL compared against artificial threshold {min_qual:g}"
        else:
            raise RuntimeError(f"Unknown QC step in Day21 source module: {step}")

        steps.append({"check": step, "passed": passed, "detail": detail})
        if not passed:
            status, reason = FAIL_REASON[step]
            break

    original_status, original_reason = check_qc(row, reference, intervals, min_qual=min_qual)
    if (status, reason) != (original_status, original_reason):
        raise AssertionError("Day22 trace differs from the tested Day21 QC decision")
    return {"variant_id": row["variant_id"],
            "qc_status": status, "reason": reason,
            "bed_anchor_0based": row["pos_1based"] - 1,
            "checked_steps": steps,
            "not_checked_steps": QC_ORDER[len(steps):],
            "source_line": row.get("source_line"),
            "qual": row["qual"], "filter": row["filter"]}


def make_probe(label, pos, ref, alt, qual=60.0, filt="PASS", chrom="chr1"):
    return {"label": label, "variant_id": f"{chrom}:{pos}:{ref}>{alt}",
            "chrom": chrom, "pos_1based": pos, "ref": ref, "alt": alt,
            "qual": qual, "filter": filt, "source_line": None}


def extra_probe_rows():
    """Explicit isolated examples that DO exercise otherwise shadowed branches."""
    return [
        make_probe("LOW_QUAL_INSIDE_BED", 105, "A", "T", qual=10),
        make_probe("REF_MISMATCH_INSIDE_BED", 106, "G", "T"),
        make_probe("LOW_QUAL_OUTSIDE_BED", 120, "A", "T", qual=10),
        make_probe("REF_MISMATCH_OUTSIDE_BED", 140, "G", "T"),
        make_probe("BED_FIRST_INCLUDED_POS", 91, "A", "C"),
        make_probe("BED_LAST_INCLUDED_POS", 115, "A", "C"),
        make_probe("BED_FIRST_EXCLUDED_POS", 116, "A", "C"),
        make_probe("SECOND_BED_FIRST_INCLUDED_POS", 146, "A", "C"),
        make_probe("INDEL_NEEDS_REVIEW", 150, "A", "AT"),
        make_probe("MISSING_QUAL_NEEDS_REVIEW", 160, "A", "T", qual=None),
    ]


def qc_report(vcf_text, reference, intervals, min_qual=20.0):
    validate_min_qual(min_qual)
    if QC_ORDER != ["FILTER", "BED", "REF", "INDEL", "MISSING_QUAL", "LOW_QUAL"]:
        raise RuntimeError("Unexpected Day21 QC order; do not silently change lesson")
    input_rows = parse_synthetic_vcf(vcf_text)
    rows = [trace_qc(row, reference, intervals, min_qual) for row in input_rows]
    probes = [{"label": row["label"], **trace_qc(row, reference, intervals, min_qual)}
              for row in extra_probe_rows()]
    counts = Counter(row["qc_status"] for row in rows)
    reasons = Counter(row["reason"] for row in rows)
    assert sum(counts.values()) == len(input_rows)
    return {
        "data_kind": EVIDENCE_CLASS,
        "assembly": ASSEMBLY,
        "mode": "OFFLINE_ONLY",
        "qc_policy": "FIRST_DECIDING_REASON_ONLY",
        "qc_order": list(QC_ORDER),
        "bed_policy": "VCF_1BASED_POS_MINUS_ONE_ANCHOR_FOR_TOY_BED",
        "min_qual_artificial_threshold": float(min_qual),
        "candidate_count": len(input_rows),
        "selected_count": counts.get("SELECTED", 0),
        "review_required_count": counts.get("REVIEW_REQUIRED", 0),
        "rejected_count": counts.get("REJECTED", 0),
        "reason_counts": dict(sorted(reasons.items())),
        "qc_rows": rows,
        "supplementary_probe_cases": probes,
        "fixture_source": "DAY21_BYTE_IDENTICAL_SYNTHETIC_FIXTURES",
        "source_parser_and_qc_module": "day21_offline_vcf_review_lab.py",
        "real_genome_used": False,
        "real_api_calls_this_run": 0,
        "atlas_query_performed": False,
        "clinical_interpretation": False,
    }


def render_html(report):
    """Single-file static HTML dashboard; every data-derived string is escaped."""
    def q(value):
        return escape("null" if value is None else str(value), quote=True)

    rows_html = "\n".join(
        "<tr>" + "".join(f"<td>{q(value)}</td>" for value in (
            row["variant_id"], row["qc_status"], row["reason"], row["bed_anchor_0based"],
            row["qual"], " → ".join(s["check"] for s in row["checked_steps"]))) + "</tr>"
        for row in report["qc_rows"]
    )
    probes_html = "\n".join(
        "<tr>" + "".join(f"<td>{q(value)}</td>" for value in (
            row["label"], row["variant_id"], row["qc_status"], row["reason"],
            " → ".join(s["check"] for s in row["checked_steps"]))) + "</tr>"
        for row in report["supplementary_probe_cases"]
    )
    details_html = "\n".join(
        "<details><summary>" + q(row["variant_id"]) + " · " + q(row["reason"]) + "</summary>" +
        "<ol>" + "".join("<li>" + q(s["check"]) + "：" +
                            ("PASS" if s["passed"] else "STOP") + " — " + q(s["detail"]) + "</li>"
                            for s in row["checked_steps"]) + "</ol>" +
        "<p>未繼續檢查：" + q(", ".join(row["not_checked_steps"]) or "無") + "</p></details>"
        for row in report["qc_rows"]
    )
    return f'''<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Day22 | Offline VCF QC Inspector</title>
<style>:root{{font:16px/1.6 system-ui,"Noto Sans TC",sans-serif;color:#263b51;background:#f1f4f8}}
body{{max-width:1140px;margin:auto;padding:32px 18px 64px}}h1{{margin:0;font-size:1.8rem}}.lead{{color:#526477}}
.notice{{background:#fff4dc;border-left:6px solid #bc833b;padding:16px;margin:20px 0;border-radius:8px}}
.stats{{display:grid;grid-template-columns:repeat(auto-fit,minmax(130px,1fr));gap:12px}}
.stats div,section{{background:white;border:1px solid #dce4ed;border-radius:14px;padding:16px;margin-top:15px}}
.number{{font-size:2rem;font-weight:800;color:#346b97}}.sub{{font-size:.85rem;color:#667789}}
table{{width:100%;border-collapse:collapse}}th,td{{padding:9px;text-align:left;border-bottom:1px solid #e9ebed;vertical-align:top}}
th{{background:#e9eff5}}section{{overflow-x:auto}}details{{padding:12px;border-bottom:1px solid #dce4ed}}
summary{{cursor:pointer;font-weight:600}}code{{background:#f2f5f8;padding:2px 4px}}.warn{{color:#8d3639}}
footer{{font-size:.84rem;color:#53657b;margin-top:24px}}</style></head><body>
<h1>Day 22｜VCF 資料前處理 — FILTER、QUAL、REF 與 BED QC</h1>
<p class="lead">從 Day21 同一批人工候選資料，觀察「第一個決定性原因」如何產生。</p>
<div class="notice"><strong>OFFLINE · SYNTHETIC EDUCATIONAL DEMO</strong><br>
參考序列是 200 個人工 A，基因組標記 TOY_NOT_GRCH38。沒有患者、真實模型、AlphaGenome API 或臨床判讀。</div>
<div class="stats"><div><span class="number">{report['candidate_count']}</span><br><span class="sub">ALT 候選</span></div>
<div><span class="number">{report['selected_count']}</span><br><span class="sub">SELECTED</span></div>
<div><span class="number">{report['review_required_count']}</span><br><span class="sub">QC 待審</span></div>
<div><span class="number">{report['rejected_count']}</span><br><span class="sub">REJECTED</span></div></div>
<section><h2>1｜QC 主表</h2><p>判定順序：{q(' → '.join(report['qc_order']))}。</p>
<p>只記錄第一個決定性原因，不表示其他尚未執行的檢查已通過。</p>
<table><thead><tr><th>變異</th><th>狀態</th><th>原因</th><th>BED anchor</th><th>QUAL</th><th>已檢查步驟</th></tr></thead><tbody>{rows_html}</tbody></table></section>
<section><h2>2｜打開 QC 決策路徑</h2>{details_html}</section>
<section><h2>3｜另外建立、位於 BED 內的人工探針</h2>
<p>預設 POS 120 的低 QUAL、POS 140 的 REF mismatch 都先被 BED 排除。用 POS 105/106 才能真正檢查 LOW_QUAL 和 REF_MISMATCH。</p>
<table><thead><tr><th>探針</th><th>變異</th><th>狀態</th><th>原因</th><th>檢查路徑</th></tr></thead><tbody>{probes_html}</tbody></table></section>
<footer>資料來源：DAY21_BYTE_IDENTICAL_SYNTHETIC_FIXTURES · {q(report['assembly'])} · atlas_query_performed=false · real_api_calls_this_run=0。<br>
僅練習教學 QC；此處用 POS−1 作 BED anchor 判定，indel 仍須正常化及更嚴謹的跨區間處理。</footer></body></html>'''


def execute_offline(vcf_path, reference_path, regions_path, output_dir, min_qual=20.0):
    paths = {"vcf": Path(vcf_path), "reference": Path(reference_path),
             "regions": Path(regions_path)}
    raw = {kind: path.read_text(encoding="utf-8") for kind, path in paths.items()}
    report = qc_report(raw["vcf"], load_toy_reference(raw["reference"]),
                       parse_toy_bed(raw["regions"]), min_qual=min_qual)
    report["input_provenance"] = {
        kind: {"filename": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        for kind, path in paths.items()
    }
    folder = Path(output_dir)
    folder.mkdir(parents=True, exist_ok=True)
    json_path = folder / "day22_qc_report.json"
    html_path = folder / "day22_qc_dashboard.html"
    json_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    html_path.write_text(render_html(report), encoding="utf-8")
    return {"status": "PASS", "mode": "OFFLINE_ONLY",
            "candidate_count": report["candidate_count"],
            "selected_count": report["selected_count"],
            "review_required_count": report["review_required_count"],
            "rejected_count": report["rejected_count"],
            "low_qual_probe": next(p["reason"] for p in report["supplementary_probe_cases"]
                                   if p["label"] == "LOW_QUAL_INSIDE_BED"),
            "ref_mismatch_probe": next(p["reason"] for p in report["supplementary_probe_cases"]
                                       if p["label"] == "REF_MISMATCH_INSIDE_BED"),
            "json_file": str(json_path), "html_file": str(html_path),
            "real_api_calls_this_run": 0}


def main(argv=None):
    p = argparse.ArgumentParser(description="Day22 offline synthetic QC inspector (never LIVE)")
    p.add_argument("--mode", choices=["offline"], required=True)
    p.add_argument("--vcf", default="day21_synthetic.vcf")
    p.add_argument("--reference", default="day21_toy_reference.json")
    p.add_argument("--regions", default="day21_synthetic_regions.bed")
    p.add_argument("--output-dir", default="day22_output")
    p.add_argument("--min-qual", type=float, default=20.0,
                   help="Artificial QC threshold for educational purposes only")
    args = p.parse_args(argv)
    summary = execute_offline(args.vcf, args.reference, args.regions,
                              args.output_dir, min_qual=args.min_qual)
    print(json.dumps(summary, ensure_ascii=False, sort_keys=True))
    return summary


if __name__ == "__main__":
    main()
