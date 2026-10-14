"""Day 21: offline-only synthetic VCF -> QC -> annotation -> review packet.

Educational example using an invented chr1 reference, BED and annotation scores.
Not a complete VCF parser, genome annotation pipeline, or medical classifier.
No SDK, HTTP, credential loading, patient data, or real API calls.
Python 3.12+ standard library only.
"""
from __future__ import annotations

import argparse
import hashlib
from html import escape
import json
import math
from pathlib import Path
import re

VCF_FIELDS = ("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
ASSEMBLY = "TOY_NOT_GRCH38"
SOURCE_FLAG = "##source=SYNTHETIC_EDUCATIONAL_DEMO"
REFERENCE_FLAG = "##reference=TOY_NOT_GRCH38"
ALLELE_RE = re.compile(r"[ACGT]+\Z")
VARIANT_RE = re.compile(r"chr[0-9XYM]+:[1-9][0-9]*:[ACGT]+>[ACGT]+\Z")
QC_ORDER = ["FILTER", "BED", "REF", "INDEL", "MISSING_QUAL", "LOW_QUAL"]


def parse_synthetic_vcf(text):
    """Parse just the documented 8-column, sample-free synthetic VCF subset."""
    if not isinstance(text, str) or SOURCE_FLAG not in text.splitlines() or REFERENCE_FLAG not in text.splitlines():
        raise ValueError("VCF must explicitly identify its synthetic source and toy reference")
    seen_header = False
    rows = []
    variant_ids = set()
    for line_number, line in enumerate(text.splitlines(), start=1):
        if not line.strip() or line.startswith("##"):
            continue
        if line.startswith("#"):
            if seen_header or tuple(line.split("\t")) != VCF_FIELDS:
                raise ValueError("Only the eight-column educational VCF header is supported")
            seen_header = True
            continue
        if not seen_header:
            raise ValueError("VCF records require the eight-column header")
        cells = line.split("\t")
        if len(cells) != len(VCF_FIELDS):
            raise ValueError(f"Line {line_number}: expected 8 VCF fields")
        chrom, pos_raw, _identifier, ref, alts, qual_raw, filt, _info = cells
        if not re.fullmatch(r"chr[0-9XYM]+", chrom):
            raise ValueError(f"Line {line_number}: invalid chromosome")
        if not pos_raw.isdecimal() or int(pos_raw) < 1:
            raise ValueError(f"Line {line_number}: invalid 1-based POS")
        pos = int(pos_raw)
        if ALLELE_RE.fullmatch(ref) is None or not alts:
            raise ValueError(f"Line {line_number}: invalid REF/ALT")
        if qual_raw == ".":
            qual = None
        else:
            try:
                qual = float(qual_raw)
            except ValueError as exc:
                raise ValueError(f"Line {line_number}: invalid QUAL") from exc
            if not math.isfinite(qual) or qual < 0:
                raise ValueError(f"Line {line_number}: QUAL must be finite and nonnegative")
        if not filt or filt == ".":
            raise ValueError(f"Line {line_number}: FILTER must be present")
        for alt in alts.split(","):
            if ALLELE_RE.fullmatch(alt) is None or alt == ref:
                raise ValueError(f"Line {line_number}: invalid ALT allele")
            vid = f"{chrom}:{pos}:{ref}>{alt}"
            if vid in variant_ids:
                raise ValueError(f"Line {line_number}: duplicate variant identifier")
            variant_ids.add(vid)
            rows.append({"variant_id": vid, "chrom": chrom, "pos_1based": pos,
                         "ref": ref, "alt": alt, "qual": qual, "filter": filt,
                         "source_line": line_number})
    if not seen_header or not rows:
        raise ValueError("Require synthetic header and at least one variant")
    return rows


def load_toy_reference(text):
    obj = json.loads(text)
    if not isinstance(obj, dict) or obj.get("data_kind") != "SYNTHETIC_FIXTURE" or obj.get("assembly") != ASSEMBLY:
        raise ValueError("Reference must be an explicit synthetic fixture")
    chromosomes = obj.get("chromosomes")
    if not isinstance(chromosomes, dict) or not chromosomes:
        raise ValueError("Reference chromosomes object is required")
    for chrom, sequence in chromosomes.items():
        if not isinstance(chrom, str) or not isinstance(sequence, str) or ALLELE_RE.fullmatch(sequence) is None:
            raise ValueError("Invalid toy reference sequence")
    return chromosomes


def parse_toy_bed(text):
    lines = [line for line in text.splitlines() if line.strip() and not line.startswith("#")]
    if not text.startswith("#SYNTHETIC_FIXTURE") or not lines:
        raise ValueError("BED must be an explicit synthetic fixture")
    intervals = []
    for line in lines:
        cells = line.split("\t")
        if len(cells) != 3 or not re.fullmatch(r"chr[0-9XYM]+", cells[0]):
            raise ValueError("BED must have chrom, start, end")
        try:
            start, end = int(cells[1]), int(cells[2])
        except ValueError as exc:
            raise ValueError("BED coordinates must be integers") from exc
        if start < 0 or end <= start:
            raise ValueError("BED requires nonempty 0-based half-open intervals")
        intervals.append((cells[0], start, end))
    return intervals


def load_toy_scores(text):
    obj = json.loads(text)
    if not isinstance(obj, dict) or obj.get("data_kind") != "SYNTHETIC_FIXTURE" or obj.get("assembly") != ASSEMBLY:
        raise ValueError("Scores must be an explicit synthetic fixture")
    scores = obj.get("scores")
    if not isinstance(scores, dict):
        raise ValueError("Score mapping must be an object")
    for variant, value in scores.items():
        if VARIANT_RE.fullmatch(variant) is None or isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError("Invalid synthetic score key/value")
    return scores


def in_bed(chrom, pos_1based, intervals):
    anchor_0based = pos_1based - 1
    return any(c == chrom and start <= anchor_0based < end
               for c, start, end in intervals)


def check_qc(row, reference, intervals, min_qual=20.0):
    """Return first deciding QC status/reason (not all possible reasons)."""
    if row["filter"] != "PASS":
        return "REJECTED", "FILTER_NOT_PASS"
    if not in_bed(row["chrom"], row["pos_1based"], intervals):
        return "REJECTED", "OUTSIDE_TEACHING_BED"
    seq = reference.get(row["chrom"], "")
    start = row["pos_1based"] - 1
    if seq[start:start + len(row["ref"])] != row["ref"]:
        return "REJECTED", "REF_MISMATCH"
    if len(row["ref"]) != 1 or len(row["alt"]) != 1:
        return "REVIEW_REQUIRED", "INDEL_NOT_NORMALIZED"
    if row["qual"] is None:
        return "REVIEW_REQUIRED", "MISSING_QUAL"
    if row["qual"] < min_qual:
        return "REJECTED", "LOW_QUAL"
    return "SELECTED", "TEACHING_QC_PASS"


def build_review_packet(vcf_text, reference, intervals, scores):
    all_candidates = parse_synthetic_vcf(vcf_text)
    qc_rows = []
    review_packet = []
    for row in all_candidates:
        status, reason = check_qc(row, reference, intervals)
        qc_rows.append({"variant_id": row["variant_id"], "qc_status": status,
                        "reason": reason, "qual": row["qual"],
                        "filter": row["filter"], "source_line": row["source_line"]})
        if status == "SELECTED":
            score = scores.get(row["variant_id"])
            review_packet.append({"variant_id": row["variant_id"],
                                  "synthetic_score": score,
                                  "score_status": "MISSING_SCORE" if score is None else "SYNTHETIC_SCORE",
                                  "review_status": "REVIEW_REQUIRED",
                                  "evidence_class": "SYNTHETIC_FIXTURE"})
    report = {
        "evidence_class": "SYNTHETIC_EDUCATIONAL_DEMO",
        "mode": "OFFLINE_ONLY",
        "assembly": ASSEMBLY,
        "candidate_count": len(all_candidates),
        "selected_count": len(review_packet),
        "review_required_qc_count": sum(row["qc_status"] == "REVIEW_REQUIRED" for row in qc_rows),
        "rejected_count": sum(row["qc_status"] == "REJECTED" for row in qc_rows),
        "missing_score_count": sum(row["synthetic_score"] is None for row in review_packet),
        "qc_order": QC_ORDER,
        "qc_rows": qc_rows,
        "review_packet": review_packet,
        "score_source": "SYNTHETIC_FIXTURE_ONLY",
        "atlas_query_performed": False,
        "clinical_interpretation": False,
        "real_api_calls_this_run": 0,
    }
    if report["selected_count"] + report["review_required_qc_count"] + report["rejected_count"] != report["candidate_count"]:
        raise RuntimeError("QC group counts do not match candidates")
    return report


def render_html(report):
    """Self-contained static, offline teaching UI; escape ALL external strings."""
    def col(value):
        if value is None:
            return '<span class="missing">null</span>'
        return escape(str(value), quote=True)

    qc = "\n".join(
        '<tr><td>' + col(r["variant_id"]) + '</td><td><span class="pill">' +
        col(r["qc_status"]) + '</span></td><td>' + col(r["reason"]) +
        '</td><td>' + col(r["qual"]) + '</td></tr>' for r in report["qc_rows"])
    packets = "\n".join(
        '<tr><td>' + col(r["variant_id"]) + '</td><td>' +
        col(r["synthetic_score"]) + '</td><td>' + col(r["score_status"]) +
        '</td><td>' + col(r["review_status"]) + '</td></tr>' for r in report["review_packet"])
    return f'''<!doctype html>
<html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Day 21 | Synthetic VCF Review Packet</title><style>
:root{{color-scheme:light;font-family:system-ui,"Noto Sans TC",sans-serif;color:#283348;background:#f7f5f2}}
body{{margin:0 auto;padding:28px 22px 60px;max-width:1100px}}
h1{{font-size:1.8rem;margin:0 0 6px}}.intro{{color:#536378}}.warn{{background:#fff2dd;border-left:5px solid #bd783b;padding:14px 18px;margin:20px 0;border-radius:8px}}
.stats{{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin:20px 0}}
.stat{{background:#fff;border:1px solid #e1e2e5;border-radius:14px;padding:18px}}.num{{font-size:2rem;font-weight:750;color:#285b87}}.label{{font-size:.88rem;color:#5b6982}}
section{{background:#fff;border:1px solid #e1e2e5;border-radius:14px;padding:18px;margin-top:20px;overflow-x:auto}}
table{{width:100%;border-collapse:collapse;font-size:.9rem}}th,td{{text-align:left;padding:11px;border-bottom:1px solid #eee}}th{{background:#eef2fa}}.pill{{font-size:.8rem;background:#e6edf4;border-radius:7px;padding:4px 7px}}.missing{{color:#ae3b41;font-weight:700}}
footer{{font-size:.83rem;color:#5b6982;margin-top:30px}}</style></head><body>
<h1>Day 21 | 從 VCF 到研究審閱封包</h1><div class="intro">QC → Fixture 註釋 → Evidence / Human Review</div>
<div class="warn"><strong>SYNTHETIC EDUCATIONAL DEMO</strong>：所有資料人工合成，非 GRCh38、非 AlphaGenome API 分數，未進行臨床判定。OFFLINE_ONLY。</div>
<div class="stats"><div class="stat"><div class="num">{report['candidate_count']}</div><div class="label">ALT 候選</div></div>
<div class="stat"><div class="num">{report['selected_count']}</div><div class="label">SELECTED</div></div>
<div class="stat"><div class="num">{report['review_required_qc_count']}</div><div class="label">QC REVIEW_REQUIRED</div></div>
<div class="stat"><div class="num">{report['rejected_count']}</div><div class="label">REJECTED</div></div></div>
<section><h2>第一站｜QC 全部候選（保留失敗理由）</h2><p>順序：{col(' → '.join(report['qc_order']))}；每筆記錄第一個決定性原因。</p>
<table><thead><tr><th>Variant</th><th>QC Status</th><th>Reason</th><th>QUAL</th></tr></thead><tbody>{qc}</tbody></table></section>
<section><h2>第二站｜可供人工審閱的候選封包</h2><p>只有 SELECTED 進入此表；但所有八筆仍保留在 QC 表。缺分數記為 null，不填零。</p>
<table><thead><tr><th>Variant</th><th>Artificial Score</th><th>Score Status</th><th>Review Status</th></tr></thead><tbody>{packets}</tbody></table></section>
<footer>Evidence: SYNTHETIC_FIXTURE_ONLY · assembly: {col(report['assembly'])} · atlas_query_performed: false · real_api_calls_this_run: 0 · 研究用教學，非醫療建議。</footer>
</body></html>'''


def hash_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def execute_offline(vcf_path, reference_path, regions_path, scores_path, output_dir):
    paths = {"vcf": Path(vcf_path), "reference": Path(reference_path),
             "regions": Path(regions_path), "scores": Path(scores_path)}
    data = {key: path.read_text(encoding="utf-8") for key, path in paths.items()}
    report = build_review_packet(data["vcf"], load_toy_reference(data["reference"]),
                                 parse_toy_bed(data["regions"]), load_toy_scores(data["scores"]))
    report["input_provenance"] = {key: {"filename": path.name, "sha256": hash_file(path)}
                                  for key, path in paths.items()}
    folder = Path(output_dir)
    folder.mkdir(parents=True, exist_ok=True)
    report_path = folder / "day21_review_packet.json"
    html_path = folder / "day21_review_panel.html"
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    html_path.write_text(render_html(report), encoding="utf-8")
    return {"status": "PASS", "mode": "OFFLINE_ONLY",
            "candidate_count": report["candidate_count"],
            "selected_count": report["selected_count"],
            "review_required_qc_count": report["review_required_qc_count"],
            "rejected_count": report["rejected_count"],
            "missing_score_count": report["missing_score_count"],
            "json_file": str(report_path), "html_file": str(html_path),
            "real_api_calls_this_run": 0}


def main(argv=None):
    parser = argparse.ArgumentParser(description="Offline synthetic VCF to research review packet")
    parser.add_argument("--mode", choices=["offline"], required=True, help="LIVE intentionally unavailable")
    parser.add_argument("--vcf", default="day21_synthetic.vcf")
    parser.add_argument("--reference", default="day21_toy_reference.json")
    parser.add_argument("--regions", default="day21_synthetic_regions.bed")
    parser.add_argument("--scores", default="day21_synthetic_scores.json")
    parser.add_argument("--output-dir", default="day21_output")
    args = parser.parse_args(argv)
    summary = execute_offline(args.vcf, args.reference, args.regions, args.scores, args.output_dir)
    print(json.dumps(summary, sort_keys=True, ensure_ascii=False))


if __name__ == "__main__":
    main()
