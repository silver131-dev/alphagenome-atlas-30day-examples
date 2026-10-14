"""Day 25 | OFFLINE genomic-locus REF/ALT visualization lesson.

This reads Day 24's SYNTHETIC report and its original SYNTHETIC track fixture.
It does not call AlphaGenome, Atlas, Gemini, a network endpoint or a genome FASTA.
Only Python 3.12+ standard library modules are required for the student lab.
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
import xml.etree.ElementTree as ET

EVIDENCE_CLASS = "SYNTHETIC_EDUCATIONAL_DEMO"
MODE = "OFFLINE_ONLY"
ASSEMBLY = "TOY_NOT_GRCH38"
INPUT_SCHEMA = "day24.synthetic.refalt_delta.v1"
TRACK_SCHEMA = "day24.synthetic.refalt_tracks.v1"
FIGURE_SCHEMA = "day25.synthetic.locus_plot.v1"
VARIANT_RE = re.compile(r"(chr(?:[1-9]|1[0-9]|2[0-2]|X|Y|M)):([1-9][0-9]*):([ACGT])>([ACGT])\Z")
SHA_RE = re.compile(r"[0-9a-f]{64}\Z")
SVG_NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG_NS)


def load_json(path: str | Path, *, max_bytes: int = 1_000_000):
    """Read a bounded local fixture; NEVER download or deserialize pickle files."""
    file = Path(path)
    if not file.is_file() or file.stat().st_size > max_bytes:
        raise ValueError("Local fixture missing or exceeds 1 MB")
    try:
        return json.loads(file.read_text(encoding="utf-8"), parse_constant=lambda s: (_ for _ in ()).throw(ValueError("Non-finite JSON is invalid")))
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError("Invalid local UTF-8 JSON") from exc


def sha256_file(path: str | Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _exact(record, key, expected):
    if key not in record or type(record[key]) is not type(expected) or record[key] != expected:
        raise ValueError(f"Incorrect synthetic input contract: {key}")


def _finite_unit_number(value):
    return type(value) in (float, int) and math.isfinite(value) and 0 <= value <= 1


def verify_inputs(report, track, *, track_sha=None):
    """Reject model, patient, malformed or mismatched data before plotting."""
    if not isinstance(report, dict) or not isinstance(track, dict):
        raise ValueError("Day24 report and track must be JSON objects")
    contract = {
        "schema_version": INPUT_SCHEMA, "evidence_class": EVIDENCE_CLASS,
        "mode": MODE, "assembly": ASSEMBLY,
        "track_label": "ARTIFICIAL_SIGNAL_PROXY",
        "track_unit": "UNIT_INTERVAL_0_TO_1_TOY",
        "biosample": "TOY_CELL_CONTEXT_NOT_REAL_TISSUE",
        "coordinate_system": "1-based positions", "delta_definition": "ALT_proxy - REF_proxy",
        "mechanism_status": "HYPOTHESIS_ONLY_NOT_VALIDATED",
        "model_prediction": False, "official_atlas_data": False,
        "reference_allele_verified": False, "clinical_interpretation": False,
        "atlas_query_performed": False, "real_api_calls_this_run": 0,
    }
    for key, expected in contract.items():
        _exact(report, key, expected)
    track_contract = {
        "schema_version": TRACK_SCHEMA, "evidence_class": EVIDENCE_CLASS,
        "mode": MODE, "assembly": ASSEMBLY,
        "track_label": "ARTIFICIAL_SIGNAL_PROXY",
        "track_unit": "UNIT_INTERVAL_0_TO_1_TOY",
        "biosample": "TOY_CELL_CONTEXT_NOT_REAL_TISSUE",
        "coordinate_system": "1-based positions", "model_prediction": False,
        "official_atlas_data": False, "real_api_calls_this_run": 0,
    }
    for key, expected in track_contract.items():
        _exact(track, key, expected)
    variant_id = report.get("variant_id")
    m = VARIANT_RE.fullmatch(variant_id) if isinstance(variant_id, str) else None
    if not m or m.group(3) == m.group(4) or track.get("variant_id") != variant_id:
        raise ValueError("Invalid or mismatched synthetic SNV")
    variant_pos = int(m.group(2))
    rows = report.get("point_rows")
    pos, refs, alts = (track.get(k) for k in ("positions_1based", "ref", "alt"))
    if not all(isinstance(v, list) for v in (rows, pos, refs, alts)):
        raise ValueError("Point arrays must be lists")
    if not (2 <= len(rows) <= 1000 and len(rows) == len(pos) == len(refs) == len(alts)):
        raise ValueError("Point arrays must be equally sized, 2-1000 points")
    if type(report.get("point_count")) is not int or report["point_count"] != len(rows):
        raise ValueError("Day24 point_count is incorrect")
    for i, (row, p, ref, alt) in enumerate(zip(rows, pos, refs, alts)):
        if not isinstance(row, dict) or type(p) is not int or p < 1:
            raise ValueError("Malformed position")
        if i and p != pos[i - 1] + 1:
            raise ValueError("Nonconsecutive 1-based positions")
        if not (_finite_unit_number(ref) and _finite_unit_number(alt)):
            raise ValueError("Invalid or non-finite proxy value")
        expected = {
            "position_1based": p, "ref_proxy": ref, "alt_proxy": alt,
            "delta_alt_minus_ref": round(alt - ref, 3),
            "variant_here": p == variant_pos,
        }
        for k, v in expected.items():
            if k not in row or type(row[k]) is not type(v) or row[k] != v:
                # int vs float is allowed for scientifically identical proxy values only
                if k in ("ref_proxy", "alt_proxy", "delta_alt_minus_ref") and _finite_unit_number(row.get(k)) and abs(row[k] - v) <= 1e-9:
                    continue
                raise ValueError(f"Day24 row {i} does not match fixture: {k}")
    if not (pos[0] <= variant_pos <= pos[-1]):
        raise ValueError("Variant outside plotted interval")
    peak_i = max(range(len(rows)), key=lambda i: abs(rows[i]["delta_alt_minus_ref"]))
    _exact(report, "largest_absolute_delta_position_1based", pos[peak_i])
    if type(report.get("largest_absolute_delta")) not in (float, int) or report["largest_absolute_delta"] != rows[peak_i]["delta_alt_minus_ref"]:
        raise ValueError("Largest signed delta incorrectly recorded")
    provenance = report.get("input_provenance")
    if not isinstance(provenance, dict) or not isinstance(provenance.get("track_fixture"), dict):
        raise ValueError("Day24 provenance required")
    reference = provenance["track_fixture"]
    if not isinstance(reference.get("filename"), str) or not SHA_RE.fullmatch(str(reference.get("sha256"))):
        raise ValueError("Malformed upstream track hash")
    if track_sha is not None and track_sha != reference["sha256"]:
        raise ValueError("Source track SHA-256 does not match Day24 report")
    # The score from Day23 is a different proxy. It is not plotted or mixed.
    link = report.get("day23_annotation_link")
    if link is not None:
        if not isinstance(link, dict) or link.get("variant_id") != variant_id or link.get("qc_status") != "SELECTED" or link.get("review_status") != "REVIEW_REQUIRED":
            raise ValueError("Inconsistent Day23 link")
        if link.get("score_status") not in ("MISSING_SCORE", "SYNTHETIC_SCORE"):
            raise ValueError("Day23 status invalid")
        if link["score_status"] == "MISSING_SCORE" and link.get("synthetic_score") is not None:
            raise ValueError("Day23 missing score must remain null")
        if link["score_status"] == "SYNTHETIC_SCORE" and (type(link.get("synthetic_score")) not in (int, float) or not math.isfinite(link["synthetic_score"])):
            raise ValueError("Day23 synthetic score must be finite")
    return report


def make_figure_summary(report, *, report_sha: str, track_sha: str):
    verify_inputs(report, {
        "schema_version": TRACK_SCHEMA, "evidence_class": EVIDENCE_CLASS,
        "mode": MODE, "assembly": ASSEMBLY,
        "track_label": report["track_label"], "track_unit": report["track_unit"],
        "biosample": report["biosample"], "coordinate_system": report["coordinate_system"],
        "variant_id": report["variant_id"],
        "positions_1based": [r["position_1based"] for r in report["point_rows"]],
        "ref": [r["ref_proxy"] for r in report["point_rows"]],
        "alt": [r["alt_proxy"] for r in report["point_rows"]],
        "model_prediction": False, "official_atlas_data": False,
        "real_api_calls_this_run": 0,
    }, track_sha=track_sha)
    return {
        "schema_version": FIGURE_SCHEMA,
        "evidence_class": EVIDENCE_CLASS,
        "mode": MODE,
        "assembly": ASSEMBLY,
        "variant_id": report["variant_id"],
        "coordinate_system": "1-based positions",
        "track_unit": "UNIT_INTERVAL_0_TO_1_TOY",
        "point_count": report["point_count"],
        "track_count": 2,
        "delta_definition": "ALT_proxy - REF_proxy",
        "largest_absolute_delta_position_1based": report["largest_absolute_delta_position_1based"],
        "largest_absolute_delta": report["largest_absolute_delta"],
        "day23_linked": report.get("day23_annotation_link") is not None,
        "day23_score_combined_with_delta": False,
        "source_day24_report_sha256": report_sha,
        "source_day24_track_sha256": track_sha,
        "upstream_day23_source_sha256": (report.get("input_provenance", {}).get("day23_packet") or {}).get("sha256"),
        "upstream_native_model_response_retained": False,
        "model_prediction": False,
        "official_atlas_data": False,
        "reference_allele_verified": False,
        "clinical_interpretation": False,
        "real_api_calls_this_run": 0,
        "figure_caption": "SYNTHETIC proxy signals. Five points; not actual AlphaGenome outputs or clinical evidence.",
    }


def draw_svg(report):
    """Render a single 1920×1080 16:9 figure with aligned signed-delta panel."""
    rows = report["point_rows"]
    variant_id = report["variant_id"]
    n = len(rows)
    variant_pos = int(VARIANT_RE.fullmatch(variant_id).group(2))
    root = ET.Element(f"{{{SVG_NS}}}svg", {"width": "1920", "height": "1080", "viewBox": "0 0 1920 1080", "role": "img", "aria-label": "Synthetic REF ALT genomic track with aligned delta; not model predictions"})

    def elem(tag, attrs=None, label=None, parent=None):
        child = ET.SubElement(parent if parent is not None else root, f"{{{SVG_NS}}}{tag}", {k: str(v) for k, v in (attrs or {}).items()})
        if label is not None:
            child.text = str(label)
        return child

    def text(x, y, label, size=24, fill="#444956", weight="normal", anchor=None):
        a = {"x": x, "y": y, "font-family": "Arial, sans-serif", "font-size": size, "fill": fill, "font-weight": weight}
        if anchor:
            a["text-anchor"] = anchor
        return elem("text", a, label)

    elem("rect", {"width": 1920, "height": 1080, "fill": "#E5DEDA"})
    text(108, 88, "Day 25  |  SYNTHETIC REF / ALT GENOMIC TRACKS", 49, "#323c4c", "bold")
    text(110, 133, "Same coordinates  ·  two proxy tracks  ·  one signed difference", 26, "#5d6570")
    elem("rect", {"x": 110, "y": 162, "width": 1700, "height": 66, "rx": 18, "fill": "#f6f0ed"})
    text(140, 205, variant_id, 25, "#4d5d78", "bold")
    text(530, 205, "TOY_NOT_GRCH38  ·  1-based  ·  artificial units 0–1", 23, "#65707e")
    text(1530, 205, "NOT MODEL DATA", 22, "#9d4655", "bold")
    elem("rect", {"x": 105, "y": 256, "width": 1710, "height": 392, "rx": 22, "fill": "#fffaf8"})
    elem("rect", {"x": 105, "y": 668, "width": 1710, "height": 218, "rx": 22, "fill": "#fffaf8"})
    text(150, 310, "REF / ALT proxy signal", 30, "#394457", "bold")
    text(150, 720, "Signed Delta = ALT − REF", 30, "#394457", "bold")
    left, right = 285, 1710
    top, bottom = 349, 567
    x = lambda i: left + (right - left) * i / (n - 1)
    y = lambda v: bottom - (bottom - top) * v
    for mark in (0, 0.25, 0.5, 0.75, 1.0):
        yy = y(mark)
        elem("line", {"x1": left, "y1": f"{yy:.2f}", "x2": right, "y2": f"{yy:.2f}", "stroke": "#e4e7e8", "stroke-width": 2})
        text(245, yy + 9, f"{mark:g}", 23, "#727b84", anchor="end")
    for name, field, color in (("REF", "ref_proxy", "#507d98"), ("ALT", "alt_proxy", "#c27486")):
        coords = " ".join(f"{x(i):.2f},{y(row[field]):.2f}" for i, row in enumerate(rows))
        elem("polyline", {"points": coords, "stroke": color, "stroke-width": 7, "stroke-linecap": "round", "stroke-linejoin": "round", "fill": "none", "data-track": name})
        for i, row in enumerate(rows):
            dot = elem("circle", {"cx": f"{x(i):.2f}", "cy": f"{y(row[field]):.2f}", "r": 9, "fill": color, "stroke": "#fff", "stroke-width": 2, "data-position": row["position_1based"], "data-track": name})
            elem("title", label=f"{name} position {row['position_1based']}: {row[field]:g}", parent=dot)
    for idx, (name, color) in enumerate((("REF", "#507d98"), ("ALT", "#c27486"))):
        start = 1480 + idx * 140
        elem("line", {"x1": start, "x2": start + 35, "y1": 309, "y2": 309, "stroke": color, "stroke-width": 7})
        text(start + 44, 319, name, 24, color, "bold")
    for i, row in enumerate(rows):
        if n <= 12 or i % max(1, n // 9) == 0 or i == n - 1:
            text(x(i), 612, str(row["position_1based"]), 24, "#505d6a", anchor="middle")
    pos_i = next(i for i, r in enumerate(rows) if r["position_1based"] == variant_pos)
    elem("line", {"x1": f"{x(pos_i):.2f}", "y1": 324, "x2": f"{x(pos_i):.2f}", "y2": 853, "stroke": "#a18458", "stroke-width": 3, "stroke-dasharray": "10,9", "data-variant-marker": str(variant_pos)})
    text(x(pos_i) + 16, 348, "Variant", 21, "#907248", "bold")
    deltas = [r["delta_alt_minus_ref"] for r in rows]
    scale = max(max(abs(v) for v in deltas), .1)
    zero = 784
    elem("line", {"x1": left, "y1": zero, "x2": right, "y2": zero, "stroke": "#9ea6a7", "stroke-width": 2.5})
    for i, (row, d) in enumerate(zip(rows, deltas)):
        e = zero - d / scale * 65
        if d:
            elem("line", {"x1": f"{x(i):.2f}", "y1": zero, "x2": f"{x(i):.2f}", "y2": f"{e:.2f}", "stroke": "#937993", "stroke-width": 13, "stroke-linecap": "round", "data-delta-pos": str(row["position_1based"])})
        text(x(i), 861, f"{d:+.2f}", 22, "#5d5472", anchor="middle")
    p = report["largest_absolute_delta_position_1based"]
    d = report["largest_absolute_delta"]
    text(110, 926, f"Largest absolute Delta at 1-based {p}: signed value {d:+.2f}", 27, "#384459", "bold")
    text(110, 968, "SYNTHETIC EDUCATIONAL DATA  |  No Atlas API, no genome reference validation, no clinical inference", 23, "#884f5a")
    # Deliberately leave a clear lower margin for future video captions.
    return ET.tostring(root, encoding="unicode")


def draw_html(report, summary, svg):
    """Create browser-readable offline lesson with SVG and fully escaped text."""
    h = lambda t: escape(str(t), quote=True)
    rows = "\n".join("<tr>" + "".join(f"<td>{h(r[f])}</td>" for f in ("position_1based", "ref_proxy", "alt_proxy", "delta_alt_minus_ref", "variant_here")) + "</tr>" for r in report["point_rows"])
    return f'''<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Day 25 | 合成基因座 REF/ALT 軌跡</title>
<style>body{{margin:0;background:#E5DEDA;color:#293545;font-family:system-ui, sans-serif}}main{{max-width:1240px;margin:auto;padding:25px}}.card{{background:#fffaf8;border-radius:20px;padding:22px;margin:18px 0}}svg{{width:100%;height:auto;border-radius:12px}}table{{border-collapse:collapse;width:100%}}td,th{{padding:10px;border-bottom:1px solid #e2d4d0;text-align:right}}td:first-child,th:first-child{{text-align:left}}.warning{{color:#873a51;font-weight:700}}.nav a{{margin-right:16px;color:#35587b}}</style></head>
<body><main><h1>Day 25 | 合成基因座軌跡圖</h1><p>REF 與 ALT 使用共同 1-based 位置；下方 Delta = ALT − REF，負號不可消失。</p>
<div class="card">{svg}</div>
<div class="card"><h2>逐點數值核對</h2><table><thead><tr><th>位置 1-based</th><th>REF proxy</th><th>ALT proxy</th><th>ALT − REF</th><th>變異位置？</th></tr></thead><tbody>{rows}</tbody></table></div>
<div class="card"><h2>證據說明</h2><p class="warning">SYNTHETIC_EDUCATIONAL_DEMO · TOY_NOT_GRCH38 · 0 real API calls</p>
<p>這是 Day 24 人工 JSON 的繪圖，不是 AlphaGenome Atlas／Model 的預測結果。Day 23 的人工分數不會拿來與此 Delta 加總。</p>
<p>Day24 report SHA-256: <code>{h(summary['source_day24_report_sha256'])}</code></p>
<p>Day24 track SHA-256: <code>{h(summary['source_day24_track_sha256'])}</code></p>
<p>模型機制：NOT_VALIDATED · 人類審閱：REQUIRED · 臨床解讀：不適用</p>
<p class="nav"><a href="day25_locus_track.csv">CSV 數據表</a><a href="day25_locus_summary.json">JSON 概要</a><a href="day25_locus_tracks.svg">SVG 圖</a></p></div></main></body></html>'''


def emit_outputs(report_path, source_track_path, output_dir):
    report = load_json(report_path)
    track = load_json(source_track_path)
    track_sha = sha256_file(source_track_path)
    verify_inputs(report, track, track_sha=track_sha)
    summary = make_figure_summary(report, report_sha=sha256_file(report_path), track_sha=track_sha)
    svg = draw_svg(report)
    folder = Path(output_dir)
    folder.mkdir(parents=True, exist_ok=True)
    json_path = folder / "day25_locus_summary.json"
    csv_path = folder / "day25_locus_track.csv"
    svg_path = folder / "day25_locus_tracks.svg"
    html_path = folder / "day25_locus_view.html"
    json_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["position_1based", "ref_proxy", "alt_proxy", "delta_alt_minus_ref", "variant_here"])
        writer.writeheader()
        writer.writerows(report["point_rows"])
    svg_path.write_text(svg, encoding="utf-8")
    html_path.write_text(draw_html(report, summary, svg), encoding="utf-8")
    return {
        "status": "PASS", "mode": MODE, "evidence_class": EVIDENCE_CLASS,
        "variant_id": summary["variant_id"],
        "point_count": summary["point_count"], "track_count": 2,
        "largest_absolute_delta_position_1based": summary["largest_absolute_delta_position_1based"],
        "largest_absolute_delta": summary["largest_absolute_delta"],
        "json_file": str(json_path), "csv_file": str(csv_path),
        "svg_file": str(svg_path), "html_file": str(html_path),
        "real_api_calls_this_run": 0,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description="Draw 16:9 synthetic REF/ALT locus tracks from Day24")
    parser.add_argument("--mode", choices=["offline"], required=True)
    parser.add_argument("--day24-report", default="day24_track_delta_report.json")
    parser.add_argument("--source-track", default="day24_synthetic_tracks.json")
    parser.add_argument("--output-dir", default="day25_output")
    args = parser.parse_args(argv)
    try:
        result = emit_outputs(args.day24_report, args.source_track, args.output_dir)
    except (ValueError, KeyError, OSError, TypeError) as exc:
        parser.error(str(exc))
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))
    return result


if __name__ == "__main__":
    main()
