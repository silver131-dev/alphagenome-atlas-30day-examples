"""Day 26: OFFLINE synthetic PPM-like DNA motif preference lesson.

A toy A>G substitution changes one position in an artificial DNA sequence.
This is NOT AlphaGenome, Atlas, TF-binding affinity, log-odds PWM or wet lab.
No network access, API key, patient data or third-party packages are needed.
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
INPUT_SCHEMA = "day26.synthetic.ppm_context.v1"
OUTPUT_SCHEMA = "day26.synthetic.ppm_comparison.v1"
DAY25_SCHEMA = "day25.synthetic.locus_plot.v1"
BASES = "ACGT"
VARIANT_RE = re.compile(r"(chr(?:[1-9]|1[0-9]|2[0-2]|X|Y|M)):([1-9][0-9]*):([ACGT])>([ACGT])\Z")
MOTIF_RE = re.compile(r"[A-Z][A-Z0-9_]{0,31}\Z")
SVG_NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG_NS)


def reject_nonfinite(value):
    raise ValueError(f"Nonstandard JSON number is forbidden: {value}")


def reject_duplicates(pairs):
    obj = {}
    for key, value in pairs:
        if key in obj:
            raise ValueError(f"Duplicate JSON key: {key}")
        obj[key] = value
    return obj


def load_json(path, *, max_bytes=200_000):
    """Load only small, local JSON; do not read pickle or remote URLs."""
    path = Path(path)
    if not path.is_file() or path.stat().st_size > max_bytes:
        raise ValueError("Local JSON file missing or exceeds size limit")
    try:
        return json.loads(path.read_text(encoding="utf-8"),
                          parse_constant=reject_nonfinite,
                          object_pairs_hook=reject_duplicates)
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError("Invalid UTF-8 JSON") from exc


def file_sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(65536), b""):
            h.update(block)
    return h.hexdigest()


def require_exact(record, key, expected):
    if key not in record or type(record[key]) is not type(expected) or record[key] != expected:
        raise ValueError(f"Input contract mismatch: {key}")


def validate_fixture(fixture):
    """Validate toy sequence, 1-based SNV and PPM probability rows."""
    if not isinstance(fixture, dict):
        raise ValueError("Fixture must be a JSON object")
    expected = {
        "schema_version": INPUT_SCHEMA, "evidence_class": EVIDENCE_CLASS,
        "mode": MODE, "assembly": ASSEMBLY, "scan_strand": "+",
        "model_prediction": False, "official_motif_atlas": False,
        "clinical_interpretation": False, "real_api_calls_this_run": 0,
    }
    for k, v in expected.items():
        require_exact(fixture, k, v)

    chrom = fixture.get("chrom")
    start = fixture.get("start_1based")
    sequence = fixture.get("ref_sequence")
    variant = fixture.get("variant_id")
    anchor = fixture.get("anchor_start_1based")
    if (not isinstance(chrom, str) or type(start) is not int or start < 1
            or not isinstance(sequence, str) or not 7 <= len(sequence) <= 21
            or any(b not in BASES for b in sequence) or type(anchor) is not int):
        raise ValueError("Invalid toy chrom/start/reference sequence/anchor")
    match = VARIANT_RE.fullmatch(variant) if isinstance(variant, str) else None
    if not match or match.group(1) != chrom or match.group(3) == match.group(4):
        raise ValueError("Invalid synthetic single-base variant")
    position, ref, alt = int(match.group(2)), match.group(3), match.group(4)
    if not start <= position < start + len(sequence):
        raise ValueError("Variant is outside synthetic sequence")
    if sequence[position - start] != ref:
        raise ValueError("REF does not match this *toy* sequence")

    high = fixture.get("preferred_probability")
    low = fixture.get("other_probability")
    if (type(high) not in (float, int) or type(low) not in (float, int)
            or not math.isfinite(high) or not math.isfinite(low)
            or not 0 < low < high < 1 or not math.isclose(high + 3 * low, 1.0, abs_tol=1e-9)):
        raise ValueError("Artificial PPM row must sum to 1 with higher preferred base")
    motifs = fixture.get("motifs")
    if not isinstance(motifs, list) or not 1 <= len(motifs) <= 2:
        raise ValueError("Expected 1-2 synthetic motifs")
    motif_ids = set()
    motif_lengths = set()
    for motif in motifs:
        if not isinstance(motif, dict):
            raise ValueError("Motif must be an object")
        motif_id, preferred = motif.get("motif_id"), motif.get("preferred_sequence")
        if (not isinstance(motif_id, str) or not MOTIF_RE.fullmatch(motif_id)
                or motif_id in motif_ids):
            raise ValueError("Invalid or duplicate motif ID")
        if (not isinstance(preferred, str) or not 3 <= len(preferred) <= 12
                or any(b not in BASES for b in preferred)):
            raise ValueError("Invalid motif preferred sequence")
        motif_ids.add(motif_id)
        motif_lengths.add(len(preferred))
    if len(motif_lengths) != 1:
        raise ValueError("Teaching motifs must have identical width")
    width = motif_lengths.pop()
    if width > len(sequence) or not start <= anchor <= start + len(sequence) - width:
        raise ValueError("Anchor is outside synthetic sequence")
    if not anchor <= position < anchor + width:
        raise ValueError("Anchor motif must overlap the selected variant")
    return {"chrom": chrom, "start": start, "ref": sequence,
            "alt": sequence[:position - start] + alt + sequence[position - start + 1:],
            "position": position, "variant": variant, "anchor": anchor,
            "width": width, "motifs": motifs, "high": float(high), "low": float(low)}


def synthetic_ppm(preferred_sequence, high=0.7, low=0.1):
    """PPM-like row probabilities, not MEME log-odds PWM scores."""
    if not isinstance(preferred_sequence, str) or not preferred_sequence or any(c not in BASES for c in preferred_sequence):
        raise ValueError("Invalid preferred DNA sequence")
    if (type(high) not in (float, int) or type(low) not in (float, int)
            or not 0 < low < high < 1 or not math.isclose(high + 3 * low, 1.0, abs_tol=1e-9)):
        raise ValueError("PPM-like probabilities must sum to one")
    return [{base: (float(high) if base == wanted else float(low)) for base in BASES}
            for wanted in preferred_sequence]


def mean_probability_proxy(sequence, ppm):
    """Uncalibrated teaching arithmetic; NOT a binding probability."""
    if not isinstance(sequence, str) or not isinstance(ppm, list) or not ppm or len(sequence) != len(ppm):
        raise ValueError("Sequence and PPM width differ")
    if any(c not in BASES for c in sequence):
        raise ValueError("Sequence must contain A/C/G/T")
    for row in ppm:
        if (not isinstance(row, dict) or set(row) != set(BASES)
                or any(type(x) not in (float, int) or not math.isfinite(x) or x < 0 for x in row.values())
                or not math.isclose(sum(row.values()), 1.0, abs_tol=1e-9)):
            raise ValueError("Invalid PPM probability row")
    return round(sum(row[base] for row, base in zip(ppm, sequence)) / len(sequence), 6)


def scan_all_windows(ctx):
    """Calculate toy mean preference on the positive strand for every window."""
    rows = []
    for motif in ctx["motifs"]:
        ppm = synthetic_ppm(motif["preferred_sequence"], ctx["high"], ctx["low"])
        for i in range(len(ctx["ref"]) - ctx["width"] + 1):
            left = ctx["start"] + i
            ref_word = ctx["ref"][i:i + ctx["width"]]
            alt_word = ctx["alt"][i:i + ctx["width"]]
            ref_score = mean_probability_proxy(ref_word, ppm)
            alt_score = mean_probability_proxy(alt_word, ppm)
            rows.append({"motif_id": motif["motif_id"],
                         "window_start_1based": left,
                         "window_end_1based_inclusive": left + ctx["width"] - 1,
                         "ref_kmer": ref_word, "alt_kmer": alt_word,
                         "ref_ppm_mean": ref_score, "alt_ppm_mean": alt_score,
                         "delta_alt_minus_ref": round(alt_score - ref_score, 6),
                         "variant_in_window": left <= ctx["position"] < left + ctx["width"]})
    return rows


def validate_day25_link(day25, ctx):
    """Link ONLY the variant identifier; never combine artificial score scales."""
    if not isinstance(day25, dict):
        raise ValueError("Day25 link must be a JSON object")
    for key, expected in {
        "schema_version": DAY25_SCHEMA, "evidence_class": EVIDENCE_CLASS,
        "mode": MODE, "assembly": ASSEMBLY, "model_prediction": False,
        "official_atlas_data": False, "clinical_interpretation": False,
        "real_api_calls_this_run": 0,
        "day23_score_combined_with_delta": False,
    }.items():
        require_exact(day25, key, expected)
    if day25.get("variant_id") != ctx["variant"]:
        raise ValueError("Day25 and Day26 variant IDs do not match")
    if day25.get("track_unit") != "UNIT_INTERVAL_0_TO_1_TOY":
        raise ValueError("Day25 track unit is not the teaching fixture")
    return True


def build_report(fixture, *, fixture_sha256, day25=None, day25_sha256=None):
    ctx = validate_fixture(fixture)
    if not re.fullmatch(r"[0-9a-f]{64}", fixture_sha256):
        raise ValueError("Source SHA-256 must be hexadecimal")
    if day25 is not None:
        validate_day25_link(day25, ctx)
        if not isinstance(day25_sha256, str) or not re.fullmatch(r"[0-9a-f]{64}", day25_sha256):
            raise ValueError("Day25 SHA-256 missing")
    window_rows = scan_all_windows(ctx)
    anchor_results = []
    for motif in ctx["motifs"]:
        row = next(r for r in window_rows
                   if r["motif_id"] == motif["motif_id"] and r["window_start_1based"] == ctx["anchor"])
        delta = row["delta_alt_minus_ref"]
        anchor_results.append({
            "motif_id": motif["motif_id"], "preferred_sequence": motif["preferred_sequence"],
            "window_start_1based": ctx["anchor"],
            "window_end_1based_inclusive": ctx["anchor"] + ctx["width"] - 1,
            "ref_kmer": row["ref_kmer"], "alt_kmer": row["alt_kmer"],
            "ref_ppm_mean": row["ref_ppm_mean"], "alt_ppm_mean": row["alt_ppm_mean"],
            "delta_alt_minus_ref": delta,
            "teaching_interpretation": ("HIGHER_TOY_PREFERENCE" if delta > 0 else
                                        "LOWER_TOY_PREFERENCE" if delta < 0 else "NO_TOY_DIFFERENCE"),
        })
    result = {
        "schema_version": OUTPUT_SCHEMA,
        "evidence_class": EVIDENCE_CLASS, "mode": MODE, "assembly": ASSEMBLY,
        "variant_id": ctx["variant"], "chrom": ctx["chrom"],
        "ref_sequence_start_1based": ctx["start"], "ref_sequence": ctx["ref"],
        "alt_sequence": ctx["alt"], "variant_position_1based": ctx["position"],
        "scan_strand": "+", "scanned_reverse_complement": False,
        "scan_window_count_per_motif": len(ctx["ref"]) - ctx["width"] + 1,
        "motif_count": len(ctx["motifs"]), "window_rows_count": len(window_rows),
        "variant_overlapping_window_rows": sum(r["variant_in_window"] for r in window_rows),
        "ppm_definition": {"preferred_probability": ctx["high"],
                           "other_probability": ctx["low"],
                           "row_sum": ctx["high"] + 3 * ctx["low"],
                           "score_definition": "arithmetic_mean_of_per_position_probabilities",
                           "log_odds_pwm_performed": False},
        "anchor_start_1based": ctx["anchor"], "anchor_results": anchor_results,
        "window_rows": window_rows,
        "input_provenance": {"fixture_filename": "day26_synthetic_motif_context.json",
                             "fixture_sha256": fixture_sha256,
                             "day25_summary_sha256": day25_sha256},
        "day25_linked": day25 is not None,
        "day25_link_scope": "VARIANT_ID_ONLY_NOT_NUMERIC_SCORE",
        "day25_delta_combined_with_ppm": False,
        "motif_binding_confirmed": False,
        "motif_sequence_scanned": True,
        "motif_scan_note": "FORWARD_STRAND_ARTIFICIAL_PPM_NOT_REAL_MOTIF_DISCOVERY",
        "official_motif_atlas": False,
        "reference_allele_verified_against_grch38": False,
        "model_prediction": False,
        "mini_gene_assay_performed": False,
        "clinical_interpretation": False,
        "review_status": "REVIEW_REQUIRED",
        "real_api_calls_this_run": 0,
    }
    return result


def draw_svg(report):
    """Generate 1920x1080 embedded-safe SVG using only stdlib XML."""
    if report.get("schema_version") != OUTPUT_SCHEMA or report.get("evidence_class") != EVIDENCE_CLASS:
        raise ValueError("Not a Day26 synthetic PPM report")
    root = ET.Element(f"{{{SVG_NS}}}svg", {"viewBox": "0 0 1920 1080", "width": "1920", "height": "1080",
                                             "role": "img", "aria-label": "Synthetic PPM comparison for one A to G variant"})
    def add(tag, attrs=None, body=None, parent=root):
        node = ET.SubElement(parent, f"{{{SVG_NS}}}{tag}", attrs or {})
        if body is not None:
            node.text = str(body)
        return node
    def text(x, y, value, size=30, color="#303044", weight="400", parent=root):
        return add("text", {"x": str(x), "y": str(y), "font-family": "Arial, sans-serif",
                            "font-size": str(size), "fill": color, "font-weight": weight}, value, parent)
    add("rect", {"x": "0", "y": "0", "width": "1920", "height": "1080", "fill": "#E5DEDA"})
    add("rect", {"x": "62", "y": "58", "width": "1796", "height": "940", "rx": "30", "fill": "#FFFDFC"})
    text(111, 129, "DAY 26  |  REF / ALT  -  SYNTHETIC PPM PREFERENCE", 38, "#283949", "700")
    text(111, 170, "Toy variant chr1:110:A>G  |  motif overlap at chr1:109-113  |  NOT AlphaGenome output", 22, "#7C6976")
    left = 145
    gap = 71
    seq_start = report["ref_sequence_start_1based"]
    # Show 17 positions; variant at 110 is emphasized.
    for label, sequence, y in (("REF", report["ref_sequence"], 268), ("ALT", report["alt_sequence"], 390)):
        text(110, y + 15, label, 28, "#293B48", "700")
        for i, base in enumerate(sequence):
            pos = seq_start + i
            highlight = pos == report["variant_position_1based"]
            x = left + 135 + i * gap
            add("rect", {"x": str(x), "y": str(y - 40), "width": "56", "height": "73", "rx": "12",
                         "fill": "#F7D9A8" if highlight else "#F1EFF1",
                         "stroke": "#BA7B45" if highlight else "#DAD6DE", "stroke-width": "2"})
            text(x + 18, y + 8, base, 30, "#513B39" if highlight else "#304052", "700")
            if label == "REF":
                text(x + 8, y - 50, str(pos), 15, "#927D86")
    text(110, 472, "A -> G at position 110 changes one letter, not an experimentally measured binding event.", 24, "#71596A")
    anchors = report["anchor_results"]
    for i, row in enumerate(anchors):
        x0 = 109 + i * 855
        y0 = 515
        add("rect", {"x": str(x0), "y": str(y0), "width": "805", "height": "348", "rx": "23",
                     "fill": "#F4F6F8", "stroke": "#E1DFE5", "stroke-width": "2"})
        text(x0 + 33, y0 + 52, row["motif_id"] + "  |  consensus " + row["preferred_sequence"], 28, "#243649", "700")
        text(x0 + 33, y0 + 87, f"REF: {row['ref_kmer']}       ALT: {row['alt_kmer']}", 24, "#796874")
        # Horizontal bars have a fixed 0..1 toy score scale.
        maxw = 540
        for label, value, offset, color in (("REF", row["ref_ppm_mean"], 134, "#598EAE"),
                                            ("ALT", row["alt_ppm_mean"], 214, "#D98D8A")):
            text(x0 + 33, y0 + offset + 22, label, 23, "#3B4653", "600")
            add("rect", {"x": str(x0 + 128), "y": str(y0 + offset), "width": str(maxw), "height": "33", "rx": "9", "fill": "#E3E5E9"})
            add("rect", {"x": str(x0 + 128), "y": str(y0 + offset), "width": str(round(value * maxw, 2)),
                         "height": "33", "rx": "9", "fill": color})
            text(x0 + 690, y0 + offset + 26, f"{value:.2f}", 25, "#38495A", "700")
        delta = row["delta_alt_minus_ref"]
        text(x0 + 33, y0 + 320, f"Delta = ALT - REF: {delta:+.2f}  |  toy preference only", 26,
             "#AD5A56" if delta < 0 else "#2D7775", "700")
    text(115, 934, "PPM-like means mean of artificial position probabilities; NOT log-odds PWM or TF binding affinity.", 25, "#655463")
    text(115, 971, "FORWARD STRAND ONLY   |   TOY_NOT_GRCH38   |   CLINICAL INTERPRETATION: NO   |   REAL API CALLS: 0", 21, "#7E6876")
    return ET.tostring(root, encoding="unicode", xml_declaration=False)


def draw_html(report, svg):
    """Static report preview, no JavaScript/CDN/network."""
    def esc(value):
        return escape(str(value), quote=True)
    table = "".join(
        "<tr>" + "".join(f"<td>{esc(row[col])}</td>" for col in (
            "motif_id", "window_start_1based", "window_end_1based_inclusive",
            "ref_kmer", "alt_kmer", "ref_ppm_mean", "alt_ppm_mean",
            "delta_alt_minus_ref", "variant_in_window")) + "</tr>"
        for row in report["window_rows"])
    headings = "".join(f"<th>{esc(n)}</th>" for n in (
        "Motif", "Start 1-based", "End inclusive", "REF", "ALT", "REF toy score", "ALT toy score", "ALT-REF", "Contains SNV"))
    return ("<!doctype html><html lang='zh-Hant'><head><meta charset='utf-8'>"
            "<meta name='viewport' content='width=device-width, initial-scale=1'>"
            "<title>Day26 | Synthetic PPM comparison</title><style>"
            "body{margin:0;background:#E5DEDA;color:#283849;font:16px system-ui,Arial,sans-serif}"
            "main{margin:auto;max-width:1300px;padding:28px}svg{width:100%;height:auto;border-radius:16px}"
            ".notice{border-left:5px solid #b46756;background:#fff;padding:16px;margin:16px 0}"
            ".table{overflow-x:auto;background:#fff;border-radius:14px;padding:16px}"
            "table{border-collapse:collapse;min-width:990px;width:100%}td,th{padding:9px;border-bottom:1px solid #ddd;text-align:left}"
            "th{background:#e9edef}tr:nth-child(even){background:#fafafa}"
            "a{color:#205b7e}</style></head><body><main>"
            f"<h1>Day26｜從序列偏好理解 REF/ALT 差異</h1><p>{esc(report['variant_id'])} | OFFLINE</p>"
            "<div class='notice'><strong>SYNTHETIC EDUCATIONAL DEMO</strong>: 人工序列及 PPM-like 分數，"
            "不是 Google AlphaGenome／Atlas 模型結果、實驗 TF 結合或臨床證據。</div>"
            + svg +
            "<div class='table'><h2>所有窗口（正股、1-based、兩端包含）</h2><table><thead><tr>"
            + headings + "</tr></thead><tbody>" + table + "</tbody></table></div>"
            "<p>一張完整圖不等於有真實生物學效度。Day25 和本日僅連結同一個人工 variant ID，分數不可相加。</p>"
            "<p><a href='day26_motif_report.json'>JSON</a> · <a href='day26_motif_scan.csv'>CSV</a> · "
            "<a href='day26_ppm_preference.svg'>SVG</a></p></main></body></html>")


def write_outputs(fixture_path, output_dir, *, day25_summary_path=None):
    fixture_path = Path(fixture_path)
    fixture = load_json(fixture_path)
    day25 = load_json(day25_summary_path) if day25_summary_path is not None else None
    report = build_report(fixture,
                          fixture_sha256=file_sha256(fixture_path),
                          day25=day25,
                          day25_sha256=file_sha256(day25_summary_path) if day25 is not None else None)
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    filenames = {
        "json_file": outdir / "day26_motif_report.json",
        "csv_file": outdir / "day26_motif_scan.csv",
        "svg_file": outdir / "day26_ppm_preference.svg",
        "html_file": outdir / "day26_motif_view.html",
    }
    filenames["json_file"].write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    columns = ["motif_id", "window_start_1based", "window_end_1based_inclusive", "ref_kmer", "alt_kmer",
               "ref_ppm_mean", "alt_ppm_mean", "delta_alt_minus_ref", "variant_in_window"]
    with filenames["csv_file"].open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(report["window_rows"])
    svg = draw_svg(report)
    filenames["svg_file"].write_text(svg + "\n", encoding="utf-8")
    filenames["html_file"].write_text(draw_html(report, svg), encoding="utf-8")
    return {
        "status": "PASS", "mode": MODE, "variant_id": report["variant_id"],
        "motif_count": report["motif_count"], "window_rows_count": report["window_rows_count"],
        "variant_overlapping_window_rows": report["variant_overlapping_window_rows"],
        "anchor_results": [{"motif_id": row["motif_id"], "ref_ppm_mean": row["ref_ppm_mean"],
                            "alt_ppm_mean": row["alt_ppm_mean"], "delta_alt_minus_ref": row["delta_alt_minus_ref"]}
                           for row in report["anchor_results"]],
        "day25_linked": report["day25_linked"], "model_prediction": False,
        "real_api_calls_this_run": 0,
        **{name: str(path) for name, path in filenames.items()},
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description="Day 26 synthetic PPM-like motif lesson, offline only")
    parser.add_argument("--mode", required=True, choices=["offline"])
    parser.add_argument("--fixture", default="day26_synthetic_motif_context.json")
    parser.add_argument("--day25-summary", default=None,
                        help="Optional verified Day25 SYNTHETIC summary; link variant ID only")
    parser.add_argument("--output-dir", default="day26_output")
    args = parser.parse_args(argv)
    result = write_outputs(args.fixture, args.output_dir, day25_summary_path=args.day25_summary)
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))
    return result


if __name__ == "__main__":
    main()
