"""Day 20: OFFLINE-only synthetic motif/interval coordinate learning lab.

All motif labels are invented BED-like intervals, not motif sequence hits.
No AlphaGenome SDK, network access, patient data, or real API calls.
Python 3.12+ standard library only.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
from html import escape
import json
from pathlib import Path
import re

REGION_PATTERN = re.compile(r"^(chr(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|Toy)):(\d+)-(\d+)$")
VARIANT_PATTERN = re.compile(
    r"^(chr(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|Toy)):(\d+):([ACGT])>([ACGT])$"
)
MOTIF_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$")
MAX_WINDOW_BP = 10_000


@dataclass(frozen=True)
class MotifInterval:
    motif_id: str
    chrom: str
    start_0based: int
    end_0based: int
    evidence_class: str = "SYNTHETIC_FIXTURE"


DEFAULT_MOTIFS = (
    MotifInterval("A", "chr1", 99, 104),
    MotifInterval("B", "chr1", 109, 114),
    MotifInterval("C", "chr1", 120, 126),
    MotifInterval("D", "chr2", 99, 103),
)


def valid_motifs(rows):
    """Reject invalid or non-synthetic intervals instead of silently accepting them."""
    rows = tuple(rows)
    if not rows or len(rows) > 1000:
        raise ValueError("Require 1..1000 synthetic motif intervals")
    ids = set()
    for row in rows:
        if not isinstance(row, MotifInterval):
            raise ValueError("Motif row type must be MotifInterval")
        if not MOTIF_ID_PATTERN.fullmatch(row.motif_id):
            raise ValueError("Motif ID format is invalid")
        if row.motif_id in ids:
            raise ValueError("Duplicate motif ID")
        ids.add(row.motif_id)
        if not REGION_PATTERN.fullmatch(f"{row.chrom}:1-1"):
            raise ValueError("Motif chromosome format is invalid")
        if (type(row.start_0based) is not int or type(row.end_0based) is not int
                or row.start_0based < 0 or row.end_0based <= row.start_0based):
            raise ValueError("Motifs must use nonempty 0-based half-open intervals")
        if row.evidence_class != "SYNTHETIC_FIXTURE":
            raise ValueError("Only expressly synthetic fixtures are accepted")
    return rows


def parse_region_1based(value):
    """1-based closed region string -> chromosome and 0-based [start,end)."""
    match = REGION_PATTERN.fullmatch(value)
    if not match:
        raise ValueError("Region must have form chr1:100-120 (1-based, closed)")
    chrom, raw_start, raw_end = match.groups()
    start, end = int(raw_start), int(raw_end)
    if start < 1 or end < start or (end - start + 1) > MAX_WINDOW_BP:
        raise ValueError("Region must be positive, ordered, and <= 10000 bases")
    return chrom, start - 1, end


def parse_variant_1based(value):
    match = VARIANT_PATTERN.fullmatch(value)
    if not match:
        raise ValueError("Teaching SNV must look like chr1:111:A>G")
    chrom, pos, ref, alt = match.groups()
    position = int(pos)
    if position < 1 or ref == alt:
        raise ValueError("Teaching SNV requires position >=1 and REF != ALT")
    return chrom, position - 1, ref, alt


def overlap_bp(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def load_synthetic_tsv(path):
    """Load user-supplied explicit synthetic fixture, without any genomic lookup."""
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        needed = ["motif_id", "chrom", "start_0based", "end_0based", "evidence_class"]
        if reader.fieldnames != needed:
            raise ValueError("Fixture TSV columns must exactly match the documented header")
        out = []
        for row in reader:
            try:
                out.append(MotifInterval(
                    row["motif_id"], row["chrom"], int(row["start_0based"]),
                    int(row["end_0based"]), row["evidence_class"]
                ))
            except (TypeError, ValueError) as exc:
                raise ValueError("Fixture coordinates must be integers") from exc
    return valid_motifs(out)


def analyze(region="chr1:100-120", variant="chr1:111:A>G", motifs=DEFAULT_MOTIFS):
    motifs = valid_motifs(motifs)
    chrom, query_start, query_end = parse_region_1based(region)
    variant_chrom, variant_zero, ref, alt = parse_variant_1based(variant) if variant is not None else (None, None, None, None)
    variant_in_query = bool(variant_chrom == chrom and query_start <= variant_zero < query_end)
    overlaps = []
    for item in motifs:
        if item.chrom != chrom:
            continue
        n = overlap_bp(query_start, query_end, item.start_0based, item.end_0based)
        if n:
            overlaps.append({**asdict(item), "overlap_bp": n,
                             "contains_variant": bool(variant_in_query and item.start_0based <= variant_zero < item.end_0based)})
    return {
        "evidence_class": "SYNTHETIC_EDUCATIONAL_DEMO",
        "mode": "OFFLINE_ONLY",
        "query_1based_closed": region,
        "query_0based_half_open": {"chrom": chrom, "start": query_start, "end": query_end},
        "variant_1based": variant,
        "variant_in_query": variant_in_query,
        "variant_overlap_ids": [r["motif_id"] for r in overlaps if r["contains_variant"]],
        "overlap_count": len(overlaps),
        "overlap_ids": [r["motif_id"] for r in overlaps],
        "region_overlaps": overlaps,
        "motif_interval_count": len(motifs),
        "motif_sequence_scanned": False,
        "reference_allele_verified": False,
        "official_motif_atlas": False,
        "atlas_query_performed": False,
        "clinical_interpretation": False,
        "real_api_calls_this_run": 0,
    }


def render_html(result):
    """Produce a static offline teaching panel with escaped user-facing strings."""
    rows = []
    for item in result["region_overlaps"]:
        cells = (item["motif_id"], item["chrom"],
                 f'[{item["start_0based"]}, {item["end_0based"]})',
                 str(item["overlap_bp"]), "YES" if item["contains_variant"] else "NO")
        rows.append("<tr>" + "".join(f"<td>{escape(x)}</td>" for x in cells) + "</tr>")
    query = escape(result["query_1based_closed"])
    bed = result["query_0based_half_open"]
    bed_display = escape(f'{bed["chrom"]}:[{bed["start"]},{bed["end"]})')
    variant = escape(result["variant_1based"] or "not provided")
    table = "\n".join(rows) if rows else '<tr><td colspan="5">No overlaps in synthetic fixture</td></tr>'
    return f"""<!doctype html>
<html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Day 20 | Offline motif intervals</title><style>
body{{font-family:system-ui,-apple-system,sans-serif;background:#E5DEDA;color:#27313b;margin:0;padding:2rem}}
main{{background:#fffdfb;border-radius:18px;max-width:900px;margin:auto;padding:2rem;box-shadow:0 8px 32px #0001}}
small{{color:#765f6b}} h1{{font-size:1.55rem}} .chip{{display:inline-block;background:#eaf4f7;padding:.35rem .7rem;border-radius:20px;margin:.2rem}}
table{{width:100%;border-collapse:collapse;margin-top:1.2rem}}th,td{{padding:.7rem;border-bottom:1px solid #ddd;text-align:left}}
th{{background:#f7f2f0}}.warn{{background:#fff2e8;padding:1rem;border-radius:10px;border-left:4px solid #bd7659}}
</style></head><body><main>
<small>SYNTHETIC EDUCATIONAL DEMO · OFFLINE_ONLY · no AlphaGenome API</small>
<h1>Day 20 · 基因座區間探索</h1>
<p><span class="chip">1-based: {query}</span> <span class="chip">0-based BED: {bed_display}</span></p>
<p>教學變異：<b>{variant}</b>；位於人工標籤中的候選：<b>{escape(", ".join(result["variant_overlap_ids"]) or "none")}</b></p>
<p>重疊人工區間共 <b>{result["overlap_count"]}</b> 筆。</p>
<table><thead><tr><th>Motif label</th><th>Chromosome</th><th>BED interval</th><th>Overlap (bp)</th><th>Contains SNV</th></tr></thead><tbody>{table}</tbody></table>
<p class="warn">這些是人工命名的 BED 標籤，不是 DNA motif 序列掃描結果，也不是蛋白質結合、AlphaGenome 官方 Atlas 資料或臨床判讀。</p>
</main></body></html>
"""


def run(output_dir="day20_output", region="chr1:100-120", variant="chr1:111:A>G", motifs=DEFAULT_MOTIFS):
    result = analyze(region=region, variant=variant, motifs=motifs)
    folder = Path(output_dir)
    folder.mkdir(parents=True, exist_ok=True)
    json_file = folder / "day20_interval_report.json"
    html_file = folder / "day20_motif_view.html"
    json_file.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    html_file.write_text(render_html(result), encoding="utf-8")
    return {"status": "PASS", "mode": "OFFLINE_ONLY", "query": region,
            "overlap_count": result["overlap_count"], "overlap_ids": result["overlap_ids"],
            "variant_overlap_ids": result["variant_overlap_ids"],
            "json_file": str(json_file), "html_file": str(html_file),
            "real_api_calls_this_run": 0}


def main(argv=None):
    parser = argparse.ArgumentParser(description="Day 20 offline interval/motif label lesson")
    parser.add_argument("--mode", choices=["offline"], default="offline")
    parser.add_argument("--region", default="chr1:100-120", help="1-based inclusive query")
    parser.add_argument("--variant", default="chr1:111:A>G", help="Synthetic SNV; not REF-verified")
    parser.add_argument("--fixture", help="Optional SYNTHETIC_FIXTURE TSV")
    parser.add_argument("--output-dir", default="day20_output")
    args = parser.parse_args(argv)
    motifs = load_synthetic_tsv(args.fixture) if args.fixture else DEFAULT_MOTIFS
    print(json.dumps(run(args.output_dir, args.region, args.variant, motifs),
                     ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
