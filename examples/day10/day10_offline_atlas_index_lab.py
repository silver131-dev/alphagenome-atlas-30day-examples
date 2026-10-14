"""Day 10 — offline, synthetic coordinate-index lesson.

This is NOT Google DeepMind's 1 PB Atlas infrastructure or an AlphaGenome SDK
client. It never accesses a network or reads credentials. Python 3.12+ stdlib.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import html
import io
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from xml.etree import ElementTree as ET

DATA_KIND = "SYNTHETIC_EDUCATIONAL_DEMO"
ASSEMBLY = "TOY_NOT_GRCH38"
CHROMS = ("chrToy", "chrDemo")
BIN_WIDTH = 10
BASES = frozenset("ACGT")
REGION_PATTERN = re.compile(r"^(chrToy|chrDemo):([1-9][0-9]*)-([1-9][0-9]*)$")
VARIANT_PATTERN = re.compile(r"^(chrToy|chrDemo):([1-9][0-9]*):([ACGT])>([ACGT])$")
SVG_NS = "http://www.w3.org/2000/svg"
INDEX_FILE = "day10_catalog_index.json"
REPORT_FILE = "day10_query_report.json"
CSV_FILE = "day10_query_rows.csv"
HTML_FILE = "day10_index_dashboard.html"
SVG_FILE = "day10_index_map.svg"
MANIFEST_FILE = "day10_sha256_manifest.json"


def canonical_json(value: object) -> bytes:
    return (json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True,
                       allow_nan=False) + "\n").encode("utf-8")


def sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def sha_file(path: Path) -> str:
    return sha256(path.read_bytes())


def parse_region(region: str) -> tuple[str, int, int]:
    if not isinstance(region, str) or (m := REGION_PATTERN.fullmatch(region)) is None:
        raise ValueError("Region must look like chrToy:100-120 (1-based closed)")
    chrom, lo, hi = m.group(1), int(m.group(2)), int(m.group(3))
    if hi < lo or hi - lo + 1 > 200:
        raise ValueError("Region must be ordered and span at most 200 teaching bases")
    return chrom, lo, hi


def parse_variant(variant: str) -> tuple[str, int, str, str]:
    if not isinstance(variant, str) or (m := VARIANT_PATTERN.fullmatch(variant)) is None:
        raise ValueError("Variant must look like chrToy:110:C>G (1-based SNV)")
    chrom, pos, ref, alt = m.group(1), int(m.group(2)), m.group(3), m.group(4)
    if ref == alt:
        raise ValueError("REF and ALT must differ")
    return chrom, pos, ref, alt


def validate_row(row: object) -> dict:
    if not isinstance(row, dict) or set(row) != {"chrom", "pos_1based", "ref", "alt", "toy_proxy_score"}:
        raise ValueError("Unexpected or missing row fields (private patient fields prohibited)")
    chrom, pos, ref, alt = row["chrom"], row["pos_1based"], row["ref"], row["alt"]
    if chrom not in CHROMS or type(pos) is not int or not 1 <= pos <= 1_000_000:
        raise ValueError("Only chrToy/chrDemo, positive teaching positions are supported")
    if ref not in BASES or alt not in BASES or ref == alt or not isinstance(ref, str) or not isinstance(alt, str):
        raise ValueError("Expected distinct single-letter REF/ALT")
    score = row["toy_proxy_score"]
    if score is not None and (type(score) not in (int, float) or not math.isfinite(score) or not 0 <= score <= 1):
        raise ValueError("Toy proxy score must be finite in [0,1] or null")
    return dict(row)


def validate_catalog(catalog: object) -> list[dict]:
    if not isinstance(catalog, dict) or set(catalog) != {"data_kind", "assembly", "rows", "description"}:
        raise ValueError("Invalid catalog schema")
    if catalog["data_kind"] != DATA_KIND or catalog["assembly"] != ASSEMBLY:
        raise ValueError("Real datasets / real reference genomes not permitted")
    if not isinstance(catalog["description"], str) or not catalog["description"]:
        raise ValueError("Fixture description required")
    rows = catalog["rows"]
    if not isinstance(rows, list) or not 1 <= len(rows) <= 100:
        raise ValueError("Expected 1–100 synthetic variants")
    validated = [validate_row(x) for x in rows]
    keys = [(r["chrom"], r["pos_1based"], r["ref"], r["alt"]) for r in validated]
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicate synthetic variant")
    refs: dict[tuple[str, int], set[str]] = defaultdict(set)
    for row in validated:
        refs[(row["chrom"], row["pos_1based"])].add(row["ref"])
    if any(len(value) != 1 for value in refs.values()):
        raise ValueError("Conflicting REF letters at a single teaching position")
    return sorted(validated, key=lambda r: (r["chrom"], r["pos_1based"], r["ref"], r["alt"]))


def bucket_id(pos: int, width: int = BIN_WIDTH) -> int:
    if type(pos) is not int or pos < 1 or type(width) is not int or width < 1:
        raise ValueError("Positive 1-based position and bin width required")
    return (pos - 1) // width


def load_fixture(path: Path) -> tuple[list[dict], str]:
    raw = path.read_bytes()
    doc = json.loads(raw)
    return validate_catalog(doc), sha256(raw)


def create_index(folder: Path, rows: list[dict], fixture_sha: str) -> dict:
    groups: dict[tuple[str, int], list[dict]] = defaultdict(list)
    for row in rows:
        groups[(row["chrom"], bucket_id(row["pos_1based"]))].append(row)
    (folder / "chunks").mkdir(exist_ok=False)
    buckets: list[dict] = []
    for (chrom, bid), bucket_rows in sorted(groups.items()):
        filename = f"chunks/{chrom}_{bid:06d}.json.gz"
        payload = {"data_kind": DATA_KIND, "assembly": ASSEMBLY,
                   "chrom": chrom, "bin_id": bid, "rows": bucket_rows}
        # fixed mtime means bitwise-repeatable within the same Python/gzip implementation
        compressed = gzip.compress(canonical_json(payload), mtime=0)
        (folder / filename).write_bytes(compressed)
        buckets.append({"chrom": chrom, "bin_id": bid, "path": filename,
                        "record_count": len(bucket_rows), "sha256": sha256(compressed)})
    index = {"data_kind": DATA_KIND, "assembly": ASSEMBLY,
             "index_kind": "USER_CREATED_SYNTHETIC_10_BP_BUCKETS_NOT_ATLAS_INTERNALS",
             "coordinate_system": "1-based positions; interval [start,end] inclusive",
             "bin_width_bp": BIN_WIDTH, "total_records": len(rows),
             "source_fixture_sha256": fixture_sha, "bucket_count": len(buckets),
             "buckets": buckets, "api_calls": 0}
    (folder / INDEX_FILE).write_bytes(canonical_json(index))
    return index


def read_bucket(folder: Path, meta: dict) -> list[dict]:
    name = meta.get("path")
    if (not isinstance(name, str) or not name.startswith("chunks/") or "\\" in name or
        ".." in Path(name).parts or Path(name).is_absolute()):
        raise ValueError("Unsafe chunk path")
    raw = (folder / name).read_bytes()
    if sha256(raw) != meta["sha256"]:
        raise ValueError("Chunk SHA mismatch: " + meta["path"])
    if len(raw) > 200_000:
        raise ValueError("Teaching gzip file is too large")
    body = gzip.decompress(raw)
    if len(body) > 500_000:
        raise ValueError("Teaching uncompressed chunk is too large")
    doc = json.loads(body)
    if (not isinstance(doc, dict) or doc.get("data_kind") != DATA_KIND or
        doc.get("assembly") != ASSEMBLY or doc.get("chrom") != meta["chrom"] or
        doc.get("bin_id") != meta["bin_id"] or not isinstance(doc.get("rows"), list)):
        raise ValueError("Chunk content / metadata mismatch")
    rows = [validate_row(x) for x in doc["rows"]]
    if len(rows) != meta["record_count"] or any(
        row["chrom"] != meta["chrom"] or bucket_id(row["pos_1based"]) != meta["bin_id"] for row in rows
    ):
        raise ValueError("Chunk rows do not match its bin metadata")
    return rows


def query_index(folder: Path, index: dict, region: str, variant: str | None) -> dict:
    chrom, low, high = parse_region(region)
    v = None if variant is None else parse_variant(variant)
    if v is not None and (v[0] != chrom or not low <= v[1] <= high):
        raise ValueError("Variant must be inside the searched region and chromosome")
    min_bin, max_bin = bucket_id(low), bucket_id(high)
    selected = [m for m in index["buckets"] if
                m["chrom"] == chrom and min_bin <= m["bin_id"] <= max_bin]
    accessed = [r for meta in selected for r in read_bucket(folder, meta)]
    matches = sorted((r for r in accessed if low <= r["pos_1based"] <= high),
                     key=lambda r: (r["pos_1based"], r["ref"], r["alt"]))
    v_matches = [r for r in matches if v is not None and
                 (r["chrom"], r["pos_1based"], r["ref"], r["alt"]) == v]
    return {
        "status": "PASS", "mode": "OFFLINE_ONLY", "evidence_class": DATA_KIND,
        "assembly": ASSEMBLY,
        "query_region_1based_closed": region,
        "query_region_0based_half_open": {"chrom": chrom, "start": low - 1, "end": high},
        "query_variant_1based": variant,
        "total_catalog_records": index["total_records"],
        "total_index_buckets": index["bucket_count"],
        "buckets_read": len(selected), "bucket_ids_read": [m["bin_id"] for m in selected],
        "rows_read_from_buckets": len(accessed), "region_hit_count": len(matches),
        "variant_hit_count": len(v_matches),
        "missing_proxy_count": sum(r["toy_proxy_score"] is None for r in matches),
        "rows": [dict(r, variant_id=f'{r["chrom"]}:{r["pos_1based"]}:{r["ref"]}>{r["alt"]}',
                      score_status="MISSING_SCORE" if r["toy_proxy_score"] is None else "SYNTHETIC_SCORE")
                 for r in matches],
        "variant_hits": [f'{r["chrom"]}:{r["pos_1based"]}:{r["ref"]}>{r["alt"]}' for r in v_matches],
        "full_scan_comparison_is_a_row_count_not_a_speed_benchmark": True,
        "google_atlas_internal_index_verified": False,
        "atlas_query_performed": False, "model_prediction_performed": False,
        "clinical_interpretation": False, "real_api_calls_this_run": 0,
    }


def write_csv(folder: Path, report: dict) -> None:
    with (folder / CSV_FILE).open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["variant_id", "chrom", "pos_1based", "ref", "alt", "toy_proxy_score", "score_status", "evidence_class"])
        for r in report["rows"]:
            writer.writerow([r["variant_id"], r["chrom"], r["pos_1based"], r["ref"], r["alt"],
                             "" if r["toy_proxy_score"] is None else r["toy_proxy_score"],
                             r["score_status"], DATA_KIND])


def render_html(index: dict, report: dict) -> str:
    esc = html.escape
    row_html = "\n".join(
        f'<tr><td>{esc(r["variant_id"])}</td><td>{esc(str(r["toy_proxy_score"])) if r["toy_proxy_score"] is not None else "<em>null</em>"}</td>'
        f'<td>{esc(r["score_status"])}</td></tr>' for r in report["rows"])
    block_html = "\n".join(
        f'<div class="block {"active" if b["bin_id"] in report["bucket_ids_read"] and b["chrom"]==parse_region(report["query_region_1based_closed"])[0] else ""}">'
        f'<strong>{esc(b["chrom"])} · BIN {b["bin_id"]}</strong><small>{b["record_count"]} rows</small></div>'
        for b in index["buckets"])
    return f'''<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Day 10 | Synthetic Atlas Index</title><style>
:root{{--bg:#E5DEDA;--ink:#263848;--rose:#B67B80;--blue:#638F9B;--card:#fffefa;}}
*{{box-sizing:border-box}}body{{background:var(--bg);color:var(--ink);font-family:Arial,"Noto Sans TC",sans-serif;margin:0;line-height:1.55}}
main{{max-width:1100px;margin:32px auto;padding:24px}}header{{background:var(--card);padding:32px;border-radius:18px;border:1px solid #cbbfbb}}
h1{{margin:4px 0;font-size:32px}}.label{{font-size:13px;letter-spacing:.08em;color:#8b5260;font-weight:700}}
.grid{{display:grid;grid-template-columns:1fr 1fr;gap:18px;margin-top:20px}}section{{background:var(--card);border-radius:16px;padding:22px}}
.stat{{font-size:40px;font-weight:700;color:var(--rose)}}.blocks{{display:grid;grid-template-columns:repeat(2,1fr);gap:8px}}
.block{{border:1px solid #c5c5c5;background:#eeeded;padding:12px;border-radius:8px}}.block.active{{background:#d4eae8;border:2px solid var(--blue)}}
.block small{{display:block}}table{{width:100%;border-collapse:collapse;font-size:14px}}th,td{{padding:10px;text-align:left;border-bottom:1px solid #ddd}}
.warning{{background:#fcf0e1;border-left:5px solid #b48065;padding:12px;margin-top:18px}}
@media(max-width:720px){{.grid{{grid-template-columns:1fr}}}}
</style></head><body><main><header><span class="label">OFFLINE · SYNTHETIC · NOT THE GOOGLE INDEX</span>
<h1>Day 10 — 座標索引與 Atlas 查詢</h1><p>1 PB 是 Google 官方發布的 Atlas 規模；這裡只有 {index['total_records']} 筆人工 SNV，索引演算法由本課程自行設計。</p>
<div class="warning">本範例不呼叫 API，不證明 Google 使用 B-Tree、Tabix 或特定分塊策略；數值不能判斷致病性。</div></header>
<div class="grid"><section><h2>Lookup by genomic coordinate</h2><div class="label">{esc(report['query_region_1based_closed'])} → 0-based [{report['query_region_0based_half_open']['start']},{report['query_region_0based_half_open']['end']})</div>
<p class="stat">{report['region_hit_count']} hits</p><p>{report['total_catalog_records']} total rows · {report['rows_read_from_buckets']} rows decoded from {report['buckets_read']}/{report['total_index_buckets']} buckets.</p><p>計算的是讀取筆數，沒有測量時間或效能加速倍數。</p></section>
<section><h2>Index buckets</h2><div class="blocks">{block_html}</div></section>
<section style="grid-column:1/-1"><h2>Retrieved artificial variants</h2><table><thead><tr><th>Variant (1-based)</th><th>Artificial proxy score</th><th>Status</th></tr></thead><tbody>{row_html}</tbody></table>
<p>Exact variant: {esc(str(report['query_variant_1based']))} · {report['variant_hit_count']} match(es).</p></section></div>
<p>Google Atlas API: <code>query_variant()</code> / <code>query_interval()</code> (NOT CALLED). Model API: <code>predict_variant()</code> (NOT CALLED).</p></main></body></html>'''


def render_svg(index: dict, report: dict) -> str:
    ET.register_namespace("", SVG_NS)
    svg = ET.Element(f"{{{SVG_NS}}}svg", {"viewBox": "0 0 1920 1080", "width": "1920", "height": "1080"})
    def rect(x, y, w, h, fill, radius=0, stroke=None):
        attrs={"x":str(x),"y":str(y),"width":str(w),"height":str(h),"fill":fill}
        if radius: attrs["rx"]=str(radius)
        if stroke: attrs["stroke"]=stroke
        ET.SubElement(svg,f"{{{SVG_NS}}}rect",attrs)
    def txt(x,y,label,size=30,fill="#233547",weight=None):
        attrs={"x":str(x),"y":str(y),"fill":fill,"font-size":str(size),"font-family":"Arial, sans-serif"}
        if weight:attrs["font-weight"]=str(weight)
        ET.SubElement(svg,f"{{{SVG_NS}}}text",attrs).text=str(label)
    rect(0,0,1920,1080,"#E5DEDA")
    txt(110,115,"ALPHAGENOME ATLAS  |  DAY 10",28,"#8A5D66","bold")
    txt(110,190,"From 1 PB Atlas to coordinate lookup",54,"#263747","bold")
    txt(110,240,"Toy chunk index  /  OFFLINE SYNTHETIC DATA",27,"#557684")
    rect(1180,88,625,84,"#f8f4f1",18)
    txt(1210,142,"REAL ATLAS: ~9 BILLION SNVs / 1 PB",27,"#8a5d66","bold")
    rect(92,300,540,574,"#FAF8F5",24)
    rect(665,300,1152,574,"#FAF8F5",24)
    txt(130,365,"FULL SCAN — toy baseline",33,"#263747","bold")
    txt(130,460,f'{report["total_catalog_records"]}',88,"#B47C83","bold")
    txt(130,505,"catalog rows to inspect",26,"#496879")
    txt(130,590,"No timing measurement",26,"#72646A")
    txt(130,633,"No performance claim",26,"#72646A")
    txt(130,720,"Small local example != Atlas internals",22,"#8A5D66")
    txt(712,365,"INDEXED WINDOW — toy chunks",33,"#263747","bold")
    txt(712,417,f'Query: {report["query_region_1based_closed"]}',28,"#557684")
    # Show selected and unselected chunks in deterministic grid of 4 across
    for i,b in enumerate(index["buckets"]):
        col=i%4;row=i//4;x=710+col*266;y=465+row*130
        on=b["bin_id"] in report["bucket_ids_read"] and b["chrom"]==parse_region(report["query_region_1based_closed"])[0]
        rect(x,y,243,110,"#D4EBE4" if on else "#EEE8E5",14,"#8DB4B0" if on else "#DDCCCB")
        txt(x+18,y+36,f'{b["chrom"]}: bin {b["bin_id"]}',22,"#284954")
        txt(x+18,y+76,f'{b["record_count"]} rows' + ('  READ' if on else ''),23,"#287A77" if on else "#7D858B","bold" if on else None)
    txt(715,778,f'{report["buckets_read"]}/{report["total_index_buckets"]} bins read  |  {report["rows_read_from_buckets"]} rows decoded  |  {report["region_hit_count"]} hits',33,"#467F87","bold")
    txt(715,826,"1-based closed interval -> 0-based half-open",24,"#625E67")
    rect(92,914,1725,102,"#EEDDD7",20)
    txt(130,961,"NOT GOOGLE INTERNALS: B-Tree / Tabix / storage latency are unverified.",26,"#804954","bold")
    txt(130,996,"Atlas API lookup != Model API prediction  |  NO LIVE API  |  NO CLINICAL CLAIM",23,"#804954")
    return ET.tostring(svg,encoding="unicode",xml_declaration=True)


def run(fixture: Path, region: str, variant: str|None, dest: Path) -> dict:
    if dest.exists():
        raise ValueError("Refusing to overwrite output directory; choose a new --output-dir")
    rows, fixture_sha=load_fixture(fixture)
    # Validate all request fields BEFORE creating output directory.
    chrom, low, high = parse_region(region)
    if variant is not None:
        vc,vp,_,_=parse_variant(variant)
        if vc != chrom or not low <= vp <= high:
            raise ValueError("Variant must be within region")
    dest.mkdir(parents=True)
    index=create_index(dest,rows,fixture_sha)
    report=query_index(dest,index,region,variant)
    (dest/REPORT_FILE).write_bytes(canonical_json(report))
    write_csv(dest,report)
    (dest/HTML_FILE).write_text(render_html(index,report),encoding="utf-8")
    (dest/SVG_FILE).write_text(render_svg(index,report),encoding="utf-8")
    entries=[]
    for path in sorted(dest.rglob("*")):
        if path.is_file():
            entries.append({"path":path.relative_to(dest).as_posix(), "sha256":sha_file(path),"size_bytes":path.stat().st_size})
    manifest={"data_kind":DATA_KIND,"assembly":ASSEMBLY,"source_fixture_sha256":fixture_sha,
              "files":entries,"file_count":len(entries),"unsigned_integrity_list":True,
              "real_api_calls_this_run":0}
    (dest/MANIFEST_FILE).write_bytes(canonical_json(manifest))
    return {"status":"PASS","mode":"OFFLINE_ONLY","total_catalog_records":index["total_records"],
            "total_index_buckets":index["bucket_count"],"buckets_read":report["buckets_read"],
            "rows_read_from_buckets":report["rows_read_from_buckets"],"region_hit_count":report["region_hit_count"],
            "variant_hit_count":report["variant_hit_count"],"missing_proxy_count":report["missing_proxy_count"],
            "output_artifact_count":len(entries),"output_dir":str(dest),"real_api_calls_this_run":0}


def verify(fixture: Path, dest: Path) -> dict:
    if not dest.is_dir():
        raise ValueError("Output directory not found")
    manifest=json.loads((dest/MANIFEST_FILE).read_text(encoding="utf-8"))
    rows, actual_sha=load_fixture(fixture)
    if (manifest.get("data_kind")!=DATA_KIND or manifest.get("source_fixture_sha256")!=actual_sha or
        manifest.get("assembly")!=ASSEMBLY or manifest.get("unsigned_integrity_list") is not True):
        raise ValueError("Manifest metadata or source fixture SHA mismatch")
    entries=manifest.get("files")
    if not isinstance(entries,list) or len(entries)!=manifest.get("file_count"):
        raise ValueError("Manifest file count invalid")
    available={p.relative_to(dest).as_posix() for p in dest.rglob("*") if p.is_file() and p.name != MANIFEST_FILE}
    registered=[e.get("path") for e in entries]
    if set(registered)!=available or len(set(registered))!=len(entries):
        raise ValueError("Extra/missing files or duplicate manifest path")
    for meta in entries:
        name=meta["path"]
        if name.startswith("/") or ".." in Path(name).parts or "\\" in name:
            raise ValueError("Unsafe file path")
        path=dest/name
        if sha_file(path)!=meta["sha256"] or path.stat().st_size!=meta["size_bytes"]:
            raise ValueError("Manifest checksum mismatch: "+name)
    index=json.loads((dest/INDEX_FILE).read_text(encoding="utf-8"))
    if (index.get("data_kind")!=DATA_KIND or index.get("assembly")!=ASSEMBLY or
        index.get("source_fixture_sha256")!=actual_sha or index.get("total_records")!=len(rows) or
        index.get("bucket_count")!=len(index.get("buckets",[]))):
        raise ValueError("Index metadata mismatch")
    all_rows=[r for b in index["buckets"] for r in read_bucket(dest,b)]
    if sorted(all_rows,key=lambda r:(r["chrom"],r["pos_1based"],r["ref"],r["alt"]))!=rows:
        raise ValueError("Index does not reconstruct original synthetic catalog")
    saved=json.loads((dest/REPORT_FILE).read_text(encoding="utf-8"))
    recomputed=query_index(dest,index,saved["query_region_1based_closed"],saved["query_variant_1based"])
    if saved!=recomputed:
        raise ValueError("Query report differs from recomputed values")
    return {"status":"PASS","checksums":"MATCH","files_verified":len(entries),
            "catalog_records_verified":len(rows),"region_hit_count":saved["region_hit_count"],
            "index_reconstructed":True,"real_api_calls_this_run":0}


def main(argv: list[str]|None=None) -> int:
    parser=argparse.ArgumentParser(description="Day10: toy local Atlas-like index, no network")
    parser.add_argument("command",choices=("run","verify"))
    parser.add_argument("--mode",choices=("offline","live"),default="offline")
    parser.add_argument("--fixture",default="day10_synthetic_atlas_catalog.json")
    parser.add_argument("--region",default="chrToy:100-120")
    parser.add_argument("--variant",default="chrToy:110:C>G")
    parser.add_argument("--output-dir",default="day10_output")
    args=parser.parse_args(argv)
    try:
        if args.mode!="offline":
            raise ValueError("LIVE_API_NOT_AVAILABLE — offline teaching only")
        if args.command=="run": result=run(Path(args.fixture),args.region,args.variant,Path(args.output_dir))
        else:result=verify(Path(args.fixture),Path(args.output_dir))
    except (ValueError, KeyError, OSError, TypeError, json.JSONDecodeError) as e:
        print(json.dumps({"status":"FAIL","error":str(e),"real_api_calls_this_run":0},ensure_ascii=False,sort_keys=True),file=sys.stderr)
        return 2
    print(json.dumps(result,ensure_ascii=False,sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
