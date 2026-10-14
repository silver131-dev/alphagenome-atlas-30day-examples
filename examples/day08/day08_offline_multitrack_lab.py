"""Day 08: local-only single-base synthetic splice / ATAC track lesson.

The data are fictional; values have no calibration or biological provenance.
No AlphaGenome SDK, API keys, networking or patient data are accessed.
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
from xml.etree import ElementTree as ET

EVIDENCE_CLASS = "SYNTHETIC_EDUCATIONAL_DEMO"
ASSEMBLY = "TOY_NOT_GRCH38"
FILES = (
    "day08_multitrack_report.json",
    "day08_basewise_tracks.csv",
    "day08_atac_display_bins.csv",
    "day08_multitrack.svg",
    "day08_multitrack_dashboard.html",
)
SVG_NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG_NS)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def render_json(obj: object) -> str:
    return json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True) + "\n"


def load_fixture(path: str | Path) -> tuple[dict, str]:
    data = Path(path).read_bytes()
    fixture = json.loads(data.decode("utf-8"))
    if not isinstance(fixture, dict):
        raise ValueError("Fixture must be a JSON object")
    if fixture.get("evidence_class") != EVIDENCE_CLASS:
        raise ValueError("Not explicitly synthetic educational data")
    if fixture.get("assembly") != ASSEMBLY:
        raise ValueError("The teaching fixture must not claim a real assembly")
    if fixture.get("model_prediction") is not False:
        raise ValueError("The teaching fixture cannot claim a model prediction")
    if fixture.get("chrom") != "chrToy":
        raise ValueError("Teaching example requires chrToy")
    start, end = fixture.get("start_1based"), fixture.get("end_1based")
    if type(start) is not int or type(end) is not int or start < 1 or end < start:
        raise ValueError("Invalid 1-based closed coordinates")
    length = end - start + 1
    if length < 6 or length > 300 or length % 6 != 0:
        raise ValueError("Synthetic span must be 6–300 bases, multiple of 6")
    sequence = fixture.get("ref_sequence")
    if not isinstance(sequence, str) or len(sequence) != length or re.fullmatch(r"[ACGT]+", sequence) is None:
        raise ValueError("REF sequence does not match the genomic interval")
    variant = fixture.get("variant")
    if not isinstance(variant, dict) or set(variant) != {"chrom", "position_1based", "ref", "alt"}:
        raise ValueError("Variant must provide four fields")
    pos = variant["position_1based"]
    if (variant["chrom"] != "chrToy" or type(pos) is not int or not start <= pos <= end or
        variant["ref"] not in "ACGT" or variant["alt"] not in "ACGT" or
        len(variant["ref"]) != 1 or len(variant["alt"]) != 1 or
        variant["ref"] == variant["alt"] or sequence[pos - start] != variant["ref"]):
        raise ValueError("Variant is not a consistent teaching SNV")
    tracks = fixture.get("tracks")
    if not isinstance(tracks, dict) or set(tracks) != {"splice_donor_plus", "atac_accessibility"}:
        raise ValueError("Exactly two named synthetic tracks are needed")
    for key in tracks:
        track = tracks[key]
        if not isinstance(track, dict) or set(track) != {"ref", "alt", "resolution_bp", "units"}:
            raise ValueError(f"Invalid track contract: {key}")
        if type(track["resolution_bp"]) is not int or track["resolution_bp"] != 1:
            raise ValueError("This lesson's SOURCE tracks must be 1 bp")
        if track["units"] != "ARTIFICIAL_0_TO_1_PROXY":
            raise ValueError("Only synthetic unitless proxy values are allowed")
        for allele in ("ref", "alt"):
            values = track[allele]
            if not isinstance(values, list) or len(values) != length:
                raise ValueError("Track values must span the entire interval")
            for value in values:
                if value is None:
                    if key == "atac_accessibility":
                        raise ValueError("ATAC proxy missing values not allowed in this fixture")
                elif (type(value) not in (int, float) or not math.isfinite(value) or
                      not 0 <= value <= 1):
                    raise ValueError("Track values must be finite proxies within [0,1] or null")
        for ref, alt in zip(track["ref"], track["alt"]):
            if (ref is None) != (alt is None):
                raise ValueError("REF/ALT missingness must align for this teaching fixture")
    if fixture.get("synthetic_biosample") != "NEURON_LIKE_MOCK":
        raise ValueError("Unknown synthetic biosample")
    return fixture, sha256_bytes(data)


def nearest(value: float | None) -> float | None:
    return round(float(value), 6) if value is not None else None


def rows_from_fixture(fixture: dict) -> list[dict]:
    start, end = fixture["start_1based"], fixture["end_1based"]
    sp = fixture["tracks"]["splice_donor_plus"]
    atac = fixture["tracks"]["atac_accessibility"]
    rows = []
    for index, pos in enumerate(range(start, end + 1)):
        sref, salt, aref, aalt = sp["ref"][index], sp["alt"][index], atac["ref"][index], atac["alt"][index]
        rows.append({
            "position_1based": pos,
            "base_ref": fixture["ref_sequence"][index],
            "splice_ref": nearest(sref),
            "splice_alt": nearest(salt),
            "splice_delta": nearest(salt - sref) if sref is not None and salt is not None else None,
            "atac_ref": nearest(aref),
            "atac_alt": nearest(aalt),
            "atac_delta": nearest(aalt - aref),
        })
    return rows


def peak(rows: list[dict], column: str) -> dict:
    valid = [r for r in rows if r[column] is not None]
    if not valid:
        raise ValueError("No observed (synthetic) track values")
    best = max(valid, key=lambda r: abs(r[column]))
    return {"position_1based": best["position_1based"], "delta_alt_minus_ref": best[column]}


def make_display_bins(rows: list[dict], width: int = 6) -> list[dict]:
    if type(width) is not int or width <= 0 or len(rows) % width:
        raise ValueError("Display bin size must evenly divide the source track")
    bins = []
    for i in range(0, len(rows), width):
        chunk = rows[i:i + width]
        ref = sum(r["atac_ref"] for r in chunk) / width
        alt = sum(r["atac_alt"] for r in chunk) / width
        bins.append({
            "start_1based": chunk[0]["position_1based"],
            "end_1based": chunk[-1]["position_1based"],
            "start_0based": chunk[0]["position_1based"] - 1,
            "end_0based_exclusive": chunk[-1]["position_1based"],
            "atac_ref_mean": nearest(ref),
            "atac_alt_mean": nearest(alt),
            "atac_delta_mean": nearest(alt - ref),
            "aggregation": "USER_CREATED_MEAN_FOR_DISPLAY_NOT_MODEL_RESOLUTION",
        })
    return bins


def build_report(fixture: dict, fixture_sha256: str) -> dict:
    rows = rows_from_fixture(fixture)
    bins = make_display_bins(rows)
    v = fixture["variant"]
    missing = [r["position_1based"] for r in rows if r["splice_ref"] is None]
    variant_bin = next(b for b in bins if b["start_1based"] <= v["position_1based"] <= b["end_1based"])
    return {
        "status": "PASS",
        "mode": "OFFLINE_ONLY",
        "evidence_class": EVIDENCE_CLASS,
        "assembly": ASSEMBLY,
        "source_sha256": fixture_sha256,
        "synthetic_biosample": fixture["synthetic_biosample"],
        "variant": f"{v['chrom']}:{v['position_1based']}:{v['ref']}>{v['alt']}",
        "interval_1based_closed": f"{fixture['chrom']}:{fixture['start_1based']}-{fixture['end_1based']}",
        "base_count": len(rows),
        "source_resolution_bp": {"splice_sites": 1, "atac": 1},
        "atac_display_bin_width_bp": 6,
        "atac_display_bin_count": len(bins),
        "splice_max_abs_delta": peak(rows, "splice_delta"),
        "atac_max_abs_delta": peak(rows, "atac_delta"),
        "variant_display_bin_atac_delta_mean": variant_bin["atac_delta_mean"],
        "missing_splice_positions_1based": missing,
        "missing_splice_count": len(missing),
        "basewise_rows": rows,
        "user_display_bins": bins,
        "model_prediction": False,
        "official_atac_binned_to_6bp": False,
        "atlas_query_performed": False,
        "real_api_calls_this_run": 0,
        "clinical_interpretation": False,
    }


def render_csv(items: list[dict], columns: tuple[str, ...]) -> str:
    out = io.StringIO(newline="")
    writer = csv.DictWriter(out, fieldnames=columns, lineterminator="\n")
    writer.writeheader()
    writer.writerows(items)
    return out.getvalue()


def render_svg(report: dict) -> str:
    root = ET.Element(f"{{{SVG_NS}}}svg", {"width": "1920", "height": "1080", "viewBox": "0 0 1920 1080", "role": "img",
                                          "aria-label": "Synthetic splice donor 1bp and ATAC 1bp REF ALT tracks, with user-created 6bp display averages"})

    def el(parent, tag, attrs=None, text=None):
        node = ET.SubElement(parent, f"{{{SVG_NS}}}{tag}", attrs or {})
        if text is not None:
            node.text = str(text)
        return node

    def txt(parent, x, y, text, size=28, color="#2e3440", weight="400"):
        return el(parent, "text", {"x": str(x), "y": str(y), "font-family": "Arial, Helvetica, sans-serif",
                                   "font-size": str(size), "font-weight": weight, "fill": color}, text)

    def rect(parent, x, y, w, h, color, rx=20, stroke=None):
        attrs = {"x":str(x), "y":str(y), "width":str(w), "height":str(h), "rx":str(rx), "fill":color}
        if stroke:
            attrs.update({"stroke":stroke,"stroke-width":"2"})
        return el(parent, "rect", attrs)

    rect(root, 0, 0, 1920, 1080, "#E5DEDA", 0)
    txt(root, 72, 80, "DAY 08  |  SINGLE-BASE MULTI-TRACK ATLAS", 43, "#283d4d", "700")
    txt(root, 72, 120, "SYNTHETIC EDUCATIONAL DEMO  ·  REF / ALT  ·  NOT AN ALPHAGENOME PREDICTION", 23, "#78636c")
    rect(root, 69, 144, 1780, 80, "#faf6f3", 22)
    txt(root, 100, 192, f"Variant: {report['variant']}     Splice +0.70 at 1012     ATAC -0.36 at 1012", 28, "#445565", "700")
    positions = [r["position_1based"] for r in report["basewise_rows"]]
    rows = report["basewise_rows"]
    left, right = 165, 1760
    start, end = positions[0], positions[-1]

    def x_for(pos):
        return left + (pos - start) * (right - left) / (end - start)

    def panel(y, h, label, description, dataname, c1, c2):
        rect(root, 69, y, 1780, h, "#faf7f5", 24)
        txt(root, 100, y+41, label, 27, "#364859", "700")
        txt(root, 650, y+41, description, 18, "#756d75")
        top=y+79;bot=y+h-46
        for value in [0,0.5,1]:
            yp=bot-value*(bot-top)
            el(root,"line",{"x1":str(left),"y1":str(yp),"x2":str(right),"y2":str(yp),"stroke":"#e4dcd9","stroke-width":"2"})
            txt(root, 100, yp+7, f"{value:g}", 19, "#857984")
        for row in rows:
            if (row["position_1based"]-start)%4 == 0 or row["position_1based"]==end:
                xp=x_for(row["position_1based"])
                txt(root, xp-18, bot+29, str(row["position_1based"]), 16, "#877880")
        vp=int(report["variant"].split(":")[1]);vx=x_for(vp)
        el(root,"line",{"x1":str(vx),"y1":str(top),"x2":str(vx),"y2":str(bot),"stroke":"#ba8c57","stroke-width":"3","stroke-dasharray":"6,7"})
        for name, color in [(f"{dataname}_ref",c1),(f"{dataname}_alt",c2)]:
            segments=[]; points=[]
            for row in rows:
                value=row[name]
                if value is None:
                    if points:segments.append(points);points=[]
                    continue
                points.append((x_for(row["position_1based"]),bot-value*(bot-top)))
            if points:segments.append(points)
            for points in segments:
                if len(points)>1:
                    el(root,"polyline",{"points":" ".join(f"{x:.1f},{y:.1f}" for x,y in points),"fill":"none","stroke":color,"stroke-width":"5","stroke-linejoin":"round"})
                for px,py in points:
                    el(root,"circle",{"cx":f"{px:.1f}","cy":f"{py:.1f}","r":"4.1","fill":color})
        txt(root, 1400, y+44, "● REF", 19, c1, "700")
        txt(root, 1530, y+44, "● ALT", 19, c2, "700")

    panel(250, 314, "SPLICE DONOR (1 bp)", "per-base synthetic site proxy  |  gap = missing value", "splice", "#5d86a3", "#c77d83")
    panel(588, 310, "ATAC ACCESSIBILITY (1 bp)", "per-base synthetic chromatin accessibility proxy", "atac", "#5d86a3", "#c77d83")
    rect(root, 69, 922, 1780, 105, "#faf7f5", 22)
    txt(root, 100, 961, "DISPLAY-ONLY MEAN (6 bp)", 24, "#364859", "700")
    labels=[f"{b['start_1based']}–{b['end_1based']}: {b['atac_delta_mean']:+.2f}" for b in report["user_display_bins"]]
    for i,label in enumerate(labels):
        xx=112+i*431
        rect(root, xx, 976, 405, 38, "#ebdddd" if i==1 else "#e5eeed", 9)
        txt(root, xx+14, 1003, label, 21, "#56606b", "700" if i==1 else "400")
    txt(root, 87, 1061, "ATAC source is 1 bp. Grouping into 6 bp here is our own plotting choice, not the official model resolution.", 20, "#6a5f65")
    return ET.tostring(root, encoding="unicode") + "\n"


def render_html(report: dict, svg: str) -> str:
    columns=("position_1based","base_ref","splice_ref","splice_alt","splice_delta","atac_ref","atac_alt","atac_delta")
    head="".join(f"<th>{html.escape(col)}</th>" for col in columns)
    rows=[]
    for row in report["basewise_rows"]:
        tds="".join(f"<td>{html.escape(str(row[col]) if row[col] is not None else 'null')}</td>" for col in columns)
        css=' class="variant"' if row["position_1based"]==int(report["variant"].split(":")[1]) else ""
        rows.append(f"<tr{css}>{tds}</tr>")
    return """<!doctype html>
<html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Day08 | Synthetic Multimodal Tracks</title>
<style>body{background:#E5DEDA;color:#253642;font-family:system-ui,'Noto Sans TC',sans-serif;margin:0 auto;max-width:1240px;padding:18px}h1{font-size:30px}article{background:#faf7f5;padding:18px;border-radius:14px;margin:18px 0}svg{width:100%;height:auto}table{width:100%;border-collapse:collapse;font-size:14px}th,td{text-align:left;padding:8px;border-bottom:1px solid #dfd8d3}.variant{background:#f4e9d8}.warn{font-weight:700;color:#983f53}</style></head><body>
<h1>Day08 | 單鹼基解析度圖譜</h1>
<p class="warn">SYNTHETIC EDUCATIONAL DEMO — NOT AlphaGenome output; NOT clinical evidence. No real API calls.</p>
<article><h2>同一組座標，不同的分子意義</h2><p>兩條來源軌跡都以 1 bp 儲存。剪接位點代理值與 ATAC 染色質代理值不可視為同一單位；下方 6 bp 平均僅是視覺化操作。</p>
""" + svg + "</article><article><h2>每一個位置都可追溯</h2><p>空白表示缺值 null，不代表 0；負的 ALT−REF 只代表人工數值下降。</p><div style='overflow-x:auto'><table><thead><tr>" + head + "</tr></thead><tbody>" + "\n".join(rows) + "</tbody></table></div></article></body></html>\n"


def materialize(fixture: dict, digest: str) -> dict[str, str]:
    report=build_report(fixture,digest)
    return {
        FILES[0]: render_json(report),
        FILES[1]: render_csv(report["basewise_rows"],tuple(report["basewise_rows"][0])),
        FILES[2]: render_csv(report["user_display_bins"],tuple(report["user_display_bins"][0])),
        FILES[3]: render_svg(report),
        FILES[4]: render_html(report,render_svg(report)),
    }


def run(fixture_path: str | Path, output_dir: str | Path) -> dict:
    fixture,digest=load_fixture(fixture_path)
    content=materialize(fixture,digest)
    out=Path(output_dir)
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise FileExistsError("Output directory is not empty. Choose a new output directory")
    out.mkdir(parents=True,exist_ok=True)
    for name, value in content.items():
        (out/name).write_text(value,encoding="utf-8")
    manifest={"input_sha256":digest,"files":{name:sha256_bytes(value.encode("utf-8")) for name,value in sorted(content.items())},
              "evidence_class":EVIDENCE_CLASS,"signed_provenance":False}
    (out/"day08_artifact_manifest.json").write_text(render_json(manifest),encoding="utf-8")
    report=build_report(fixture,digest)
    return {"status":"PASS","mode":"OFFLINE_ONLY","base_count":report["base_count"],
            "splice_max_delta":report["splice_max_abs_delta"]["delta_alt_minus_ref"],
            "atac_max_delta":report["atac_max_abs_delta"]["delta_alt_minus_ref"],
            "atac_6bp_display_mean_delta":report["variant_display_bin_atac_delta_mean"],
            "missing_splice_count":report["missing_splice_count"],"artifact_count":len(content),
            "real_api_calls_this_run":0}


def verify(fixture_path: str | Path, output_dir: str | Path) -> dict:
    fixture,digest=load_fixture(fixture_path)
    out=Path(output_dir)
    manifest=json.loads((out/"day08_artifact_manifest.json").read_text(encoding="utf-8"))
    if (not isinstance(manifest,dict) or manifest.get("input_sha256")!=digest or
        manifest.get("evidence_class")!=EVIDENCE_CLASS or manifest.get("signed_provenance") is not False):
        raise ValueError("Input or provenance mismatch")
    expected=materialize(fixture,digest)
    if not isinstance(manifest.get("files"),dict) or set(manifest["files"])!=set(expected):
        raise ValueError("Output manifest file list mismatch")
    for name,text in expected.items():
        actual=(out/name).read_bytes()
        expected_sha=sha256_bytes(text.encode("utf-8"))
        if actual != text.encode("utf-8") or manifest["files"][name]!=expected_sha:
            raise ValueError(f"Source-replay/sha mismatch: {name}")
    return {"status":"PASS","checksums":"MATCH","files_verified":len(expected),"source_sha256":digest,
            "model_prediction":False,"real_api_calls_this_run":0}


def main(argv: list[str] | None = None) -> int:
    p=argparse.ArgumentParser(description="Day08 offline-only multi-track teaching Lab")
    p.add_argument("command",choices=["run","verify"])
    p.add_argument("--mode",required=True,choices=["offline","live"])
    p.add_argument("--fixture",default="day08_synthetic_multitrack.json")
    p.add_argument("--output-dir",default="day08_output")
    a=p.parse_args(argv)
    if a.mode!="offline":
        p.error("LIVE is intentionally unsupported: no real API calls are available")
    try:
        result=run(a.fixture,a.output_dir) if a.command=="run" else verify(a.fixture,a.output_dir)
    except (ValueError,FileNotFoundError,FileExistsError,KeyError,TypeError,UnicodeError,json.JSONDecodeError) as exc:
        p.error(str(exc))
    print(json.dumps(result,ensure_ascii=False,sort_keys=True))
    return 0


if __name__=="__main__":
    sys.exit(main())
