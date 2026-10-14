"""Day 04: synthetic AI Studio MVP -> Python data-contract teaching exercise.

No Google AI Studio interaction, no AlphaGenome API, no patient data, no network.
Uses Python 3.12+ standard library only. Do not use with real patient records.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import re
from pathlib import Path
from xml.etree import ElementTree as ET

EVIDENCE = "SYNTHETIC_EDUCATIONAL_DEMO"
SCHEMA = "DAY04_TOY_MVP_V1"
TRACK_IDS = ("SPLICING_PROXY", "ACCESSIBILITY_PROXY")
ARTIFACTS = (
    "day04_mvp_contract.json",
    "day04_locus_tracks.csv",
    "day04_mvp_bridge_map.md",
    "day04_mvp_two_views.svg",
    "day04_mvp_dashboard.html",
)
VARIANT_RX = re.compile(r"^chrToy:([1-9][0-9]*):([ACGT])>([ACGT])$")
REGION_RX = re.compile(r"^chrToy:([1-9][0-9]*)-([1-9][0-9]*)$")
NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", NS)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def read_spec(path: str | Path) -> dict:
    file = Path(path)
    if not file.is_file() or file.stat().st_size > 100_000:
        raise ValueError("Missing or oversized synthetic fixture")
    data = json.loads(file.read_text(encoding="utf-8"))
    validate_spec(data)
    return data


def parse_variant(value: str) -> tuple[int, str, str]:
    if not isinstance(value, str):
        raise ValueError("Variant must be a string")
    match = VARIANT_RX.fullmatch(value)
    if not match or match.group(2) == match.group(3):
        raise ValueError("Only chrToy:position:REF>ALT single-base substitutions are supported")
    return int(match.group(1)), match.group(2), match.group(3)


def parse_region(value: str) -> tuple[int, int]:
    if not isinstance(value, str):
        raise ValueError("Region must be a string")
    match = REGION_RX.fullmatch(value)
    if not match:
        raise ValueError("Expected chrToy:start-end region")
    start, end = map(int, match.groups())
    if start > end or end - start > 100_000:
        raise ValueError("Invalid or oversized toy interval")
    return start, end


def validate_spec(data: object) -> None:
    if not isinstance(data, dict) or set(data) != {
        "schema_version", "evidence_class", "assembly", "variant_query", "region_query",
        "biosample_label", "variant_view", "locus_view"
    }:
        raise ValueError("Unexpected fixture schema")
    if (data["schema_version"], data["evidence_class"], data["assembly"]) != (
        SCHEMA, EVIDENCE, "TOY_NOT_GRCH38"
    ):
        raise ValueError("This lesson requires explicitly synthetic data")
    if data["biosample_label"] != "SYNTHETIC_NEURON_LIKE_CONTEXT":
        raise ValueError("Biosample must be an unmistakably synthetic label")
    pos, ref, alt = parse_variant(data["variant_query"])
    start, end = parse_region(data["region_query"])
    if not start <= pos <= end:
        raise ValueError("Toy variant must lie within the toy interval")
    vv = data["variant_view"]
    if not isinstance(vv, dict) or set(vv) != {
        "avi_phred", "score_status", "feature_importance_available"
    } or vv["avi_phred"] is not None or vv["score_status"] != "NOT_QUERIED" or vv["feature_importance_available"] is not False:
        raise ValueError("Artificial MVP cannot contain a fake official AVI score")
    lv = data["locus_view"]
    if not isinstance(lv, dict) or set(lv) != {"positions_1based", "tracks"}:
        raise ValueError("Invalid locus view structure")
    points = lv["positions_1based"]
    if not isinstance(points, list) or not 3 <= len(points) <= 100 or any(
        type(p) is not int or not start <= p <= end for p in points
    ) or sorted(set(points)) != points or pos not in points:
        raise ValueError("Invalid, out-of-order or missing toy track coordinates")
    tracks = lv["tracks"]
    if not isinstance(tracks, list) or len(tracks) != len(TRACK_IDS) or tuple(
        t.get("track_id") if isinstance(t, dict) else None for t in tracks
    ) != TRACK_IDS:
        raise ValueError("Expected ordered toy splicing/accessibility tracks")
    for track in tracks:
        if set(track) != {"track_id", "values"} or not isinstance(track["values"], list) or len(track["values"]) != len(points):
            raise ValueError("Invalid track shape")
        for value in track["values"]:
            if value is not None and (
                type(value) not in (int, float) or not math.isfinite(value) or not 0 <= value <= 1
            ):
                raise ValueError("Track values must be finite toy values in [0,1] or null")


def build_contract(data: dict, source_sha256: str) -> dict:
    validate_spec(data)
    start, end = parse_region(data["region_query"])
    return {
        "schema_version": SCHEMA,
        "mode": "OFFLINE_ONLY",
        "evidence_class": EVIDENCE,
        "assembly": "TOY_NOT_GRCH38",
        "variant_view": {
            "query": data["variant_query"],
            "avi_phred": None,
            "score_status": "NOT_QUERIED",
            "feature_importance_available": False,
        },
        "locus_view": {
            "query": data["region_query"],
            "region_0based_half_open": {"start": start - 1, "end": end},
            "biosample": data["biosample_label"],
            "positions_1based": list(data["locus_view"]["positions_1based"]),
            "tracks": data["locus_view"]["tracks"],
        },
        "reverse_engineering": {
            "first_prototype": "GOOGLE_AI_STUDIO_DESIGN_TASK_NOT_EXECUTED_HERE",
            "web_stack": "REACT_UI_AND_POSSIBLE_NODE_BACKEND",
            "python_port": "SEPARATE_LOCAL_TEACHING_IMPLEMENTATION",
            "real_atlas_sdk_compatible": "NOT_TESTED",
        },
        "source_fixture_sha256": source_sha256,
        "reference_allele_verified": False,
        "real_api_calls_this_run": 0,
        "real_model_predictions": False,
        "clinical_interpretation": False,
        "external_transfer_allowed": False,
        "status": "PASS",
    }


def locus_csv(data: dict) -> str:
    points = data["locus_view"]["positions_1based"]
    rows = ["position_1based,splicing_proxy,accessibility_proxy"]
    a, b = (track["values"] for track in data["locus_view"]["tracks"])
    for x, left, right in zip(points, a, b):
        def fmt(v):
            return "" if v is None else format(v, ".6g")
        rows.append(f"{x},{fmt(left)},{fmt(right)}")
    return "\n".join(rows) + "\n"


def reverse_mapping_md(data: dict) -> str:
    return f"""# Day 04 | Reverse-engineering bridge (offline teaching)

**Evidence:** {EVIDENCE}. No Google AI Studio app was generated by this script.

| What the student sees in AI Studio MVP | Data-contract concept | Python function |
|---|---|---|
| Variant View tab | One `chrToy:position:REF>ALT` query | `parse_variant()` |
| Locus View tab | One 1-based closed interval | `parse_region()` |
| Region display | 0-based half-open for BED-style conversion | `build_contract()` |
| AVI score placeholder | Not queried: `null`, never fabricated | `validate_spec()` / `build_contract()` |
| Two synthetic track layers | Positions, values, `null` vs. true zero | `validate_spec()` / `locus_csv()` |
| Card + genome browser | JSON -> HTML/SVG rendering | `build_html()` / `build_svg()` |
| Download / audit | SHA-256 verifies local byte consistency | `run_lab()` / `verify_lab()` |

## From prototype to real analysis

1. Describe the screen in Google AI Studio Build Mode. Use only fictitious `chrToy` inputs and **no external fetches**.
2. Inspect generated React (and possible Node.js) files under **Code**; optionally export a ZIP. AI Studio does not automatically convert it into Python.
3. Independently implement and test the Python data contract. The included HTML is a local reconstruction, **not** an actual exported AI Studio app.
4. Day 05 compares AlphaFold and AlphaGenome; Day 16 teaches official SDK readiness; Day 17 introduces Atlas query examples. Separate authorization is needed before any live API call.

## Observed contract

- Variant: `{data['variant_query']}` — `avi_phred=null`, `NOT_QUERIED`.
- Region: `{data['region_query']}` — five artificial locus points (see fixture).
- Synthetic biosample: `{data['biosample_label']}`.
- `reference_allele_verified=false`; `real_api_calls_this_run=0`; `real_model_predictions=false`.

The demonstration checks software/data handling only: **not** Google DeepMind architecture, actual variant effects, clinical significance or security certification.
"""


def _svg_elem(parent, tag: str, attrs: dict, txt: str | None = None):
    el = ET.SubElement(parent, f"{{{NS}}}{tag}", {k: str(v) for k, v in attrs.items()})
    if txt is not None:
        el.text = txt
    return el


def build_svg(data: dict) -> str:
    validate_spec(data)
    root = ET.Element(f"{{{NS}}}svg", {"viewBox": "0 0 1920 1080", "width": "1920", "height": "1080", "role": "img", "aria-label": "Synthetically illustrated Variant View and Locus View for Day 04"})
    _svg_elem(root, "rect", {"x": 0, "y": 0, "width": 1920, "height": 1080, "fill": "#E5DEDA"})
    _svg_elem(root, "rect", {"x": 58, "y": 54, "width": 1804, "height": 966, "rx": 34, "fill": "#FDFBF9"})
    _svg_elem(root, "text", {"x": 110, "y": 139, "font-family": "Arial, sans-serif", "font-size": 43, "font-weight": "bold", "fill": "#403D45"}, "Day 04  |  Prototype first, explain with Python")
    _svg_elem(root, "text", {"x": 110, "y": 192, "font-family": "Arial, sans-serif", "font-size": 24, "fill": "#766D79"}, "Google AI Studio MVP  →  Code inspection  →  Python data contract")
    # Two distinct views
    _svg_elem(root, "rect", {"x": 100, "y": 260, "width": 810, "height": 587, "rx": 25, "fill": "#F5EEEF", "stroke": "#DCC8D1", "stroke-width": 2})
    _svg_elem(root, "rect", {"x": 938, "y": 260, "width": 880, "height": 587, "rx": 25, "fill": "#EEF4F6", "stroke": "#C9DCE2", "stroke-width": 2})
    _svg_elem(root, "text", {"x": 144, "y": 337, "font-family": "Arial, sans-serif", "font-size": 39, "font-weight": "bold", "fill": "#514052"}, "Variant View")
    _svg_elem(root, "text", {"x": 985, "y": 337, "font-family": "Arial, sans-serif", "font-size": 39, "font-weight": "bold", "fill": "#315464"}, "Locus View")
    _svg_elem(root, "text", {"x": 144, "y": 400, "font-family": "monospace", "font-size": 28, "fill": "#5E4A54"}, data["variant_query"])
    _svg_elem(root, "text", {"x": 985, "y": 400, "font-family": "monospace", "font-size": 28, "fill": "#315464"}, data["region_query"])
    _svg_elem(root, "text", {"x": 144, "y": 475, "font-family": "Arial, sans-serif", "font-size": 27, "fill": "#766371"}, "AVI PHRED")
    _svg_elem(root, "text", {"x": 144, "y": 540, "font-family": "Arial, sans-serif", "font-size": 48, "font-weight": "bold", "fill": "#AA667B"}, "NOT QUERIED")
    _svg_elem(root, "rect", {"x": 144, "y": 585, "width": 672, "height": 85, "rx": 17, "fill": "#EEDDE3"})
    _svg_elem(root, "text", {"x": 163, "y": 637, "font-family": "Arial, sans-serif", "font-size": 24, "fill": "#6E5361"}, "Score = null; no API or clinical interpretation")
    _svg_elem(root, "text", {"x": 144, "y": 749, "font-family": "Arial, sans-serif", "font-size": 26, "fill": "#6E5361"}, "One REF → ALT at one coordinate")
    # Locus graph
    _svg_elem(root, "text", {"x": 988, "y": 471, "font-family": "Arial, sans-serif", "font-size": 24, "fill": "#587581"}, "SYNTHETIC splicing / accessibility proxy")
    plot_x, plot_y, plot_w, plot_h = 1010, 510, 690, 228
    for line_i in range(4):
        y = plot_y + line_i * (plot_h / 3)
        _svg_elem(root, "line", {"x1": plot_x, "y1": y, "x2": plot_x + plot_w, "y2": y, "stroke": "#D3E1E6", "stroke-width": 2})
    points = data["locus_view"]["positions_1based"]
    colors = ("#A96F85", "#5CA5A1")
    for track, color in zip(data["locus_view"]["tracks"], colors):
        segments = []
        current = []
        for index, val in enumerate(track["values"]):
            x = plot_x + (index / (len(points) - 1)) * plot_w
            if val is None:
                if current:
                    segments.append(current)
                current = []
            else:
                y = plot_y + plot_h - val * plot_h * 1.8
                current.append((x,y))
                _svg_elem(root, "circle", {"cx": x, "cy": y, "r": 7, "fill": color})
        if current: segments.append(current)
        for segment in segments:
            if len(segment)>1:
                _svg_elem(root, "polyline", {"points": " ".join(f"{x:.1f},{y:.1f}" for x,y in segment), "stroke":color, "stroke-width":6, "fill":"none", "stroke-linecap":"round", "stroke-linejoin":"round"})
    for i,pos in enumerate(points):
        x=plot_x+(i/(len(points)-1))*plot_w
        _svg_elem(root,"text",{"x":x,"y":plot_y+plot_h+45,"text-anchor":"middle","font-family":"Arial, sans-serif","font-size":21,"fill":"#657883"},str(pos))
    _svg_elem(root, "text", {"x": 985, "y": 811, "font-family": "Arial, sans-serif", "font-size": 22, "fill": "#587581"}, "Missing = gap; zero = plotted at baseline")
    _svg_elem(root, "rect", {"x":100,"y":885,"width":1718,"height":80,"rx":20,"fill":"#EEE8E7"})
    _svg_elem(root, "text", {"x":125,"y":934,"font-family":"Arial, sans-serif","font-size":27,"font-weight":"bold","fill":"#695D69"}, "SYNTHETIC EDUCATIONAL DEMO   •   MODEL PREDICTIONS: 0   •   API CALLS: 0")
    return ET.tostring(root, encoding="unicode") + "\n"


def build_html(data: dict, report: dict) -> str:
    validate_spec(data)
    sample = html.escape(data["biosample_label"], quote=True)
    variant = html.escape(data["variant_query"], quote=True)
    region = html.escape(data["region_query"], quote=True)
    points = data["locus_view"]["positions_1based"]
    tracks = data["locus_view"]["tracks"]
    rows = "".join(
        "<tr><td>" + str(pos) + "</td>" + "".join(
            "<td>" + ("MISSING" if track["values"][i] is None else f"{track['values'][i]:.3f}") + "</td>" for track in tracks
        ) + "</tr>" for i,pos in enumerate(points)
    )
    return f'''<!doctype html>
<html lang="zh-Hant"><head><meta charset="utf-8" />
<meta http-equiv="Content-Security-Policy" content="default-src 'none'; style-src 'unsafe-inline'; script-src 'unsafe-inline'; base-uri 'none'; form-action 'none'" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>Day 04 | Synthetic two-view MVP</title>
<style>
:root{{color-scheme:light}}*{{box-sizing:border-box}}body{{margin:0;font:16px system-ui,'Noto Sans TC',sans-serif;background:#E5DEDA;color:#383740}}
main{{max-width:1100px;margin:2.5rem auto;background:#fefcfb;border-radius:24px;padding:30px;box-shadow:0 12px 36px #4d3b4320}}
h1{{margin-bottom:4px}}.sub{{color:#6e6470}}.tags{{margin:16px 0;padding:12px 16px;background:#f4eeef;border-radius:12px}}
nav{{display:flex;gap:12px;margin:22px 0}}button{{padding:11px 22px;border-radius:18px;border:1px solid #a4a0a9;background:white;cursor:pointer;font:inherit}}button[aria-selected="true"]{{background:#4e505e;color:white}}
section{{border:1px solid #d6d1d5;border-radius:18px;padding:22px;background:#fbf9f8}}.panel-title{{font-size:1.4rem;font-weight:750}}
.notqueried{{font-size:2.5rem;font-weight:800;color:#a56c82}}.disclaimer{{margin-top:24px;background:#fff3ee;border-left:4px solid #a57e6c;padding:16px}}
table{{width:100%;border-collapse:collapse;margin-top:15px}}td,th{{text-align:left;padding:12px;border-bottom:1px solid #e7e3e6}}
small{{color:#675e69}}code{{font-family:ui-monospace,Consolas,monospace}}
</style></head><body><main>
<h1>Day 04｜先看見 MVP，再拆解 Python</h1>
<p class="sub">本機重建的教學畫面；不是 Google AI Studio 實際匯出結果，也不是官方 AlphaGenome Atlas。</p>
<div class="tags">OFFLINE ONLY ・ TOY_NOT_GRCH38 ・ NO REAL API ・ NO PATIENT DATA</div>
<nav aria-label="View selector"><button id="tab-variant" aria-selected="true" aria-controls="view-variant" type="button">Variant View</button>
<button id="tab-locus" aria-selected="false" aria-controls="view-locus" type="button">Locus View</button></nav>
<section id="view-variant"><p class="panel-title">Variant View</p><p>查詢：<code>{variant}</code></p>
<p>AVI PHRED：<strong class="notqueried">NOT QUERIED</strong></p>
<p><small>feature importance: unavailable; score = null。這不是預測模型的輸出。</small></p>
<p>一個基因座位置 + REF → ALT，適合設計單一變異詳情卡片。</p></section>
<section id="view-locus" hidden><p class="panel-title">Locus View</p><p>區間：<code>{region}</code></p>
<p>人工生物樣本標籤：<code>{sample}</code></p>
<table><thead><tr><th>1-based 位置</th><th>Splicing proxy</th><th>Accessibility proxy</th></tr></thead><tbody>{rows}</tbody></table>
<p><small>MISSING (null) 與真正的 0.000 不同。這些是人工教學數字，不是實驗或 AlphaGenome 分數。</small></p></section>
<div class="disclaimer"><strong>Evidence gate:</strong> {EVIDENCE}. `reference_allele_verified=false`； `real_api_calls_this_run=0`；不作臨床解讀。</div>
</main><script>
(function(){{
 const tabs=[document.getElementById('tab-variant'),document.getElementById('tab-locus')];
 const views=[document.getElementById('view-variant'),document.getElementById('view-locus')];
 function select(index){{tabs.forEach((tab,i)=>{{tab.setAttribute('aria-selected',String(i===index));views[i].hidden=(i!==index)}})}}
 tabs.forEach((tab,index)=>tab.addEventListener('click',()=>select(index)));
}})();
</script></body></html>
'''


def run_lab(fixture_path: str | Path, output_dir: str | Path) -> dict:
    data = read_spec(fixture_path)
    source_hash = sha256_file(Path(fixture_path))
    report = build_contract(data, source_hash)
    folder = Path(output_dir)
    if folder.exists():
        raise FileExistsError("Existing output directory will not be overwritten")
    folder.mkdir(parents=True)
    artifacts = {
        "day04_mvp_contract.json": json.dumps(report, ensure_ascii=False, sort_keys=True, indent=2) + "\n",
        "day04_locus_tracks.csv": locus_csv(data),
        "day04_mvp_bridge_map.md": reverse_mapping_md(data),
        "day04_mvp_two_views.svg": build_svg(data),
        "day04_mvp_dashboard.html": build_html(data, report),
    }
    for name in ARTIFACTS:
        (folder/name).write_text(artifacts[name],encoding="utf-8")
    manifest = {
        "schema_version": SCHEMA,
        "mode": "OFFLINE_ONLY",
        "source_fixture_sha256": source_hash,
        "files": {name:sha256_file(folder/name) for name in ARTIFACTS},
        "real_api_calls_this_run":0,
    }
    (folder/"day04_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    return {"status":"PASS","mode":"OFFLINE_ONLY","variant":data["variant_query"],"region":data["region_query"],"locus_position_count":len(data["locus_view"]["positions_1based"]),"artifact_count":len(ARTIFACTS),"real_api_calls_this_run":0,"official_atlas_query_performed":False}


def verify_lab(fixture_path: str | Path, output_dir: str | Path) -> dict:
    data = read_spec(fixture_path)
    folder = Path(output_dir)
    manifest = json.loads((folder/"day04_manifest.json").read_text(encoding="utf-8"))
    if not isinstance(manifest,dict) or manifest.get("schema_version")!=SCHEMA or manifest.get("mode")!="OFFLINE_ONLY" or manifest.get("real_api_calls_this_run")!=0:
        raise ValueError("Invalid manifest schema/mode")
    if manifest.get("source_fixture_sha256")!=sha256_file(Path(fixture_path)):
        raise ValueError("Fixture SHA-256 mismatch")
    files = manifest.get("files")
    if not isinstance(files,dict) or set(files)!=set(ARTIFACTS):
        raise ValueError("Manifest file list mismatch")
    for name in ARTIFACTS:
        if not isinstance(files[name],str) or files[name] != sha256_file(folder/name):
            raise ValueError(f"Output checksum mismatch: {name}")
    report=json.loads((folder/"day04_mvp_contract.json").read_text(encoding="utf-8"))
    if report!=build_contract(data,sha256_file(Path(fixture_path))):
        raise ValueError("Data contract mismatch")
    if (folder/"day04_locus_tracks.csv").read_text(encoding="utf-8")!=locus_csv(data):
        raise ValueError("CSV source mismatch")
    return {"status":"PASS","checksums":"MATCH","files_verified":len(ARTIFACTS),"source_fixture_verified":True,"real_api_calls_this_run":0}


def main(argv: list[str] | None = None) -> None:
    ap=argparse.ArgumentParser(description="Day 04: offline AI Studio MVP/Python bridge")
    ap.add_argument("command", choices=("run","verify"))
    ap.add_argument("--mode",choices=("offline",),required=True)
    ap.add_argument("--fixture",default="day04_synthetic_mvp.json")
    ap.add_argument("--output-dir",default="day04_output")
    args=ap.parse_args(argv)
    result=(run_lab if args.command=="run" else verify_lab)(args.fixture,args.output_dir)
    print(json.dumps(result,ensure_ascii=False,sort_keys=True))


if __name__=="__main__":
    main()
