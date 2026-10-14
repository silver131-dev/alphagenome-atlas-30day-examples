#!/usr/bin/env python3
"""Day 09: synthetic REF/ALT track difference vs coordinate shift artifact.

Standard-library-only educational experiment. NO AlphaGenome or other API calls.
Only fictional chrToy values are accepted. Not a reimplementation of model
shift augmentation, indel alignment, or an experimental/clinical assay.
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
from pathlib import Path
from typing import Any
import xml.etree.ElementTree as ET

KIND = "SYNTHETIC_EDUCATIONAL_DEMO"
ASSEMBLY = "TOY_NOT_GRCH38"
VARIANT_RE = re.compile(r"chrToy:([1-9][0-9]*):([ACGT])>([ACGT])\Z")
ARTIFACT_NAMES = (
    "day09_ref_alt_report.json",
    "day09_effect_delta.csv",
    "day09_shift_alignment.csv",
    "day09_delta_shift_figure.svg",
    "day09_delta_shift_dashboard.html",
    "day09_evidence_notes.md",
)
SVG_NS = "http://www.w3.org/2000/svg"


def canonical_json(value: Any) -> bytes:
    return (json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2,
                       allow_nan=False) + "\n").encode("utf-8")


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_fixture(path: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    record = json.loads(raw)
    if not isinstance(record, dict):
        raise ValueError("Fixture must be a JSON object")
    return record


def _positions(obj: Any, name: str) -> list[int]:
    if not isinstance(obj, list) or len(obj) < 2:
        raise ValueError(f"{name} needs at least two 1-based positions")
    if any(type(x) is not int or x < 1 for x in obj):
        raise ValueError(f"{name} must contain positive integer positions")
    if any(b != a + 1 for a, b in zip(obj, obj[1:])):
        raise ValueError(f"{name} must be strictly consecutive")
    return obj


def _signals(obj: Any, name: str, expected: int) -> list[float]:
    if not isinstance(obj, list) or len(obj) != expected:
        raise ValueError(f"{name} length mismatch")
    if any(type(x) not in (float, int) or not math.isfinite(x) or not 0 <= x <= 1
           for x in obj):
        raise ValueError(f"{name} requires finite synthetic values from 0 to 1")
    return [float(x) for x in obj]


def validated_fixture(record: dict[str, Any]) -> dict[str, Any]:
    if record.get("data_kind") != KIND or record.get("assembly") != ASSEMBLY:
        raise ValueError("Explicit synthetic chrToy fixture required")
    variant = record.get("variant")
    match = VARIANT_RE.fullmatch(variant) if isinstance(variant, str) else None
    if not match or match.group(2) == match.group(3):
        raise ValueError("Only a fictional chrToy 1-based SNV REF>ALT is accepted")
    vp, ref_base, _ = int(match.group(1)), match.group(2), match.group(3)
    coords = _positions(record.get("positions_1based"), "positions_1based")
    seq = record.get("reference_sequence")
    if not isinstance(seq, str) or not re.fullmatch("[ACGT]+", seq) or len(seq) != len(coords):
        raise ValueError("Synthetic reference sequence length/bases invalid")
    if vp not in coords or seq[coords.index(vp)] != ref_base:
        raise ValueError("Variant REF does not match synthetic reference string")
    ref = _signals(record.get("reference_track"), "reference_track", len(coords))
    alt = _signals(record.get("alternate_track"), "alternate_track", len(coords))
    sw = record.get("shift_only_windows")
    if not isinstance(sw, dict):
        raise ValueError("Missing separate shift-only windows")
    rpos = _positions(sw.get("ref_positions_1based"), "ref window positions")
    apos = _positions(sw.get("alt_positions_1based"), "alt window positions")
    if len(rpos) != len(apos) or len(rpos) != len(coords) - 1:
        raise ValueError("Shift probe must compare equal-length one-base-shifted views")
    if rpos != coords[:-1] or apos != coords[1:]:
        raise ValueError("Shift probe must shift one coordinate, not change biology")
    rv = _signals(sw.get("ref_values"), "shift ref values", len(rpos))
    av = _signals(sw.get("alt_values"), "shift alt values", len(apos))
    ref_map = dict(zip(coords, ref))
    if any(abs(val - ref_map[p]) > 1e-12 for p, val in zip(rpos, rv)) or any(
            abs(val - ref_map[p]) > 1e-12 for p, val in zip(apos, av)):
        raise ValueError("Shift-only experiment must use the SAME unchanged signal")
    if not isinstance(record.get("label"), str) or not record["label"].strip():
        raise ValueError("Synthetic label required")
    return record


def build_report(record: dict[str, Any], fixture_sha256: str) -> dict[str, Any]:
    record = validated_fixture(record)
    coords = record["positions_1based"]
    ref = record["reference_track"]
    alt = record["alternate_track"]
    delta = [round(b - a, 6) for a, b in zip(ref, alt)]
    idx = max(range(len(delta)), key=lambda i: abs(delta[i]))
    sw = record["shift_only_windows"]
    rp, ap = sw["ref_positions_1based"], sw["alt_positions_1based"]
    rv, av = sw["ref_values"], sw["alt_values"]
    naive = [round(a - r, 6) for r, a in zip(rv, av)]
    rm, am = dict(zip(rp, rv)), dict(zip(ap, av))
    common = sorted(rm.keys() & am.keys())
    aligned = [round(am[p] - rm[p], 6) for p in common]
    return {
        "evidence_class": KIND,
        "assembly": ASSEMBLY,
        "variant": record["variant"],
        "input_provenance": {"fixture_name": "day09_synthetic_refalt_shift.json",
                              "fixture_sha256": fixture_sha256},
        "coordinate_system": "1-based closed genomic positions",
        "unit": "synthetic proxy 0..1, NOT a probability or experimental measurement",
        "effect_comparison": {
            "point_count": len(coords),
            "delta_definition": "ALT_MINUS_REF at identical genomic positions",
            "largest_absolute_delta_position_1based": coords[idx],
            "largest_absolute_delta": delta[idx],
            "rows": [
                {"position_1based": p, "ref_proxy": r, "alt_proxy": a, "delta": d}
                for p, r, a, d in zip(coords, ref, alt, delta)
            ],
        },
        "shift_only_probe": {
            "biological_signal_changed": False,
            "ref_window_start_1based": rp[0],
            "alt_window_start_1based": ap[0],
            "naive_comparison": "ALT[index] minus REF[index]; WRONG when windows differ",
            "naive_max_absolute_delta": max(abs(v) for v in naive),
            "naive_rows": [
                {"ref_position_1based": r, "alt_position_1based": a,
                 "naive_delta": d} for r, a, d in zip(rp, ap, naive)
            ],
            "coordinate_aligned_shared_position_count": len(common),
            "aligned_max_absolute_delta": max(abs(v) for v in aligned),
            "aligned_rows": [{"position_1based": p, "aligned_delta": d}
                             for p, d in zip(common, aligned)],
            "unmatched_ref_positions": sorted(set(rp) - set(ap)),
            "unmatched_alt_positions": sorted(set(ap) - set(rp)),
        },
        "scientific_limits": [
            "Synthetic signals do not establish DNA effects or pathogenicity.",
            "Shift-only probe is coordinate post-processing, NOT AlphaGenome Shift Augmentation.",
            "No indel alignment, liftover, experimental measurement, or clinical inference.",
            "Output window and assay / ontology must match before real track comparison.",
        ],
        "model_prediction": False,
        "genome_reference_verified": False,
        "clinical_interpretation": False,
        "real_api_calls_this_run": 0,
    }


def csv_text(headers: tuple[str, ...], rows: list[dict[str, Any]]) -> str:
    f = io.StringIO(newline="")
    writer = csv.DictWriter(f, fieldnames=headers, lineterminator="\n", extrasaction="raise")
    writer.writeheader()
    writer.writerows({k: r.get(k, "") for k in headers} for r in rows)
    return f.getvalue()


def svg_text(report: dict[str, Any]) -> str:
    effect = report["effect_comparison"]
    shift = report["shift_only_probe"]
    rows = effect["rows"]
    positions = [r["position_1based"] for r in rows]
    start, end = positions[0], positions[-1]
    width, height = 1920, 1080
    # Inline SVG is made only from validated synthetic positions and numbers.
    def e(parent: ET.Element, tag: str, attrs: dict[str, str], label: str | None = None) -> ET.Element:
        node = ET.SubElement(parent, f"{{{SVG_NS}}}{tag}", attrs)
        if label is not None:
            node.text = label
        return node

    root = ET.Element(f"{{{SVG_NS}}}svg", {
        "xmlns": SVG_NS, "width": str(width), "height": str(height),
        "viewBox": f"0 0 {width} {height}", "role": "img",
        "aria-label": "Synthetic REF ALT comparison and false delta from a shifted coordinate window",
    })
    e(root, "rect", {"width": "1920", "height": "1080", "fill": "#E5DEDA"})
    e(root, "rect", {"x": "45", "y": "42", "width": "1830", "height": "990",
                     "rx": "30", "fill": "#FFF9F7", "stroke": "#D2C0BD", "stroke-width": "2"})
    def t(x: int, y: int, txt: str, size: int = 26, fill: str = "#3D3B48", weight: str = "400") -> None:
        e(root, "text", {"x": str(x), "y": str(y), "font-family": "Noto Sans CJK TC, Arial, sans-serif",
                         "font-size": str(size), "font-weight": weight, "fill": fill}, txt)
    t(90, 112, "DAY 09  |  REF/ALT Difference & Coordinate Alignment", 43, "#28273D", "700")
    t(93, 163, "SYNTHETIC EDUCATIONAL DEMO     •     chrToy:104:A>G     •     NO REAL API", 21, "#766770")
    e(root, "rect", {"x": "80", "y": "198", "width": "1120", "height": "690", "rx": "24", "fill": "#FFFFFF", "stroke": "#E8DADA"})
    e(root, "rect", {"x": "1230", "y": "198", "width": "600", "height": "690", "rx": "24", "fill": "#F7F3FD", "stroke": "#DBD2E5"})
    t(110, 253, "01   同一座標的 REF / ALT 差分", 31, "#403D5B", "700")
    t(110, 293, "ALT - REF, per genomic position", 22, "#747184")
    t(1260, 253, "02   視窗平移會造出假差異", 29, "#403D5B", "700")
    t(1260, 293, "Same signal; ALT view starts +1 bp", 22, "#747184")

    lx, rx, top, bottom = 162, 1150, 345, 645
    def xx(pos: int) -> float: return lx + (pos - start) * (rx - lx) / (end - start)
    def yy(v: float) -> float: return bottom - v * (bottom - top)
    for level in (0.0, .25, .5, .75, 1.0):
        y = yy(level)
        e(root, "line", {"x1": str(lx), "x2": str(rx), "y1": f"{y:.2f}", "y2": f"{y:.2f}",
                         "stroke": "#EEE7EA", "stroke-width": "2"})
        t(105, int(y)+7, f"{level:.2f}", 18, "#77717D")
    for pos in positions:
        x=xx(pos)
        e(root, "line", {"x1": f"{x:.2f}", "x2": f"{x:.2f}", "y1": str(top), "y2": "826",
                         "stroke": "#F1E9EC", "stroke-width": "1"})
        t(int(x)-20, 864, str(pos), 21, "#61596A")
    variant_pos = int(report["variant"].split(":")[1])
    e(root, "line", {"x1": f"{xx(variant_pos):.2f}", "x2": f"{xx(variant_pos):.2f}",
                     "y1": str(top), "y2": "825", "stroke": "#B78AB2", "stroke-width": "3",
                     "stroke-dasharray": "8,6"})

    for key, color, label in (("ref_proxy", "#4B789E", "REF"), ("alt_proxy", "#C77D88", "ALT")):
        pts=" ".join(f"{xx(r['position_1based']):.2f},{yy(r[key]):.2f}" for r in rows)
        e(root, "polyline", {"points": pts, "fill": "none", "stroke": color,
                              "stroke-width": "6", "stroke-linecap": "round", "stroke-linejoin": "round",
                              "data-track": label})
        for r in rows:
            e(root, "circle", {"cx": f"{xx(r['position_1based']):.2f}", "cy": f"{yy(r[key]):.2f}",
                               "r": "7", "fill": color, "data-track": label})
    t(880, 332, "REF", 24, "#4B789E", "700")
    t(990, 332, "ALT", 24, "#C77D88", "700")
    t(110, 696, "逐點 Delta (ALT - REF)", 22, "#55495D", "700")
    dybase, dh=766, 73
    e(root, "line", {"x1": str(lx), "x2": str(rx), "y1": str(dybase), "y2": str(dybase), "stroke": "#98919A", "stroke-width": "2"})
    for r in rows:
        d=r["delta"];x=xx(r["position_1based"]);yy2=dybase-d*160
        e(root,"line",{"x1": f"{x:.2f}", "x2": f"{x:.2f}", "y1": str(dybase),
                        "y2": f"{yy2:.2f}", "stroke": "#D18B9F" if d < 0 else "#649FA2",
                        "stroke-width":"18", "stroke-linecap":"round"})
    t(110, 819, "最大絕對差分：104 位點 = −0.30 (人工數值)", 22, "#9C6077", "700")

    naive=shift["naive_rows"]
    points=[r["naive_delta"] for r in naive]
    t(1260, 358, f"按陣列索引相減   max |Delta| = {shift['naive_max_absolute_delta']:.2f}", 23, "#A66177", "700")
    ax0, ax1=1270,1792
    mid=530
    e(root,"line",{"x1":str(ax0),"x2":str(ax1),"y1":str(mid),"y2":str(mid),"stroke":"#A49EAC","stroke-width":"2"})
    bw=(ax1-ax0)/len(points)
    for i,d in enumerate(points):
        xx0=ax0+i*bw+15; endy=mid-d*230
        e(root,"rect",{"x":f"{xx0:.1f}","y":f"{min(mid,endy):.1f}","width":f"{bw-30:.1f}",
                       "height":f"{max(2,abs(endy-mid)):.1f}","rx":"4","fill":"#C77D88" if d<0 else "#649FA2"})
        t(int(xx0), 597, str(i), 16, "#746D80")
    t(1260, 645, "同一個生物訊號，只是視窗起點不同！", 21, "#686074")
    e(root,"rect",{"x":"1252","y":"675","width":"558","height":"147","rx":"20","fill":"#E4F3EF"})
    t(1275, 717, "依共同基因組座標對齊", 25, "#27675F", "700")
    t(1275, 763, f"max |Delta| = {shift['aligned_max_absolute_delta']:.1f}   •   shared = {shift['coordinate_aligned_shared_position_count']} bp", 29, "#26785E", "700")
    t(90, 944, "重點：比對 REF/ALT 時，位置、模態、組織與輸出單位都必須一致。", 26, "#4C4657", "700")
    t(92, 989, "這是座標後處理的人工教學，不是官方 Shift Augmentation、Indel 對齊或模型預測。", 21, "#826F7D")
    return ET.tostring(root, encoding="unicode") + "\n"


def html_text(report: dict[str, Any], svg: str) -> str:
    erows = report["effect_comparison"]["rows"]
    nrows = report["shift_only_probe"]["naive_rows"]
    def tr(cells: list[Any]) -> str:
        return "<tr>"+"".join(f"<td>{html.escape(str(c))}</td>" for c in cells)+"</tr>"
    left="\n".join(tr([r["position_1based"],r["ref_proxy"],r["alt_proxy"],r["delta"]]) for r in erows)
    right="\n".join(tr([r["ref_position_1based"],r["alt_position_1based"],r["naive_delta"]]) for r in nrows)
    return f'''<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Day 09 — REF/ALT 差分與座標對齊</title>
<style>body{{margin:0;background:#E5DEDA;color:#292735;font:17px/1.7 system-ui,Arial,sans-serif}}main{{max-width:1160px;margin:30px auto;padding:24px;background:#fffbfa;border-radius:20px}}svg{{width:100%;height:auto}}table{{width:100%;border-collapse:collapse;background:white}}th,td{{padding:9px 12px;border-bottom:1px solid #eadedb;text-align:left}}.columns{{display:grid;grid-template-columns:1fr 1fr;gap:20px}}.note{{background:#f2e9ee;border-left:5px solid #9b7899;padding:18px}}@media(max-width:780px){{.columns{{grid-template-columns:1fr}}}}</style></head><body><main>
<h1>Day 09｜REF/ALT 差分與視窗位移</h1><p><strong>SYNTHETIC ONLY / NO API</strong> — chrToy:104:A&gt;G</p>{svg}
<p class="note">真正的人工 REF/ALT 訊號差異在位置 104 為 −0.30；同一訊號只改了視窗起點，若錯誤按索引相減會出現非零假差異，正確依座標對齊則歸零。</p>
<div class="columns"><section><h2>同一座標 REF/ALT</h2><table><thead><tr><th>1-based pos</th><th>REF</th><th>ALT</th><th>Delta</th></tr></thead><tbody>{left}</tbody></table></section>
<section><h2>錯誤：不同視窗按索引相減</h2><table><thead><tr><th>REF pos</th><th>ALT pos</th><th>假 Delta</th></tr></thead><tbody>{right}</tbody></table></section></div>
<h2>科學證據界線</h2><p>這不是 AlphaGenome 的 Shift Augmentation 實作；也不涉及插入缺失座標對齊或真實分子預測。依座標對齊是畫圖及差分前的基本資料處理規則。</p>
<p>來源 SHA-256：<code>{html.escape(report['input_provenance']['fixture_sha256'])}</code></p>
</main></body></html>\n'''


def note_text(report: dict[str, Any]) -> str:
    e=report["effect_comparison"];s=report["shift_only_probe"]
    return f"""# Day 09 — Synthetic REF / ALT Evidence Notes

- Evidence: `{KIND}`; assembly `{ASSEMBLY}`; variant `{report['variant']}`.
- Fixture SHA-256: `{report['input_provenance']['fixture_sha256']}`.
- Aligned same-position REF/ALT max delta: `{e['largest_absolute_delta']:+.2f}` at 1-based `{e['largest_absolute_delta_position_1based']}`.
- **No biological change** in the shift-only comparison, but **index-only subtraction** gives max absolute delta `{s['naive_max_absolute_delta']:.2f}`.
- Match by genomic coordinate: `{s['coordinate_aligned_shared_position_count']}` shared positions, max absolute delta `{s['aligned_max_absolute_delta']:.1f}`.
- `real_api_calls_this_run=0`. This is not a model window-augmentation test.
- No claim of genotype pathogenicity, experimentally validated effect, genome build, tissue calibration, or real model signal.

**Interpretation:** in real model outputs, compare the same strand, biosample, assay, units, and genomic coordinate before calculating ALT−REF. The toy single-base offset does not reproduce the actual model's shift augmentation method.\n"""


def run(fixture: Path, output_dir: Path) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError("Output directory is not empty; choose a new directory")
    data_bytes=fixture.read_bytes()
    report=build_report(json.loads(data_bytes),digest(data_bytes))
    outputs={
        ARTIFACT_NAMES[0]:canonical_json(report),
        ARTIFACT_NAMES[1]:csv_text(("position_1based","ref_proxy","alt_proxy","delta"),
                                   report["effect_comparison"]["rows"]).encode("utf-8"),
        ARTIFACT_NAMES[2]:csv_text(("ref_position_1based","alt_position_1based","naive_delta"),
                                   report["shift_only_probe"]["naive_rows"]).encode("utf-8"),
    }
    svg=svg_text(report)
    outputs[ARTIFACT_NAMES[3]]=svg.encode("utf-8")
    outputs[ARTIFACT_NAMES[4]]=html_text(report,svg).encode("utf-8")
    outputs[ARTIFACT_NAMES[5]]=note_text(report).encode("utf-8")
    output_dir.mkdir(parents=True,exist_ok=True)
    for name, data in outputs.items():
        (output_dir/name).write_bytes(data)
    manifest={"evidence_class":KIND,"source_fixture_name":fixture.name,
              "source_fixture_sha256":digest(data_bytes),
              "files":{name:digest(blob) for name,blob in sorted(outputs.items())}}
    (output_dir/"day09_sha256_manifest.json").write_bytes(canonical_json(manifest))
    return {"status":"PASS","mode":"OFFLINE_ONLY","variant":report["variant"],
            "effect_max_delta":report["effect_comparison"]["largest_absolute_delta"],
            "naive_shift_max_abs_delta":report["shift_only_probe"]["naive_max_absolute_delta"],
            "aligned_shift_max_abs_delta":report["shift_only_probe"]["aligned_max_absolute_delta"],
            "artifact_count":len(outputs),"real_api_calls_this_run":0,
            "output_dir":str(output_dir)}


def verify(fixture: Path, report_dir: Path) -> dict[str, Any]:
    data_bytes=fixture.read_bytes()
    record=validated_fixture(json.loads(data_bytes))
    manifest=json.loads((report_dir/"day09_sha256_manifest.json").read_text(encoding="utf-8"))
    if manifest.get("evidence_class")!=KIND or manifest.get("source_fixture_name") != fixture.name:
        raise ValueError("Evidence class/source filename mismatch")
    if manifest.get("source_fixture_sha256") != digest(data_bytes):
        raise ValueError("Source fixture content changed")
    if set(manifest.get("files",{})) != set(ARTIFACT_NAMES):
        raise ValueError("Manifest file list differs from declared artifacts")
    mismatched=[name for name in ARTIFACT_NAMES if not (report_dir/name).is_file() or
                digest((report_dir/name).read_bytes()) != manifest["files"][name]]
    if mismatched:
        raise ValueError("Artifact checksums do not match: "+", ".join(mismatched))
    expected=build_report(record,digest(data_bytes))
    found=json.loads((report_dir/ARTIFACT_NAMES[0]).read_text(encoding="utf-8"))
    if found != expected:
        raise ValueError("Report contract differs from current validated fixture")
    return {"status":"PASS","checksums":"MATCH","artifacts_verified":len(ARTIFACT_NAMES),
            "fixture_verified":True,"real_api_calls_this_run":0}


def main(argv: list[str] | None = None) -> None:
    ap=argparse.ArgumentParser(description="Day 09: local-only synthetic REF/ALT shift contrast")
    ap.add_argument("command",choices=("run","verify"))
    ap.add_argument("--mode",choices=("offline",),default="offline")
    ap.add_argument("--fixture",default="day09_synthetic_refalt_shift.json")
    ap.add_argument("--output-dir",default="day09_output")
    args=ap.parse_args(argv)
    if args.mode != "offline":
        ap.error("Only offline mode is supported")
    result=run(Path(args.fixture),Path(args.output_dir)) if args.command=="run" else verify(
        Path(args.fixture),Path(args.output_dir))
    print(json.dumps(result,ensure_ascii=False,sort_keys=True))

if __name__=="__main__":
    main()
