"""Day 06: synthetic 100-kb versus 1-Mb context-window teaching lab.

Only Python standard library; no AlphaGenome SDK, network, API key,
reference FASTA, chromatin-contact measurements, or patient data.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import io
import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

DATA_KIND = "SYNTHETIC_EDUCATIONAL_DEMO"
GENOME_BUILD = "TOY_NOT_GRCH38"
LOCAL_WIDTH = 100_000  # Arbitrary educational comparison window, not an SDK contract.
LONG_WIDTH = 1_048_576  # Exactly 2**20 bp (official 1-Mb-size sequence length).
ARTIFACTS = (
    "day06_context_report.json",
    "day06_context_features.csv",
    "day06_context_map.svg",
    "day06_context_dashboard.html",
)
MANIFEST_NAME = "day06_sha256_manifest.json"
FEATURE_ID_RE = re.compile(r"^[A-Z][A-Z0-9_]{1,39}$")
SVG_NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG_NS)


def integer(value: Any, name: str, *, lo: int = 1, hi: int | None = None) -> int:
    if type(value) is not int or value < lo or (hi is not None and value > hi):
        raise ValueError(f"Invalid {name}: expected integer in allowed range")
    return value


def validate_fixture(data: Any) -> dict[str, Any]:
    """Reject malformed, ambiguous, or non-synthetic inputs."""
    if not isinstance(data, dict):
        raise ValueError("Fixture must be a JSON object")
    if data.get("evidence_class") != DATA_KIND or data.get("assembly") != GENOME_BUILD:
        raise ValueError("Synthetic fixture and TOY_NOT_GRCH38 labels required")
    if data.get("chromosome") != "chrToy":
        raise ValueError("Only chrToy is supported")
    genome_length = integer(data.get("toy_chromosome_length_bp"), "toy chromosome length", lo=LONG_WIDTH)
    anchor = integer(data.get("anchor_1based"), "anchor", hi=genome_length)
    windows = data.get("windows_bp")
    if not isinstance(windows, dict) or windows != {"local": LOCAL_WIDTH, "expanded": LONG_WIDTH}:
        raise ValueError("Expected 100000 and 1048576 bp comparison windows")
    if data.get("contact_map_measured") is not False or data.get("experimental_validation_performed") is not False:
        raise ValueError("No contact or experimental validation is permitted in this lesson")
    features = data.get("features")
    if not isinstance(features, list) or len(features) < 3 or len(features) > 100:
        raise ValueError("Provide between 3 and 100 synthetic features")
    seen_ids: set[str] = set()
    seen_positions: set[int] = set()
    target_count = 0
    for feature in features:
        if not isinstance(feature, dict):
            raise ValueError("Feature must be an object")
        fid = feature.get("id")
        if not isinstance(fid, str) or not FEATURE_ID_RE.fullmatch(fid) or fid in seen_ids:
            raise ValueError("Invalid or duplicated feature id")
        seen_ids.add(fid)
        pos = integer(feature.get("position_1based"), f"{fid} position", hi=genome_length)
        if pos in seen_positions:
            raise ValueError("Duplicated 1-based feature position")
        seen_positions.add(pos)
        kind = feature.get("kind")
        if kind not in ("TARGET_GENE", "ENHANCER_CANDIDATE"):
            raise ValueError("Unsupported feature type")
        if kind == "TARGET_GENE":
            target_count += 1
            if pos != anchor:
                raise ValueError("The target must be at the anchor position")
        if feature.get("evidence_class") != "SYNTHETIC_FIXTURE":
            raise ValueError("Features must be explicitly synthetic")
    if target_count != 1:
        raise ValueError("Exactly one TARGET_GENE is required")
    hypotheses = data.get("hypothetical_links")
    if not isinstance(hypotheses, list) or len(hypotheses) > 20:
        raise ValueError("Hypotheses must be a list of at most twenty entries")
    for h in hypotheses:
        if not isinstance(h, dict) or h.get("candidate_id") not in seen_ids or h.get("target_id") not in seen_ids:
            raise ValueError("Hypothesis references missing feature")
        f = next(f for f in features if f["id"] == h["candidate_id"])
        t = next(f for f in features if f["id"] == h["target_id"])
        if (f["kind"] != "ENHANCER_CANDIDATE" or t["kind"] != "TARGET_GENE"
                or h.get("status") != "HYPOTHESIS_NOT_MEASURED"):
            raise ValueError("Links must be unmeasured enhancer-to-target hypotheses")
    return data


def centered_window(anchor_1based: int, width: int, chromosome_length: int) -> dict[str, int]:
    """Compute exact-width, clipped 0-based half-open interval on a toy chromosome.

    This is an educational coordinate calculation, not the SDK's resize algorithm.
    """
    chromosome_length = integer(chromosome_length, "chromosome length")
    width = integer(width, "width", hi=chromosome_length)
    anchor_1based = integer(anchor_1based, "anchor", hi=chromosome_length)
    anchor0 = anchor_1based - 1
    start = max(0, min(anchor0 - width // 2, chromosome_length - width))
    end = start + width
    if not (start <= anchor0 < end):
        raise ValueError("Anchor not in computed window")
    return {"start_0based": start, "end_0based_exclusive": end,
            "start_1based_inclusive": start + 1, "end_1based_inclusive": end,
            "width_bp": width}


def in_window(pos_1based: int, window: dict[str, int]) -> bool:
    integer(pos_1based, "position")
    return window["start_0based"] <= pos_1based - 1 < window["end_0based_exclusive"]


def build_report(fixture: dict[str, Any]) -> dict[str, Any]:
    data = validate_fixture(fixture)
    anchor = data["anchor_1based"]
    length = data["toy_chromosome_length_bp"]
    local = centered_window(anchor, LOCAL_WIDTH, length)
    expanded = centered_window(anchor, LONG_WIDTH, length)
    rows = []
    for f in sorted(data["features"], key=lambda x: (x["position_1based"], x["id"])):
        pos = f["position_1based"]
        rows.append({"id": f["id"], "kind": f["kind"], "position_1based": pos,
                     "distance_to_target_bp": abs(pos - anchor),
                     "local_100kb": in_window(pos, local),
                     "expanded_2pow20bp": in_window(pos, expanded),
                     "evidence_class": "SYNTHETIC_FIXTURE"})
    local_ids = [f["id"] for f in rows if f["local_100kb"]]
    expanded_ids = [f["id"] for f in rows if f["expanded_2pow20bp"]]
    newly_seen = [f["id"] for f in rows if f["expanded_2pow20bp"] and not f["local_100kb"]]
    return {
        "status": "PASS", "mode": "OFFLINE_ONLY",
        "evidence_class": DATA_KIND, "assembly": GENOME_BUILD,
        "chromosome": data["chromosome"], "anchor_1based": anchor,
        "local_window": local, "expanded_window": expanded,
        "local_count": len(local_ids), "expanded_count": len(expanded_ids),
        "newly_seen_count": len(newly_seen),
        "local_feature_ids": local_ids, "expanded_feature_ids": expanded_ids,
        "newly_seen_feature_ids": newly_seen,
        "hypothetical_links": data["hypothetical_links"],
        "features": rows,
        "contact_map_measured": False, "chromatin_loop_inferred": False,
        "gene_regulation_established": False, "experiment_performed": False,
        "official_alphagenome_prediction": False, "clinical_interpretation": False,
        "real_api_calls_this_run": 0,
    }


def json_text(obj: Any) -> str:
    return json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n"


def csv_text(report: dict[str, Any]) -> str:
    columns = ("id", "kind", "position_1based", "distance_to_target_bp", "local_100kb",
               "expanded_2pow20bp", "evidence_class")
    buf = io.StringIO(newline="")
    writer = csv.DictWriter(buf, fieldnames=columns, lineterminator="\n")
    writer.writeheader()
    for row in report["features"]:
        writer.writerow({**row, "local_100kb": str(row["local_100kb"]).lower(),
                         "expanded_2pow20bp": str(row["expanded_2pow20bp"]).lower()})
    return buf.getvalue()


def svg_text(report: dict[str, Any]) -> str:
    def el(parent: ET.Element, name: str, attrs: dict[str, str] | None = None, text: str | None = None) -> ET.Element:
        node = ET.SubElement(parent, f"{{{SVG_NS}}}{name}", attrs or {})
        if text is not None:
            node.text = text
        return node

    root = ET.Element(f"{{{SVG_NS}}}svg", {"width":"1920", "height":"1080", "viewBox":"0 0 1920 1080",
                                                 "role":"img", "aria-label":"Synthetic DNA context window diagram comparing 100kb and 1048576bp, no 3D contacts measured"})
    el(root,"rect",{"x":"0","y":"0","width":"1920","height":"1080","fill":"#E5DEDA"})
    el(root,"text",{"x":"92","y":"105","font-family":"Arial,sans-serif","font-size":"52","font-weight":"700","fill":"#24364B"},"DAY 06  /  LONG-RANGE CONTEXT")
    el(root,"text",{"x":"94","y":"152","font-family":"Arial,sans-serif","font-size":"26","fill":"#5A6270"},"What a 100 kb window misses, a 1,048,576 bp window can include")
    el(root,"rect",{"x":"75","y":"190","width":"1770","height":"425","rx":"22","fill":"#F9F7F5","fill-opacity":"0.90"})
    el(root,"rect",{"x":"75","y":"642","width":"1770","height":"226","rx":"22","fill":"#F9F7F5","fill-opacity":"0.90"})
    el(root,"text",{"x":"112","y":"246","font-family":"Arial,sans-serif","font-size":"32","font-weight":"700","fill":"#24364B"},"EXPANDED / 1,048,576 bp (2^20)")
    el(root,"text",{"x":"112","y":"681","font-family":"Arial,sans-serif","font-size":"32","font-weight":"700","fill":"#24364B"},"LOCAL / 100,000 bp (toy comparison)")
    el(root,"text",{"x":"112","y":"284","font-family":"Arial,sans-serif","font-size":"18","fill":"#657586"},"1-based inclusive: "+str(report["expanded_window"]["start_1based_inclusive"])+" - "+str(report["expanded_window"]["end_1based_inclusive"]))
    expanded = report["expanded_window"]
    local = report["local_window"]
    def px(pos:int,window:dict[str,int])->float:
        return 160 + (pos-window["start_1based_inclusive"]) / (window["width_bp"]-1) * 1580
    el(root,"line",{"x1":"160","y1":"418","x2":"1740","y2":"418","stroke":"#8997A7","stroke-width":"7","stroke-linecap":"round"})
    target = next(f for f in report["features"] if f["kind"]=="TARGET_GENE")
    tx=px(target["position_1based"],expanded)
    for row in report["features"]:
        if not row["expanded_2pow20bp"]:continue
        x=px(row["position_1based"],expanded)
        is_target=row["kind"]=="TARGET_GENE"
        color="#355C7D" if is_target else "#B97860"
        el(root,"circle",{"cx":f"{x:.2f}","cy":"418","r":"13" if is_target else "9","fill":color})
        dy=348 if is_target else (487 if row["id"]=="NEAR_ENHANCER" else (365 if row["id"]=="DISTAL_ENHANCER" else 501))
        el(root,"line",{"x1":f"{x:.2f}","y1":"418","x2":f"{x:.2f}","y2":str(dy+(14 if dy<418 else -18)),"stroke":color,"stroke-width":"2"})
        if row["id"]=="DISTAL_ENHANCER":
            el(root,"text",{"x":f"{x-6:.2f}","y":str(dy-9),"font-family":"Arial,sans-serif","font-size":"23","font-weight":"700","fill":color},"DISTAL? (+400 kb)")
        elif row["id"]=="TARGET":
            el(root,"text",{"x":f"{x-60:.2f}","y":str(dy-12),"font-family":"Arial,sans-serif","font-size":"23","font-weight":"700","fill":color},"TARGET")
        elif row["id"]=="NEAR_ENHANCER":
            el(root,"text",{"x":f"{x+8:.2f}","y":str(dy+9),"font-family":"Arial,sans-serif","font-size":"23","fill":color},"NEAR (+10 kb)")
        elif row["id"].startswith("LEFT_") or row["id"].startswith("RIGHT_"):
            el(root,"text",{"x":f"{x+9 if row["id"].startswith("LEFT") else x-147:.2f}","y":"531","font-family":"Arial,sans-serif","font-size":"18","fill":color},row["id"].replace("_", " "))
    # Dashed arc is explicitly hypothetical and cannot be read as a measured 3D contact.
    distal=next(f for f in report["features"] if f["id"]=="DISTAL_ENHANCER")
    dx=px(distal["position_1based"],expanded)
    path=f"M {tx:.2f} 403 Q {(tx+dx)/2:.2f} 260 {dx:.2f} 403"
    el(root,"path",{"d":path,"fill":"none","stroke":"#B89565","stroke-width":"5","stroke-dasharray":"14 11"})
    el(root,"text",{"x":f"{(tx+dx)/2-164:.2f}","y":"315","font-family":"Arial,sans-serif","font-size":"19","font-weight":"700","fill":"#966F3F"},"HYPOTHETICAL LINK ONLY")
    el(root,"text",{"x":"115","y":"575","font-family":"Arial,sans-serif","font-size":"21","fill":"#4E6474"},f"Visible: {report['expanded_count']} elements   |   Added vs 100 kb: {report['newly_seen_count']}   |   Dotted arc: NOT OBSERVED")
    el(root,"line",{"x1":"160","y1":"759","x2":"1740","y2":"759","stroke":"#A5AFB8","stroke-width":"8","stroke-linecap":"round"})
    for row in report["features"]:
        if not row["local_100kb"]:continue
        x=px(row["position_1based"],local)
        is_target=row["kind"]=="TARGET_GENE"
        color="#355C7D" if is_target else "#B97860"
        el(root,"circle",{"cx":f"{x:.2f}","cy":"759","r":"13" if is_target else "10","fill":color})
        el(root,"text",{"x":f"{x-51:.2f}","y":"813","font-family":"Arial,sans-serif","font-size":"21","font-weight":"700","fill":color},"TARGET" if is_target else "NEAR")
    el(root,"text",{"x":"110","y":"846","font-family":"Arial,sans-serif","font-size":"17","fill":"#667585"},"Only the target and nearby candidate fit into this shorter window")
    el(root,"text",{"x":"94","y":"931","font-family":"Arial,sans-serif","font-size":"28","font-weight":"700","fill":"#2D4056"},"SYNTHETIC COORDINATES | NOT A PHYSICAL CONTACT MAP | NO MODEL CALL")
    el(root,"text",{"x":"94","y":"966","font-family":"Arial,sans-serif","font-size":"20","fill":"#657181"},"Increasing sequence context makes distant candidates visible; it does not establish regulatory causality.")
    return ET.tostring(root,encoding="unicode",xml_declaration=False)+"\n"


def html_text(report: dict[str, Any], svg: str) -> str:
    columns=("id", "kind", "position_1based", "distance_to_target_bp", "local_100kb", "expanded_2pow20bp")
    head=''.join('<th>'+html.escape(k)+'</th>' for k in columns)
    body=''.join('<tr>'+''.join('<td>'+html.escape(str(row[k]))+'</td>' for k in columns)+'</tr>' for row in report["features"])
    return ("<!doctype html>\n<html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">"
            "<title>Day 06 synthetic context windows</title><style>body{font:16px system-ui;margin:0 auto;max-width:1250px;padding:24px;background:#E5DEDA;color:#253650}"
            "svg{width:100%;height:auto}table{width:100%;border-collapse:collapse;background:white}th,td{border:1px solid #bbb;padding:8px;text-align:left}"
            "p{line-height:1.6}</style></head><body><h1>Day 06 / Context Window Lab</h1>"
            "<p><strong>SYNTHETIC EDUCATIONAL DATA.</strong> This exercise does not infer 3D contact, regulatory activity, clinical risk, or real AlphaGenome results.</p>"
            +svg+"<h2>Window comparison</h2><p>100 kb contains " +str(report["local_count"])+" records; 1,048,576 bp contains "
            +str(report["expanded_count"])+" records; newly visible: "+str(report["newly_seen_count"])+".</p>"
            +"<table><thead><tr>"+head+"</tr></thead><tbody>"+body+"</tbody></table>"
            "<p>No network request, no patient data, no chromatin contact measurement, no real model prediction. No clinical inference.</p>"
            "</body></html>\n")


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_fixture(path: Path) -> tuple[dict[str, Any], bytes]:
    raw=path.read_bytes()
    try:
        obj=json.loads(raw)
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise ValueError("Invalid UTF-8 JSON fixture") from exc
    return validate_fixture(obj),raw


def make_outputs(fixture: dict[str, Any]) -> tuple[dict[str, Any],dict[str, bytes]]:
    report=build_report(fixture)
    svg=svg_text(report)
    content={
       ARTIFACTS[0]:json_text(report).encode("utf-8"),
       ARTIFACTS[1]:csv_text(report).encode("utf-8"),
       ARTIFACTS[2]:svg.encode("utf-8"),
       ARTIFACTS[3]:html_text(report,svg).encode("utf-8"),
    }
    return report, content


def manifest_text(fixture_path:Path,raw:bytes,files:dict[str,bytes])->str:
    return json_text({"evidence_class":DATA_KIND,"source_filename":fixture_path.name,
                      "source_sha256":sha256(raw),"artifacts_sha256":{name:sha256(b) for name,b in files.items()},
                      "manifest_is_digital_signature":False,"real_api_calls_this_run":0})


def run(fixture_path:Path,output_dir:Path)->dict[str,Any]:
    fixture,raw=load_fixture(fixture_path)
    report,files=make_outputs(fixture)
    manifest=manifest_text(fixture_path,raw,files).encode("utf-8")
    if output_dir.exists():
        raise FileExistsError("Refusing to overwrite existing output directory")
    output_dir.mkdir(parents=True,exist_ok=False)
    for name,b in files.items():
        (output_dir/name).write_bytes(b)
    (output_dir/MANIFEST_NAME).write_bytes(manifest)
    return {"status":"PASS","mode":"OFFLINE_ONLY","local_count":report["local_count"],
            "expanded_count":report["expanded_count"],"newly_seen_count":report["newly_seen_count"],
            "expanded_width_bp":LONG_WIDTH,"artifact_count":len(files),
            "real_api_calls_this_run":0}


def verify(fixture_path:Path,output_dir:Path)->dict[str,Any]:
    fixture,raw=load_fixture(fixture_path)
    expected_report,expected_files=make_outputs(fixture)
    path=output_dir/MANIFEST_NAME
    manifest=json.loads(path.read_text(encoding="utf-8"))
    expected_manifest=json.loads(manifest_text(fixture_path,raw,expected_files))
    if manifest!=expected_manifest:
        raise ValueError("Source/artifact manifest mismatch")
    for filename,expected in expected_files.items():
        if (output_dir/filename).read_bytes()!=expected:
            raise ValueError("Artifact mismatch: "+filename)
    return {"status":"PASS","checksums":"MATCH", "files_verified":len(expected_files),
            "source_sha256":sha256(raw),"local_count":expected_report["local_count"],
            "expanded_count":expected_report["expanded_count"],"real_api_calls_this_run":0}


def main(argv:list[str]|None=None)->int:
    parser=argparse.ArgumentParser(description="Synthetic Day06 sequence context teaching demo")
    parser.add_argument("command",choices=["run","verify"])
    parser.add_argument("--mode",default="offline",choices=["offline","live"])
    parser.add_argument("--fixture",default="day06_synthetic_context.json")
    parser.add_argument("--output-dir",default="day06_output")
    args=parser.parse_args(argv)
    if args.mode!="offline":
        parser.error("LIVE mode is deliberately unavailable; no API access")
    try:
        if args.command=="run":
            result=run(Path(args.fixture),Path(args.output_dir))
        else:
            result=verify(Path(args.fixture),Path(args.output_dir))
    except (ValueError,OSError,TypeError,KeyError,json.JSONDecodeError) as exc:
        parser.exit(status=2,message="DAY06_ERROR: "+str(exc)+"\n")
    print(json.dumps(result,ensure_ascii=False,sort_keys=True))
    return 0

if __name__=="__main__":
    sys.exit(main())
