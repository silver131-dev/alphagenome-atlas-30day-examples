"""Day 27 | A small, OFFLINE AlphaVariant-style teaching CLI.

Orchestrates the *existing* Day21-Day26 synthetic educational modules.
No AlphaGenome SDK, keys, HTTP, patient files, real reference genome or LIVE mode.
Not the private GitHub v0.1.0 `alpha-variant` executable.
Python 3.12+ standard library only.
"""
from __future__ import annotations

import argparse
import csv
from html import escape
import importlib
import json
import hashlib
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import sys
import tempfile
import xml.etree.ElementTree as ET

EVIDENCE = "SYNTHETIC_EDUCATIONAL_DEMO"
MODE = "OFFLINE_ONLY"
ASSEMBLY = "TOY_NOT_GRCH38"
SCHEMA = "day27.synthetic.pipeline.v1"
STAGES = (21, 22, 23, 24, 25, 26)
SOURCE_MODULES = {
    21: "day21_offline_vcf_review_lab",
    22: "day22_offline_qc_explorer_lab",
    23: "day23_offline_annotation_lab",
    24: "day24_offline_refalt_track_lab",
    25: "day25_offline_locus_plot_lab",
    26: "day26_offline_ppm_motif_lab",
}
RUN_FILES = {
    21: "day21/day21_review_packet.json",
    22: "day22/day22_qc_report.json",
    23: "day23/day23_annotation_report.json",
    24: "day24/day24_track_delta_report.json",
    25: "day25/day25_locus_summary.json",
    26: "day26/day26_motif_report.json",
}
SHA_RE = re.compile(r"[0-9a-f]{64}\Z")
SVG_NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG_NS)


def file_hash(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def load_json(path: Path, *, max_bytes=2_000_000) -> dict:
    path = Path(path)
    if not path.is_file() or path.stat().st_size > max_bytes:
        raise ValueError(f"Local JSON missing or too large: {path.name}")
    def reject_duplicates(pairs):
        value = {}
        for k,v in pairs:
            if k in value:
                raise ValueError("Duplicate JSON key")
            value[k] = v
        return value
    obj = json.loads(path.read_text(encoding="utf-8"),
                     parse_constant=lambda x: (_ for _ in ()).throw(ValueError("Non-finite JSON")),
                     object_pairs_hook=reject_duplicates)
    if not isinstance(obj, dict):
        raise ValueError("Expected a JSON object")
    return obj


def verify_bundle(root: Path) -> dict:
    """Check for accidental drift against bundle manifest (not a signed security attestation)."""
    root = Path(root)
    manifest_path = root / "DAY27_SOURCE_MANIFEST.json"
    manifest = load_json(manifest_path)
    if manifest.get("schema_version") != "day27.input_bundle.v1":
        raise ValueError("Day27 fixture/module manifest has unexpected schema")
    files = manifest.get("files")
    expected_files = {f"{name}.py" for name in SOURCE_MODULES.values()}
    expected_files |= {
        "fixtures/day21_synthetic.vcf", "fixtures/day21_toy_reference.json",
        "fixtures/day21_synthetic_regions.bed", "fixtures/day21_synthetic_scores.json",
        "fixtures/day24_synthetic_tracks.json", "fixtures/day26_synthetic_motif_context.json",
    }
    if not isinstance(files,dict) or set(files)!= expected_files:
        raise ValueError("Source manifest file list differs from approved Day 21–26 bundle")
    for relative,record in files.items():
        path = root / relative
        if (not isinstance(record,dict) or not SHA_RE.fullmatch(str(record.get("sha256")))
                or not path.is_file() or path.is_symlink()
                or path.stat().st_size != record.get("bytes")
                or file_hash(path) != record["sha256"]):
            raise ValueError(f"Bundle drift: {relative}")
    return {"manifest_sha256":file_hash(manifest_path), "files_checked":len(files)}


def stage_modules() -> dict:
    """Import only local, already hash-checked Python modules; never install/use SDK."""
    return {day:importlib.import_module(name) for day,name in SOURCE_MODULES.items()}


def check_stage_contract(stage_results: dict, root: Path) -> dict:
    """Validate small cross-day contracts without mixing independently invented scores."""
    reports = {day:load_json(root / RUN_FILES[day]) for day in STAGES}
    d21,d22,d23,d24,d25,d26 = (reports[day] for day in STAGES)
    for day,report in reports.items():
        if (report.get("assembly") != ASSEMBLY or report.get("real_api_calls_this_run") != 0
                or report.get("clinical_interpretation") is not False):
            # Day22 uses a data_kind flag rather than evidence_class; still synthetic.
            raise ValueError(f"Day {day} unexpectedly describes live, real or clinical data")
    if d21.get("evidence_class") != EVIDENCE or d21.get("mode") != MODE:
        raise ValueError("Day21 not a verified synthetic report")
    if d22.get("data_kind") != EVIDENCE or d22.get("mode") != MODE:
        raise ValueError("Day22 not a verified synthetic report")
    for day in (23,24,25,26):
        if reports[day].get("evidence_class") != EVIDENCE or reports[day].get("mode") != MODE:
            raise ValueError(f"Day {day} not explicitly synthetic/offline")
    counts = (d21.get("candidate_count"),d21.get("selected_count"),d21.get("review_required_qc_count"),d21.get("rejected_count"))
    if counts != (8,3,2,3) or (d22.get("candidate_count"),d22.get("selected_count"),d22.get("review_required_count"),d22.get("rejected_count")) != counts:
        raise ValueError("Day21 and Day22 QC counts are inconsistent")
    if (d23.get("candidate_count"),d23.get("selected_count"),d23.get("annotated_count"),d23.get("missing_score_count")) != (8,3,2,1):
        raise ValueError("Day23 annotation contract mismatch")
    if (d23.get("qc_rows") != d21.get("qc_rows") or
            [r["variant_id"] for r in d23.get("annotations",[])] !=
            [r["variant_id"] for r in d21["review_packet"]]):
        raise ValueError("Day23 annotations do not match the Day21 selected variants")
    for row in d23["annotations"]:
        if (row.get("review_status") != "REVIEW_REQUIRED" or
                row.get("score_status") not in ("SYNTHETIC_SCORE","MISSING_SCORE")):
            raise ValueError("Day23 annotation review semantics differ")
    probes = {r.get("label"):r.get("reason") for r in d22.get("supplementary_probe_cases",[])}
    if (probes.get("LOW_QUAL_INSIDE_BED"),probes.get("REF_MISMATCH_INSIDE_BED")) != ("LOW_QUAL","REF_MISMATCH"):
        raise ValueError("Day22 two essential independent QC probes missing")
    variant=d24.get("variant_id")
    if variant != "chr1:110:A>G" or variant != d25.get("variant_id") or variant != d26.get("variant_id"):
        raise ValueError("Pipeline stages do not describe the same toy variant ID")
    if (d24.get("day23_annotation_link") is None or d24.get("model_prediction") is not False or
            d24.get("official_atlas_data") is not False or
            d24.get("largest_absolute_delta_position_1based") != 110 or
            d24.get("largest_absolute_delta") != -0.3):
        raise ValueError("Day24 toy-track link / delta mismatch")
    if (d25.get("model_prediction") is not False or
            d25.get("day23_score_combined_with_delta") is not False or
            d25.get("largest_absolute_delta") != -0.3):
        raise ValueError("Day25 visualization did not preserve the original toy delta")
    if (d26.get("official_motif_atlas") is not False or
            d26.get("day25_linked") is not True or
            d26.get("day25_delta_combined_with_ppm") is not False or
            d26.get("ppm_definition",{}).get("log_odds_pwm_performed") is not False):
        raise ValueError("Day26 PPM source and separate-scale contract mismatch")
    anchor={r["motif_id"]:r["delta_alt_minus_ref"] for r in d26.get("anchor_results",[])}
    if anchor != {"MOTIF_A":-0.12,"MOTIF_B":0.12}:
        raise ValueError("Day26 anchor motif preferences have drifted")
    for day,res in stage_results.items():
        if res.get("status") != "PASS" or res.get("real_api_calls_this_run") != 0:
            raise ValueError(f"Day {day} did not report offline success")
    return {"variant_id":variant,"candidate_count":8,"selected_count":3,"qc_review_count":2,
            "qc_rejected_count":3,"annotated_count":2,"missing_score_count":1,
            "largest_absolute_delta_position_1based":110,"largest_absolute_delta":-0.3,
            "ppm_anchor_deltas":anchor,"biological_mechanism_status":"NOT_VALIDATED",
            "score_combination_performed":False}


def draw_pipeline_svg() -> str:
    """1920x1080 schematic of the six steps, not biological results."""
    r=ET.Element(f"{{{SVG_NS}}}svg",{"viewBox":"0 0 1920 1080","width":"1920","height":"1080",
                                      "role":"img","aria-label":"Six synthetic offline analysis stages, Day 21 to Day 26"})
    def add(name,attrs,text=None):
        el=ET.SubElement(r,f"{{{SVG_NS}}}{name}",attrs)
        if text is not None: el.text=text
    add("rect",{"width":"1920","height":"1080","fill":"#E5DEDA"})
    add("text",{"x":"100","y":"126","font-size":"58","font-weight":"700","fill":"#343946"},"Day 27  |  AlphaVariant Teaching CLI")
    add("text",{"x":"101","y":"185","font-size":"30","fill":"#7A5865"},"Google AI Studio prototype   /   Python modules   /   offline research review packet")
    cards=[(21,"VCF Import","8 artificial ALT candidates"),(22,"QC Trace","3 SELECTED · 2 REVIEW · 3 REJECTED"),
           (23,"Annotation","2 toy scores · 1 missing"),(24,"REF / ALT Delta","5 points · max delta -0.30"),
           (25,"Locus SVG","Traceable toy visualization"),(26,"Motif PPM-like","Motif A -0.12 · B +0.12")]
    for i,(day,title,sub) in enumerate(cards):
        x=95+(i%3)*587; y=270+(i//3)*300
        add("rect",{"x":str(x),"y":str(y),"width":"532","height":"228","rx":"26","fill":"#FFFBF9","stroke":"#D5C3C5","stroke-width":"3"})
        add("circle",{"cx":str(x+65),"cy":str(y+59),"r":"33","fill":"#D6B8B8"})
        add("text",{"x":str(x+41),"y":str(y+71),"font-size":"30","font-weight":"700","fill":"#28384B"},str(day))
        add("text",{"x":str(x+108),"y":str(y+72),"font-size":"37","font-weight":"700","fill":"#343946"},title)
        add("text",{"x":str(x+37),"y":str(y+150),"font-size":"25","fill":"#655A66"},sub)
        add("text",{"x":str(x+37),"y":str(y+195),"font-size":"21","fill":"#357A71"},"OFFLINE  ·  SYNTHETIC  ·  PASS")
        if i in (0,1,3,4):
            add("line",{"x1":str(x+537),"y1":str(y+116),"x2":str(x+560),"y2":str(y+116),"stroke":"#A67E89","stroke-width":"5"})
            add("polygon",{"points":f"{x+558},{y+103} {x+577},{y+116} {x+558},{y+129}","fill":"#A67E89"})
    add("text",{"x":"103","y":"920","font-size":"34","font-weight":"700","fill":"#854D58"},
        "Evidence is NOT clinical: different toy score scales are NEVER added together.")
    add("text",{"x":"103","y":"975","font-size":"27","fill":"#7A6871"},
        "TOY_NOT_GRCH38   |   No AlphaGenome API calls   |   No patient data   |   Manual review required")
    return ET.tostring(r,encoding="unicode")+"\n"


def render_html(summary:dict)->str:
    h=lambda x:escape(str(x),quote=True)
    entries="\n".join(
        f"<tr><td>Day {day}</td><td>{h(SOURCE_MODULES[day])}</td>"
        f"<td><a href='{h(RUN_FILES[day])}'>JSON</a></td><td>PASS</td></tr>" for day in STAGES)
    return f'''<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>Day27 | AlphaVariant OFFLINE CLI</title>
<style>body{{margin:0;background:#E5DEDA;color:#334155;font:16px/1.6 system-ui,"Noto Sans TC",sans-serif}}
main{{max-width:1250px;margin:auto;padding:30px 22px 70px}}
.hero{{background:#fffaf8;padding:24px;border-radius:18px}}h1{{margin:0}}
.warn{{border-left:5px solid #a76565;padding:14px;background:#fff1e7;margin:20px 0;border-radius:9px}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px}}
.card,section{{padding:20px;background:#fffaf8;border-radius:16px;margin:14px 0}}
.num{{font-size:2rem;font-weight:800;color:#356681}}.sub{{font-size:.8rem;color:#64748b}}
table{{width:100%;border-collapse:collapse}}td,th{{text-align:left;padding:10px;border-bottom:1px solid #e7d7d3}}
th{{background:#f4e8e3}}a{{color:#315b83}}img{{width:100%;max-height:540px;object-fit:contain}}
footer{{font-size:.85rem;color:#697386}}</style></head><body><main>
<div class="hero"><h1>Day 27｜封裝分析流程 — 打造自己的 AlphaVariant CLI</h1>
<p>一條指令串起 Day 21–26 的六個獨立、合成的分析模組。</p></div>
<div class="warn"><strong>OFFLINE ONLY · SYNTHETIC · TOY_NOT_GRCH38</strong><br>
沒有 AlphaGenome Atlas／Model LIVE API、沒有病人資料、沒有醫療判斷。
這是教學 CLI，不是私人 GitHub v0.1.0 的 alpha-variant 正式程式。</div>
<div class="cards"><div class="card"><div class="num">{h(summary['candidate_count'])}</div><div class="sub">ALT 候選</div></div>
<div class="card"><div class="num">{h(summary['selected_count'])}</div><div class="sub">通過教學 QC</div></div>
<div class="card"><div class="num">{h(summary['missing_score_count'])}</div><div class="sub">缺人工分數</div></div>
<div class="card"><div class="num">6</div><div class="sub">獨立模組</div></div></div>
<section><h2>依序執行的六個階段</h2><table><thead><tr><th>階段</th><th>重用模組</th><th>證據</th><th>結果</th></tr></thead><tbody>{entries}</tbody></table></section>
<section><h2>CLI Pipeline Overview</h2><img src="day27_pipeline_flow.svg" alt="Synthetic pipeline stages Day 21 through 26"></section>
<section><h2>Day 25 曲線與 Day 26 序列偏好</h2>
<p>兩份不同的合成資料，只以人工變異識別碼串聯，不合併分數，也不推論真實生物機制。</p>
<p><a href="day25/day25_locus_view.html">Day25 軌跡檢視</a> · <a href="day26/day26_motif_view.html">Day26 PPM 檢視</a></p>
<p><a href="day27_pipeline_summary.json">整合摘要 JSON</a> ·
<a href="day27_run_manifest.json">SHA-256 檔案紀錄</a> · <a href="day27_stages.csv">分階段 CSV</a></p></section>
<footer>Evidence: {h(EVIDENCE)} · {h(ASSEMBLY)} · Real API calls: 0 · Review: REQUIRED</footer>
</main></body></html>'''


def relative_hashes(stage_root:Path)->dict:
    result={}
    for path in sorted(stage_root.rglob("*")):
        if not path.is_file() or path.is_symlink():
            continue
        relative=path.relative_to(stage_root).as_posix()
        if relative == "day27_run_manifest.json":
            continue
        result[relative]={"sha256":file_hash(path),"bytes":path.stat().st_size}
    return result


def run_pipeline(output_dir:Path,*,bundle_root:Path|None=None)->dict:
    """Orchestrate verified bundled modules; stage to temp dir before rename."""
    bundle_root=Path(bundle_root) if bundle_root is not None else Path(__file__).resolve().parent
    verified=verify_bundle(bundle_root)
    target=Path(output_dir).expanduser().absolute()
    if target.exists() or target.is_symlink():
        raise ValueError("Output directory already exists; choose a fresh --output-dir (never overwrite)")
    target.parent.mkdir(parents=True,exist_ok=True)
    modules=stage_modules()
    fixtures=bundle_root / "fixtures"
    vcf=fixtures / "day21_synthetic.vcf"
    reference=fixtures / "day21_toy_reference.json"
    regions=fixtures / "day21_synthetic_regions.bed"
    scores=fixtures / "day21_synthetic_scores.json"
    track=fixtures / "day24_synthetic_tracks.json"
    motif=fixtures / "day26_synthetic_motif_context.json"
    temp_root=Path(tempfile.mkdtemp(prefix=".day27-working-",dir=str(target.parent)))
    try:
        results={}
        results[21]=modules[21].execute_offline(vcf,reference,regions,scores,temp_root/"day21")
        results[22]=modules[22].execute_offline(vcf,reference,regions,temp_root/"day22")
        results[23]=modules[23].execute_offline(vcf,reference,regions,scores,temp_root/"day23",
                                                 input_day21=temp_root/RUN_FILES[21])
        results[24]=modules[24].emit_artifacts(track,temp_root/RUN_FILES[23],temp_root/"day24")
        results[25]=modules[25].emit_outputs(temp_root/RUN_FILES[24],track,temp_root/"day25")
        results[26]=modules[26].write_outputs(motif,temp_root/"day26",
                                              day25_summary_path=temp_root/RUN_FILES[25])
        checked=check_stage_contract(results,temp_root)
        stages=[{"day":day,"module":SOURCE_MODULES[day]+".py","status":"PASS",
                 "json_file":RUN_FILES[day],"json_sha256":file_hash(temp_root/RUN_FILES[day])}
                for day in STAGES]
        summary={"schema_version":SCHEMA,"evidence_class":EVIDENCE,"mode":MODE,
                 "assembly":ASSEMBLY,"status":"PASS","pipeline_stage_count":6,
                 "modules_reused_unmodified":True,"source_bundle_file_count":verified["files_checked"],
                 "source_manifest_sha256":verified["manifest_sha256"],
                 "stages":stages,"artifact_manifest":"day27_run_manifest.json",
                 **checked,"score_source":"SYNTHETIC_FIXTURE_ONLY",
                 "model_prediction":False,"official_atlas_data":False,
                 "atlas_query_performed":False,"clinical_interpretation":False,
                 "manual_review":"REQUIRED","real_api_calls_this_run":0}
        svg=draw_pipeline_svg()
        (temp_root/"day27_pipeline_flow.svg").write_text(svg,encoding="utf-8")
        (temp_root/"day27_review_dashboard.html").write_text(render_html(summary),encoding="utf-8")
        with (temp_root/"day27_stages.csv").open("w",encoding="utf-8",newline="") as f:
            w=csv.DictWriter(f,fieldnames=["day","module","status","json_file","json_sha256"],lineterminator="\n")
            w.writeheader();w.writerows(stages)
        (temp_root/"day27_pipeline_summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
        artifact_list=relative_hashes(temp_root)
        manifest={"schema_version":"day27.synthetic.run_manifest.v1","evidence_class":EVIDENCE,
                  "mode":MODE,"source_manifest_sha256":verified["manifest_sha256"],
                  "artifact_count":len(artifact_list),"artifacts":artifact_list,
                  "integrity_scope":"all output artifacts excluding only this manifest",
                  "not_a_digital_signature":True,"real_api_calls_this_run":0}
        (temp_root/"day27_run_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
        if target.exists():
            raise ValueError("Output directory appeared during execution; refusing overwrite")
        temp_root.rename(target)
    finally:
        if temp_root.exists(): shutil.rmtree(temp_root)
    return {"status":"PASS","mode":MODE,"variant_id":checked["variant_id"],
            "pipeline_stage_count":6,"candidate_count":checked["candidate_count"],
            "selected_count":checked["selected_count"],
            "missing_score_count":checked["missing_score_count"],
            "largest_absolute_delta":checked["largest_absolute_delta"],
            "motif_a_delta":checked["ppm_anchor_deltas"]["MOTIF_A"],
            "motif_b_delta":checked["ppm_anchor_deltas"]["MOTIF_B"],
            "summary_file":str(target/"day27_pipeline_summary.json"),
            "html_file":str(target/"day27_review_dashboard.html"),
            "manifest_file":str(target/"day27_run_manifest.json"),
            "real_api_calls_this_run":0}


def verify_run(run_dir:Path)->dict:
    folder=Path(run_dir)
    summary=load_json(folder/"day27_pipeline_summary.json")
    manifest=load_json(folder/"day27_run_manifest.json")
    if (summary.get("schema_version") != SCHEMA or summary.get("status")!="PASS" or
            summary.get("mode") != MODE or summary.get("real_api_calls_this_run")!=0 or
            summary.get("clinical_interpretation") is not False or
            summary.get("pipeline_stage_count")!=6 or
            manifest.get("schema_version")!="day27.synthetic.run_manifest.v1" or
            manifest.get("evidence_class")!=EVIDENCE or manifest.get("mode")!=MODE or
            manifest.get("source_manifest_sha256") != summary.get("source_manifest_sha256")):
        raise ValueError("Run metadata is not a valid Day27 synthetic record")
    files=manifest.get("artifacts")
    if not isinstance(files,dict) or not files or manifest.get("artifact_count") != len(files):
        raise ValueError("Missing artifact manifest entries")
    base=folder.resolve()
    for relative,expected in files.items():
        p=PurePosixPath(relative)
        if (p.is_absolute() or not p.parts or any(x in ("..","") for x in p.parts)
                or relative == "day27_run_manifest.json"
                or not isinstance(expected,dict) or not SHA_RE.fullmatch(str(expected.get("sha256")))):
            raise ValueError("Manifest path/checksum is invalid")
        path=folder.joinpath(*p.parts)
        if (not path.is_file() or path.is_symlink()
                or not path.resolve().is_relative_to(base)
                or path.stat().st_size != expected.get("bytes")
                or file_hash(path) != expected["sha256"]):
            raise ValueError(f"Output file checksum mismatch: {relative}")
    stages=summary.get("stages")
    if not isinstance(stages,list) or [entry.get("day") for entry in stages] != list(STAGES):
        raise ValueError("Run stages missing/out of order")
    for stage in stages:
        relative=stage.get("json_file")
        if relative != RUN_FILES[stage["day"]] or relative not in files:
            raise ValueError("Day stage artifact not in manifest")
        if stage.get("json_sha256") != files[relative]["sha256"]:
            raise ValueError("Stage JSON hash differs from artifact manifest")
    verified_inputs=verify_bundle(Path(__file__).resolve().parent)
    if verified_inputs["manifest_sha256"] != summary.get("source_manifest_sha256"):
        raise ValueError("CLI source manifest differs from run source provenance")
    return {"status":"PASS","mode":MODE,"stage_count":len(stages),
            "artifact_count":len(files),"checksums":"MATCH", "real_api_calls_this_run":0}


def main(argv=None)->int:
    parser=argparse.ArgumentParser(description="AlphaVariant-style EDUCATIONAL CLI: 100% offline synthetic Data (not GitHub alpha-variant)")
    sub=parser.add_subparsers(dest="command",required=True)
    run=sub.add_parser("run",help="Orchestrate Day21-Day26 synthetic stages")
    run.add_argument("--mode",required=True,choices=["demo","offline"],help="LIVE is not implemented and cannot be selected")
    run.add_argument("--output-dir",default="day27_output")
    verify=sub.add_parser("verify",help="Check run artifact hashes and stage metadata")
    verify.add_argument("--run-dir",default="day27_output")
    args=parser.parse_args(argv)
    try:
        if args.command=="run":
            result=run_pipeline(Path(args.output_dir))
        else:
            result=verify_run(Path(args.run_dir))
    except (ValueError, OSError, KeyError, TypeError, RuntimeError, AssertionError) as exc:
        parser.exit(2,f"DAY27_OFFLINE_ERROR: {exc}\n")
    print(json.dumps(result,ensure_ascii=False,sort_keys=True,allow_nan=False))
    return 0


if __name__=="__main__":
    sys.exit(main())
