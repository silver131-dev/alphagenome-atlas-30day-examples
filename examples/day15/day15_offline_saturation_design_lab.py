"""Day 15: synthetic saturation *design* exercise, Python 3.12+ stdlib only.

Four toy DNA positions, three alternatives each, planned repeats and controls.
No physical experiments, no real biology, no clinical claims, no APIs/network.
This tool is NOT a CRISPR design tool, protocol generator or patient analysis.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import io
import json
import random
from pathlib import Path
from typing import Any

BASES = "ACGT"
DATA_KIND = "SYNTHETIC_EDUCATIONAL_DEMO"
CONTROLS = ("WT_REFERENCE", "MOCK_PROCESS_CONTROL", "PERTURBATION_PLACEHOLDER")
ARTIFACT_NAMES = (
    "day15_design_plan.json",
    "day15_snv_candidates.csv",
    "day15_planned_slots.csv",
    "day15_design_notes.md",
    "day15_saturation_map.svg",
    "day15_design_dashboard.html",
)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def validate_spec(spec: Any) -> dict[str, Any]:
    """Accept only a fixed-size fictional sequence and planned design metadata."""
    if not isinstance(spec, dict):
        raise ValueError("Design must be a JSON object")
    expected = {
        "schema_version", "data_kind", "assembly", "chromosome", "start_1based",
        "reference_sequence", "planned_repeats", "randomization_seed",
        "hypothesis_label", "assay_route", "clinical_interpretation",
        "experiment_performed",
    }
    if set(spec) != expected:
        raise ValueError("Unexpected or missing fields (no patient fields allowed)")
    fixed = {
        "schema_version": "DAY15_SYNTHETIC_DESIGN_V1",
        "data_kind": DATA_KIND,
        "assembly": "TOY_NOT_GRCH38",
        "chromosome": "chrToy",
        "hypothesis_label": "SYNTHETIC_SPLICE_QUESTION_ONLY",
        "assay_route": "SPLICING_REPORTER_DESIGN_ONLY",
        "clinical_interpretation": False,
        "experiment_performed": False,
    }
    for k, wanted in fixed.items():
        if type(spec.get(k)) is not type(wanted) or spec[k] != wanted:
            raise ValueError(f"{k} must be the explicit synthetic teaching value")
    sequence = spec.get("reference_sequence")
    if not isinstance(sequence, str) or not (1 <= len(sequence) <= 12):
        raise ValueError("Provide 1–12 fictional reference bases")
    if any(base not in BASES for base in sequence):
        raise ValueError("Reference must contain uppercase A/C/G/T only")
    pos = spec.get("start_1based")
    if type(pos) is not int or not (1 <= pos <= 1000000):
        raise ValueError("start_1based must be a positive toy integer")
    repeats = spec.get("planned_repeats")
    if type(repeats) is not int or not (1 <= repeats <= 6):
        raise ValueError("planned_repeats must be from 1 to 6")
    seed = spec.get("randomization_seed")
    if type(seed) is not int or not (0 <= seed <= 99999999):
        raise ValueError("randomization_seed must be 0–99999999")
    return dict(spec)


def read_spec(path: str | Path) -> tuple[dict[str, Any], str]:
    path = Path(path)
    if not path.is_file() or path.is_symlink():
        raise ValueError("Expected a real local JSON file")
    if path.stat().st_size > 100000:
        raise ValueError("Synthetic design file too large")
    content = path.read_bytes()
    try:
        decoded = json.loads(content.decode("utf-8"))
    except (ValueError, UnicodeError) as exc:
        raise ValueError("Invalid UTF-8 JSON fixture") from exc
    return validate_spec(decoded), sha256_bytes(content)


def enumerate_snvs(spec: dict[str, Any]) -> list[dict[str, Any]]:
    spec = validate_spec(spec)
    candidates = []
    for i, ref in enumerate(spec["reference_sequence"]):
        position = spec["start_1based"] + i
        for alt in BASES:
            if alt == ref:
                continue
            candidates.append({
                "candidate_id": f"SNV_{i + 1:02d}_{ref}_{alt}",
                "variant_id": f"chrToy:{position}:{ref}>{alt}",
                "position_1based": position,
                "ref": ref,
                "alt": alt,
                "design_class": "SYNTHETIC_SNVS_ONLY",
                "model_prediction": "NOT_RUN",
                "observed_effect": None,
                "measurement_status": "NOT_MEASURED",
            })
    return candidates


def allocate_slots(spec: dict[str, Any], candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Make pre-declared teaching slots, not actual independent biological samples."""
    spec = validate_spec(spec)
    if candidates != enumerate_snvs(spec):
        raise ValueError("Candidate rows must exactly match the synthetic substitution space")
    items = [(c["candidate_id"], "VARIANT", c["variant_id"]) for c in candidates]
    items.extend((code, "CONTROL_PLACEHOLDER", None) for code in CONTROLS)
    rows = []
    for repeat in range(1, spec["planned_repeats"] + 1):
        shuffled = list(items)
        random.Random(spec["randomization_seed"] + repeat).shuffle(shuffled)
        for j, (candidate_id, kind, variant_id) in enumerate(shuffled, 1):
            rows.append({
                "slot_id": f"RUN{repeat}_SLOT{j:02d}",
                "planned_run": f"R{repeat}",
                "candidate_id": candidate_id,
                "item_type": kind,
                "variant_id": variant_id,
                "observed_readout": None,
                "measurement_status": "NOT_MEASURED",
                "biological_independence_verified": False,
            })
    return rows


def build_plan(spec: dict[str, Any], fixture_sha256: str) -> dict[str, Any]:
    spec = validate_spec(spec)
    if not isinstance(fixture_sha256, str) or len(fixture_sha256) != 64 or any(c not in "0123456789abcdef" for c in fixture_sha256):
        raise ValueError("Expected a lowercase SHA-256 digest")
    candidates = enumerate_snvs(spec)
    slots = allocate_slots(spec, candidates)
    expected_count = len(spec["reference_sequence"]) * 3
    if len(candidates) != expected_count or len(slots) != (expected_count + len(CONTROLS)) * spec["planned_repeats"]:
        raise AssertionError("Toy saturation coverage contract failed")
    return {
        "schema_version": "DAY15_SYNTHETIC_PLAN_V1",
        "evidence_class": DATA_KIND,
        "assembly": spec["assembly"],
        "hypothesis_label": spec["hypothesis_label"],
        "assay_route": spec["assay_route"],
        "reference": {
            "chromosome": "chrToy",
            "start_1based": spec["start_1based"],
            "sequence": spec["reference_sequence"],
        },
        "expected_possible_snvs": expected_count,
        "candidate_count": len(candidates),
        "toy_snv_coverage_fraction": 1.0,
        "control_placeholder_count": len(CONTROLS),
        "control_placeholder_labels": list(CONTROLS),
        "planned_repeats": spec["planned_repeats"],
        "planned_slots_per_run": len(candidates) + len(CONTROLS),
        "planned_slot_count": len(slots),
        "randomization_seed": spec["randomization_seed"],
        "candidate_variants": candidates,
        "planned_slots": slots,
        "fixture_provenance": {"source_kind": "LOCAL_SYNTHETIC_FIXTURE", "sha256": fixture_sha256},
        "readout_available": False,
        "independent_biological_replicates_confirmed": False,
        "experiment_performed": False,
        "saturation_genome_editing_performed": False,
        "minigene_assay_performed": False,
        "reporter_assay_performed": False,
        "clinical_interpretation": False,
        "research_ethics_approval": "NOT_ASSESSED",
        "real_api_calls_this_run": 0,
        "interpretation_status": "DESIGN_ONLY_NOT_BIOLOGICAL_EVIDENCE",
    }


def to_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n"


def to_csv(rows: list[dict[str, Any]], columns: tuple[str, ...]) -> str:
    out = io.StringIO(newline="")
    writer = csv.DictWriter(out, fieldnames=list(columns), lineterminator="\n", extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow(row)
    return out.getvalue()


def render_svg(plan: dict[str, Any]) -> str:
    """Static 1920x1080 overview, only displaying fictional design metadata."""
    seq = plan["reference"]["sequence"]
    positions = [plan["reference"]["start_1based"] + i for i in range(len(seq))]
    if len(seq) > 12:
        raise ValueError("Too many toy positions for visualization")
    def esc(x: Any) -> str:
        return html.escape(str(x), quote=True)
    nodes = []
    nodes.append('<svg xmlns="http://www.w3.org/2000/svg" width="1920" height="1080" viewBox="0 0 1920 1080" role="img" aria-label="Synthetic Day15 saturation design overview">')
    nodes.append('<rect width="1920" height="1080" fill="#E5DEDA"/>')
    nodes.append('<rect x="68" y="62" width="1784" height="952" rx="34" fill="#F9F7F5" stroke="#D0C5BF" stroke-width="2"/>')
    nodes.append('<text x="130" y="148" font-family="Arial,sans-serif" font-size="56" fill="#26384A" font-weight="700">Day 15 | From prediction to testable design</text>')
    nodes.append('<text x="132" y="202" font-family="Arial,sans-serif" font-size="25" fill="#576A78">SYNTHETIC / IN SILICO ONLY  ·  Not a wet-lab result</text>')
    boxes = [
        (132, 276, 500, 190, "01 | Hypothesis", "Model signal is a question", "not a confirmed mechanism", "#E7EBF1"),
        (710, 276, 500, 190, "02 | Variant space", f'{len(seq)} toy bases × 3 alternatives', f'{plan["candidate_count"]} planned SNVs', "#F2E4E0"),
        (1288, 276, 500, 190, "03 | Controls", "WT / mock / perturbation", "all are unvalidated placeholders", "#E5F0EC"),
    ]
    for x,y,w,h,heading,main,sub,color in boxes:
        nodes.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="26" fill="{color}"/>')
        nodes.append(f'<text x="{x+28}" y="{y+52}" font-family="Arial,sans-serif" font-size="30" font-weight="700" fill="#24354A">{esc(heading)}</text>')
        nodes.append(f'<text x="{x+28}" y="{y+105}" font-family="Arial,sans-serif" font-size="29" fill="#344961">{esc(main)}</text>')
        nodes.append(f'<text x="{x+28}" y="{y+146}" font-family="Arial,sans-serif" font-size="23" fill="#687987">{esc(sub)}</text>')
    nodes.append('<text x="134" y="559" font-family="Arial,sans-serif" font-size="34" font-weight="700" fill="#26384A">Four toy bases, three alternatives per position</text>')
    if len(seq) <= 6:
        cell_w, gap, x0 = 365, 20, 140
    else:
        cell_w, gap, x0 = min(250, int(1560/len(seq))), 12, 138
    for i, (position, ref) in enumerate(zip(positions, seq)):
        x = x0+i*(cell_w+gap)
        if x+cell_w>1800: break
        nodes.append(f'<rect x="{x}" y="595" width="{cell_w}" height="170" rx="20" fill="#FFFFFF" stroke="#D9DBE1" stroke-width="2"/>')
        nodes.append(f'<text x="{x+20}" y="639" font-family="Arial,sans-serif" font-size="24" fill="#57677A">Toy pos {position}</text>')
        nodes.append(f'<text x="{x+20}" y="695" font-family="Arial,sans-serif" font-size="42" font-weight="700" fill="#2B5B77">{esc(ref)}</text>')
        alts=' / '.join(b for b in BASES if b!=ref)
        nodes.append(f'<text x="{x+95}" y="695" font-family="Arial,sans-serif" font-size="27" fill="#C37160">→ {esc(alts)}</text>')
    nodes.append('<rect x="132" y="812" width="1656" height="116" rx="22" fill="#293B50"/>')
    nodes.append(f'<text x="174" y="861" font-family="Arial,sans-serif" font-size="30" fill="#F8F7F6">{plan["planned_repeats"]} planned runs × {plan["planned_slots_per_run"]} slots = {plan["planned_slot_count"]} assignments</text>')
    nodes.append('<text x="174" y="903" font-family="Arial,sans-serif" font-size="22" fill="#D7DCE4">NO RNA measurements · NO cell culture · NO editing · NO API calls</text>')
    nodes.append('<text x="134" y="982" font-family="Arial,sans-serif" font-size="22" fill="#8B646A">Research question → enumeration → controls → allocation → future independent validation</text>')
    nodes.append('</svg>')
    return '\n'.join(nodes)+'\n'


def render_notes(plan: dict[str, Any]) -> str:
    return (f"# Day15 Synthetic Study Design — Not Experimental Evidence\n\n"
            f"- Reference: `chrToy` `{plan['reference']['sequence']}` at 1-based {plan['reference']['start_1based']}\n"
            f"- Intended substitution coverage: **{plan['candidate_count']}/{plan['expected_possible_snvs']}** toy SNVs\n"
            f"- Control labels: {', '.join(CONTROLS)} (all **unvalidated placeholders**)\n"
            f"- Planned runs: **{plan['planned_repeats']}**; proposed slots: **{plan['planned_slot_count']}**\n"
            f"- Source fixture SHA-256: `{plan['fixture_provenance']['sha256']}`\n"
            f"- Readouts: **NOT_MEASURED**; experiment: **NOT_PERFORMED**\n"
            f"- Repeat labels are design bookkeeping, **not verified biological replicates**.\n"
            f"- Saturation genome editing, mini-gene and reporter assays: **NOT_PERFORMED**.\n"
            f"- No validated assay-positive control, no phenotype effect, no clinical inference.\n"
            f"- SHA-256 shows file identity, not scientific validity or source authenticity.\n"
            f"- Future real study requires separate laboratory design, ethics, biosafety, controls, readout quality and independent validation.\n")


def render_html(plan: dict[str, Any], svg: str) -> str:
    rows = ''.join('<tr>' + ''.join(f'<td>{html.escape(str(c[k]))}</td>' for k in ('candidate_id', 'variant_id', 'measurement_status'))+'</tr>' for c in plan['candidate_variants'])
    return ('<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">'
            '<title>Day15 Synthetic Experiment Design</title><style>'
            'body{margin:0;background:#E5DEDA;color:#273648;font:17px/1.7 system-ui,sans-serif}main{max-width:1050px;margin:40px auto;padding:36px;background:#fff;border-radius:22px}'
            'h1{margin:0 0 8px}table{width:100%;border-collapse:collapse}th,td{border-bottom:1px solid #ddd;padding:10px;text-align:left}th{background:#f3eeeb}'
            'svg{width:100%;height:auto}small{color:#755a5a}.pill{display:inline-block;background:#f4ece8;border-radius:30px;padding:5px 15px}</style></head><body><main>'
            '<p class="pill">SYNTHETIC EDUCATIONAL DEMO · OFFLINE ONLY</p><h1>Day15 | 從模型假說到實驗設計</h1>'
            f'<p>這是一份人工的排程草案：{plan["candidate_count"]} 種 SNV，{plan["control_placeholder_count"]} 個對照佔位，{plan["planned_repeats"]} 次預定執行，共 {plan["planned_slot_count"]} 個欄位。</p>'
            '<p><strong>NOT_MEASURED / NOT_PERFORMED</strong>：沒有執行任何基因編輯、細胞培養、RNA 測量或 AlphaGenome API。</p>'
            f'{svg}<h2>人工候選一覽</h2><table><thead><tr><th>候選</th><th>位置與替換</th><th>觀測狀態</th></tr></thead><tbody>{rows}</tbody></table>'
            '<p><small>每項「planned run」只是規劃的重複標籤，不表示已取得獨立生物複本。此文件不是實驗方案批准或臨床報告。</small></p></main></body></html>\n')


def artifact_bytes(plan: dict[str, Any]) -> dict[str, bytes]:
    candidates = plan['candidate_variants']
    slots = plan['planned_slots']
    svg = render_svg(plan)
    artifacts = {
        "day15_design_plan.json": to_json(plan),
        "day15_snv_candidates.csv": to_csv(candidates, ('candidate_id','variant_id','position_1based','ref','alt','design_class','model_prediction','observed_effect','measurement_status')),
        "day15_planned_slots.csv": to_csv(slots, ('slot_id','planned_run','candidate_id','item_type','variant_id','observed_readout','measurement_status','biological_independence_verified')),
        "day15_design_notes.md": render_notes(plan),
        "day15_saturation_map.svg": svg,
        "day15_design_dashboard.html": render_html(plan,svg),
    }
    if set(artifacts) != set(ARTIFACT_NAMES):
        raise AssertionError('Artifact list incomplete')
    return {name: value.encode('utf-8') for name, value in artifacts.items()}


def run(fixture: str | Path, output_dir: str | Path) -> dict[str, Any]:
    spec, digest = read_spec(fixture)
    plan = build_plan(spec, digest)
    folder = Path(output_dir)
    if folder.is_symlink() or folder.exists():
        raise FileExistsError('Refusing to overwrite an existing output directory')
    if not folder.parent.is_dir():
        raise ValueError('Output parent must already exist')
    artifacts = artifact_bytes(plan)
    folder.mkdir()
    for name, content in artifacts.items():
        (folder/name).write_bytes(content)
    manifest = {
        "schema_version": "DAY15_MANIFEST_V1",
        "evidence_class": DATA_KIND,
        "fixture_sha256": digest,
        "file_sha256": {name:sha256_bytes(data) for name,data in sorted(artifacts.items())},
        "published_or_validated": False,
        "real_api_calls_this_run": 0,
    }
    (folder/'day15_manifest.json').write_text(to_json(manifest),encoding='utf-8')
    return {"status":"PASS", "mode":"OFFLINE_ONLY", "candidate_count":plan['candidate_count'],
            "control_placeholder_count":len(CONTROLS),"planned_repeats":plan['planned_repeats'],
            "planned_slot_count":plan['planned_slot_count'],"artifact_count":len(artifacts)+1,
            "experiment_performed":False,"real_api_calls_this_run":0}


def verify(fixture: str | Path, output_dir: str | Path) -> dict[str, Any]:
    _, fixture_sha = read_spec(fixture)
    folder = Path(output_dir)
    manifest_path = folder/'day15_manifest.json'
    if not manifest_path.is_file() or manifest_path.is_symlink():
        raise ValueError('Manifest missing or not a regular file')
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    if manifest.get('schema_version') != 'DAY15_MANIFEST_V1' or manifest.get('fixture_sha256')!=fixture_sha:
        raise ValueError('Fixture or manifest mismatch')
    hashes = manifest.get('file_sha256')
    if not isinstance(hashes,dict) or set(hashes)!=set(ARTIFACT_NAMES):
        raise ValueError('Artifact inventory mismatch')
    for name in ARTIFACT_NAMES:
        path=folder/name
        if not path.is_file() or path.is_symlink():
            raise ValueError(f'Missing or unsupported artifact: {name}')
        if sha256_bytes(path.read_bytes())!=hashes[name]:
            raise ValueError(f'Artifact content changed: {name}')
    plan=json.loads((folder/'day15_design_plan.json').read_text(encoding='utf-8'))
    if plan.get('experiment_performed') is not False or plan.get('clinical_interpretation') is not False or plan.get('real_api_calls_this_run')!=0:
        raise ValueError('Evidence status contract mismatch')
    if len(plan.get('candidate_variants',[]))!=plan.get('candidate_count') or len(plan.get('planned_slots',[]))!=plan.get('planned_slot_count'):
        raise ValueError('Row count mismatch')
    return {"status":"PASS","checksums":"MATCH","artifact_count":len(ARTIFACT_NAMES),
            "evidence_class":DATA_KIND,"experiment_performed":False,"real_api_calls_this_run":0}


def main(argv: list[str] | None = None) -> int:
    parser=argparse.ArgumentParser(description='Day15 local synthetic saturation design (NOT a wet-lab protocol)')
    parser.add_argument('command', choices=['run','verify'])
    parser.add_argument('--mode',choices=['offline'],default='offline')
    parser.add_argument('--fixture',default='day15_synthetic_design.json')
    parser.add_argument('--output-dir',default='day15_output')
    args=parser.parse_args(argv)
    result=run(args.fixture,args.output_dir) if args.command=='run' else verify(args.fixture,args.output_dir)
    print(json.dumps(result,ensure_ascii=False,sort_keys=True))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
