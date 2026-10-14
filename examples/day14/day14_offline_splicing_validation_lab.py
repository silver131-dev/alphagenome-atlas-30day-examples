"""Day 14 — standalone, local-only synthetic splice-validation DESIGN lesson.

This script neither predicts DNM1 splicing nor runs a wet-lab experiment.
It never uses a patient VCF, credentials, or a network API.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import io
import json
from pathlib import Path
from typing import Any

DATA_KIND = "SYNTHETIC_EDUCATIONAL_DESIGN"
SCENARIO_ID = "DNM1_CRYPTIC_SPLICE_TEACHING"
CONSTRUCTS = {
    "WT": "REFERENCE_CONTROL",
    "ALT": "CANDIDATE_VARIANT",
    "REVERT": "SEQUENCE_REVERSION_CONTROL",
    "EMPTY": "VECTOR_BACKGROUND_CONTROL",
}
CONTEXTS = {
    "NEURONAL_LIKE": "HYPOTHESIS_RELEVANT_MODEL",
    "COMPARATOR": "CONTEXT_COMPARATOR",
}
ASSAYS = ["RT_PCR", "JUNCTION_SEQUENCING"]
FIELDNAMES = [
    "slot_id", "cell_context", "construct", "planned_biological_replicate",
    "expected_pattern", "rt_pcr_result", "junction_sequence_result",
    "observed_sample", "experiment_status",
]
OUTPUT_NAMES = [
    "day14_splicing_design.json", "day14_planned_matrix.csv",
    "day14_hypothesis_review.md", "day14_design_flow.svg",
    "day14_design_dashboard.html",
]


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def json_bytes(obj: Any) -> bytes:
    return (json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode("utf-8")


def check_fixture(fixture: Any) -> dict[str, Any]:
    """Reject extra fields and allow only an explicitly fixed *fictional* design."""
    if not isinstance(fixture, dict):
        raise ValueError("The teaching fixture must be a JSON object")
    expected = {
        "data_kind", "scenario_id", "gene_symbol", "reference_allele",
        "alternative_allele", "genomic_context", "putative_extension_nt",
        "constructs", "cell_contexts", "planned_biological_replicates",
        "planned_assays",
    }
    if set(fixture) != expected:
        raise ValueError("Fixture fields do not match the synthetic design schema")
    required_literals = {
        "data_kind": DATA_KIND,
        "scenario_id": SCENARIO_ID,
        "gene_symbol": "DNM1",
        "reference_allele": "G",
        "alternative_allele": "A",
        "genomic_context": "TOY_NO_COORDINATES",
        "planned_assays": ASSAYS,
    }
    for name, literal in required_literals.items():
        if fixture[name] != literal:
            raise ValueError(f"Unexpected fictional fixture value: {name}")
    if (not isinstance(fixture["putative_extension_nt"], int)
            or isinstance(fixture["putative_extension_nt"], bool)
            or not 1 <= fixture["putative_extension_nt"] <= 1000):
        raise ValueError("putative_extension_nt must be a positive teaching integer")
    replicates = fixture["planned_biological_replicates"]
    if not isinstance(replicates, int) or isinstance(replicates, bool) or not 2 <= replicates <= 6:
        raise ValueError("Expected 2–6 PLANNED biological replicate labels")
    for key, expected_roles in (("constructs", CONSTRUCTS), ("cell_contexts", CONTEXTS)):
        raw = fixture[key]
        if not isinstance(raw, list) or len(raw) != len(expected_roles):
            raise ValueError(f"Incorrect {key} count")
        if any(not isinstance(p, dict) or set(p) != {"id", "role"} for p in raw):
            raise ValueError(f"Unexpected {key} entry")
        roles = {p["id"]: p["role"] for p in raw if isinstance(p["id"], str)}
        if len(roles) != len(raw) or roles != expected_roles:
            raise ValueError(f"Incorrect or duplicated {key} identifiers and roles")
    return fixture


def proposed_pattern(construct: str, context: str) -> str:
    """Qualitative *hypothesis*, never an observed or predicted probability."""
    if construct == "EMPTY":
        return "ASSAY_BACKGROUND_ONLY"
    if construct == "WT":
        return "BASELINE_TO_MEASURE"
    if construct == "REVERT":
        return "RESTORATION_TO_TEST"
    if context == "NEURONAL_LIKE":
        return "HYPOTHESIS_ALT_MORE_CRYPTIC_JUNCTION_THAN_WT"
    return "CONTEXT_EFFECT_UNCERTAIN"


def build_plan(fixture: dict[str, Any], source_name: str, source_sha: str) -> dict[str, Any]:
    check_fixture(fixture)
    slots = []
    n = fixture["planned_biological_replicates"]
    for context in fixture["cell_contexts"]:
        for construct in fixture["constructs"]:
            for replicate in range(1, n + 1):
                slots.append({
                    "slot_id": f"PLAN_{len(slots) + 1:03d}",
                    "cell_context": context["id"],
                    "construct": construct["id"],
                    "planned_biological_replicate": replicate,
                    "expected_pattern": proposed_pattern(construct["id"], context["id"]),
                    "rt_pcr_result": None,
                    "junction_sequence_result": None,
                    "observed_sample": False,
                    "experiment_status": "NOT_PERFORMED",
                })
    extension_nt = fixture["putative_extension_nt"]
    return {
        "report_kind": DATA_KIND,
        "scenario_id": SCENARIO_ID,
        "gene_symbol": "DNM1",
        "variant_notation": "G>A (NO_COORDINATE_TOY)",
        "assembly": "TOY_NOT_GRCH38",
        "evidence_class": "SYNTHETIC_EXPERIMENT_DESIGN_ONLY",
        "source_provenance": {"filename": source_name, "sha256": source_sha},
        "construct_count": len(CONSTRUCTS),
        "cell_context_count": len(CONTEXTS),
        "planned_biological_replicates": n,
        "planned_slot_count": len(slots),
        "observed_sample_count": 0,
        "putative_extension_nt": extension_nt,
        "in_frame_if_in_coding_region": extension_nt % 3 == 0,
        "potential_added_codons_if_translated": extension_nt // 3 if extension_nt % 3 == 0 else None,
        "protein_function_verified": False,
        "splicing_measured": False,
        "experiment_performed": False,
        "live_api_calls": 0,
        "clinical_interpretation": False,
        "planned_assays": list(ASSAYS),
        "orthogonal_evidence_needed": [
            "WT versus ALT RT-PCR band pattern comparison",
            "Sequence-confirmed splice junction and exact insertion length",
            "Sequence-reversion comparison to probe specificity",
            "Consider context-dependent endogenous RNA evidence if available",
            "Evaluate protein consequence separately; RNA result alone is insufficient",
        ],
        "hypothesis_rejection_examples": [
            "No reproducible ALT-versus-WT splice difference in the planned relevant model",
            "Same unexpected band appears in empty vector or WT control",
            "Sequenced junction does not match the proposed cryptic acceptor",
            "Reversion fails to restore the baseline pattern under comparable conditions",
        ],
        "limitations": [
            "All conditions and replicates are PLANNED, no cells or patient samples measured",
            "The 39 nt is a teaching assumption; 13 codons is conditional arithmetic, not protein evidence",
            "Mini-gene constructs may omit endogenous chromatin, transcription and tissue context",
            "Neither molecular causality nor clinical pathogenicity has been established",
        ],
        "governance_review_status": "MANUAL_REVIEW_REQUIRED",
        "external_transfer_allowed": False,
        "data_release_allowed": False,
        "slots": slots,
    }


def csv_text(plan: dict[str, Any]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    for slot in plan["slots"]:
        writer.writerow({k: "" if slot[k] is None else slot[k] for k in FIELDNAMES})
    return output.getvalue()


def summary_markdown(plan: dict[str, Any]) -> str:
    checks = "\n".join("- " + q for q in plan["hypothesis_rejection_examples"])
    return (
        "# Day 14 | Synthetic DNM1 Splicing Validation Design\n\n"
        "**Design only:** no DNM1 RNA was measured; no patients or cell experiments.\n\n"
        f"- Planned constructs: {plan['construct_count']} (WT / ALT / REVERT / EMPTY)\n"
        f"- Planned contexts: {plan['cell_context_count']} (NEURONAL_LIKE / COMPARATOR)\n"
        f"- Planned biological replicate labels per combination: {plan['planned_biological_replicates']}\n"
        f"- Planned slots: {plan['planned_slot_count']}; observed samples: 0\n"
        f"- Hypothetical inclusion: {plan['putative_extension_nt']} nucleotides, possibly "
        f"{plan['potential_added_codons_if_translated']} extra codons IF in frame and translated.\n"
        "- Planned readouts: RT-PCR and sequencing of junctions.\n"
        "- Observed splice fraction / protein function / clinical significance: **NOT MEASURED**.\n"
        "- External transfer: denied; human review: required.\n\n"
        "## What would challenge the hypothesis?\n\n" + checks + "\n\n"
        "A study report is not an experiment result. This exercise is not an anonymisation or clinical compliance system.\n"
    )


def svg_text(plan: dict[str, Any]) -> str:
    # All displayed strings are from the fixed synthetic schema.
    title = html.escape("Day 14 | DNM1 splice hypothesis: design before conclusions")
    bg = "#E5DEDA"
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="1920" height="1080" viewBox="0 0 1920 1080">',
        '<defs><linearGradient id="panel" x1="0" x2="1"><stop stop-color="#fff9f6"/><stop offset="1" stop-color="#f7f0f2"/></linearGradient></defs>',
        f'<rect width="1920" height="1080" fill="{bg}"/>',
        '<g font-family="DejaVu Sans,Arial,sans-serif" fill="#39414b">',
        f'<text x="110" y="125" font-size="48" font-weight="700">{title}</text>',
        '<text x="110" y="177" font-size="26">Synthetic educational design — no experiment, no patient data, no API call</text>',
        '<rect x="92" y="219" width="1736" height="582" rx="27" fill="url(#panel)" stroke="#d0bec1" stroke-width="2"/>',
        '<text x="142" y="280" font-size="31" font-weight="700">From molecular hypothesis to independent verification</text>',
        '<text x="144" y="341" font-size="25">DNA G&gt;A (toy)  →  putative cryptic acceptor  →  possible 39 nt extension</text>',
        '<text x="146" y="381" font-size="21" fill="#745c66">39 nt = 13 codons only IF coding frame is retained. Protein effect is NOT measured.</text>',
        '<text x="143" y="456" font-size="25" font-weight="700">Planned constructs</text>',
    ]
    cards = [
        (142, "WT", "baseline"), (555, "ALT", "candidate"),
        (968, "REVERT", "specificity"), (1381, "EMPTY", "background"),
    ]
    for x, name, sub in cards:
        fill = "#e8d6d7" if name == "ALT" else "#edf2f1"
        lines += [
            f'<rect x="{x}" y="483" width="350" height="162" rx="21" fill="{fill}" stroke="#c3b5b6"/>',
            f'<text x="{x+30}" y="545" font-size="33" font-weight="700">{name}</text>',
            f'<text x="{x+30}" y="591" font-size="23">{sub}</text>',
        ]
    lines += [
        '<text x="143" y="705" font-size="25">2 cell contexts × 4 constructs × 3 planned replicates = 24 planned slots</text>',
        '<text x="143" y="747" font-size="24">Independent checks: RT-PCR → junction sequencing → revert / controls</text>',
        '<text x="109" y="871" font-size="31" font-weight="700">What would falsify the hypothesis?</text>',
        '<text x="109" y="919" font-size="25">No reproducible ALT shift  •  unexpected WT / empty band  •  wrong junction  •  no rescue</text>',
        '<text x="109" y="985" font-size="21" fill="#855c63">NOT_PERFORMED   |   OBSERVED_SAMPLES=0   |   MANUAL_REVIEW_REQUIRED</text>',
        '</g></svg>\n',
    ]
    return "\n".join(lines)


def html_text(plan: dict[str, Any], svg: str) -> str:
    heads = "".join("<th>" + html.escape(k) + "</th>" for k in FIELDNAMES)
    rows = []
    for slot in plan["slots"]:
        cells = []
        for key in FIELDNAMES:
            value = slot[key]
            cells.append("<td>" + html.escape("NOT_MEASURED" if value is None else str(value)) + "</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    return (
        '<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8">'
        '<title>Day 14 | Synthetic Splicing Design</title>'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        '<style>body{margin:0;background:#E5DEDA;color:#344052;font:16px system-ui,sans-serif}'
        'main{max-width:1200px;margin:auto;padding:22px}'
        'svg{width:100%;height:auto;background:#E5DEDA}'
        '.warn{padding:15px;border:1px solid #9d7375;background:#fff9f5;border-radius:12px}'
        'table{border-collapse:collapse;width:100%;font-size:12px;background:#fff}'
        'td,th{padding:9px;border:1px solid #d4c9cb;text-align:left}'
        '.tablewrap{overflow-x:auto}h1{font-size:26px}</style></head><body><main>'
        '<h1>Day 14｜DNM1 剪接假說：獨立驗證設計</h1>'
        '<p class="warn">SYNTHETIC EDUCATIONAL DESIGN — 無病人、無模型查詢、無濕實驗；所有配置均為規劃。</p>'
        + svg + '<h2>Planned experiment matrix（不是實測結果）</h2><div class="tablewrap"><table><thead><tr>'
        + heads + '</tr></thead><tbody>' + "".join(rows) + '</tbody></table></div>'
        '<p>RNA / Protein / Clinical result: NOT MEASURED. Data release: BLOCKED. Human review: REQUIRED.</p>'
        '</main></body></html>\n'
    )


def render_artifacts(plan: dict[str, Any]) -> dict[str, bytes]:
    svg = svg_text(plan)
    return {
        OUTPUT_NAMES[0]: json_bytes(plan),
        OUTPUT_NAMES[1]: csv_text(plan).encode("utf-8"),
        OUTPUT_NAMES[2]: summary_markdown(plan).encode("utf-8"),
        OUTPUT_NAMES[3]: svg.encode("utf-8"),
        OUTPUT_NAMES[4]: html_text(plan, svg).encode("utf-8"),
    }


def source_plan(fixture_path: Path) -> tuple[dict[str, Any], bytes]:
    data = fixture_path.read_bytes()
    fixture = json.loads(data)
    check_fixture(fixture)
    return build_plan(fixture, fixture_path.name, digest(data)), data


def run(fixture_path: Path, output_dir: Path) -> dict[str, Any]:
    plan, source = source_plan(fixture_path)
    if output_dir.exists() and (not output_dir.is_dir() or any(output_dir.iterdir())):
        raise FileExistsError("Choose a NEW or empty output directory; do not overwrite prior evidence")
    rendered = render_artifacts(plan)
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, data in rendered.items():
        (output_dir / name).write_bytes(data)
    manifest = {
        "manifest_kind": "SYNTHETIC_DAY14_MANIFEST_NOT_SIGNED",
        "source_fixture": {"filename": fixture_path.name, "sha256": digest(source)},
        "artifact_hashes": {name: digest(data) for name, data in rendered.items()},
        "signature_present": False,
        "experiment_performed": False,
    }
    (output_dir / "day14_artifact_manifest.json").write_bytes(json_bytes(manifest))
    return {
        "status": "PASS", "mode": "OFFLINE_ONLY",
        "planned_slot_count": plan["planned_slot_count"],
        "observed_sample_count": 0, "construct_count": plan["construct_count"],
        "cell_context_count": plan["cell_context_count"],
        "artifact_count": len(rendered),
        "experiment_performed": False,
        "real_api_calls_this_run": 0,
        "output_dir": str(output_dir),
    }


def verify(fixture_path: Path, output_dir: Path) -> dict[str, Any]:
    plan, source = source_plan(fixture_path)
    manifest = json.loads((output_dir / "day14_artifact_manifest.json").read_text(encoding="utf-8"))
    expected_hashes = {name: digest(data) for name, data in render_artifacts(plan).items()}
    if (set(manifest) != {"manifest_kind", "source_fixture", "artifact_hashes", "signature_present", "experiment_performed"}
            or manifest["manifest_kind"] != "SYNTHETIC_DAY14_MANIFEST_NOT_SIGNED"
            or manifest["signature_present"] is not False
            or manifest["experiment_performed"] is not False
            or manifest["source_fixture"] != {"filename": fixture_path.name, "sha256": digest(source)}
            or manifest["artifact_hashes"] != expected_hashes):
        raise ValueError("Manifest or source fixture does not match the expected synthetic run")
    for filename, sha in expected_hashes.items():
        if digest((output_dir / filename).read_bytes()) != sha:
            raise ValueError(f"Artifact checksum mismatch: {filename}")
    return {
        "status": "PASS", "checksums": "MATCH",
        "verified_artifact_count": len(expected_hashes),
        "planned_slot_count": plan["planned_slot_count"],
        "experiment_performed": False,
        "real_api_calls_this_run": 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="DNM1 synthetic minigene validation DESIGN lesson")
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("run", "verify"):
        p = sub.add_parser(name)
        p.add_argument("--fixture", default="day14_synthetic_splicing_design.json")
        p.add_argument("--output-dir", default="day14_output")
        if name == "run":
            p.add_argument("--mode", choices=["offline"], required=True)
    args = parser.parse_args()
    fixture = Path(args.fixture)
    folder = Path(args.output_dir)
    result = run(fixture, folder) if args.command == "run" else verify(fixture, folder)
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
