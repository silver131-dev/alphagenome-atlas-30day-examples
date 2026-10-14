"""Day 05 | Choose an evidence route before choosing an AI model.

SYNTHETIC EDUCATIONAL DEMO ONLY.
No AlphaFold / AlphaGenome execution; no API, network, patients or diagnosis.
Python 3.12+ standard library.
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
from pathlib import Path
from typing import Any

DATA_KIND = "SYNTHETIC_EDUCATIONAL_DEMO"
ASSEMBLY = "TOY_NOT_GRCH38"
SCHEMA_VERSION = "1.0"
QUESTIONS = frozenset({
    "PROTEIN_STRUCTURE", "MOLECULAR_COMPLEX", "RNA_SPLICING",
    "GENE_EXPRESSION", "CHROMATIN_ACCESSIBILITY", "UNSPECIFIED",
})
CONSEQUENCES = frozenset({"MISSENSE", "PROMOTER", "INTRONIC", "SPLICE_REGION", "UNKNOWN"})
CASE_KEYS = frozenset({"case_id", "variant", "consequence", "research_question"})
TOP_KEYS = frozenset({"data_kind", "assembly", "schema_version", "cases"})
TOY_ID_RE = re.compile(r"TOY_[A-Z0-9_]{2,40}\Z")
TOY_VARIANT_RE = re.compile(r"chrToy:([1-9][0-9]*):([ACGT])>([ACGT])\Z")
ARTIFACTS = (
    "day05_routing_report.json", "day05_routing_table.csv",
    "day05_model_map.svg", "day05_routing_dashboard.html",
    "day05_learning_card.md",
)


def json_text(obj: Any) -> str:
    return json.dumps(obj, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False) + "\n"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def read_fixture(path: Path) -> dict[str, Any]:
    """Accept only a tiny, explicit synthetic fixture schema."""
    if not path.is_file() or path.is_symlink() or path.stat().st_size > 65536:
        raise ValueError("Fixture missing, symlinked, or too large")
    obj = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(obj, dict) or frozenset(obj) != TOP_KEYS:
        raise ValueError("Unexpected top-level fixture fields")
    if (obj["data_kind"] != DATA_KIND or obj["assembly"] != ASSEMBLY
            or obj["schema_version"] != SCHEMA_VERSION):
        raise ValueError("Only explicit toy educational fixtures are accepted")
    cases = obj["cases"]
    if not isinstance(cases, list) or not 1 <= len(cases) <= 20:
        raise ValueError("Expected 1–20 synthetic cases")
    seen: set[str] = set()
    for case in cases:
        if not isinstance(case, dict) or frozenset(case) != CASE_KEYS:
            raise ValueError("Unexpected case fields, including possible private identifiers")
        case_id = case["case_id"]
        variant = case["variant"]
        if not isinstance(case_id, str) or not TOY_ID_RE.fullmatch(case_id):
            raise ValueError("Invalid toy case identifier")
        if case_id in seen:
            raise ValueError("Duplicate case ID")
        seen.add(case_id)
        if not isinstance(variant, str):
            raise ValueError("Variant must be a toy string")
        match = TOY_VARIANT_RE.fullmatch(variant)
        if not match or match.group(2) == match.group(3):
            raise ValueError("Only valid chrToy SNVs with REF != ALT are allowed")
        if not isinstance(case["consequence"], str) or case["consequence"] not in CONSEQUENCES:
            raise ValueError("Unknown synthetic consequence label")
        if not isinstance(case["research_question"], str) or case["research_question"] not in QUESTIONS:
            raise ValueError("Unknown research question category")
    return obj


def route_case(case: dict[str, str]) -> dict[str, Any]:
    """Route a *question* to a teaching workflow, never predict effects."""
    question = case["research_question"]
    consequence = case["consequence"]
    if question in {"PROTEIN_STRUCTURE", "MOLECULAR_COMPLEX"} and consequence == "MISSENSE":
        route = "STRUCTURE_EVIDENCE_START"
        family = "ALPHAFOLD_FAMILY"
        next_evidence = ("EXPERIMENTAL_COMPLEX_AND_FUNCTION" if question == "MOLECULAR_COMPLEX"
                         else "EXPERIMENTAL_STRUCTURE_AND_FUNCTION")
    elif question in {"RNA_SPLICING", "GENE_EXPRESSION", "CHROMATIN_ACCESSIBILITY"}:
        route = "REGULATORY_EVIDENCE_START"
        family = "ALPHAGENOME_MODEL"
        next_evidence = {
            "RNA_SPLICING": "TISSUE_RELEVANT_RNA_OR_SPLICE_ASSAY",
            "GENE_EXPRESSION": "TISSUE_RELEVANT_EXPRESSION_ASSAY",
            "CHROMATIN_ACCESSIBILITY": "TISSUE_RELEVANT_CHROMATIN_ASSAY",
        }[question]
    else:
        route = "QUESTION_REVIEW_REQUIRED"
        family = "NONE_SELECTED"
        next_evidence = "CLARIFY_CONSEQUENCE_AND_RESEARCH_QUESTION"
    return {
        "case_id": case["case_id"],
        "variant": case["variant"],
        "consequence": consequence,
        "research_question": question,
        "candidate_route": route,
        "model_family_to_study": family,
        "next_evidence": next_evidence,
        "analysis_status": "PLANNED_ONLY",
        "predicted_effect": None,
        "clinical_risk": None,
        "model_called": False,
        "manual_review": "REQUIRED",
        "evidence_class": "SYNTHETIC_FIXTURE",
    }


def build_report(fixture: dict[str, Any], source_hash: str) -> dict[str, Any]:
    rows = [route_case(case) for case in fixture["cases"]]
    counts = {
        "structure_route_count": sum(r["candidate_route"] == "STRUCTURE_EVIDENCE_START" for r in rows),
        "regulatory_route_count": sum(r["candidate_route"] == "REGULATORY_EVIDENCE_START" for r in rows),
        "question_review_count": sum(r["candidate_route"] == "QUESTION_REVIEW_REQUIRED" for r in rows),
    }
    return {
        "status": "PASS",
        "mode": "OFFLINE_ONLY",
        "evidence_class": DATA_KIND,
        "assembly": ASSEMBLY,
        "input_sha256": source_hash,
        "case_count": len(rows),
        **counts,
        "model_predictions_made": 0,
        "real_api_calls_this_run": 0,
        "clinical_interpretation": False,
        "notices": [
            "This is a question-to-evidence teaching router, not AlphaFold or AlphaGenome.",
            "AlphaFold 3 can model molecular complexes, not just single-protein structures.",
            "AlphaGenome regulatory tracks are not a complete time-course simulation of development.",
            "Sequence annotation or model suggestions are not causal or clinical evidence.",
        ],
        "rows": rows,
    }


def render_csv(rows: list[dict[str, Any]]) -> str:
    fields = (
        "case_id", "variant", "consequence", "research_question",
        "candidate_route", "model_family_to_study", "next_evidence",
        "analysis_status", "predicted_effect", "clinical_risk", "manual_review",
        "evidence_class",
    )
    buff = io.StringIO(newline="")
    w = csv.DictWriter(buff, fieldnames=fields, lineterminator="\n", extrasaction="ignore")
    w.writeheader()
    for row in rows:
        w.writerow({key: ("" if row[key] is None else row[key]) for key in fields})
    return buff.getvalue()


def render_svg(report: dict[str, Any]) -> str:
    """A static, self-contained and editable 16:9 educational diagram."""
    a, b, c = (report[k] for k in (
        "structure_route_count", "regulatory_route_count", "question_review_count"))
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="1920" height="1080" viewBox="0 0 1920 1080" role="img" aria-label="Day 05 educational overview: AlphaFold protein-structure questions, AlphaGenome gene-regulation questions, and manual review. Synthetic example only.">
<rect width="1920" height="1080" fill="#E5DEDA"/>
<defs><linearGradient id="surface" x1="0" y1="0" x2="1" y2="1"><stop stop-color="#ffffff" stop-opacity=".88"/><stop offset="1" stop-color="#f5eefa" stop-opacity=".72"/></linearGradient></defs>
<text x="110" y="95" font-size="39" font-family="DejaVu Sans,Arial,sans-serif" fill="#453b44" font-weight="bold">AlphaGenome Atlas | Day 05</text>
<text x="110" y="157" font-size="45" font-family="DejaVu Sans,Arial,sans-serif" fill="#302c36" font-weight="bold">Two questions, two kinds of evidence</text>
<text x="110" y="201" font-size="24" font-family="DejaVu Sans,Arial,sans-serif" fill="#675f68">Google AI Studio MVP → Python question router → evidence review</text>
<rect x="95" y="250" width="834" height="505" rx="34" fill="url(#surface)" stroke="#bc9b9b" stroke-width="3"/>
<rect x="990" y="250" width="834" height="505" rx="34" fill="url(#surface)" stroke="#8aa2bb" stroke-width="3"/>
<rect x="129" y="288" width="250" height="50" rx="20" fill="#eed6d4"/>
<text x="154" y="323" font-size="29" font-family="DejaVu Sans,Arial,sans-serif" font-weight="bold" fill="#76565b">AlphaFold family</text>
<rect x="1024" y="288" width="302" height="50" rx="20" fill="#dceaf1"/>
<text x="1050" y="323" font-size="29" font-family="DejaVu Sans,Arial,sans-serif" font-weight="bold" fill="#3c657e">AlphaGenome model</text>
<text x="143" y="408" font-size="35" font-family="DejaVu Sans,Arial,sans-serif" fill="#443942" font-weight="bold">Protein structure</text>
<text x="143" y="458" font-size="27" font-family="DejaVu Sans,Arial,sans-serif" fill="#665760">Amino-acid sequence</text>
<path d="M 145 483 L 815 483" stroke="#bc989e" stroke-width="3" stroke-dasharray="12 12"/>
<path d="M 230 580 C 290 486 332 620 401 545 S 510 658 597 547 S 680 641 752 535" fill="none" stroke="#bc8288" stroke-width="12" stroke-linecap="round"/>
<text x="143" y="685" font-size="26" font-family="DejaVu Sans,Arial,sans-serif" fill="#544950">Question: shape / molecular complex?</text>
<text x="1038" y="408" font-size="35" font-family="DejaVu Sans,Arial,sans-serif" fill="#384958" font-weight="bold">Gene regulation</text>
<text x="1038" y="458" font-size="27" font-family="DejaVu Sans,Arial,sans-serif" fill="#526c80">DNA sequence and molecular readouts</text>
<path d="M 1040 483 L 1708 483" stroke="#8aa6b8" stroke-width="3" stroke-dasharray="12 12"/>
<path d="M 1080 565 C 1140 480 1210 642 1270 565 S 1400 483 1470 565 S 1590 645 1706 530" fill="none" stroke="#7ca6b4" stroke-width="12" stroke-linecap="round"/>
<path d="M 1080 620 C 1140 548 1210 690 1270 610 S 1400 538 1470 614 S 1600 700 1706 593" fill="none" stroke="#c6a0c2" stroke-width="8" stroke-linecap="round"/>
<text x="1038" y="685" font-size="26" font-family="DejaVu Sans,Arial,sans-serif" fill="#435b6b">Question: RNA / splicing / chromatin?</text>
<rect x="95" y="803" width="1729" height="142" rx="31" fill="#f5f0ef" stroke="#cbbebc" stroke-width="2"/>
<text x="132" y="859" font-size="29" font-family="DejaVu Sans,Arial,sans-serif" fill="#4f4149" font-weight="bold">Seven synthetic question cards</text>
<text x="132" y="909" font-size="29" font-family="DejaVu Sans,Arial,sans-serif" fill="#654d54">Structure: {a}   |   Regulation: {b}   |   Review: {c}</text>
<text x="111" y="1008" font-size="24" font-family="DejaVu Sans,Arial,sans-serif" fill="#645961">SYNTHETIC EDUCATIONAL DEMO  •  Model calls: 0  •  No clinical predictions  •  Human review required</text>
</svg>'''


def render_html(report: dict[str, Any]) -> str:
    rows = []
    for r in report["rows"]:
        cells = (r[key] for key in (
            "case_id", "variant", "consequence", "research_question", "candidate_route",
            "model_family_to_study", "next_evidence", "manual_review"))
        rows.append("<tr>" + "".join(f"<td>{html.escape(str(v))}</td>" for v in cells) + "</tr>")
    headers = ("Toy case", "Toy variant", "Consequence", "Research question", "Candidate route",
               "Model family to study", "Next evidence", "Review")
    th = "".join(f"<th>{html.escape(x)}</th>" for x in headers)
    return f'''<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>Day 05 | Protein structure vs. genome regulation</title>
<style>body{{font-family:system-ui,Arial,sans-serif;background:#E5DEDA;color:#322e35;margin:0;padding:35px}}
main{{max-width:1300px;margin:auto;background:#fffdfa;border-radius:20px;padding:35px}}
table{{border-collapse:collapse;width:100%;font-size:14px}}th,td{{padding:11px;border-bottom:1px solid #ded8d6;text-align:left;vertical-align:top}}
th{{background:#f0e7e4}}.note{{padding:14px;background:#eeeaf0;border-left:5px solid #947fa2}}
.cards{{display:flex;flex-wrap:wrap;gap:15px;margin:18px 0}}.card{{background:#f4eded;padding:22px;border-radius:15px}}
small{{color:#59565d}}@media(max-width:700px){{main{{padding:18px;overflow:auto}}}}</style></head>
<body><main><p>AlphaGenome Atlas | Day 05</p><h1>Protein structure vs. genomic regulation</h1>
<p class="note"><strong>SYNTHETIC EDUCATIONAL DEMO.</strong> No AlphaFold or AlphaGenome model was executed. No API, no patient data and no clinical interpretation. All seven cases are fictional.</p>
<p>Step 1: Google AI Studio mock dashboard → Step 2: Python evidence router → Step 3: review missing evidence.</p>
<div class="cards"><div class="card"><strong>Structure route</strong><h2>{report['structure_route_count']}</h2></div>
<div class="card"><strong>Regulatory route</strong><h2>{report['regulatory_route_count']}</h2></div>
<div class="card"><strong>Question review</strong><h2>{report['question_review_count']}</h2></div></div>
<table><thead><tr>{th}</tr></thead><tbody>{''.join(rows)}</tbody></table>
<h2>Evidence boundary</h2><p>AlphaFold 3 also studies molecular complexes. AlphaGenome predicts DNA-related molecular readouts, not a complete time-course simulation. The route is a teaching suggestion only.</p>
<small>Input SHA-256: {html.escape(report['input_sha256'])} | real_api_calls_this_run=0 | REVIEW_REQUIRED</small>
</main></body></html>'''


def render_markdown(report: dict[str, Any]) -> str:
    lines = [
        "# Day 05 | Question-to-evidence learning card", "",
        "**SYNTHETIC EDUCATIONAL DEMO — NOT REAL MODEL OUTPUT**", "",
        "AlphaFold-family tools help investigate protein or complex structures; AlphaGenome helps investigate DNA regulatory readouts.",
        "Neither is a direct clinical diagnosis. Model execution and API calls: **0**.", "",
        "| Toy case | Question | Candidate route | Evidence to seek |",
        "|---|---|---|---|",
    ]
    for row in report["rows"]:
        lines.append("| " + " | ".join(str(row[k]) for k in (
            "case_id", "research_question", "candidate_route", "next_evidence")) + " |")
    lines += ["", "## Reader check", "", 
              "- Does the input require DNA sequence or amino-acid sequence?",
              "- Is the research question about structure, expression, splicing or chromatin?",
              "- Is the evidence a prediction, actual experiment, or only an intended next step?",
              "- Was any model run? **NO**. Is manual review needed? **YES**.", ""]
    return "\n".join(lines)


def output_payloads(report: dict[str, Any]) -> dict[str, bytes]:
    return {
        "day05_routing_report.json": json_text(report).encode("utf-8"),
        "day05_routing_table.csv": render_csv(report["rows"]).encode("utf-8"),
        "day05_model_map.svg": render_svg(report).encode("utf-8"),
        "day05_routing_dashboard.html": render_html(report).encode("utf-8"),
        "day05_learning_card.md": render_markdown(report).encode("utf-8"),
    }


def run_lab(fixture_path: Path, output_dir: Path) -> dict[str, Any]:
    fixture = read_fixture(fixture_path)
    if output_dir.exists() or output_dir.is_symlink():
        raise FileExistsError("Output exists; use a fresh --output-dir")
    if output_dir.parent.is_symlink():
        raise ValueError("Output parent may not be a symlink")
    report = build_report(fixture, sha256_file(fixture_path))
    payloads = output_payloads(report)
    output_dir.mkdir(parents=True, exist_ok=False)
    for name, content in payloads.items():
        (output_dir / name).write_bytes(content)
    manifest = {
        "status": "PASS", "evidence_class": DATA_KIND,
        "source_fixture_name": fixture_path.name, "source_fixture_sha256": sha256_file(fixture_path),
        "files": {name: sha256_bytes(payloads[name]) for name in ARTIFACTS},
        "real_api_calls_this_run": 0,
    }
    (output_dir / "day05_sha256_manifest.json").write_text(json_text(manifest), encoding="utf-8")
    return {
        "status": "PASS", "mode": "OFFLINE_ONLY",
        "case_count": report["case_count"],
        "structure_route_count": report["structure_route_count"],
        "regulatory_route_count": report["regulatory_route_count"],
        "question_review_count": report["question_review_count"],
        "artifact_count": len(ARTIFACTS), "output_dir": str(output_dir),
        "model_predictions_made": 0,
        "real_api_calls_this_run": 0,
    }


def verify_lab(fixture_path: Path, output_dir: Path) -> dict[str, Any]:
    fixture = read_fixture(fixture_path)
    manifest_path = output_dir / "day05_sha256_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict) or frozenset(manifest) != frozenset({
        "status", "evidence_class", "source_fixture_name", "source_fixture_sha256",
        "files", "real_api_calls_this_run",
    }):
        raise ValueError("Invalid manifest schema")
    expected_source_sha = sha256_file(fixture_path)
    if manifest["source_fixture_name"] != fixture_path.name or manifest["source_fixture_sha256"] != expected_source_sha:
        raise ValueError("Fixture path/name or SHA-256 does not match manifest")
    if (manifest["status"] != "PASS" or manifest["evidence_class"] != DATA_KIND
            or manifest["real_api_calls_this_run"] != 0):
        raise ValueError("Invalid manifest provenance fields")
    if not isinstance(manifest["files"], dict) or set(manifest["files"]) != set(ARTIFACTS):
        raise ValueError("Artifact allowlist mismatch")
    expected_payloads = output_payloads(build_report(fixture, expected_source_sha))
    for name in ARTIFACTS:
        path = output_dir / name
        if not path.is_file() or path.is_symlink():
            raise ValueError(f"Missing / symlinked artifact: {name}")
        expected_bytes = expected_payloads[name]
        actual_bytes = path.read_bytes()
        if (sha256_bytes(actual_bytes) != manifest["files"][name]
                or actual_bytes != expected_bytes):
            raise ValueError(f"Artifact integrity failure: {name}")
    if set(p.name for p in output_dir.iterdir()) != set(ARTIFACTS) | {"day05_sha256_manifest.json"}:
        raise ValueError("Unexpected output artifact")
    return {
        "status": "PASS", "mode": "OFFLINE_ONLY", "checksums": "MATCH",
        "files_verified": len(ARTIFACTS), "case_count": len(fixture["cases"]),
        "model_predictions_made": 0, "real_api_calls_this_run": 0,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Day 05 synthetic model-question router; offline only")
    parser.add_argument("command", choices=("run", "verify"))
    parser.add_argument("--mode", required=True, choices=("offline",))
    parser.add_argument("--fixture", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        result = (run_lab(args.fixture, args.output_dir) if args.command == "run"
                  else verify_lab(args.fixture, args.output_dir))
    except (ValueError, FileNotFoundError, FileExistsError, OSError, json.JSONDecodeError) as exc:
        print(f"Day 05 validation error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
