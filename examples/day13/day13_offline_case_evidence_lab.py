"""Day 13 — local-only fictional DNM1-inspired evidence interpretation lab.

Strictly synthetic teaching input, no real patient, API, variant lookup, scoring,
clinical classification, disease diagnosis or functional experiment.
Only Python standard library modules are used.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import io
import json
import re
from pathlib import Path
from typing import Any

SCHEMA = "DAY13_SYNTHETIC_EVIDENCE_V1"
KIND = "SYNTHETIC_EDUCATIONAL_CASE"
STAGES = ("PANEL", "WES", "WGS", "BLOOD_RNA", "MODEL", "VALIDATION")
EVIDENCE = {
    "PANEL": ("GENE_PANEL", "INCONCLUSIVE"),
    "WES": ("EXOME_REVIEW", "INCONCLUSIVE"),
    "WGS": ("GENOME_ANNOTATION", "TOY_VARIANT_PRESENT"),
    "BLOOD_RNA": ("RNA_CONTEXT", "INCONCLUSIVE"),
    "MODEL": ("PREDICTION_HYPOTHESIS", "HYPOTHESIS_ONLY"),
    "VALIDATION": ("INDEPENDENT_RNA_VALIDATION", "NOT_PERFORMED"),
}
FIELDS = ("step_id", "evidence_type", "source_id", "status", "observed_value", "review_question")
OUTPUTS = (
    "day13_evidence_packet.json",
    "day13_evidence_table.csv",
    "day13_hypothesis_questions.md",
    "day13_evidence_map.svg",
    "day13_review_dashboard.html",
)
ALLOWED_KEYS = {"schema_version", "data_kind", "scenario_id", "gene_symbol", "assembly",
                "variant_id", "phenotype_label", "clinical_interpretation",
                "real_api_calls", "experiment_performed", "evidence_rows"}
SAFE_ID = re.compile(r"[A-Z][A-Z0-9_]{2,48}\Z", re.ASCII)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_fixture(path: str | Path) -> tuple[dict[str, Any], str]:
    p = Path(path)
    if not p.is_file() or p.is_symlink() or p.stat().st_size > 50000:
        raise ValueError("Expected a small ordinary synthetic JSON fixture")
    data = p.read_bytes()
    try:
        record = json.loads(data)
    except (ValueError, UnicodeDecodeError) as exc:
        raise ValueError("Not valid JSON") from exc
    validate(record)
    return record, sha(data)


def validate(r: Any) -> None:
    if not isinstance(r, dict) or set(r) != ALLOWED_KEYS:
        raise ValueError("Incorrect synthetic fixture fields (patient data forbidden)")
    fixed = {
        "schema_version": SCHEMA,
        "data_kind": KIND,
        "scenario_id": "DNM1_INSPIRED_TOY_EVIDENCE",
        "gene_symbol": "DNM1",
        "assembly": "TOY_NOT_GRCH38",
        "variant_id": "chrToy:101:G>A",
        "phenotype_label": "SYNTHETIC_EARLY_SEIZURES",
        "clinical_interpretation": False,
        "real_api_calls": 0,
        "experiment_performed": False,
    }
    for k, v in fixed.items():
        if type(r.get(k)) is not type(v) or r[k] != v:
            raise ValueError(f"{k}: non-synthetic or unsupported value")
    rows = r.get("evidence_rows")
    if not isinstance(rows, list) or len(rows) != len(STAGES):
        raise ValueError("Expected precisely six evidence stages")
    if [row.get("step_id") if isinstance(row, dict) else None for row in rows] != list(STAGES):
        raise ValueError("Missing, duplicate or reordered evidence stages")
    for row in rows:
        if not isinstance(row, dict) or set(row) != set(FIELDS):
            raise ValueError("Invalid evidence row")
        step = row["step_id"]
        if row["evidence_type"] != EVIDENCE[step][0] or row["status"] != EVIDENCE[step][1]:
            raise ValueError("Unexpected evidence type or status")
        if not isinstance(row["source_id"], str) or not SAFE_ID.fullmatch(row["source_id"]):
            raise ValueError("Synthetic source ID must be safe")
        if not isinstance(row["review_question"], str) or not (8 <= len(row["review_question"]) <= 128) or "\n" in row["review_question"]:
            raise ValueError("Invalid review question")
        value = row["observed_value"]
        if step == "WGS":
            if value != "TOY_INTRONIC_CANDIDATE":
                raise ValueError("Expected synthetic WGS annotation label")
        elif value is not None:
            raise ValueError("Unmeasured or modeled values must remain JSON null")


def packet(r: dict[str, Any], fixture_sha256: str) -> dict[str, Any]:
    validate(r)
    counts = {label: sum(x["status"] == label for x in r["evidence_rows"])
              for label in ("INCONCLUSIVE", "TOY_VARIANT_PRESENT", "HYPOTHESIS_ONLY", "NOT_PERFORMED")}
    return {
        "schema_version": "DAY13_EVIDENCE_PACKET_V1", "evidence_class": KIND,
        "source_fixture": "day13_synthetic_evidence_case.json", "source_fixture_sha256": fixture_sha256,
        "gene_symbol": r["gene_symbol"], "assembly": r["assembly"],
        "variant_id": r["variant_id"], "phenotype_label": r["phenotype_label"],
        "step_count": len(STAGES), "status_counts": counts,
        "observed_values_present": 1, "model_score": None,
        "evidence_rows": r["evidence_rows"],
        "clinical_classification": None, "review_status": "REVIEW_REQUIRED",
        "experiment_performed": False, "atlas_query_performed": False,
        "real_api_calls_this_run": 0,
        "conclusion": "RESEARCH_HYPOTHESIS_ONLY",
    }


def json_bytes(obj: Any) -> bytes:
    return (json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode("utf-8")


def table_bytes(report: dict[str, Any]) -> bytes:
    mem = io.StringIO(newline="")
    writer = csv.DictWriter(mem, fieldnames=FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(report["evidence_rows"])
    return mem.getvalue().encode("utf-8")


def questions(report: dict[str, Any]) -> bytes:
    out = [
        "# Day 13｜合成 DNM1 證據閱讀卡",
        "",
        "**這不是病人紀錄或 AlphaGenome 回應。** 來源為固定 toy JSON。",
        "",
        "- 模型提示只是研究假說；沒有取得剪接模型分數。",
        "- 血液 RNA 標記 INCONCLUSIVE，不代表陰性證明或已測得零。",
        "- WGS 欄位是人工註解，既無 GRCh38 座標，也未驗證參考鹼基。",
        "- 功能檢驗 NOT_PERFORMED；沒有任何真實樣本或生物複本。",
        "",
        "## 六個關鍵審閱問題",
        "",
    ]
    for row in report["evidence_rows"]:
        out.append(f"- **{row['step_id']}**：{row['review_question']}")
    out += ["", "## 下一步", "", "Day 14 規劃 WT、ALT 與回復對照，先寫下能推翻剪接假說的觀測。", "",
            "**Clinical classification:** null；**Review:** REQUIRED；**Real API calls:** 0。", ""]
    return "\n".join(out).encode("utf-8")


def svg_bytes(report: dict[str, Any]) -> bytes:
    cols = {"INCONCLUSIVE": ("#FFF3E9", "#8F5942"),
            "TOY_VARIANT_PRESENT": ("#DAEBED", "#2E6975"),
            "HYPOTHESIS_ONLY": ("#EADFFC", "#725A9F"),
            "NOT_PERFORMED": ("#F4E3EA", "#91556A")}
    rows: list[str] = []
    for i, item in enumerate(report["evidence_rows"]):
        x, y = 92 + (i % 3) * 590, 325 + (i // 3) * 270
        bg, fg = cols[item["status"]]
        rows += [
            f'<rect x="{x}" y="{y}" width="548" height="222" rx="26" fill="{bg}" stroke="{fg}" stroke-opacity=".32"/>',
            f'<text x="{x+28}" y="{y+52}" class="stage">{html.escape(item["step_id"])}</text>',
            f'<text x="{x+28}" y="{y+98}" class="state" fill="{fg}">{html.escape(item["status"])}</text>',
            f'<text x="{x+28}" y="{y+148}" class="detail">{html.escape(item["evidence_type"])}</text>',
            f'<text x="{x+28}" y="{y+186}" class="detail">來源：人工合成資料</text>',
        ]
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="1920" height="1080" viewBox="0 0 1920 1080" role="img" aria-labelledby="title desc">
    <title id="title">Day 13｜模型線索走向可追溯證據</title>
    <desc id="desc">六個人工證據步驟：PANEL、WES、WGS、BLOOD RNA、MODEL、VALIDATION。全部是合成教學資料，未實驗也未呼叫 API。</desc>
    <style>text{{font-family:'Noto Sans CJK TC','Microsoft JhengHei',Arial,sans-serif;}}
    .h{{font-size:60px;font-weight:700;fill:#473D45;}} .sub{{font-size:27px;fill:#685D62;}}
    .stage{{font-size:38px;font-weight:700;fill:#433B47;}} .state{{font-size:25px;font-weight:700;}}
    .detail{{font-size:22px;fill:#5C5962;}}</style>
    <rect width="1920" height="1080" fill="#E5DEDA"/>
    <text class="h" x="92" y="125">Day 13｜DNM1 個案導讀</text>
    <text class="sub" x="92" y="184">從模型線索走到可追溯證據・Python 離線教學・TOY_NOT_GRCH38</text>
    <text class="sub" x="92" y="240">不是病人資料｜不含正式 AVI 分數｜沒有生物實驗</text>
    {''.join(rows)}
    <text class="sub" x="92" y="955">INCONCLUSIVE ≠ ZERO　／　HYPOTHESIS_ONLY ≠ DIAGNOSIS</text>
    <text class="sub" x="92" y="1015">Review required · Clinical classification: null · Real API calls: 0</text>
    </svg>'''
    return svg.encode("utf-8")


def html_bytes(report: dict[str, Any]) -> bytes:
    trs = "".join(
        "<tr>" + "".join(f"<td>{html.escape(str(row[k]) if row[k] is not None else 'null')}</td>"
                           for k in ("step_id", "evidence_type", "status", "observed_value", "review_question")) + "</tr>"
        for row in report["evidence_rows"]
    )
    body = f'''<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"/>
    <title>Day 13｜合成證據閱讀卡</title><style>
    body{{font-family:system-ui,'Microsoft JhengHei',sans-serif;color:#443D47;background:#E5DEDA;margin:0;padding:45px;}}
    main{{max-width:1150px;margin:auto;background:#fffaf7;border-radius:22px;padding:32px;}}
    table{{border-collapse:collapse;width:100%;font-size:14px;}}th,td{{text-align:left;padding:13px;border-bottom:1px solid #dacdd0;}}
    .warning{{background:#F7E9E9;border-left:5px solid #AC7782;padding:15px;}}
    </style></head><body><main><h1>Day 13｜DNM1 個案導讀</h1>
    <p>變異：{html.escape(report['variant_id'])}（純人工座標）／六階段證據審閱</p>
    <p class="warning">SYNTHETIC · TOY_NOT_GRCH38 · 模型與實驗資料未取得 · 不能判定致病性</p>
    <table><thead><tr><th>階段</th><th>證據種類</th><th>目前狀態</th><th>觀測值</th><th>下一個研究問題</th></tr></thead>
    <tbody>{trs}</tbody></table><p>Review status: {html.escape(report['review_status'])} ／ API calls: 0</p>
    <p>本頁為可重現的本地教育範例。不能視為個案資料或醫療建議。</p></main></body></html>'''
    return (body + "\n").encode("utf-8")


def collect(r: dict[str, Any], fixture_sha256: str) -> dict[str, bytes]:
    p = packet(r, fixture_sha256)
    return {OUTPUTS[0]: json_bytes(p), OUTPUTS[1]: table_bytes(p), OUTPUTS[2]: questions(p),
            OUTPUTS[3]: svg_bytes(p), OUTPUTS[4]: html_bytes(p)}


def run(fixture: str, output_dir: str) -> dict[str, Any]:
    r, checksum = read_fixture(fixture)
    target = Path(output_dir)
    if target.is_symlink() or target.exists():
        raise ValueError("Output directory already exists; use a fresh location")
    if not target.parent.is_dir():
        raise ValueError("Output parent must already exist")
    items = collect(r, checksum)
    target.mkdir()
    for name, content in items.items():
        (target / name).write_bytes(content)
    manifest = {"schema": "DAY13_MANIFEST_V1", "fixture_name": Path(fixture).name,
                "fixture_sha256": checksum, "artifacts": {k: sha(v) for k, v in sorted(items.items())},
                "real_api_calls_this_run": 0}
    (target / "day13_manifest.json").write_bytes(json_bytes(manifest))
    return {"status": "PASS", "mode": "OFFLINE_ONLY", "step_count": len(STAGES),
            "inconclusive_count": 3, "hypothesis_only_count": 1, "validation_performed": False,
            "artifact_count": len(items), "manifest_file": str(target / "day13_manifest.json"),
            "output_dir": str(target), "real_api_calls_this_run": 0}


def verify(fixture: str, output_dir: str) -> dict[str, Any]:
    r, checksum = read_fixture(fixture)
    folder = Path(output_dir)
    manifest_file = folder / "day13_manifest.json"
    if folder.is_symlink() or not folder.is_dir() or not manifest_file.is_file() or manifest_file.is_symlink():
        raise ValueError("Missing ordinary output directory and manifest")
    try:
        manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
    except (ValueError, UnicodeDecodeError) as exc:
        raise ValueError("Invalid manifest") from exc
    expected = collect(r, checksum)
    if not isinstance(manifest, dict) or set(manifest) != {"schema", "fixture_name", "fixture_sha256", "artifacts", "real_api_calls_this_run"}:
        raise ValueError("Invalid manifest contract")
    if manifest["schema"] != "DAY13_MANIFEST_V1" or manifest["fixture_sha256"] != checksum or manifest["fixture_name"] != Path(fixture).name:
        raise ValueError("Fixture provenance mismatch")
    if type(manifest["real_api_calls_this_run"]) is not int or manifest["real_api_calls_this_run"] != 0:
        raise ValueError("Real API count invalid")
    if set(manifest["artifacts"]) != set(OUTPUTS):
        raise ValueError("Missing or extra manifest entries")
    for name, generated in expected.items():
        p = folder / name
        if p.is_symlink() or not p.is_file():
            raise ValueError("Missing ordinary artifact")
        actual = p.read_bytes()
        if sha(actual) != manifest["artifacts"][name] or actual != generated:
            raise ValueError(f"Artifact content differs: {name}")
    if {p.name for p in folder.iterdir()} != set(OUTPUTS) | {"day13_manifest.json"}:
        raise ValueError("Unexpected output artifacts")
    return {"status": "PASS", "checksums": "MATCH", "source_fixture_sha256": checksum,
            "artifact_count": len(expected), "review_status": "REVIEW_REQUIRED", "real_api_calls_this_run": 0}


def main() -> None:
    parser = argparse.ArgumentParser(description="Synthetic-only DNM1 evidence lesson; no API")
    parser.add_argument("command", choices=("run", "verify"))
    parser.add_argument("--mode", default="offline", choices=("offline", "live"))
    parser.add_argument("--fixture", default="day13_synthetic_evidence_case.json")
    parser.add_argument("--output-dir", default="day13_output")
    args = parser.parse_args()
    if args.mode != "offline":
        parser.error("LIVE is not implemented or authorized")
    try:
        result = run(args.fixture, args.output_dir) if args.command == "run" else verify(args.fixture, args.output_dir)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
