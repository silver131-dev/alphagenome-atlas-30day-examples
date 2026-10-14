"""Day 16: inspect local AlphaGenome SDK readiness without requesting any API.

Two safe modes:
- offline: stdlib-only checks; does not import the AlphaGenome SDK.
- inspect-sdk: imports the installed SDK to inspect exported interfaces, but
  never creates a client, loads secret values, or performs a model request.

Use entirely synthetic/local information. NOT an authentication test.
"""
from __future__ import annotations

import argparse
import html
import importlib
from importlib import metadata
import json
import os
from pathlib import Path
import sys

EXPECTED_SDK = "0.10.0"
KEY_VARIABLES = ("ALPHA_GENOME_API_KEY", "ALPHAGENOME_API_KEY")
REPORT_NAME = "day16_readiness_report.json"
HTML_NAME = "day16_readiness_dashboard.html"


def sdk_distribution_version() -> str | None:
    """Read local package metadata; never import SDK or request the network."""
    try:
        return metadata.version("alphagenome")
    except metadata.PackageNotFoundError:
        return None


def api_interface_snapshot() -> dict[str, bool | str]:
    """Inspect importable symbols only: NO client factory invocation."""
    try:
        atlas = importlib.import_module("alphagenome.atlas.atlas")
        model = importlib.import_module("alphagenome.models.dna_client")
        genome = importlib.import_module("alphagenome.data.genome")
    except (ImportError, AttributeError, OSError) as exc:
        return {"inspection": "UNAVAILABLE", "error_class": type(exc).__name__}

    variant = getattr(genome, "Variant", None)
    output_type = getattr(model, "OutputType", None)
    result = {
        "inspection": "LOCAL_IMPORT_ONLY",
        "atlas_create_available": callable(getattr(atlas, "create", None)),
        "model_create_available": callable(getattr(model, "create", None)),
        "variant_from_str_available": callable(getattr(variant, "from_str", None)),
        "model_rna_seq_enum_available": hasattr(output_type, "RNA_SEQ"),
    }
    result["all_checked_exports_available"] = all(
        value for key, value in result.items() if key.endswith("_available")
    )
    return result


def build_report(
    *, mode: str = "offline", expected_sdk: str = EXPECTED_SDK,
    installed_sdk: str | None = None,
    key_names_present: set[str] | None = None,
    python_version: tuple[int, int, int] | None = None,
    in_virtualenv: bool | None = None,
    inspected_api: dict | None = None,
) -> dict:
    """Assemble a teaching-only report, never touching secret values."""
    if mode not in ("offline", "inspect-sdk"):
        raise ValueError("Only local offline and inspect-sdk modes are permitted")
    if not expected_sdk or not isinstance(expected_sdk, str):
        raise ValueError("Expected SDK version must be a nonempty string")
    version = python_version or tuple(sys.version_info[:3])
    if len(version) != 3 or any(not isinstance(n, int) or n < 0 for n in version):
        raise ValueError("Invalid Python version")
    active = (sys.prefix != sys.base_prefix) if in_virtualenv is None else in_virtualenv
    installed = sdk_distribution_version() if installed_sdk is None else installed_sdk
    # Only check key NAMES. Never evaluate, print, store or hash key VALUES.
    names = set(os.environ) if key_names_present is None else key_names_present
    key_states = {name: ("PRESENT_NOT_VALIDATED" if name in names else "NOT_SET_IN_PROCESS")
                  for name in KEY_VARIABLES}
    sdk_status = ("NOT_INSTALLED" if installed is None
                  else "VERSION_MATCH" if installed == expected_sdk else "VERSION_DIFFERS")
    report = {
        "report_kind": "LOCAL_EDUCATIONAL_ENVIRONMENT_CHECK",
        "mode": "OFFLINE_ONLY" if mode == "offline" else "LOCAL_SDK_IMPORT_INSPECTION",
        "status": "PASS" if version[:2] >= (3, 10) else "REVIEW_REQUIRED",
        "python_version": ".".join(str(n) for n in version),
        "python_supported": version[:2] >= (3, 10),
        "virtualenv_active": bool(active),
        "sdk": {"distribution": "alphagenome", "installed_version": installed,
                "expected_version": expected_sdk, "status": sdk_status},
        "api_key_presence_only": key_states,
        "api_client_created": False,
        "api_key_validated": False,
        "model_or_atlas_response_obtained": False,
        "real_api_calls_this_run": 0,
        "current_live_access": "NOT_TESTED",
        "clinical_use": False,
        "git_visibility_check": "NOT_PERFORMED",
        "publication_status": "PUBLICATION_ON_HOLD",
        "notice": "SDK presence and environment-variable presence are NOT live access proof.",
    }
    if mode == "inspect-sdk":
        check = inspected_api if inspected_api is not None else api_interface_snapshot()
        report["sdk_interface_check"] = check
        if sdk_status != "VERSION_MATCH" or check.get("all_checked_exports_available") is not True:
            report["status"] = "REVIEW_REQUIRED"
    return report


def render_html(report: dict) -> str:
    """Make a small standalone local HTML teaching card; escape dynamic text."""
    sdk = report["sdk"]
    def safe(value: object) -> str:
        return html.escape(str(value), quote=True)
    rows = (
        ("Python", report["python_version"]),
        ("Virtual environment", "ACTIVE" if report["virtualenv_active"] else "NOT_ACTIVE"),
        ("Installed SDK", sdk["installed_version"] or "NOT_INSTALLED"),
        ("SDK expected", sdk["expected_version"]),
        ("SDK status", sdk["status"]),
        ("Atlas / Model API calls", report["real_api_calls_this_run"]),
        ("LIVE access", report["current_live_access"]),
        ("Publication", report["publication_status"]),
    )
    table_rows = "\n".join(f"<tr><th>{safe(k)}</th><td>{safe(v)}</td></tr>" for k, v in rows)
    key_rows = "\n".join(
        f"<li>{safe(k)}: {safe(v)}</li>"
        for k, v in report["api_key_presence_only"].items()
    )
    return f"""<!doctype html>
<html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Day 16 · Local API Readiness</title>
<style>
body{{margin:0;font-family:system-ui,sans-serif;color:#29333b;background:#E5DEDA}}
main{{max-width:850px;margin:7vh auto;padding:32px;background:#fffdfb;border-radius:20px}}
h1{{margin-top:0}} table{{width:100%;border-collapse:collapse}}
th,td{{padding:12px;border-bottom:1px solid #ded7d3;text-align:left}}
th{{width:52%}} small{{color:#525d66}} .pill{{font-weight:700;color:#864a34}}
</style></head><body><main>
<small>ALPHAGENOME ATLAS · DAY 16 · LOCAL ONLY</small>
<h1>API 環境準備檢查</h1>
<p>這不是 API 登入成功證明，也沒有進行模型預測。</p>
<table>{table_rows}</table>
<h2>金鑰環境變數（僅記錄是否存在）</h2><ul>{key_rows}</ul>
<p class="pill">研究用途教學 · OFFLINE · LIVE NOT TESTED</p>
<small>{safe(report['notice'])}</small>
</main></body></html>"""


def write_outputs(output_dir: Path, report: dict) -> dict[str, str]:
    """Create a new, empty output directory; never overwrite existing files."""
    if output_dir.exists():
        raise FileExistsError(f"Output exists: {output_dir}. Choose a new --output-dir")
    output_dir.mkdir(parents=True, exist_ok=False)
    json_path = output_dir / REPORT_NAME
    html_path = output_dir / HTML_NAME
    with json_path.open("x", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=2, sort_keys=True)
        f.write("\n")
    with html_path.open("x", encoding="utf-8") as f:
        f.write(render_html(report))
    return {"json_file": str(json_path), "html_file": str(html_path)}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Offline-only AlphaGenome environment readiness")
    parser.add_argument("--mode", choices=["offline", "inspect-sdk"], default="offline")
    parser.add_argument("--expected-sdk", default=EXPECTED_SDK)
    parser.add_argument("--output-dir", default="day16_output")
    args = parser.parse_args(argv)
    report = build_report(mode=args.mode, expected_sdk=args.expected_sdk)
    try:
        paths = write_outputs(Path(args.output_dir), report)
    except (FileExistsError, OSError) as e:
        parser.error(str(e))
    print(json.dumps({
        "status": report["status"],
        "mode": report["mode"],
        "sdk_status": report["sdk"]["status"],
        "live_access": report["current_live_access"],
        "real_api_calls_this_run": 0,
        **paths,
    }, ensure_ascii=False, sort_keys=True))
    return 0 if report["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
