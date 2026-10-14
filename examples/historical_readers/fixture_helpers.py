"""Strict, offline-only fixture provenance and checksum validation.

The historic_reader filenames are retained ONLY for compatibility with the
proposed public tutorial file tree. No historical AlphaGenome values are used.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

SCHEMA = "alpha-atlas-synthetic-fixture-v1"
EXPECTED_FILES = frozenset({
    "day17_synthetic_avi_fixture.json",
    "day18_synthetic_splicing_fixture.json",
})
ROOT = Path(__file__).resolve().parent


class FixtureError(ValueError):
    """Fail-closed validation for a missing/modified/non-synthetic fixture."""


def _read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise FixtureError(f"Fixture JSON unreadable: {path.name}: {exc}") from exc


def finite_number(value: Any, name: str, lower: float | None = None,
                  upper: float | None = None) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise FixtureError(f"{name} must be numeric")
    number = float(value)
    if not math.isfinite(number) or (lower is not None and number < lower) or (
        upper is not None and number > upper
    ):
        raise FixtureError(f"{name} out of range")
    return number


def validate_source(source: Any) -> None:
    if not isinstance(source, dict) or source.get("kind") != "SYNTHETIC_FIXTURE":
        raise FixtureError("Not a synthetic teaching fixture")
    expected = {
        "created_for": "offline_education",
        "origin": "hand_authored",
        "model_called": False,
        "api_calls": 0,
        "real_model_output": False,
        "patient_data": False,
    }
    if set(source) != ({"kind"} | set(expected)):
        raise FixtureError("Unexpected synthetic source metadata keys")
    for field, value in expected.items():
        if type(source.get(field)) is not type(value) or source.get(field) != value:
            raise FixtureError(f"Source declaration invalid: {field}")


def validate_variant(variant: Any) -> None:
    if not isinstance(variant, dict):
        raise FixtureError("Missing synthetic variant")
    if set(variant) != {"chromosome", "position_1based", "ref", "alt", "label"}:
        raise FixtureError("Unexpected variant identity keys")
    if variant.get("chromosome") != "chrToy" or variant.get("label") != "FICTIONAL_VARIANT":
        raise FixtureError("Real or unidentified chromosome/variant not allowed")
    pos = variant.get("position_1based")
    if type(pos) is not int or pos <= 0:
        raise FixtureError("Position must be a positive 1-based integer")
    if variant.get("ref") not in ("A", "C", "G", "T") or variant.get("alt") not in ("A", "C", "G", "T"):
        raise FixtureError("Invalid REF/ALT base")
    if variant["ref"] == variant["alt"]:
        raise FixtureError("REF and ALT must differ")


def load_verified_fixture(file_name: str, expected_day: int,
                          root: Path = ROOT) -> tuple[dict[str, Any], str]:
    """Hash and source-type checks: integrity of fiction, NOT API authenticity."""
    if file_name not in EXPECTED_FILES:
        raise FixtureError("File not on fixture allowlist")
    manifest = _read_json(root / "SYNTHETIC_FIXTURE_MANIFEST.json")
    if not isinstance(manifest, dict) or manifest.get("schema") != SCHEMA:
        raise FixtureError("Manifest schema invalid")
    if manifest.get("origin") != "HAND_AUTHORED_SYNTHETIC_ONLY":
        raise FixtureError("Manifest cannot establish synthetic provenance")
    fixture_map = manifest.get("files")
    if not isinstance(fixture_map, dict) or set(fixture_map) != EXPECTED_FILES:
        raise FixtureError("Manifest fixture list incomplete or unexpected")
    hash_value = fixture_map.get(file_name)
    if not isinstance(hash_value, str) or len(hash_value) != 64:
        raise FixtureError("Invalid manifest SHA-256")
    try:
        raw = (root / file_name).read_bytes()
    except OSError as exc:
        raise FixtureError("Fixture missing") from exc
    computed = hashlib.sha256(raw).hexdigest()
    if computed != hash_value:
        raise FixtureError(f"FIXTURE_SHA256_MISMATCH: {file_name}")
    fixture = _read_json(root / file_name)
    if not isinstance(fixture, dict) or fixture.get("schema") != SCHEMA:
        raise FixtureError("Fixture schema invalid")
    if type(fixture.get("day")) is not int or fixture["day"] != expected_day:
        raise FixtureError("Fixture day mismatch")
    validate_source(fixture.get("source"))
    validate_variant(fixture.get("variant"))
    return fixture, computed


def provenance(digest: str) -> dict[str, Any]:
    return {
        "source_class": "SYNTHETIC_FIXTURE",
        "source_origin": "HAND_AUTHORED_SYNTHETIC_ONLY",
        "integrity_status": "FIXTURE_SHA256_VERIFIED",
        "fixture_sha256": digest,
        "real_api_calls_this_run": 0,
        "real_model_predictions": False,
        "reference_data_authenticity": "NOT_CLAIMED",
        "usage": "RESEARCH_EDUCATION_ONLY_NOT_FOR_CLINICAL_USE",
    }


def save_report(outdir: Path | None, stem: str, payload: dict[str, Any],
                markdown: str) -> None:
    if outdir is None:
        print(json.dumps(payload, indent=2, ensure_ascii=False))
        return
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / f"{stem}.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (outdir / f"{stem}.md").write_text(markdown + "\n", encoding="utf-8")
    print(f"{stem}=PASS | output_dir={outdir}")
