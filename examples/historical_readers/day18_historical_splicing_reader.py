"""Day 18 — synthetic splicing REF/ALT reader. Filename retained for migration.

Completely offline. Artificial donor/acceptor and RNA-like teaching signals.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from fixture_helpers import FixtureError, finite_number, load_verified_fixture, provenance, save_report


FIXTURE = "day18_synthetic_splicing_fixture.json"


def build_report(root: Path | None = None) -> dict:
    if root is None:
        fixture, digest = load_verified_fixture(FIXTURE, 18)
        sibling, _ = load_verified_fixture("day17_synthetic_avi_fixture.json", 17)
    else:
        fixture, digest = load_verified_fixture(FIXTURE, 18, root)
        sibling, _ = load_verified_fixture("day17_synthetic_avi_fixture.json", 17, root)
    if fixture["variant"] != sibling["variant"]:
        raise FixtureError("Day17 and Day18 teaching variant identities differ")
    context = fixture.get("context")
    if not isinstance(context, dict) or context.get("label") != "FICTIONAL_NEURAL_LIKE_CONTEXT":
        raise FixtureError("Toy tissue context missing or non-synthetic")
    if set(context) != {"label", "strand", "coordinate_system"}:
        raise FixtureError("Unexpected contextual metadata fields")
    if context.get("strand") not in ("+", "-") or context.get("coordinate_system") != "1-based":
        raise FixtureError("Unclear strand or coordinate system")
    sites = fixture.get("splice_sites")
    rna = fixture.get("rna_seq_toy_signals")
    if not isinstance(sites, list) or not sites or not isinstance(rna, list) or not rna:
        raise FixtureError("Splice sites and RNA-like signals required")
    differences = []
    for entry in sites:
        if not isinstance(entry, dict) or set(entry) != {"position_1based", "site_type", "ref", "alt"} or entry.get("site_type") not in ("donor", "acceptor"):
            raise FixtureError("Invalid synthetic splice site")
        pos = entry.get("position_1based")
        if type(pos) is not int or pos <= 0:
            raise FixtureError("Invalid 1-based splice site position")
        ref = finite_number(entry.get("ref"), "splice REF probability-like", 0, 1)
        alt = finite_number(entry.get("alt"), "splice ALT probability-like", 0, 1)
        delta = round(alt - ref, 8)
        differences.append({
            "position_1based": pos, "site_type": entry["site_type"],
            "ref": ref, "alt": alt, "delta_alt_minus_ref": delta,
            "toy_direction": "decrease" if delta < 0 else "increase" if delta > 0 else "unchanged",
        })
    toy_rna = []
    for entry in rna:
        if not isinstance(entry, dict) or set(entry) != {"position_1based", "ref", "alt"} or type(entry.get("position_1based")) is not int or entry["position_1based"] <= 0:
            raise FixtureError("Invalid toy RNA-signal coordinate")
        ref = finite_number(entry.get("ref"), "toy RNA REF", 0)
        alt = finite_number(entry.get("alt"), "toy RNA ALT", 0)
        toy_rna.append({
            "position_1based": entry["position_1based"],
            "ref": ref, "alt": alt,
            "delta_alt_minus_ref": round(alt - ref, 8),
            "units": "arbitrary_fictional_units",
        })
    return {
        "day": 18,
        "status": "SYNTHETIC_EDUCATIONAL_DEMO",
        "variant": fixture["variant"], "context": context,
        "splice_site_ref_alt_comparison": differences,
        "rna_like_ref_alt_comparison": toy_rna,
        "interpretation": {
            "delta_definition": "ALT minus REF at the same 1-based coordinate and strand",
            "negative_donor_delta": "Toy ALT signal is lower; not a validated prediction of exon skipping.",
            "caution": "Artificial probability-like splicing values and RNA-like signals are neither measured RNA-seq nor official AlphaGenome predictions. Do not make clinical claims.",
        },
        "provenance": provenance(digest),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Day18 SYNTHETIC-ONLY splicing REF/ALT reader")
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    try:
        report = build_report()
    except FixtureError as exc:
        parser.exit(2, f"DAY18_FIXTURE_QA=FAIL: {exc}\n")
    if args.verify_only:
        print("DAY18_FIXTURE_QA=PASS | SOURCE=SYNTHETIC_FIXTURE | API_CALLS=0")
        return 0
    markdown = [
        "# Day 18 | Synthetic RNA splicing REF/ALT comparison (offline)", "",
        "**SYNTHETIC EDUCATIONAL FIXTURE — NOT REAL ALPHAGENOME API OUTPUT.**", "",
        "All positions are fictional 1-based chrToy positions; signal changes are ALT minus REF.", "",
        "| Kind | Toy coordinate | REF | ALT | Delta (ALT−REF) |", "|---|---:|---:|---:|---:|",
    ]
    for row in report["splice_site_ref_alt_comparison"]:
        markdown.append(f"| {row['site_type']} | {row['position_1based']} | {row['ref']} | {row['alt']} | {row['delta_alt_minus_ref']} |")
    markdown.extend([
        "", "Synthetic RNA-like levels use arbitrary units and are not measured expression.",
        "A decreasing donor site-like signal cannot establish clinical or molecular consequences.",
        "No API, real organism reference sequence, patient record, or clinical classification.",
        f"Fixture SHA-256: {report['provenance']['fixture_sha256']}",
    ])
    save_report(args.output_dir, "day18_synthetic_splicing_report", report, "\n".join(markdown))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
