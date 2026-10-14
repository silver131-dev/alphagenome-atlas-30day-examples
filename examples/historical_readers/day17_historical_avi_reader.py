"""Day 17 — synthetic AVI score reader. Filename retained for migration.

Completely offline. Fictional educational values, NO AlphaGenome output.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

from fixture_helpers import FixtureError, finite_number, load_verified_fixture, provenance, save_report


FIXTURE = "day17_synthetic_avi_fixture.json"


def build_report(root: Path | None = None) -> dict:
    if root is None:
        fixture, digest = load_verified_fixture(FIXTURE, 17)
    else:
        fixture, digest = load_verified_fixture(FIXTURE, 17, root)
    avi = fixture.get("avi")
    if not isinstance(avi, dict) or set(avi) != {"raw_score", "quantile", "phred_like", "feature_contributions"}:
        raise FixtureError("AVI block missing or contains unreviewed fields")
    score = finite_number(avi.get("raw_score"), "raw_score")
    q = finite_number(avi.get("quantile"), "quantile", 0.0, 1.0)
    if q <= 0 or q >= 1:
        raise FixtureError("quantile must be strictly between zero and one")
    phred = finite_number(avi.get("phred_like"), "phred_like", 0.0)
    expected_phred = -10 * math.log10(1 - q)
    if not math.isclose(phred, expected_phred, abs_tol=0.02):
        raise FixtureError("Illustrative PHRED-like score inconsistent with quantile")
    feature_items = avi.get("feature_contributions")
    if not isinstance(feature_items, list) or not feature_items:
        raise FixtureError("Missing feature importance examples")
    features = []
    for item in feature_items:
        if not isinstance(item, dict) or set(item) != {"feature", "value"} or not str(item.get("feature", "")).startswith("TOY_"):
            raise FixtureError("Feature must be transparently synthetic")
        features.append({
            "feature": item["feature"],
            "fictional_contribution": finite_number(item.get("value"), "contribution"),
        })
    contribution_sum = round(sum(item["fictional_contribution"] for item in features), 8)
    residual = round(score - contribution_sum, 8)
    return {
        "day": 17,
        "status": "SYNTHETIC_EDUCATIONAL_DEMO",
        "variant": fixture["variant"],
        "avi_interpretation": {
            "fictional_raw_score": score,
            "fictional_quantile": q,
            "percentile_percent": round(q * 100, 5),
            "fictional_phred_like": phred,
            "phred_like_check": "CONSISTENT_WITH_-10log10(1-q)",
            "fictional_feature_contributions": features,
            "illustrative_feature_sum": contribution_sum,
            "illustrative_score_minus_feature_sum": residual,
            "toy_additivity_review": "REVIEW" if abs(residual) >= 1e-3 else "WITHIN_TOY_TOLERANCE",
            "limits": [
                "A score is not an ACMG classification or clinical risk estimate.",
                "Toy feature contributions do not establish causal mechanism or additivity.",
                "The quantile/PHRED-like teaching conversion does not validate a real model.",
            ],
        },
        "provenance": provenance(digest),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Day17 SYNTHETIC-ONLY AVI teaching reader")
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    try:
        report = build_report()
    except FixtureError as exc:
        parser.exit(2, f"DAY17_FIXTURE_QA=FAIL: {exc}\n")
    if args.verify_only:
        print("DAY17_FIXTURE_QA=PASS | SOURCE=SYNTHETIC_FIXTURE | API_CALLS=0")
        return 0
    avi = report["avi_interpretation"]
    markdown = "\n".join([
        "# Day 17 | Synthetic AVI interpretation (offline)", "",
        "**SYNTHETIC EDUCATIONAL FIXTURE — NOT REAL ALPHAGENOME API OUTPUT.**", "",
        f"Toy raw AVI score: {avi['fictional_raw_score']}",
        f"Toy quantile: {avi['fictional_quantile']} ({avi['percentile_percent']}th percentile)",
        f"Toy PHRED-like transformed value: {avi['fictional_phred_like']}",
        f"Toy feature-contribution sum: {avi['illustrative_feature_sum']}",
        f"Toy score-minus-feature residual: {avi['illustrative_score_minus_feature_sum']}",
        f"Toy consistency flag: {avi['toy_additivity_review']} (NOT model additivity validation)",
        "", "A high toy rank is not evidence for variant pathogenicity.",
        "Feature contributions are fictional and have no validated causal meaning.",
        "No API calls, human subject data, or clinical classification.",
        f"Fixture SHA-256: {report['provenance']['fixture_sha256']}",
    ])
    save_report(args.output_dir, "day17_synthetic_avi_report", report, markdown)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
