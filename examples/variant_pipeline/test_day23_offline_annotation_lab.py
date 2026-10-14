"""Day23 tests for the original Day21 synthetic -> structured annotation bridge."""
from __future__ import annotations

import copy
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import day23_offline_annotation_lab as d23
from day21_offline_vcf_review_lab import (
    build_review_packet, check_qc, load_toy_reference,
    load_toy_scores, parse_synthetic_vcf, parse_toy_bed,
)

BASE = Path(__file__).resolve().parent
FIXTURE = {
    "vcf": BASE / "day21_synthetic.vcf",
    "reference": BASE / "day21_toy_reference.json",
    "regions": BASE / "day21_synthetic_regions.bed",
    "scores": BASE / "day21_synthetic_scores.json",
}
UPSTREAM_SHA = {
    "vcf": "adfc480a659e581a4d69f6b7b01890016ddc6ad0d6de6e709521fc6c2804800d",
    "reference": "43e8cda00040afb123073f628b1249826f37d4fb4732aed578bf64ab28f9cd5a",
    "regions": "219e3264fc7b05d0d3c5ee8af57e07d6cb80df4e7618bcf11be4af2b7640a910",
    "scores": "77f8ed4d642a12f1c04df271abeb5e54a844a0812279ebb26fe3c9f19fc0ee75",
}


def source_packet():
    return build_review_packet(
        FIXTURE["vcf"].read_text(encoding="utf-8"),
        load_toy_reference(FIXTURE["reference"].read_text(encoding="utf-8")),
        parse_toy_bed(FIXTURE["regions"].read_text(encoding="utf-8")),
        load_toy_scores(FIXTURE["scores"].read_text(encoding="utf-8")),
    )


def input_provenance():
    return {key: {"filename": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            for key, path in FIXTURE.items()}


class Day23ScoreAndCoordinateTests(unittest.TestCase):
    def test_none_is_missing(self):
        self.assertEqual(d23.score_state(None), "MISSING_SCORE")

    def test_zero_is_not_missing(self):
        self.assertEqual(d23.score_state(0.0), "SYNTHETIC_SCORE")

    def test_negative_number_is_still_present(self):
        self.assertEqual(d23.score_state(-0.17), "SYNTHETIC_SCORE")

    def test_bool_is_not_numeric_score(self):
        with self.assertRaises(ValueError):
            d23.score_state(True)

    def test_non_numeric_score_rejected(self):
        with self.assertRaises(ValueError):
            d23.score_state("0.42")

    def test_nonfinite_scores_rejected(self):
        for value in (float("nan"), float("inf"), -float("inf")):
            with self.subTest(value=str(value)), self.assertRaises(ValueError):
                d23.score_state(value)

    def test_split_variant(self):
        self.assertEqual(d23.split_variant("chr1:110:A>C"), ("chr1", 110, "A", "C"))

    def test_bad_variant_rejected(self):
        for value in ("chrToy:100:A>G", "chr1:0:A>G", "chr1:100:A>A", "chr1:100:a>G", "=HYPERLINK(1)"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                d23.split_variant(value)


class Day23UpstreamParityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = source_packet()
        cls.reference = load_toy_reference(FIXTURE["reference"].read_text(encoding="utf-8"))
        cls.bed = parse_toy_bed(FIXTURE["regions"].read_text(encoding="utf-8"))

    def test_original_fixture_hashes_unchanged(self):
        for label, digest in UPSTREAM_SHA.items():
            self.assertEqual(hashlib.sha256(FIXTURE[label].read_bytes()).hexdigest(), digest)

    def test_seven_vcf_lines_expand_to_eight(self):
        rows = parse_synthetic_vcf(FIXTURE["vcf"].read_text(encoding="utf-8"))
        self.assertEqual(len(rows), 8)
        self.assertEqual(sum(row["pos_1based"] == 110 for row in rows), 2)

    def test_source_qc_counts(self):
        self.assertEqual([self.source[k] for k in
            ("candidate_count", "selected_count", "review_required_qc_count", "rejected_count", "missing_score_count")],
            [8, 3, 2, 3, 1])

    def test_qc_order_unchanged(self):
        self.assertEqual(self.source["qc_order"], ["FILTER", "BED", "REF", "INDEL", "MISSING_QUAL", "LOW_QUAL"])

    def test_baseline_reason_first_deciding_only(self):
        reason = {r["variant_id"]: r["reason"] for r in self.source["qc_rows"]}
        self.assertEqual(reason["chr1:120:A>T"], "OUTSIDE_TEACHING_BED")
        self.assertEqual(reason["chr1:140:G>T"], "OUTSIDE_TEACHING_BED")

    def test_matches_day21_qc_for_every_candidate(self):
        rows = parse_synthetic_vcf(FIXTURE["vcf"].read_text(encoding="utf-8"))
        by_id = {r["variant_id"]: r for r in self.source["qc_rows"]}
        for row in rows:
            with self.subTest(variant=row["variant_id"]):
                self.assertEqual(check_qc(row, self.reference, self.bed),
                                 (by_id[row["variant_id"]]["qc_status"], by_id[row["variant_id"]]["reason"]))


class Day23ContractRejectionTests(unittest.TestCase):
    def setUp(self):
        self.report = source_packet()

    def test_baseline_synthetic_packet_accepted(self):
        self.assertEqual(len(d23.validate_day21_packet(self.report)), 8)

    def test_not_dict_rejected(self):
        with self.assertRaises(ValueError):
            d23.validate_day21_packet([])

    def test_not_synthetic_rejected(self):
        self.report["evidence_class"] = "HISTORICAL_LIVE_OUTPUT"
        with self.assertRaises(ValueError):
            d23.validate_day21_packet(self.report)

    def test_wrong_assembly_rejected(self):
        self.report["assembly"] = "GRCh38"
        with self.assertRaises(ValueError):
            d23.validate_day21_packet(self.report)

    def test_live_mode_rejected(self):
        self.report["mode"] = "LIVE"
        with self.assertRaises(ValueError):
            d23.validate_day21_packet(self.report)

    def test_api_called_flag_rejected(self):
        self.report["atlas_query_performed"] = True
        with self.assertRaises(ValueError):
            d23.validate_day21_packet(self.report)

    def test_clinical_interpretation_flag_rejected(self):
        self.report["clinical_interpretation"] = True
        with self.assertRaises(ValueError):
            d23.validate_day21_packet(self.report)

    def test_real_api_count_nonzero_rejected(self):
        self.report["real_api_calls_this_run"] = 1
        with self.assertRaises(ValueError):
            d23.validate_day21_packet(self.report)

    def test_not_fixture_score_source_rejected(self):
        self.report["score_source"] = "OFFICIAL_AVI"
        with self.assertRaises(ValueError):
            d23.validate_day21_packet(self.report)

    def test_wrong_qc_order_rejected(self):
        self.report["qc_order"] = list(reversed(self.report["qc_order"]))
        with self.assertRaises(ValueError):
            d23.validate_day21_packet(self.report)

    def test_duplicate_qc_variant_rejected(self):
        self.report["qc_rows"].append(copy.deepcopy(self.report["qc_rows"][0]))
        with self.assertRaises(ValueError):
            d23.validate_day21_packet(self.report)

    def test_duplicate_annotation_rejected(self):
        self.report["review_packet"].append(copy.deepcopy(self.report["review_packet"][0]))
        with self.assertRaises(ValueError):
            d23.validate_day21_packet(self.report)

    def test_annotation_of_rejected_variant_rejected(self):
        self.report["review_packet"][0]["variant_id"] = "chr1:120:A>T"
        with self.assertRaises(ValueError):
            d23.validate_day21_packet(self.report)

    def test_missing_selected_annotation_rejected(self):
        self.report["review_packet"].pop()
        with self.assertRaises(ValueError):
            d23.validate_day21_packet(self.report)

    def test_malformed_count_rejected(self):
        self.report["selected_count"] = 4
        with self.assertRaises(ValueError):
            d23.validate_day21_packet(self.report)

    def test_score_status_mismatch_rejected(self):
        self.report["review_packet"][2]["score_status"] = "SYNTHETIC_SCORE"
        with self.assertRaises(ValueError):
            d23.validate_day21_packet(self.report)

    def test_missing_score_count_mismatch_rejected(self):
        self.report["missing_score_count"] = 0
        with self.assertRaises(ValueError):
            d23.validate_day21_packet(self.report)

    def test_invalid_score_value_rejected(self):
        self.report["review_packet"][0]["synthetic_score"] = float("inf")
        with self.assertRaises(ValueError):
            d23.validate_day21_packet(self.report)

    def test_review_status_cannot_be_approved(self):
        self.report["review_packet"][0]["review_status"] = "APPROVED"
        with self.assertRaises(ValueError):
            d23.validate_day21_packet(self.report)

    def test_qc_reason_status_mismatch_rejected(self):
        self.report["qc_rows"][0]["reason"] = "=1+1"
        with self.assertRaises(ValueError):
            d23.validate_day21_packet(self.report)

    def test_bad_source_line_rejected(self):
        self.report["qc_rows"][0]["source_line"] = 0
        with self.assertRaises(ValueError):
            d23.validate_day21_packet(self.report)

    def test_bad_qual_rejected(self):
        self.report["qc_rows"][0]["qual"] = float("nan")
        with self.assertRaises(ValueError):
            d23.validate_day21_packet(self.report)

    def test_provenance_requires_sha(self):
        with self.assertRaises(ValueError):
            d23.annotation_schema(self.report, {"vcf": {"filename": "day21_synthetic.vcf", "sha256": "BAD"}})

    def test_provenance_rejects_path_injection(self):
        with self.assertRaises(ValueError):
            d23.annotation_schema(self.report, {"vcf": {"filename": "../x", "sha256": "0"*64}})


class Day23AnnotationOutputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.schema = d23.annotation_schema(source_packet(), input_provenance())

    def test_schema_columns_exact(self):
        self.assertEqual(cls_fields := self.schema["annotation_columns"], list(d23.ANNOTATION_FIELDS))
        self.assertIn("score_status", cls_fields)

    def test_schema_counts(self):
        self.assertEqual({x: self.schema[x] for x in (
            "candidate_count", "selected_count", "annotated_count", "missing_score_count",
            "review_required_qc_count", "rejected_count")},
            {"candidate_count": 8, "selected_count": 3, "annotated_count": 2, "missing_score_count": 1,
             "review_required_qc_count": 2, "rejected_count": 3})

    def test_three_selected_rows_include_null(self):
        self.assertEqual([r["synthetic_score"] for r in self.schema["annotations"]], [0.42, -0.17, None])
        self.assertEqual([r["score_status"] for r in self.schema["annotations"]],
                         ["SYNTHETIC_SCORE", "SYNTHETIC_SCORE", "MISSING_SCORE"])

    def test_all_require_human_review(self):
        self.assertEqual({r["review_status"] for r in self.schema["annotations"]}, {"REVIEW_REQUIRED"})

    def test_no_real_avi_or_api(self):
        self.assertIs(self.schema["official_avi"], False)
        self.assertIs(self.schema["atlas_query_performed"], False)
        self.assertEqual(self.schema["real_api_calls_this_run"], 0)

    def test_valid_zero_value_keeps_nonmissing(self):
        packet = source_packet()
        packet["review_packet"][2]["synthetic_score"] = 0.0
        packet["review_packet"][2]["score_status"] = "SYNTHETIC_SCORE"
        packet["missing_score_count"] = 0
        report = d23.annotation_schema(packet, input_provenance())
        self.assertEqual(report["missing_score_count"], 0)
        self.assertEqual(report["annotated_count"], 3)
        self.assertEqual(report["annotations"][2]["synthetic_score"], 0.0)

    def test_html_escape_untrusted_content(self):
        report = copy.deepcopy(self.schema)
        report["qc_rows"][0]["variant_id"] = '<script>alert("x")</script>'
        html = d23.render_html(report)
        self.assertNotIn('<script>alert("x")</script>', html)
        self.assertIn('&lt;script&gt;', html)

    def test_json_csv_html_actual_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary = d23.emit_artifacts(self.schema, tmp)
            self.assertEqual(summary["status"], "PASS")
            self.assertEqual(summary["real_api_calls_this_run"], 0)
            for label in ("json_file", "csv_file", "html_file"):
                self.assertTrue(Path(summary[label]).is_file())
            report = json.loads(Path(summary["json_file"]).read_text(encoding="utf-8"))
            self.assertEqual(report["annotations"][2]["synthetic_score"], None)
            with Path(summary["csv_file"]).open(encoding="utf-8", newline="") as f:
                rows = list(csv.DictReader(f))
            self.assertEqual(len(rows), 3)
            self.assertEqual(rows[2]["synthetic_score"], "")
            self.assertEqual(rows[2]["score_status"], "MISSING_SCORE")
            self.assertIn("SYNTHETIC", Path(summary["html_file"]).read_text(encoding="utf-8"))

    def test_deterministic_json_and_html(self):
        with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
            p1, p2 = d23.emit_artifacts(self.schema, a), d23.emit_artifacts(self.schema, b)
            for key in ("json_file", "csv_file", "html_file"):
                self.assertEqual(Path(p1[key]).read_bytes(), Path(p2[key]).read_bytes())

    def test_default_cli_reuses_day21_fixtures(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = d23.execute_offline(FIXTURE["vcf"], FIXTURE["reference"],
                FIXTURE["regions"], FIXTURE["scores"], tmp)
            report = json.loads(Path(result["json_file"]).read_text(encoding="utf-8"))
            self.assertEqual(report["pipeline_source"], "REUSED_UNMODIFIED_DAY21_SYNTHETIC_QC")
            self.assertEqual(report["input_provenance"]["vcf"]["sha256"], UPSTREAM_SHA["vcf"])

    def test_prior_day21_json_import(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            prior = root / "day21_review_packet.json"
            prior.write_text(json.dumps(source_packet()), encoding="utf-8")
            result = d23.execute_offline(None, None, None, None, root / "output", input_day21=prior)
            report = json.loads(Path(result["json_file"]).read_text(encoding="utf-8"))
            self.assertEqual(report["input_provenance"]["day21_review_packet"]["filename"], prior.name)
            self.assertEqual(report["selected_count"], 3)

    def test_non_synthetic_prior_packet_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            prior = Path(tmp) / "input.json"
            p = source_packet()
            p["assembly"] = "GRCh38"
            prior.write_text(json.dumps(p), encoding="utf-8")
            with self.assertRaises(ValueError):
                d23.execute_offline(None, None, None, None, Path(tmp) / "out", input_day21=prior)

    def test_cli_live_is_rejected(self):
        p = subprocess.run([sys.executable, str(BASE / "day23_offline_annotation_lab.py"), "--mode", "live"],
                           text=True, capture_output=True, cwd=BASE)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("invalid choice", p.stderr)

    def test_cli_requires_mode(self):
        p = subprocess.run([sys.executable, str(BASE / "day23_offline_annotation_lab.py")],
                           text=True, capture_output=True, cwd=BASE)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("required", p.stderr)

    def test_cli_stdout_json_contract(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = subprocess.run([sys.executable, str(BASE / "day23_offline_annotation_lab.py"),
                "--mode", "offline", "--output-dir", tmp],
                text=True, capture_output=True, cwd=BASE)
            self.assertEqual(p.returncode, 0, p.stderr)
            payload = json.loads(p.stdout)
            self.assertEqual((payload["candidate_count"], payload["selected_count"],
                              payload["annotated_count"], payload["missing_score_count"]), (8, 3, 2, 1))
            self.assertEqual(payload["real_api_calls_this_run"], 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
