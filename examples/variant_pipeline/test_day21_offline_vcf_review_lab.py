"""Day 21 offline-only synthetic fixtures and report contract tests."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest

import day21_offline_vcf_review_lab as lab

ROOT = Path(__file__).resolve().parent


def fixtures():
    vcf = (ROOT / "day21_synthetic.vcf").read_text(encoding="utf-8")
    reference = lab.load_toy_reference((ROOT / "day21_toy_reference.json").read_text(encoding="utf-8"))
    regions = lab.parse_toy_bed((ROOT / "day21_synthetic_regions.bed").read_text(encoding="utf-8"))
    scores = lab.load_toy_scores((ROOT / "day21_synthetic_scores.json").read_text(encoding="utf-8"))
    return vcf, reference, regions, scores


class Day21ParserTests(unittest.TestCase):
    def test_seven_vcf_lines_expand_to_eight_candidates(self):
        self.assertEqual(len(lab.parse_synthetic_vcf(fixtures()[0])), 8)

    def test_multiallelic_110_has_g_and_c(self):
        v = lab.parse_synthetic_vcf(fixtures()[0])
        self.assertEqual([r["variant_id"] for r in v if r["pos_1based"] == 110], ["chr1:110:A>G", "chr1:110:A>C"])

    def test_missing_synthetic_marker_is_rejected(self):
        text = fixtures()[0].replace(lab.SOURCE_FLAG, "##source=HUMAN")
        with self.assertRaises(ValueError):
            lab.parse_synthetic_vcf(text)

    def test_missing_toy_reference_flag_is_rejected(self):
        text = fixtures()[0].replace(lab.REFERENCE_FLAG, "##reference=GRCh38")
        with self.assertRaises(ValueError):
            lab.parse_synthetic_vcf(text)

    def test_nine_column_vcf_row_rejected(self):
        text = fixtures()[0].replace("chr1\t100\t.\tA\tC\t60\tPASS\t.", "chr1\t100\t.\tA\tC\t60\tPASS\t.\tPATIENT")
        with self.assertRaises(ValueError):
            lab.parse_synthetic_vcf(text)

    def test_missing_vcf_header_rejected(self):
        text = fixtures()[0].replace("\t".join(lab.VCF_FIELDS) + "\n", "")
        with self.assertRaises(ValueError):
            lab.parse_synthetic_vcf(text)

    def test_duplicate_variant_rejected(self):
        text = fixtures()[0] + "chr1\t100\t.\tA\tC\t60\tPASS\t.\n"
        with self.assertRaises(ValueError):
            lab.parse_synthetic_vcf(text)

    def test_invalid_position_rejected(self):
        text = fixtures()[0].replace("chr1\t100\t", "chr1\t0\t")
        with self.assertRaises(ValueError):
            lab.parse_synthetic_vcf(text)

    def test_invalid_allele_rejected(self):
        text = fixtures()[0].replace("chr1\t100\t.\tA\tC\t", "chr1\t100\t.\tA\tN\t")
        with self.assertRaises(ValueError):
            lab.parse_synthetic_vcf(text)

    def test_ref_equals_alt_rejected(self):
        text = fixtures()[0].replace("chr1\t100\t.\tA\tC\t", "chr1\t100\t.\tA\tA\t")
        with self.assertRaises(ValueError):
            lab.parse_synthetic_vcf(text)

    def test_missing_qual_is_none_not_zero(self):
        rows = lab.parse_synthetic_vcf(fixtures()[0])
        self.assertIsNone(rows[-1]["qual"])

    def test_negative_qual_rejected(self):
        text = fixtures()[0].replace("chr1\t100\t.\tA\tC\t60\t", "chr1\t100\t.\tA\tC\t-1\t")
        with self.assertRaises(ValueError):
            lab.parse_synthetic_vcf(text)

    def test_nonfinite_qual_rejected(self):
        text = fixtures()[0].replace("chr1\t100\t.\tA\tC\t60\t", "chr1\t100\t.\tA\tC\tNaN\t")
        with self.assertRaises(ValueError):
            lab.parse_synthetic_vcf(text)

    def test_sample_data_column_disallowed(self):
        text = fixtures()[0].replace("\t".join(lab.VCF_FIELDS), "\t".join(lab.VCF_FIELDS) + "\tFORMAT\tSAMPLE")
        with self.assertRaises(ValueError):
            lab.parse_synthetic_vcf(text)


class Day21SourceTests(unittest.TestCase):
    def test_reference_is_only_toy(self):
        vcf, reference, regions, scores = fixtures()
        self.assertEqual(len(reference["chr1"]), 200)
        self.assertEqual(reference["chr1"], "A" * 200)

    def test_non_synthetic_reference_rejected(self):
        bad = '{"data_kind":"SYNTHETIC_FIXTURE","assembly":"GRCh38","chromosomes":{"chr1":"A"}}'
        with self.assertRaises(ValueError):
            lab.load_toy_reference(bad)

    def test_bed_zero_based_end_exclusive(self):
        _, _, regions, _ = fixtures()
        self.assertTrue(lab.in_bed("chr1", 115, regions))
        self.assertFalse(lab.in_bed("chr1", 116, regions))
        self.assertTrue(lab.in_bed("chr1", 91, regions))
        self.assertFalse(lab.in_bed("chr1", 90, regions))

    def test_bed_rejects_bad_range(self):
        with self.assertRaises(ValueError):
            lab.parse_toy_bed("#SYNTHETIC_FIXTURE\nchr1\t100\t100\n")

    def test_bed_rejects_missing_synthetic_marker(self):
        with self.assertRaises(ValueError):
            lab.parse_toy_bed("chr1\t90\t115\n")

    def test_scores_reject_boolean(self):
        bad = '{"data_kind":"SYNTHETIC_FIXTURE","assembly":"TOY_NOT_GRCH38","scores":{"chr1:100:A>C":true}}'
        with self.assertRaises(ValueError):
            lab.load_toy_scores(bad)

    def test_scores_reject_non_synthetic_assembly(self):
        bad = '{"data_kind":"SYNTHETIC_FIXTURE","assembly":"GRCh38","scores":{}}'
        with self.assertRaises(ValueError):
            lab.load_toy_scores(bad)


class Day21QcTests(unittest.TestCase):
    def setUp(self):
        self.vcf, self.reference, self.regions, self.scores = fixtures()
        self.rows = lab.parse_synthetic_vcf(self.vcf)

    def test_default_counts_and_first_reason(self):
        report = lab.build_review_packet(self.vcf, self.reference, self.regions, self.scores)
        self.assertEqual(tuple(report[k] for k in ("candidate_count", "selected_count", "review_required_qc_count", "rejected_count", "missing_score_count")), (8, 3, 2, 3, 1))
        self.assertEqual([r["reason"] for r in report["qc_rows"]], ["TEACHING_QC_PASS"] * 3 + ["OUTSIDE_TEACHING_BED", "FILTER_NOT_PASS", "OUTSIDE_TEACHING_BED", "INDEL_NOT_NORMALIZED", "MISSING_QUAL"])

    def test_filter_precedes_bed(self):
        self.assertEqual(lab.check_qc(self.rows[4], self.reference, self.regions), ("REJECTED", "FILTER_NOT_PASS"))

    def test_low_qual_outside_bed_is_reported_as_outside(self):
        self.assertEqual(lab.check_qc(self.rows[3], self.reference, self.regions), ("REJECTED", "OUTSIDE_TEACHING_BED"))

    def test_ref_mismatch_outside_bed_is_reported_as_outside(self):
        self.assertEqual(lab.check_qc(self.rows[5], self.reference, self.regions), ("REJECTED", "OUTSIDE_TEACHING_BED"))

    def test_independent_low_qual_inside_bed(self):
        row = dict(self.rows[0], pos_1based=105, qual=10.0)
        self.assertEqual(lab.check_qc(row, self.reference, self.regions), ("REJECTED", "LOW_QUAL"))

    def test_independent_ref_mismatch_inside_bed(self):
        row = dict(self.rows[0], pos_1based=105, ref="G")
        self.assertEqual(lab.check_qc(row, self.reference, self.regions), ("REJECTED", "REF_MISMATCH"))

    def test_indel_awaits_normalization(self):
        self.assertEqual(lab.check_qc(self.rows[6], self.reference, self.regions), ("REVIEW_REQUIRED", "INDEL_NOT_NORMALIZED"))

    def test_missing_qual_awaits_review(self):
        self.assertEqual(lab.check_qc(self.rows[7], self.reference, self.regions), ("REVIEW_REQUIRED", "MISSING_QUAL"))

    def test_threshold_20_is_selected(self):
        row = dict(self.rows[0], qual=20.0)
        self.assertEqual(lab.check_qc(row, self.reference, self.regions), ("SELECTED", "TEACHING_QC_PASS"))

    def test_scores_are_artificial_and_missing_not_zero(self):
        report = lab.build_review_packet(self.vcf, self.reference, self.regions, self.scores)
        self.assertEqual([p["synthetic_score"] for p in report["review_packet"]], [0.42, -0.17, None])
        self.assertEqual([p["score_status"] for p in report["review_packet"]], ["SYNTHETIC_SCORE", "SYNTHETIC_SCORE", "MISSING_SCORE"])

    def test_zero_is_not_missing_score(self):
        scores = dict(self.scores, **{"chr1:110:A>C": 0.0})
        report = lab.build_review_packet(self.vcf, self.reference, self.regions, scores)
        self.assertEqual(report["missing_score_count"], 0)
        self.assertEqual(report["review_packet"][2]["synthetic_score"], 0)
        self.assertEqual(report["review_packet"][2]["score_status"], "SYNTHETIC_SCORE")

    def test_selected_not_clinically_cleared(self):
        report = lab.build_review_packet(self.vcf, self.reference, self.regions, self.scores)
        self.assertTrue(all(p["review_status"] == "REVIEW_REQUIRED" for p in report["review_packet"]))
        self.assertFalse(report["clinical_interpretation"])

    def test_qc_review_rows_preserved_not_in_annotation_packet(self):
        report = lab.build_review_packet(self.vcf, self.reference, self.regions, self.scores)
        self.assertEqual(len(report["qc_rows"]), 8)
        self.assertEqual(len(report["review_packet"]), 3)
        self.assertNotIn("chr1:150:A>AT", {r["variant_id"] for r in report["review_packet"]})

    def test_no_real_model_or_atlas_claim(self):
        report = lab.build_review_packet(self.vcf, self.reference, self.regions, self.scores)
        self.assertEqual(report["real_api_calls_this_run"], 0)
        self.assertFalse(report["atlas_query_performed"])
        self.assertEqual(report["assembly"], lab.ASSEMBLY)


class Day21OutputTests(unittest.TestCase):
    def test_render_html_escapes_untrusted_text(self):
        report = lab.build_review_packet(*fixtures())
        report["qc_rows"][0]["variant_id"] = '<script>alert("x")</script>'
        html = lab.render_html(report)
        self.assertNotIn('<script>alert("x")</script>', html)
        self.assertIn('&lt;script&gt;', html)

    def test_html_contains_synthetic_disclaimer_and_missing_null(self):
        report = lab.build_review_packet(*fixtures())
        html = lab.render_html(report)
        self.assertIn('SYNTHETIC EDUCATIONAL DEMO', html)
        self.assertIn('MISSING_SCORE', html)
        self.assertIn('null', html)

    def test_report_json_and_html_written(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary = lab.execute_offline(ROOT / "day21_synthetic.vcf", ROOT / "day21_toy_reference.json", ROOT / "day21_synthetic_regions.bed", ROOT / "day21_synthetic_scores.json", Path(tmp) / "out")
            self.assertEqual(summary["status"], "PASS")
            parsed = json.loads(Path(summary["json_file"]).read_text(encoding="utf-8"))
            self.assertEqual(parsed["candidate_count"], 8)
            self.assertEqual(parsed["review_packet"][2]["synthetic_score"], None)
            self.assertTrue(Path(summary["html_file"]).is_file())

    def test_provenance_hashes_match_fixture_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary = lab.execute_offline(ROOT / "day21_synthetic.vcf", ROOT / "day21_toy_reference.json", ROOT / "day21_synthetic_regions.bed", ROOT / "day21_synthetic_scores.json", Path(tmp) / "out")
            report = json.loads(Path(summary["json_file"]).read_text(encoding="utf-8"))
            for name, meta in report["input_provenance"].items():
                self.assertEqual(meta["sha256"], lab.hash_file(ROOT / meta["filename"]))

    def test_live_mode_explicitly_rejected(self):
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as exc:
                lab.main(["--mode", "live"])
        self.assertEqual(exc.exception.code, 2)

    def test_cli_stdout_matches_contract(self):
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()) as buf:
            lab.main(["--mode", "offline", "--vcf", str(ROOT / "day21_synthetic.vcf"),
                      "--reference", str(ROOT / "day21_toy_reference.json"),
                      "--regions", str(ROOT / "day21_synthetic_regions.bed"),
                      "--scores", str(ROOT / "day21_synthetic_scores.json"),
                      "--output-dir", str(Path(tmp) / "out")])
        result = json.loads(buf.getvalue())
        self.assertEqual(result["selected_count"], 3)
        self.assertEqual(result["missing_score_count"], 1)
        self.assertEqual(result["real_api_calls_this_run"], 0)


if __name__ == "__main__":
    unittest.main()
