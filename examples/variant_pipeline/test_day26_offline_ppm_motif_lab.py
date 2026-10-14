"""Day26 synthetic PPM-like motif lab; unittest, no network or SDK."""
import copy
import csv
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

import day26_offline_ppm_motif_lab as m

ROOT = Path(__file__).resolve().parent
FIXTURE = ROOT / "day26_synthetic_motif_context.json"
DAY25 = ROOT / "day25_locus_summary.json"


def input_data():
    return json.loads(FIXTURE.read_text(encoding="utf-8"))


def altered(**kwargs):
    obj = input_data()
    obj.update(kwargs)
    return obj


def produce():
    return m.build_report(input_data(), fixture_sha256=m.file_sha256(FIXTURE),
                          day25=m.load_json(DAY25), day25_sha256=m.file_sha256(DAY25))


class FixtureContractTests(unittest.TestCase):
    def test_original_toy_sequence_is_valid(self):
        self.assertEqual(m.validate_fixture(input_data())["variant"], "chr1:110:A>G")

    def test_variant_ref_is_A_at_position_110(self):
        ctx = m.validate_fixture(input_data())
        self.assertEqual(ctx["ref"][110 - ctx["start"]], "A")

    def test_alt_sequence_is_A_to_G(self):
        ctx = m.validate_fixture(input_data())
        self.assertEqual(ctx["alt"][110 - ctx["start"]], "G")
        self.assertEqual(sum(a != b for a, b in zip(ctx["ref"], ctx["alt"])), 1)

    def test_reject_missing_evidence_class(self):
        obj = input_data(); obj.pop("evidence_class")
        with self.assertRaises(ValueError): m.validate_fixture(obj)

    def test_reject_live_mode(self):
        with self.assertRaises(ValueError): m.validate_fixture(altered(mode="LIVE"))

    def test_reject_grch38_as_fake_assembly(self):
        with self.assertRaises(ValueError): m.validate_fixture(altered(assembly="GRCh38"))

    def test_reject_fake_model_prediction_true(self):
        with self.assertRaises(ValueError): m.validate_fixture(altered(model_prediction=True))

    def test_reject_official_atlas_true(self):
        with self.assertRaises(ValueError): m.validate_fixture(altered(official_motif_atlas=True))

    def test_reject_clinical_interpretation_true(self):
        with self.assertRaises(ValueError): m.validate_fixture(altered(clinical_interpretation=True))

    def test_reject_api_calls_positive(self):
        with self.assertRaises(ValueError): m.validate_fixture(altered(real_api_calls_this_run=1))

    def test_reject_non_positive_start(self):
        with self.assertRaises(ValueError): m.validate_fixture(altered(start_1based=0))

    def test_reject_wrong_chromosome(self):
        with self.assertRaises(ValueError): m.validate_fixture(altered(chrom="chr2"))

    def test_reject_ref_mismatch_in_toy_sequence(self):
        with self.assertRaises(ValueError): m.validate_fixture(altered(variant_id="chr1:110:C>G"))

    def test_reject_variant_outside_sequence(self):
        with self.assertRaises(ValueError): m.validate_fixture(altered(variant_id="chr1:130:A>G"))

    def test_reject_non_snv(self):
        with self.assertRaises(ValueError): m.validate_fixture(altered(variant_id="chr1:110:A>GG"))

    def test_reject_same_ref_and_alt(self):
        with self.assertRaises(ValueError): m.validate_fixture(altered(variant_id="chr1:110:A>A"))

    def test_reject_invalid_dna_character(self):
        with self.assertRaises(ValueError): m.validate_fixture(altered(ref_sequence="TACCNGACATGGTCAAC"))

    def test_reject_mixed_strand(self):
        with self.assertRaises(ValueError): m.validate_fixture(altered(scan_strand="+/-"))

    def test_reject_probability_rows_not_sum_one(self):
        with self.assertRaises(ValueError): m.validate_fixture(altered(preferred_probability=0.6))

    def test_reject_bool_probability(self):
        with self.assertRaises(ValueError): m.validate_fixture(altered(preferred_probability=True))

    def test_reject_zero_probability(self):
        with self.assertRaises(ValueError): m.validate_fixture(altered(other_probability=0))

    def test_reject_duplicate_motif_ids(self):
        motifs = input_data()["motifs"]
        with self.assertRaises(ValueError): m.validate_fixture(altered(motifs=[motifs[0], motifs[0]]))

    def test_reject_empty_motifs(self):
        with self.assertRaises(ValueError): m.validate_fixture(altered(motifs=[]))

    def test_reject_invalid_motif_id(self):
        with self.assertRaises(ValueError): m.validate_fixture(altered(motifs=[{"motif_id":"Motif A", "preferred_sequence":"GACAT"}]))

    def test_reject_motif_with_different_width(self):
        with self.assertRaises(ValueError): m.validate_fixture(altered(motifs=[
            {"motif_id":"MOTIF_A", "preferred_sequence":"GACAT"},
            {"motif_id":"MOTIF_B", "preferred_sequence":"GGCATA"}]))

    def test_reject_anchor_not_overlapping_variant(self):
        with self.assertRaises(ValueError): m.validate_fixture(altered(anchor_start_1based=113))

    def test_reject_anchor_outside_sequence(self):
        with self.assertRaises(ValueError): m.validate_fixture(altered(anchor_start_1based=200))

    def test_object_must_be_dict(self):
        with self.assertRaises(ValueError): m.validate_fixture([])


class PpmMathTests(unittest.TestCase):
    def test_ppm_probability_rows_sum_one(self):
        ppm = m.synthetic_ppm("GACAT")
        for row in ppm: self.assertAlmostEqual(sum(row.values()), 1.0)

    def test_ppm_ref_A_matches_consensus(self):
        self.assertEqual(m.mean_probability_proxy("GACAT", m.synthetic_ppm("GACAT")), 0.7)

    def test_ppm_alt_A_one_mismatch(self):
        self.assertEqual(m.mean_probability_proxy("GGCAT", m.synthetic_ppm("GACAT")), 0.58)

    def test_ppm_B_ref_score(self):
        self.assertEqual(m.mean_probability_proxy("GACAT", m.synthetic_ppm("GGCAT")), 0.58)

    def test_ppm_B_alt_score(self):
        self.assertEqual(m.mean_probability_proxy("GGCAT", m.synthetic_ppm("GGCAT")), 0.7)

    def test_empty_consensus_is_invalid(self):
        with self.assertRaises(ValueError): m.synthetic_ppm("")

    def test_unknown_consensus_base_invalid(self):
        with self.assertRaises(ValueError): m.synthetic_ppm("GNCAT")

    def test_width_mismatch_invalid(self):
        with self.assertRaises(ValueError): m.mean_probability_proxy("GAC", m.synthetic_ppm("GACAT"))

    def test_sequence_with_N_invalid(self):
        with self.assertRaises(ValueError): m.mean_probability_proxy("GANAT", m.synthetic_ppm("GACAT"))

    def test_ppm_invalid_row_sum_rejected(self):
        ppm = m.synthetic_ppm("GACAT"); ppm[0]["G"] = 0.6
        with self.assertRaises(ValueError): m.mean_probability_proxy("GACAT", ppm)

    def test_mismatched_windows_have_zero_delta(self):
        rows = m.scan_all_windows(m.validate_fixture(input_data()))
        self.assertTrue(all(row["delta_alt_minus_ref"] == 0
                            for row in rows if not row["variant_in_window"]))

    def test_scan_rows_are_26(self):
        self.assertEqual(len(m.scan_all_windows(m.validate_fixture(input_data()))), 26)

    def test_ten_rows_include_snv(self):
        rows = m.scan_all_windows(m.validate_fixture(input_data()))
        self.assertEqual(sum(row["variant_in_window"] for row in rows), 10)

    def test_anchor_ppm_A_falls_and_B_rises(self):
        rows = produce()["anchor_results"]
        self.assertEqual([x["delta_alt_minus_ref"] for x in rows], [-0.12, 0.12])

    def test_1based_anchor_coordinates(self):
        rows = produce()["anchor_results"]
        self.assertTrue(all(r["window_start_1based"] == 109 and
                            r["window_end_1based_inclusive"] == 113 for r in rows))


class ProvenanceAndEvidenceTests(unittest.TestCase):
    def test_report_marks_synthetic(self):
        r = produce()
        self.assertEqual(r["evidence_class"], "SYNTHETIC_EDUCATIONAL_DEMO")
        self.assertFalse(r["official_motif_atlas"])

    def test_report_never_claims_binding(self):
        r = produce()
        self.assertFalse(r["motif_binding_confirmed"])
        self.assertFalse(r["mini_gene_assay_performed"])

    def test_report_explicitly_not_pwm(self):
        self.assertFalse(produce()["ppm_definition"]["log_odds_pwm_performed"])

    def test_no_reverse_strand_scan(self):
        r = produce()
        self.assertEqual(r["scan_strand"], "+")
        self.assertFalse(r["scanned_reverse_complement"])

    def test_not_real_grch38_ref_validation(self):
        self.assertFalse(produce()["reference_allele_verified_against_grch38"])

    def test_day25_variants_match(self):
        self.assertTrue(m.validate_day25_link(m.load_json(DAY25), m.validate_fixture(input_data())))

    def test_day25_wrong_variant_rejected(self):
        d = m.load_json(DAY25); d["variant_id"] = "chr1:108:A>G"
        with self.assertRaises(ValueError): m.validate_day25_link(d, m.validate_fixture(input_data()))

    def test_day25_live_data_rejected(self):
        d = m.load_json(DAY25); d["mode"] = "LIVE"
        with self.assertRaises(ValueError): m.validate_day25_link(d, m.validate_fixture(input_data()))

    def test_day25_model_output_rejected(self):
        d = m.load_json(DAY25); d["model_prediction"] = True
        with self.assertRaises(ValueError): m.validate_day25_link(d, m.validate_fixture(input_data()))

    def test_day25_score_does_not_combine(self):
        r = produce()
        self.assertTrue(r["day25_linked"])
        self.assertFalse(r["day25_delta_combined_with_ppm"])
        self.assertEqual(r["day25_link_scope"], "VARIANT_ID_ONLY_NOT_NUMERIC_SCORE")

    def test_no_day25_is_standalone(self):
        r = m.build_report(input_data(), fixture_sha256=m.file_sha256(FIXTURE))
        self.assertFalse(r["day25_linked"])
        self.assertIsNone(r["input_provenance"]["day25_summary_sha256"])

    def test_bad_fixture_sha_is_rejected(self):
        with self.assertRaises(ValueError): m.build_report(input_data(), fixture_sha256="abc")

    def test_bad_day25_sha_is_rejected(self):
        with self.assertRaises(ValueError): m.build_report(input_data(), fixture_sha256=m.file_sha256(FIXTURE),
                                                          day25=m.load_json(DAY25), day25_sha256="broken")

    def test_input_sha_matches_bytes(self):
        self.assertEqual(produce()["input_provenance"]["fixture_sha256"], m.file_sha256(FIXTURE))


class FileOutputTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.out = self.root / "day26_output"
        self.summary = m.write_outputs(FIXTURE, self.out, day25_summary_path=DAY25)

    def test_json_csv_svg_html_written(self):
        self.assertTrue(all(Path(self.summary[k]).is_file() for k in
                            ("json_file", "csv_file", "svg_file", "html_file")))

    def test_stdout_contract_counts(self):
        self.assertEqual(self.summary["motif_count"], 2)
        self.assertEqual(self.summary["window_rows_count"], 26)
        self.assertEqual(self.summary["variant_overlapping_window_rows"], 10)

    def test_csv_has_26_rows(self):
        with Path(self.summary["csv_file"]).open(encoding="utf-8", newline="") as h:
            self.assertEqual(len(list(csv.DictReader(h))), 26)

    def test_csv_anchor_A_B_are_correct(self):
        with Path(self.summary["csv_file"]).open(encoding="utf-8", newline="") as h:
            rows = [r for r in csv.DictReader(h) if r["window_start_1based"] == "109"]
        self.assertEqual([(r["motif_id"], r["ref_ppm_mean"], r["alt_ppm_mean"]) for r in rows],
                         [("MOTIF_A", "0.7", "0.58"), ("MOTIF_B", "0.58", "0.7")])

    def test_json_carries_full_provenance(self):
        r = json.loads(Path(self.summary["json_file"]).read_text(encoding="utf-8"))
        self.assertEqual(r["input_provenance"]["day25_summary_sha256"], m.file_sha256(DAY25))

    def test_svg_valid_xml_1920_by_1080(self):
        root = ET.parse(self.summary["svg_file"]).getroot()
        self.assertEqual(root.get("width"), "1920")
        self.assertEqual(root.get("height"), "1080")

    def test_svg_shows_both_score_directions(self):
        t = Path(self.summary["svg_file"]).read_text(encoding="utf-8")
        self.assertIn("-0.12", t)
        self.assertIn("+0.12", t)

    def test_svg_labels_synthetic(self):
        t = Path(self.summary["svg_file"]).read_text(encoding="utf-8")
        self.assertIn("SYNTHETIC", t)
        self.assertIn("TOY_NOT_GRCH38", t)

    def test_html_renders_local_links(self):
        h = Path(self.summary["html_file"]).read_text(encoding="utf-8")
        self.assertIn("day26_motif_report.json", h)
        self.assertIn("day26_ppm_preference.svg", h)

    def test_html_contains_no_javascript(self):
        h = Path(self.summary["html_file"]).read_text(encoding="utf-8")
        self.assertNotIn("<script", h.lower())
        self.assertNotIn("https://", h.lower())

    def test_html_table_is_escape_safe(self):
        r = produce(); r["window_rows"][0]["motif_id"] = "<script>alert(1)</script>"
        html = m.draw_html(r, m.draw_svg(produce()))
        self.assertIn("&lt;script&gt;", html)
        self.assertNotIn("<script>alert(1)</script>", html)

    def test_read_json_duplicate_keys_rejected(self):
        p = self.root / "duplicate.json"; p.write_text('{"x":1,"x":2}', encoding="utf-8")
        with self.assertRaises(ValueError): m.load_json(p)

    def test_read_json_nan_rejected(self):
        p = self.root / "bad.json"; p.write_text('{"x":NaN}', encoding="utf-8")
        with self.assertRaises(ValueError): m.load_json(p)

    def test_read_json_oversized_rejected(self):
        p = self.root / "large.json"; p.write_text(" " * 500, encoding="utf-8")
        with self.assertRaises(ValueError): m.load_json(p, max_bytes=100)

    def test_read_json_nonexistent_rejected(self):
        with self.assertRaises(ValueError): m.load_json(self.root / "missing.json")


class CommandLineTests(unittest.TestCase):
    def test_explicit_live_mode_rejected_by_parser(self):
        p = subprocess.run([sys.executable, str(ROOT / "day26_offline_ppm_motif_lab.py"),
                            "--mode", "live"], cwd=ROOT, capture_output=True, text=True)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("invalid choice", p.stderr)

    def test_cli_actual_stdout_json(self):
        with tempfile.TemporaryDirectory() as d:
            p = subprocess.run([sys.executable, str(ROOT / "day26_offline_ppm_motif_lab.py"),
                                "--mode", "offline", "--fixture", str(FIXTURE),
                                "--day25-summary", str(DAY25),
                                "--output-dir", str(Path(d) / "out")], cwd=ROOT,
                               capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertEqual(json.loads(p.stdout)["status"], "PASS")

    def test_cli_without_day25_optional(self):
        with tempfile.TemporaryDirectory() as d:
            p = subprocess.run([sys.executable, str(ROOT / "day26_offline_ppm_motif_lab.py"),
                                "--mode", "offline", "--fixture", str(FIXTURE),
                                "--output-dir", str(Path(d) / "out")], cwd=ROOT,
                               capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertFalse(json.loads(p.stdout)["day25_linked"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
