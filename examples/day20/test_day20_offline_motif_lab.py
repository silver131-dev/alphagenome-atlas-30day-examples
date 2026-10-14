"""Day 20 offline coordinate, fixture, CLI and HTML contract tests."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest

from day20_offline_motif_lab import (
    DEFAULT_MOTIFS, MotifInterval, analyze, load_synthetic_tsv,
    main, overlap_bp, parse_region_1based, parse_variant_1based,
    render_html, run, valid_motifs,
)


class Day20CoordinateTests(unittest.TestCase):
    def test_default_overlap_ids(self):
        report = analyze()
        self.assertEqual(report["overlap_ids"], ["A", "B"])
        self.assertEqual(report["overlap_count"], 2)
        self.assertEqual([x["overlap_bp"] for x in report["region_overlaps"]], [5, 5])

    def test_1_based_inclusive_conversion(self):
        self.assertEqual(parse_region_1based("chr1:100-120"), ("chr1", 99, 120))

    def test_right_endpoint_not_included(self):
        self.assertEqual(overlap_bp(99, 120, 120, 126), 0)
        self.assertNotIn("C", analyze()["overlap_ids"])

    def test_next_single_base_hits_c(self):
        result = analyze(region="chr1:121-121", variant="chr1:121:A>G")
        self.assertEqual(result["overlap_ids"], ["C"])
        self.assertEqual(result["variant_overlap_ids"], ["C"])
        self.assertEqual(result["region_overlaps"][0]["overlap_bp"], 1)

    def test_chromosome_must_match(self):
        self.assertEqual(analyze(region="chr2:100-101", variant="chr2:100:A>G")["overlap_ids"], ["D"])
        self.assertNotIn("D", analyze()["overlap_ids"])

    def test_first_base_coordinate(self):
        self.assertEqual(parse_region_1based("chr1:1-1"), ("chr1", 0, 1))

    def test_invalid_region_coordinates(self):
        for value in ("chr1:0-1", "chr1:4-3", "chr1:-2-5", "chr1:5", "chr23:1-4", "chr1:1-10002"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                parse_region_1based(value)

    def test_variant_in_b_only(self):
        result = analyze()
        self.assertEqual(result["variant_overlap_ids"], ["B"])
        self.assertTrue(result["variant_in_query"])
        self.assertEqual(result["region_overlaps"][0]["contains_variant"], False)
        self.assertEqual(result["region_overlaps"][1]["contains_variant"], True)

    def test_variant_outside_region_is_not_annotated(self):
        result = analyze(region="chr1:121-121")
        self.assertEqual(result["overlap_ids"], ["C"])
        self.assertFalse(result["variant_in_query"])
        self.assertEqual(result["variant_overlap_ids"], [])

    def test_optional_variant_is_null(self):
        result = analyze(variant=None)
        self.assertIsNone(result["variant_1based"])
        self.assertEqual(result["variant_overlap_ids"], [])

    def test_bad_snv_rejected(self):
        for value in ("chr1:0:A>G", "chr1:11:A>A", "chr1:11:N>G", "chr1:11:A>GG"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                parse_variant_1based(value)

    def test_query_never_claims_motif_or_api(self):
        result = analyze()
        for key in ("motif_sequence_scanned", "reference_allele_verified", "official_motif_atlas",
                    "atlas_query_performed", "clinical_interpretation"):
            self.assertFalse(result[key])
        self.assertEqual(result["real_api_calls_this_run"], 0)
        self.assertEqual(result["evidence_class"], "SYNTHETIC_EDUCATIONAL_DEMO")


class Day20FixtureTests(unittest.TestCase):
    def test_documented_tsv_matches_builtin(self):
        path = Path(__file__).with_name("day20_synthetic_motifs.tsv")
        self.assertEqual(load_synthetic_tsv(path), DEFAULT_MOTIFS)

    def test_duplicate_motif_rejected(self):
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            valid_motifs([DEFAULT_MOTIFS[0], DEFAULT_MOTIFS[0]])

    def test_non_synthetic_motif_rejected(self):
        with self.assertRaisesRegex(ValueError, "synthetic"):
            valid_motifs([MotifInterval("A", "chr1", 0, 5, "REAL_MOTIF")])

    def test_zero_length_motif_rejected(self):
        with self.assertRaises(ValueError):
            valid_motifs([MotifInterval("E", "chr1", 2, 2)])

    def test_fixture_bad_header_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "bad.tsv"
            path.write_text("motif_id\tchrom\tstart\tend\nA\tchr1\t1\t2\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                load_synthetic_tsv(path)

    def test_fixture_bad_coordinate_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "bad.tsv"
            path.write_text("motif_id\tchrom\tstart_0based\tend_0based\tevidence_class\n"
                            "A\tchr1\tfoo\t104\tSYNTHETIC_FIXTURE\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                load_synthetic_tsv(path)


class Day20ArtifactTests(unittest.TestCase):
    def test_html_renders_synthetic_disclaimer(self):
        html = render_html(analyze())
        self.assertIn("SYNTHETIC EDUCATIONAL DEMO", html)
        self.assertIn("<td>A</td>", html)
        self.assertIn("<td>B</td>", html)
        self.assertIn("0-based", html)
        self.assertNotIn("javascript:", html)

    def test_cli_creates_json_and_html(self):
        with tempfile.TemporaryDirectory() as folder:
            out = run(folder)
            report = json.loads(Path(out["json_file"]).read_text(encoding="utf-8"))
            self.assertEqual(out["overlap_ids"], ["A", "B"])
            self.assertEqual(report["variant_overlap_ids"], ["B"])
            self.assertTrue(Path(out["html_file"]).is_file())

    def test_cli_stdout_contract(self):
        with tempfile.TemporaryDirectory() as folder:
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                main(["--mode", "offline", "--output-dir", folder])
            data = json.loads(out.getvalue())
            self.assertEqual(data["mode"], "OFFLINE_ONLY")
            self.assertEqual(data["real_api_calls_this_run"], 0)

    def test_live_mode_explicitly_rejected(self):
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as exc:
                main(["--mode", "live"])
        self.assertEqual(exc.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
