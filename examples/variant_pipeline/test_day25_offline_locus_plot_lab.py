"""Day 25 | offline-only integrity, coordinate, visualization and CLI tests."""
from __future__ import annotations

import copy
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

import day25_offline_locus_plot_lab as lab

HERE = Path(__file__).resolve().parent
REPORT = HERE / "day24_track_delta_report.json"
TRACK = HERE / "day24_synthetic_tracks.json"
SVG = {"s": lab.SVG_NS}


def inputs():
    return copy.deepcopy(lab.load_json(REPORT)), copy.deepcopy(lab.load_json(TRACK))


class Day25Inputs(unittest.TestCase):
    def setUp(self):
        self.report, self.track = inputs()

    def ok(self):
        return lab.verify_inputs(self.report, self.track, track_sha=lab.sha256_file(TRACK))

    def test_day24_fixture_contract_passes(self):
        self.assertIs(self.ok(), self.report)

    def test_provenance_sha_matches_actual_track(self):
        self.assertEqual(lab.sha256_file(TRACK), self.report["input_provenance"]["track_fixture"]["sha256"])

    def test_expected_rows_and_signed_delta(self):
        self.ok()
        self.assertEqual([r["delta_alt_minus_ref"] for r in self.report["point_rows"]], [0, .25, -.3, .1, 0])
        self.assertEqual([r["position_1based"] for r in self.report["point_rows"]], [108, 109, 110, 111, 112])

    def test_at_variant_position_110(self):
        self.assertEqual([r["position_1based"] for r in self.report["point_rows"] if r["variant_here"]], [110])

    def test_reject_wrong_assembly(self):
        self.report["assembly"] = "GRCh38"
        with self.assertRaises(ValueError): self.ok()

    def test_reject_true_model_prediction(self):
        self.report["model_prediction"] = True
        with self.assertRaises(ValueError): self.ok()

    def test_reject_true_atlas_query(self):
        self.report["atlas_query_performed"] = True
        with self.assertRaises(ValueError): self.ok()

    def test_reject_nonzero_api_calls(self):
        self.report["real_api_calls_this_run"] = 1
        with self.assertRaises(ValueError): self.ok()

    def test_reject_clinical_interpretation_true(self):
        self.report["clinical_interpretation"] = True
        with self.assertRaises(ValueError): self.ok()

    def test_reject_wrong_schema_version(self):
        self.report["schema_version"] = "unknown"
        with self.assertRaises(ValueError): self.ok()

    def test_reject_altered_variant_id(self):
        self.track["variant_id"] = "chr1:110:A>T"
        with self.assertRaises(ValueError): self.ok()

    def test_reject_equal_alleles(self):
        self.track["variant_id"] = self.report["variant_id"] = "chr1:110:A>A"
        with self.assertRaises(ValueError): self.ok()

    def test_reject_wrong_point_count(self):
        self.report["point_count"] = 3
        with self.assertRaises(ValueError): self.ok()

    def test_reject_missing_point_rows(self):
        self.report.pop("point_rows")
        with self.assertRaises(ValueError): self.ok()

    def test_reject_empty_point_rows(self):
        self.report["point_rows"] = []
        with self.assertRaises(ValueError): self.ok()

    def test_reject_position_shift(self):
        self.track["positions_1based"][2] = 116
        with self.assertRaises(ValueError): self.ok()

    def test_reject_variant_outside_plot(self):
        self.track["variant_id"] = self.report["variant_id"] = "chr1:200:A>G"
        with self.assertRaises(ValueError): self.ok()

    def test_reject_bool_in_coordinates(self):
        self.track["positions_1based"][0] = True
        with self.assertRaises(ValueError): self.ok()

    def test_reject_bool_ref_value(self):
        self.track["ref"][0] = True
        with self.assertRaises(ValueError): self.ok()

    def test_reject_nan_in_track(self):
        self.track["alt"][0] = float("nan")
        with self.assertRaises(ValueError): self.ok()

    def test_reject_out_of_range_proxy(self):
        self.track["alt"][0] = 1.001
        with self.assertRaises(ValueError): self.ok()

    def test_reject_mismatched_delta(self):
        self.report["point_rows"][2]["delta_alt_minus_ref"] = +0.3
        with self.assertRaises(ValueError): self.ok()

    def test_reject_variant_location_flag_change(self):
        self.report["point_rows"][2]["variant_here"] = False
        with self.assertRaises(ValueError): self.ok()

    def test_reject_peak_sign_loss(self):
        self.report["largest_absolute_delta"] = +0.3
        with self.assertRaises(ValueError): self.ok()

    def test_reject_peak_position_mismatch(self):
        self.report["largest_absolute_delta_position_1based"] = 111
        with self.assertRaises(ValueError): self.ok()

    def test_reject_wrong_track_sha(self):
        with self.assertRaises(ValueError):
            lab.verify_inputs(self.report, self.track, track_sha="0" * 64)

    def test_reject_missing_provenance(self):
        self.report.pop("input_provenance")
        with self.assertRaises(ValueError): self.ok()

    def test_reject_bad_day23_link(self):
        self.report["day23_annotation_link"]["score_status"] = "OFFICIAL_AVI"
        with self.assertRaises(ValueError): self.ok()

    def test_no_day23_link_can_be_plotted(self):
        self.report["day23_annotation_link"] = None
        self.ok()

    def test_nonfinite_json_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            file = Path(d)/"bad.json"
            file.write_text('{"score":NaN}', encoding="utf-8")
            with self.assertRaises(ValueError): lab.load_json(file)

    def test_invalid_json_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            file=Path(d)/"bad.json"
            file.write_text('not json', encoding="utf-8")
            with self.assertRaises(ValueError): lab.load_json(file)


class Day25Outputs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.outdir = Path(cls.temp.name)/"out"
        cls.result = lab.emit_outputs(REPORT, TRACK, cls.outdir)
        cls.summary = json.loads((cls.outdir/"day25_locus_summary.json").read_text(encoding="utf-8"))
        cls.svg = (cls.outdir/"day25_locus_tracks.svg").read_text(encoding="utf-8")
        cls.root = ET.fromstring(cls.svg)
        cls.html = (cls.outdir/"day25_locus_view.html").read_text(encoding="utf-8")

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_output_file_count(self):
        self.assertEqual(len(list(self.outdir.iterdir())), 4)

    def test_stdout_summary_contract(self):
        self.assertEqual(self.result["point_count"], 5)
        self.assertEqual(self.result["track_count"], 2)
        self.assertEqual(self.result["largest_absolute_delta"], -.3)
        self.assertEqual(self.result["real_api_calls_this_run"], 0)

    def test_json_version_and_source_hashes(self):
        self.assertEqual(self.summary["schema_version"], lab.FIGURE_SCHEMA)
        self.assertEqual(self.summary["source_day24_report_sha256"], lab.sha256_file(REPORT))
        self.assertEqual(self.summary["source_day24_track_sha256"], lab.sha256_file(TRACK))

    def test_json_lacks_official_model_prediction(self):
        self.assertFalse(self.summary["model_prediction"])
        self.assertFalse(self.summary["official_atlas_data"])
        self.assertFalse(self.summary["clinical_interpretation"])

    def test_json_keeps_day23_separate(self):
        self.assertTrue(self.summary["day23_linked"])
        self.assertFalse(self.summary["day23_score_combined_with_delta"])
        self.assertNotIn("synthetic_score", self.summary)

    def test_csv_header_and_rows(self):
        with (self.outdir/"day25_locus_track.csv").open(newline="", encoding="utf-8") as h:
            rows = list(csv.DictReader(h))
        self.assertEqual(len(rows), 5)
        self.assertEqual(list(rows[0]), ["position_1based", "ref_proxy", "alt_proxy", "delta_alt_minus_ref", "variant_here"])
        self.assertEqual(rows[2]["delta_alt_minus_ref"], "-0.3")

    def test_svg_canvas_is_16_by_9(self):
        self.assertEqual(self.root.attrib["viewBox"], "0 0 1920 1080")
        self.assertEqual(self.root.attrib["width"], "1920")
        self.assertEqual(self.root.attrib["height"], "1080")

    def test_svg_background_is_required_color(self):
        self.assertTrue(any(e.attrib.get("fill") == "#E5DEDA" for e in self.root.findall("s:rect", SVG)))

    def test_svg_has_ref_and_alt_tracks(self):
        tracks = self.root.findall(".//s:polyline", SVG)
        self.assertEqual({t.attrib.get("data-track") for t in tracks}, {"REF", "ALT"})
        self.assertEqual(len(tracks), 2)

    def test_svg_has_data_points_with_titles(self):
        dots = [e for e in self.root.findall(".//s:circle", SVG) if e.attrib.get("data-track") in ("REF", "ALT")]
        self.assertEqual(len(dots), 10)
        self.assertEqual(len([d for d in dots if d.find("s:title",SVG) is not None]), 10)

    def test_svg_marks_variant_at_110(self):
        lines = [e for e in self.root.findall(".//s:line",SVG) if e.attrib.get("data-variant-marker")]
        self.assertEqual([e.attrib["data-variant-marker"] for e in lines], ["110"])

    def test_svg_shows_negative_delta(self):
        all_text=" ".join(e.text or "" for e in self.root.findall(".//s:text", SVG))
        self.assertIn("-0.30", all_text)
        self.assertIn("ALT − REF", all_text)

    def test_html_contains_embedded_svg(self):
        self.assertIn("<svg", self.html)
        self.assertIn("<table>", self.html)

    def test_html_contains_warning(self):
        self.assertIn("SYNTHETIC_EDUCATIONAL_DEMO", self.html)
        self.assertIn("TOY_NOT_GRCH38", self.html)
        self.assertIn("Day 23", self.html)

    def test_html_links_outputs(self):
        for f in ("day25_locus_track.csv", "day25_locus_summary.json", "day25_locus_tracks.svg"):
            self.assertIn(f, self.html)

    def test_no_external_js_or_remote_assets(self):
        self.assertNotIn("<script", self.html)
        self.assertNotIn("src=\"http", self.html)

    def test_svg_and_html_do_not_contain_day23_score(self):
        self.assertNotIn("-0.17", self.svg)
        self.assertNotIn("-0.17", self.html)

    def test_modified_fixture_hash_rejected(self):
        track=lab.load_json(TRACK)
        track["alt"][0]=.3
        with tempfile.TemporaryDirectory() as d:
            copyfile=Path(d)/"changed.json"
            copyfile.write_text(json.dumps(track), encoding="utf-8")
            out=Path(d)/"out"
            with self.assertRaises(ValueError): lab.emit_outputs(REPORT, copyfile, out)
            self.assertFalse(out.exists())


class Day25CommandLine(unittest.TestCase):
    def test_cli_standalone(self):
        with tempfile.TemporaryDirectory() as d:
            cmd=[sys.executable, str(HERE/"day25_offline_locus_plot_lab.py"), "--mode","offline",
                 "--day24-report",str(REPORT), "--source-track",str(TRACK), "--output-dir",str(Path(d)/"out")]
            p=subprocess.run(cmd,capture_output=True,text=True,check=True)
            data=json.loads(p.stdout)
            self.assertEqual(data["status"],"PASS")
            self.assertEqual(data["real_api_calls_this_run"],0)

    def test_live_is_explicitly_rejected(self):
        p=subprocess.run([sys.executable,str(HERE/"day25_offline_locus_plot_lab.py"),"--mode","live"],capture_output=True,text=True)
        self.assertNotEqual(p.returncode,0)
        self.assertIn("invalid choice",p.stderr)

    def test_missing_report_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=subprocess.run([sys.executable,str(HERE/"day25_offline_locus_plot_lab.py"),"--mode","offline",
             "--day24-report",str(Path(d)/"missing.json"),"--source-track",str(TRACK), "--output-dir",str(Path(d)/"out")],capture_output=True,text=True)
            self.assertNotEqual(p.returncode,0)
            self.assertFalse((Path(d)/"out").exists())


if __name__ == "__main__":
    unittest.main()
