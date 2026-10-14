"""Unittest suite for Day 08 synthetic offline multimodal track Lab."""
from __future__ import annotations

import ast
import copy
import csv
import hashlib
import io
import json
import math
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from xml.etree import ElementTree as ET

import day08_offline_multitrack_lab as lab

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "day08_offline_multitrack_lab.py"
FIXTURE = ROOT / "day08_synthetic_multitrack.json"


class FixtureContracts(unittest.TestCase):
    def setUp(self):
        self.valid, self.digest = lab.load_fixture(FIXTURE)
        self.x = copy.deepcopy(self.valid)

    def bad(self):
        with tempfile.TemporaryDirectory() as t:
            p = Path(t) / "bad.json"
            p.write_text(json.dumps(self.x, allow_nan=True), encoding="utf-8")
            with self.assertRaises((ValueError, TypeError, KeyError)):
                lab.load_fixture(p)

    def test_good_fixture(self):
        self.assertEqual(self.valid["assembly"], "TOY_NOT_GRCH38")

    def test_fixture_digest(self):
        self.assertEqual(self.digest, hashlib.sha256(FIXTURE.read_bytes()).hexdigest())

    def test_ref_allele_matches_sequence(self):
        v = self.x["variant"]
        self.assertEqual(self.x["ref_sequence"][v["position_1based"]-self.x["start_1based"]], "A")

    def test_synthetic_provenance_enforced(self):
        self.x["evidence_class"] = "HISTORICAL_LIVE_OUTPUT"; self.bad()

    def test_real_assembly_rejected(self):
        self.x["assembly"] = "GRCh38"; self.bad()

    def test_model_prediction_true_rejected(self):
        self.x["model_prediction"] = True; self.bad()

    def test_missing_model_prediction_rejected(self):
        self.x.pop("model_prediction"); self.bad()

    def test_real_chrom_rejected(self):
        self.x["chrom"] = "chr1"; self.bad()

    def test_nonint_start_rejected(self):
        self.x["start_1based"] = "1001"; self.bad()

    def test_bool_start_rejected(self):
        self.x["start_1based"] = True; self.bad()

    def test_zero_start_rejected(self):
        self.x["start_1based"] = 0; self.bad()

    def test_negative_span_rejected(self):
        self.x["end_1based"] = 1000; self.bad()

    def test_bad_span_length_rejected(self):
        self.x["end_1based"] = 1023; self.bad()

    def test_non_dna_ref_sequence_rejected(self):
        self.x["ref_sequence"] = self.x["ref_sequence"].replace("A", "N"); self.bad()

    def test_short_ref_sequence_rejected(self):
        self.x["ref_sequence"] = self.x["ref_sequence"][:-1]; self.bad()

    def test_bad_variant_fields_rejected(self):
        self.x["variant"]["source"]="dummy"; self.bad()

    def test_variant_wrong_chrom_rejected(self):
        self.x["variant"]["chrom"]="chr2"; self.bad()

    def test_variant_outside_range_rejected(self):
        self.x["variant"]["position_1based"] = 1030; self.bad()

    def test_variant_bool_position_rejected(self):
        self.x["variant"]["position_1based"] = True; self.bad()

    def test_variant_ref_mismatch_rejected(self):
        self.x["variant"]["ref"]="C"; self.bad()

    def test_variant_same_alt_rejected(self):
        self.x["variant"]["alt"]="A"; self.bad()

    def test_variant_indel_rejected(self):
        self.x["variant"]["alt"]="GT"; self.bad()

    def test_noninteger_track_resolution_rejected(self):
        self.x["tracks"]["atac_accessibility"]["resolution_bp"] = 6; self.bad()

    def test_wrong_units_rejected(self):
        self.x["tracks"]["atac_accessibility"]["units"] = "OFFICIAL_READ_COUNT"; self.bad()

    def test_extra_track_rejected(self):
        self.x["tracks"]["official_model"] = {}; self.bad()

    def test_short_track_rejected(self):
        self.x["tracks"]["splice_donor_plus"]["ref"] = self.x["tracks"]["splice_donor_plus"]["ref"][:-1]; self.bad()

    def test_negative_proxy_rejected(self):
        self.x["tracks"]["atac_accessibility"]["ref"][3] = -0.5; self.bad()

    def test_large_proxy_rejected(self):
        self.x["tracks"]["splice_donor_plus"]["alt"][3] = 1.1; self.bad()

    def test_nan_rejected(self):
        self.x["tracks"]["splice_donor_plus"]["ref"][3] = float("nan"); self.bad()

    def test_infinity_rejected(self):
        self.x["tracks"]["splice_donor_plus"]["ref"][3] = float("inf"); self.bad()

    def test_boolean_proxy_rejected(self):
        self.x["tracks"]["atac_accessibility"]["ref"][3] = True; self.bad()

    def test_one_sided_null_rejected(self):
        self.x["tracks"]["splice_donor_plus"]["ref"][0] = None; self.bad()

    def test_atac_null_rejected(self):
        self.x["tracks"]["atac_accessibility"]["ref"][0] = None
        self.x["tracks"]["atac_accessibility"]["alt"][0] = None
        self.bad()

    def test_unknown_mock_biosample_rejected(self):
        self.x["synthetic_biosample"] = "REAL_PATIENT_BLOOD"; self.bad()

    def test_json_array_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)/"bad.json"; p.write_text("[]", encoding="utf-8")
            with self.assertRaises(ValueError):lab.load_fixture(p)


class TrackComputations(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture, cls.sha = lab.load_fixture(FIXTURE)
        cls.rows = lab.rows_from_fixture(cls.fixture)
        cls.report = lab.build_report(cls.fixture, cls.sha)

    def test_base_count_24(self):self.assertEqual(len(self.rows),24)
    def test_first_base(self):self.assertEqual(self.rows[0]["position_1based"],1001)
    def test_last_base(self):self.assertEqual(self.rows[-1]["position_1based"],1024)
    def test_variant_coordinate(self):self.assertEqual(self.report["variant"],"chrToy:1012:A>G")
    def test_splice_delta_positive(self):self.assertEqual(self.report["splice_max_abs_delta"]["delta_alt_minus_ref"],0.7)
    def test_atac_delta_negative(self):self.assertEqual(self.report["atac_max_abs_delta"]["delta_alt_minus_ref"],-0.36)
    def test_splice_peak_at_variant(self):self.assertEqual(self.report["splice_max_abs_delta"]["position_1based"],1012)
    def test_atac_peak_at_variant(self):self.assertEqual(self.report["atac_max_abs_delta"]["position_1based"],1012)
    def test_splice_at_variant_values(self):
        r=self.rows[11];self.assertEqual((r["splice_ref"],r["splice_alt"]),(0.1,0.8))
    def test_atac_at_variant_values(self):
        r=self.rows[11];self.assertEqual((r["atac_ref"],r["atac_alt"]),(0.62,0.26))
    def test_zero_not_null(self):
        r=self.rows[4];self.assertEqual(r["atac_ref"],0.0);self.assertEqual(r["atac_delta"],0.0)
    def test_missing_is_null(self):
        self.assertEqual(self.rows[18]["splice_ref"],None)
        self.assertEqual(self.rows[18]["splice_delta"],None)
    def test_missing_position_is_1019(self):self.assertEqual(self.report["missing_splice_positions_1based"],[1019])
    def test_two_1bp_source_resolutions(self):self.assertEqual(self.report["source_resolution_bp"],{"splice_sites":1,"atac":1})
    def test_display_width_is_6(self):self.assertEqual(self.report["atac_display_bin_width_bp"],6)
    def test_four_6bp_bins(self):self.assertEqual(self.report["atac_display_bin_count"],4)
    def test_variant_in_second_bin(self):self.assertEqual((self.report["user_display_bins"][1]["start_1based"],self.report["user_display_bins"][1]["end_1based"]),(1007,1012))
    def test_second_bin_0based_conversion(self):self.assertEqual((self.report["user_display_bins"][1]["start_0based"],self.report["user_display_bins"][1]["end_0based_exclusive"]),(1006,1012))
    def test_bin_delta_dilution(self):self.assertEqual(self.report["variant_display_bin_atac_delta_mean"],-0.06)
    def test_binned_not_official(self):self.assertFalse(self.report["official_atac_binned_to_6bp"])
    def test_model_prediction_false(self):self.assertFalse(self.report["model_prediction"])
    def test_no_clinical(self):self.assertFalse(self.report["clinical_interpretation"])
    def test_no_live_api(self):self.assertEqual(self.report["real_api_calls_this_run"],0)
    def test_no_atlas_query(self):self.assertFalse(self.report["atlas_query_performed"])
    def test_report_source_hash(self):self.assertEqual(self.report["source_sha256"],self.sha)
    def test_bad_bin_divisor_rejected(self):
        with self.assertRaises(ValueError):lab.make_display_bins(self.rows,5)
    def test_bad_zero_bin_rejected(self):
        with self.assertRaises(ValueError):lab.make_display_bins(self.rows,0)
    def test_boolean_bin_size_rejected(self):
        with self.assertRaises(ValueError):lab.make_display_bins(self.rows,True)
    def test_peak_all_missing_rejected(self):
        with self.assertRaises(ValueError):lab.peak([{"splice_delta":None}],"splice_delta")
    def test_bin_range_no_gap(self):
        bins=self.report["user_display_bins"]
        for a,b in zip(bins,bins[1:]):self.assertEqual(a["end_0based_exclusive"],b["start_0based"])
    def test_bin_label_is_display_only(self):
        self.assertEqual(self.report["user_display_bins"][0]["aggregation"],"USER_CREATED_MEAN_FOR_DISPLAY_NOT_MODEL_RESOLUTION")


class FormatOutputs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        fixture,sha=lab.load_fixture(FIXTURE)
        cls.report=lab.build_report(fixture,sha)
        cls.contents=lab.materialize(fixture,sha)

    def test_five_output_files(self):self.assertEqual(set(self.contents),set(lab.FILES))
    def test_json_parses(self):self.assertEqual(json.loads(self.contents[lab.FILES[0]])["status"],"PASS")
    def test_base_csv_has_24_rows(self):
        self.assertEqual(len(list(csv.DictReader(io.StringIO(self.contents[lab.FILES[1]])))),24)
    def test_base_csv_missing_is_blank(self):
        rr=list(csv.DictReader(io.StringIO(self.contents[lab.FILES[1]])))
        self.assertEqual(rr[18]["splice_ref"],"")
    def test_base_csv_zero_is_zero(self):
        rr=list(csv.DictReader(io.StringIO(self.contents[lab.FILES[1]])))
        self.assertEqual(rr[4]["atac_ref"],"0.0")
    def test_binned_csv_has_4_rows(self):
        self.assertEqual(len(list(csv.DictReader(io.StringIO(self.contents[lab.FILES[2]])))),4)
    def test_svg_parses(self):
        s=ET.fromstring(self.contents[lab.FILES[3]])
        self.assertTrue(s.tag.endswith("svg"))
    def test_svg_16by9(self):
        s=ET.fromstring(self.contents[lab.FILES[3]])
        self.assertEqual((s.attrib["width"],s.attrib["height"]),("1920","1080"))
    def test_svg_background(self):self.assertIn("#E5DEDA",self.contents[lab.FILES[3]])
    def test_svg_labels_modalities(self):
        s=self.contents[lab.FILES[3]]
        self.assertIn("SPLICE DONOR (1 bp)",s)
        self.assertIn("ATAC ACCESSIBILITY (1 bp)",s)
    def test_svg_legend_y_at_header(self):
        s=ET.fromstring(self.contents[lab.FILES[3]])
        svgtexts=list(s.iter("{"+lab.SVG_NS+"}text"))
        refs=[x for x in svgtexts if x.text=="● REF"]
        self.assertEqual([x.attrib["y"] for x in refs],["294","632"])
    def test_svg_missing_track_disconnected(self):
        s=self.contents[lab.FILES[3]]
        self.assertIn("gap = missing value",s)
        self.assertIn("NOT AN ALPHAGENOME PREDICTION",s)
    def test_html_is_offline(self):
        self.assertNotIn("<script",self.contents[lab.FILES[4]])
        self.assertNotIn("https://",self.contents[lab.FILES[4]])
    def test_html_synthetic_disclaimer(self):self.assertIn("NOT AlphaGenome output",self.contents[lab.FILES[4]])
    def test_html_table(self):self.assertIn("<table>",self.contents[lab.FILES[4]])
    def test_html_null_field(self):self.assertIn("<td>null</td>",self.contents[lab.FILES[4]])
    def test_html_variant(self):self.assertIn("chrToy:1012:A&gt;G",self.contents[lab.FILES[4]])
    def test_stdlib_only_source(self):
        tree=ast.parse(SOURCE.read_text(encoding="utf-8"))
        names=set()
        for x in ast.walk(tree):
            if isinstance(x,ast.Import):names.update(n.name.split(".")[0] for n in x.names)
            if isinstance(x,ast.ImportFrom):names.add(x.module.split(".")[0])
        self.assertEqual(names,{"__future__","argparse","csv","hashlib","html","io","json","math","re","sys","pathlib","xml"})
    def test_expected_file_names(self):self.assertIn("day08_multitrack.svg",self.contents)


class CLIAcceptance(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.tmp=Path(self.temp.name)
        self.out=self.tmp/"results"

    def command(self,*args,env=None,cwd=None):
        return subprocess.run([sys.executable,str(SOURCE),*args],cwd=cwd,env=env,text=True,capture_output=True)

    def run_tool(self):return lab.run(FIXTURE,self.out)

    def test_run(self):self.assertEqual(self.run_tool()["status"],"PASS")
    def test_run_summary_count(self):self.assertEqual(self.run_tool()["artifact_count"],5)
    def test_run_no_api(self):self.assertEqual(self.run_tool()["real_api_calls_this_run"],0)
    def test_run_and_verify(self):
        self.run_tool();self.assertEqual(lab.verify(FIXTURE,self.out)["checksums"],"MATCH")
    def test_five_output_files_and_manifest(self):
        self.run_tool();self.assertEqual(len(list(self.out.iterdir())),6)
    def test_manifest_has_5_hashes(self):
        self.run_tool();d=json.loads((self.out/"day08_artifact_manifest.json").read_text())
        self.assertEqual(len(d["files"]),5)
    def test_manifest_not_signed(self):
        self.run_tool();d=json.loads((self.out/"day08_artifact_manifest.json").read_text())
        self.assertFalse(d["signed_provenance"])
    def test_overwrite_refused(self):
        self.run_tool()
        with self.assertRaises(FileExistsError):lab.run(FIXTURE,self.out)
    def test_existing_empty_dir_ok(self):
        self.out.mkdir();self.assertEqual(self.run_tool()["status"],"PASS")
    def test_tampered_svg_rejected(self):
        self.run_tool();p=self.out/"day08_multitrack.svg";p.write_text(p.read_text()+"TAMPER")
        with self.assertRaises(ValueError):lab.verify(FIXTURE,self.out)
    def test_tampered_report_rejected(self):
        self.run_tool();p=self.out/"day08_multitrack_report.json";r=json.loads(p.read_text());r["real_api_calls_this_run"]=1;p.write_text(lab.render_json(r))
        with self.assertRaises(ValueError):lab.verify(FIXTURE,self.out)
    def test_tampered_manifest_rejected(self):
        self.run_tool();p=self.out/"day08_artifact_manifest.json";r=json.loads(p.read_text());r["files"].pop("day08_multitrack.svg");p.write_text(lab.render_json(r))
        with self.assertRaises(ValueError):lab.verify(FIXTURE,self.out)
    def test_tampered_source_rejected(self):
        self.run_tool();f=self.tmp/"source_copy.json";f.write_bytes(FIXTURE.read_bytes()+b" ")
        with self.assertRaises(ValueError):lab.verify(f,self.out)
    def test_malicious_rehashed_report_rejected(self):
        self.run_tool();p=self.out/"day08_multitrack_report.json";r=json.loads(p.read_text());r["clinical_interpretation"]=True; p.write_text(lab.render_json(r))
        manifest_path=self.out/"day08_artifact_manifest.json";manifest=json.loads(manifest_path.read_text())
        manifest["files"][p.name]=lab.sha256_bytes(p.read_bytes())
        manifest_path.write_text(lab.render_json(manifest))
        with self.assertRaises(ValueError):lab.verify(FIXTURE,self.out)
    def test_cli_run_from_different_cwd(self):
        r=self.command("run","--mode","offline","--fixture",str(FIXTURE),"--output-dir",str(self.out),cwd=self.tmp)
        self.assertEqual(r.returncode,0,r.stderr)
        self.assertEqual(json.loads(r.stdout)["base_count"],24)
    def test_cli_verify_from_different_cwd(self):
        self.run_tool();r=self.command("verify","--mode","offline","--fixture",str(FIXTURE),"--output-dir",str(self.out),cwd=self.tmp)
        self.assertEqual(r.returncode,0,r.stderr)
        self.assertEqual(json.loads(r.stdout)["files_verified"],5)
    def test_live_rejected_with_api_key_set(self):
        env=dict(os.environ,ALPHAGENOME_API_KEY="dummy-secret-marker")
        r=self.command("run","--mode","live","--fixture",str(FIXTURE),"--output-dir",str(self.out),env=env)
        self.assertNotEqual(r.returncode,0)
        self.assertFalse(self.out.exists())
        self.assertNotIn("dummy-secret-marker",r.stderr)
    def test_offline_with_api_key_set(self):
        env=dict(os.environ,ALPHAGENOME_API_KEY="dummy-secret-marker")
        r=self.command("run","--mode","offline","--fixture",str(FIXTURE),"--output-dir",str(self.out),env=env)
        self.assertEqual(r.returncode,0,r.stderr)
        self.assertEqual(json.loads(r.stdout)["real_api_calls_this_run"],0)
        self.assertNotIn("dummy-secret-marker",r.stdout)
    def test_invalid_mode_rejected(self):
        r=self.command("run","--mode","realtime","--fixture",str(FIXTURE),"--output-dir",str(self.out))
        self.assertNotEqual(r.returncode,0)
    def test_missing_fixture_rejected(self):
        r=self.command("run","--mode","offline","--fixture","/missing/d08.json","--output-dir",str(self.out))
        self.assertNotEqual(r.returncode,0)
    def test_sha_reproducible_in_two_locations(self):
        x=self.tmp/"a";y=self.tmp/"b"
        lab.run(FIXTURE,x);lab.run(FIXTURE,y)
        self.assertEqual({z.name:z.read_bytes() for z in x.iterdir()}, {z.name:z.read_bytes() for z in y.iterdir()})
    def test_bad_json_rejected(self):
        f=self.tmp/"broken.json";f.write_text("{oops")
        with self.assertRaises(ValueError):lab.run(f,self.out)


if __name__ == "__main__":
    unittest.main(verbosity=2)
