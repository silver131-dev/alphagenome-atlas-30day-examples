"""Unit tests: isolated synthetic context-window exercise; no network or API."""
from __future__ import annotations

import ast
import copy
import csv
import hashlib
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

import day06_offline_context_window_lab as lab

HERE=Path(__file__).resolve().parent
FIXTURE=HERE/"day06_synthetic_context.json"


def fixture():
    return json.loads(FIXTURE.read_text(encoding="utf-8"))


class Day06CoordinateTests(unittest.TestCase):
    def test_exact_model_width(self):self.assertEqual(lab.LONG_WIDTH,2**20)
    def test_toy_local_width(self):self.assertEqual(lab.LOCAL_WIDTH,100_000)
    def test_expanded_start(self):self.assertEqual(lab.centered_window(1_000_000,lab.LONG_WIDTH,2_000_000)["start_0based"],475_711)
    def test_expanded_end(self):self.assertEqual(lab.centered_window(1_000_000,lab.LONG_WIDTH,2_000_000)["end_0based_exclusive"],1_524_287)
    def test_expanded_1based_start(self):self.assertEqual(lab.centered_window(1_000_000,lab.LONG_WIDTH,2_000_000)["start_1based_inclusive"],475_712)
    def test_expanded_1based_last(self):self.assertEqual(lab.centered_window(1_000_000,lab.LONG_WIDTH,2_000_000)["end_1based_inclusive"],1_524_287)
    def test_local_start(self):self.assertEqual(lab.centered_window(1_000_000,lab.LOCAL_WIDTH,2_000_000)["start_1based_inclusive"],950_000)
    def test_local_end(self):self.assertEqual(lab.centered_window(1_000_000,lab.LOCAL_WIDTH,2_000_000)["end_1based_inclusive"],1_049_999)
    def test_first_expanded_base_included(self):self.assertTrue(lab.in_window(475_712,lab.centered_window(1_000_000,lab.LONG_WIDTH,2_000_000)))
    def test_one_before_expanded_excluded(self):self.assertFalse(lab.in_window(475_711,lab.centered_window(1_000_000,lab.LONG_WIDTH,2_000_000)))
    def test_last_expanded_base_included(self):self.assertTrue(lab.in_window(1_524_287,lab.centered_window(1_000_000,lab.LONG_WIDTH,2_000_000)))
    def test_one_after_expanded_excluded(self):self.assertFalse(lab.in_window(1_524_288,lab.centered_window(1_000_000,lab.LONG_WIDTH,2_000_000)))
    def test_local_boundary_left(self):self.assertTrue(lab.in_window(950_000,lab.centered_window(1_000_000,lab.LOCAL_WIDTH,2_000_000)))
    def test_local_before_left_excluded(self):self.assertFalse(lab.in_window(949_999,lab.centered_window(1_000_000,lab.LOCAL_WIDTH,2_000_000)))
    def test_local_boundary_right(self):self.assertTrue(lab.in_window(1_049_999,lab.centered_window(1_000_000,lab.LOCAL_WIDTH,2_000_000)))
    def test_local_after_right_excluded(self):self.assertFalse(lab.in_window(1_050_000,lab.centered_window(1_000_000,lab.LOCAL_WIDTH,2_000_000)))
    def test_anchor_at_first_bp(self):
        w=lab.centered_window(1,lab.LONG_WIDTH,2_000_000)
        self.assertEqual((w["start_0based"],w["end_0based_exclusive"]),(0,lab.LONG_WIDTH))
    def test_anchor_at_last_bp(self):
        w=lab.centered_window(2_000_000,lab.LONG_WIDTH,2_000_000)
        self.assertEqual((w["start_0based"],w["end_0based_exclusive"]),(2_000_000-lab.LONG_WIDTH,2_000_000))
    def test_window_exact_length_at_first(self):
        w=lab.centered_window(1,lab.LONG_WIDTH,2_000_000);self.assertEqual(w["end_0based_exclusive"]-w["start_0based"],lab.LONG_WIDTH)
    def test_window_exact_length_at_last(self):
        w=lab.centered_window(2_000_000,lab.LONG_WIDTH,2_000_000);self.assertEqual(w["end_0based_exclusive"]-w["start_0based"],lab.LONG_WIDTH)
    def test_width_too_large_fails(self):
        with self.assertRaises(ValueError):lab.centered_window(4,1000,999)
    def test_anchor_zero_fails(self):
        with self.assertRaises(ValueError):lab.centered_window(0,100,1000)
    def test_anchor_too_big_fails(self):
        with self.assertRaises(ValueError):lab.centered_window(1001,100,1000)
    def test_bool_is_not_integer_anchor(self):
        with self.assertRaises(ValueError):lab.centered_window(True,100,1000)
    def test_zero_width_fails(self):
        with self.assertRaises(ValueError):lab.centered_window(40,0,1000)
    def test_negative_position_fails(self):
        with self.assertRaises(ValueError):lab.in_window(-1,lab.centered_window(10,6,100))


class Day06FixtureTests(unittest.TestCase):
    def test_good_fixture(self): self.assertIsInstance(lab.validate_fixture(fixture()),dict)
    def test_non_dict(self):
        with self.assertRaises(ValueError):lab.validate_fixture([])
    def test_wrong_data_kind(self):
        d=fixture();d["evidence_class"]="HISTORICAL_LIVE_OUTPUT"
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_wrong_genome_build(self):
        d=fixture();d["assembly"]="GRCh38"
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_wrong_chromosome(self):
        d=fixture();d["chromosome"]="chr9"
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_small_chromosome_length(self):
        d=fixture();d["toy_chromosome_length_bp"]=100000
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_anchor_missing(self):
        d=fixture();d.pop("anchor_1based")
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_anchor_float(self):
        d=fixture();d["anchor_1based"]=1_000_000.0
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_anchor_past_end(self):
        d=fixture();d["anchor_1based"]=2_000_001
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_wrong_local_width(self):
        d=fixture();d["windows_bp"]["local"]=99_999
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_wrong_long_width(self):
        d=fixture();d["windows_bp"]["expanded"]=1_000_000
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_missing_features(self):
        d=fixture();d.pop("features")
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_empty_features(self):
        d=fixture();d["features"]=[]
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_feature_nonobject(self):
        d=fixture();d["features"][0]="bad"
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_bad_feature_id(self):
        d=fixture();d["features"][0]["id"]="<script>"
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_duplicate_feature_id(self):
        d=fixture();d["features"][1]["id"]="TARGET"
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_duplicate_feature_position(self):
        d=fixture();d["features"][1]["position_1based"]=1_000_000
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_bad_position_zero(self):
        d=fixture();d["features"][1]["position_1based"]=0
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_bad_position_beyond_chromosome(self):
        d=fixture();d["features"][1]["position_1based"]=2_000_001
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_bad_feature_type(self):
        d=fixture();d["features"][1]["kind"]="MEASURED_3D_CONTACT"
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_bad_feature_evidence_class(self):
        d=fixture();d["features"][1]["evidence_class"]="VALIDATED"
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_no_target_fails(self):
        d=fixture();d["features"][0]["kind"]="ENHANCER_CANDIDATE"
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_multiple_target_fails(self):
        d=fixture();d["features"][1]["kind"]="TARGET_GENE"
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_target_anchor_mismatch(self):
        d=fixture();d["anchor_1based"]+=1
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_contact_measured_true_rejected(self):
        d=fixture();d["contact_map_measured"]=True
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_contact_missing_rejected(self):
        d=fixture();d.pop("contact_map_measured")
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_experiment_performed_true_rejected(self):
        d=fixture();d["experimental_validation_performed"]=True
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_hypothesis_missing_rejected(self):
        d=fixture();d.pop("hypothetical_links")
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_hypothesis_unknown_candidate(self):
        d=fixture();d["hypothetical_links"][0]["candidate_id"]="UNKNOWN"
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_hypothesis_wrong_target_kind_rejected(self):
        d=fixture();d["hypothetical_links"][0]["target_id"]="NEAR_ENHANCER"
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_hypothesis_measured_status_rejected(self):
        d=fixture();d["hypothetical_links"][0]["status"]="MEASURED"
        with self.assertRaises(ValueError):lab.validate_fixture(d)
    def test_hypothesis_target_as_candidate_rejected(self):
        d=fixture();d["hypothetical_links"][0]["candidate_id"]="TARGET"
        with self.assertRaises(ValueError):lab.validate_fixture(d)


class Day06ReportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.report=lab.build_report(fixture());cls.data=fixture()
    def test_local_count(self):self.assertEqual(self.report["local_count"],2)
    def test_expanded_count(self):self.assertEqual(self.report["expanded_count"],5)
    def test_newly_seen_count(self):self.assertEqual(self.report["newly_seen_count"],3)
    def test_near_visible_in_both(self):self.assertIn("NEAR_ENHANCER",self.report["local_feature_ids"])
    def test_distal_not_local(self):self.assertNotIn("DISTAL_ENHANCER",self.report["local_feature_ids"])
    def test_distal_visible_expanded(self):self.assertIn("DISTAL_ENHANCER",self.report["expanded_feature_ids"])
    def test_outside_left_not_visible_expanded(self):self.assertNotIn("LEFT_EDGE_OUT",self.report["expanded_feature_ids"])
    def test_inside_left_visible_expanded(self):self.assertIn("LEFT_EDGE_IN",self.report["expanded_feature_ids"])
    def test_outside_right_not_visible_expanded(self):self.assertNotIn("RIGHT_EDGE_OUT",self.report["expanded_feature_ids"])
    def test_inside_right_visible_expanded(self):self.assertIn("RIGHT_EDGE_IN",self.report["expanded_feature_ids"])
    def test_feature_count(self):self.assertEqual(len(self.report["features"]),7)
    def test_no_model_prediction(self):self.assertIs(self.report["official_alphagenome_prediction"],False)
    def test_no_contact_measurement(self):self.assertIs(self.report["contact_map_measured"],False)
    def test_no_contact_inference(self):self.assertIs(self.report["chromatin_loop_inferred"],False)
    def test_no_gene_causality(self):self.assertIs(self.report["gene_regulation_established"],False)
    def test_no_experiment(self):self.assertIs(self.report["experiment_performed"],False)
    def test_no_clinical_interpretation(self):self.assertIs(self.report["clinical_interpretation"],False)
    def test_no_api(self):self.assertEqual(self.report["real_api_calls_this_run"],0)
    def test_hypothesis_label(self):self.assertEqual(self.report["hypothetical_links"][0]["status"],"HYPOTHESIS_NOT_MEASURED")
    def test_report_does_not_mutate_fixture(self):
        before=copy.deepcopy(self.data);lab.build_report(self.data);self.assertEqual(self.data,before)
    def test_distal_distance_400000(self):
        f=next(x for x in self.report["features"] if x["id"]=="DISTAL_ENHANCER")
        self.assertEqual(f["distance_to_target_bp"],400000)
    def test_csv_rows(self):self.assertEqual(len(list(csv.DictReader(io.StringIO(lab.csv_text(self.report))))),7)
    def test_json_null_and_booleans(self):
        x=json.loads(lab.json_text(self.report));self.assertIs(x["contact_map_measured"],False)
    def test_svg_parses(self):
        x=ET.fromstring(lab.svg_text(self.report));self.assertEqual(x.attrib["viewBox"],"0 0 1920 1080")
    def test_svg_synthetic_label(self):self.assertIn("NOT A PHYSICAL CONTACT MAP",lab.svg_text(self.report))
    def test_svg_background(self):self.assertIn("#E5DEDA",lab.svg_text(self.report))
    def test_html_has_table(self):
        x=lab.html_text(self.report,lab.svg_text(self.report));self.assertIn("<table>",x)
    def test_html_marks_synthetic(self):
        x=lab.html_text(self.report,lab.svg_text(self.report));self.assertIn("SYNTHETIC EDUCATIONAL DATA",x)
    def test_outputs_names_fixed(self):self.assertEqual(len(lab.make_outputs(self.data)[1]),4)
    def test_output_hash_deterministic(self):
        a=lab.make_outputs(self.data)[1];b=lab.make_outputs(self.data)[1]
        self.assertEqual({k:hashlib.sha256(v).hexdigest() for k,v in a.items()},{k:hashlib.sha256(v).hexdigest() for k,v in b.items()})


class Day06ArtifactTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix="day06-test-")
        self.root=Path(self.temp.name)
        self.source=self.root/FIXTURE.name
        self.source.write_bytes(FIXTURE.read_bytes())
        self.out=self.root/"output"
    def tearDown(self):self.temp.cleanup()
    def populate(self):return lab.run(self.source,self.out)
    def test_run_writes_five_files(self):
        self.populate();self.assertEqual(len(list(self.out.iterdir())),5)
    def test_verify_after_run(self):
        self.populate();self.assertEqual(lab.verify(self.source,self.out)["checksums"],"MATCH")
    def test_run_refuses_overwrite(self):
        self.populate()
        with self.assertRaises(FileExistsError):self.populate()
    def test_change_report_rejected(self):
        self.populate();(self.out/lab.ARTIFACTS[0]).write_text('altered',encoding='utf-8')
        with self.assertRaises(ValueError):lab.verify(self.source,self.out)
    def test_change_csv_rejected(self):
        self.populate();(self.out/lab.ARTIFACTS[1]).write_text('altered',encoding='utf-8')
        with self.assertRaises(ValueError):lab.verify(self.source,self.out)
    def test_change_svg_rejected(self):
        self.populate();(self.out/lab.ARTIFACTS[2]).write_text('altered',encoding='utf-8')
        with self.assertRaises(ValueError):lab.verify(self.source,self.out)
    def test_change_html_rejected(self):
        self.populate();(self.out/lab.ARTIFACTS[3]).write_text('altered',encoding='utf-8')
        with self.assertRaises(ValueError):lab.verify(self.source,self.out)
    def test_change_manifest_rejected(self):
        self.populate();(self.out/lab.MANIFEST_NAME).write_text('"tampered"',encoding='utf-8')
        with self.assertRaises(ValueError):lab.verify(self.source,self.out)
    def test_source_change_rejected(self):
        self.populate();self.source.write_text(self.source.read_text(encoding='utf-8')+'\n',encoding='utf-8')
        with self.assertRaises(ValueError):lab.verify(self.source,self.out)
    def test_manifest_is_not_signature(self):
        self.populate();m=json.loads((self.out/lab.MANIFEST_NAME).read_text(encoding="utf-8"))
        self.assertIs(m["manifest_is_digital_signature"],False)
    def test_manifest_source_hash(self):
        self.populate();m=json.loads((self.out/lab.MANIFEST_NAME).read_text(encoding="utf-8"))
        self.assertEqual(m["source_sha256"],hashlib.sha256(self.source.read_bytes()).hexdigest())
    def test_manifest_artifact_count(self):
        self.populate();m=json.loads((self.out/lab.MANIFEST_NAME).read_text(encoding="utf-8"))
        self.assertEqual(set(m["artifacts_sha256"]),set(lab.ARTIFACTS))
    def test_identical_separate_runs(self):
        self.populate();other=self.root/"other";lab.run(self.source,other)
        for filename in (*lab.ARTIFACTS,lab.MANIFEST_NAME):self.assertEqual((self.out/filename).read_bytes(),(other/filename).read_bytes())
    def test_run_returns_zero_api(self):self.assertEqual(self.populate()["real_api_calls_this_run"],0)
    def test_cli_offline_run(self):
        p=subprocess.run([sys.executable,str(HERE/"day06_offline_context_window_lab.py"),"run","--mode","offline","--fixture",str(self.source),"--output-dir",str(self.out)],capture_output=True,text=True,cwd="/tmp",check=True)
        self.assertEqual(json.loads(p.stdout)["status"],"PASS")
    def test_cli_offline_verify(self):
        self.populate();p=subprocess.run([sys.executable,str(HERE/"day06_offline_context_window_lab.py"),"verify","--mode","offline","--fixture",str(self.source),"--output-dir",str(self.out)],capture_output=True,text=True,check=True)
        self.assertEqual(json.loads(p.stdout)["files_verified"],4)
    def test_cli_live_rejected_with_key_set(self):
        env=os.environ.copy();env["ALPHAGENOME_API_KEY"]="DO_NOT_USE_FAKE_TEST_VALUE"
        p=subprocess.run([sys.executable,str(HERE/"day06_offline_context_window_lab.py"),"run","--mode","live","--fixture",str(self.source),"--output-dir",str(self.out)],env=env,capture_output=True,text=True)
        self.assertNotEqual(p.returncode,0)
        self.assertFalse(self.out.exists())
        self.assertNotIn(env["ALPHAGENOME_API_KEY"],p.stderr+p.stdout)
    def test_cli_missing_fixture(self):
        p=subprocess.run([sys.executable,str(HERE/"day06_offline_context_window_lab.py"),"run","--fixture",str(self.root/"missing.json"),"--output-dir",str(self.out)],capture_output=True,text=True)
        self.assertNotEqual(p.returncode,0)
        self.assertFalse(self.out.exists())
    def test_invalid_json_source(self):
        self.source.write_text('{oops',encoding="utf-8")
        with self.assertRaises(ValueError):lab.run(self.source,self.out)
    def test_does_not_import_network_or_alphagenome(self):
        root=ast.parse((HERE/"day06_offline_context_window_lab.py").read_text(encoding="utf-8"))
        imports=set()
        for n in ast.walk(root):
            if isinstance(n,ast.Import):imports.update(a.name.split('.')[0] for a in n.names)
            elif isinstance(n,ast.ImportFrom) and n.module:imports.add(n.module.split('.')[0])
        self.assertFalse(imports & {"alphagenome","requests","httpx","urllib","socket","aiohttp","pysam"})
    def test_report_contains_no_raw_patient_identifiers(self):
        self.populate();t=(self.out/lab.ARTIFACTS[0]).read_text(encoding="utf-8")
        self.assertNotIn("patient_id",t);self.assertNotIn("full_name",t)
    def test_cli_is_standalone(self):
        p=subprocess.run([sys.executable,str(HERE/"day06_offline_context_window_lab.py"),"--help"],capture_output=True,text=True,cwd="/tmp")
        self.assertEqual(p.returncode,0)


if __name__=="__main__":unittest.main()
