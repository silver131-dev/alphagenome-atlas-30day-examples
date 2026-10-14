"""Day 03 offline tests; independent of AlphaGenome/Gemini SDK or network."""
from __future__ import annotations
import contextlib
import copy
import csv
import io
import json
import math
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import day03_offline_two_lenses_lab as d

HERE=Path(__file__).resolve().parent
FIXTURE=HERE/"day03_synthetic_two_lenses.json"


def fixture():
    return json.loads(FIXTURE.read_text(encoding="utf-8"))


def report():
    return d.build_report(fixture(),fixture_name=FIXTURE.name,fixture_sha256=d.sha256_file(FIXTURE))


class FixtureContractTests(unittest.TestCase):
    def test_valid_fixture(self):self.assertEqual(len(d.validate_fixture(fixture())["variants"]),7)
    def test_not_object(self):
        with self.assertRaises(ValueError):d.validate_fixture([])
    def test_only_synthetic_kind(self):
        data=fixture();data["data_kind"]="PATIENT"
        with self.assertRaises(PermissionError):d.validate_fixture(data)
    def test_toy_build_required(self):
        data=fixture();data["genome_build"]="GRCh38"
        with self.assertRaises(PermissionError):d.validate_fixture(data)
    def test_dataset_id_required(self):
        data=fixture();data["dataset_id"]="UNKNOWN"
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_disclaimer_required(self):
        data=fixture();data["dataset_note"]="clinical score"
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_missing_top_level_key(self):
        data=fixture();del data["teaching_thresholds"]
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_unexpected_top_level_key(self):
        data=fixture();data["patient_id"]="P001"
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_empty_cases(self):
        data=fixture();data["variants"]=[]
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_too_many_cases(self):
        data=fixture();data["variants"]=[data["variants"][0]]*51
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_variants_must_be_list(self):
        data=fixture();data["variants"]={}
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_thresholds_must_be_object(self):
        data=fixture();data["teaching_thresholds"]=[2,-2,.25]
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_thresholds_exact_keys(self):
        data=fixture();data["teaching_thresholds"]["clinical_risk"]=99
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_threshold_positive_zero_rejected(self):
        data=fixture();data["teaching_thresholds"]["conservation_positive_min"]=0
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_threshold_negative_positive_rejected(self):
        data=fixture();data["teaching_thresholds"]["acceleration_negative_max"]=3
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_threshold_abs_delta_gt_one(self):
        data=fixture();data["teaching_thresholds"]["molecular_abs_delta_min"]=2
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_threshold_none_rejected(self):
        data=fixture();data["teaching_thresholds"]["molecular_abs_delta_min"]=None
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_threshold_bool_rejected(self):
        data=fixture();data["teaching_thresholds"]["conservation_positive_min"]=True
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_threshold_nan_rejected(self):
        data=fixture();data["teaching_thresholds"]["conservation_positive_min"]=math.nan
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_threshold_inf_rejected(self):
        data=fixture();data["teaching_thresholds"]["acceleration_negative_max"]=-math.inf
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_invalid_case_extra_key(self):
        data=fixture();data["variants"][0]["diagnosis"]="yes"
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_invalid_case_missing_key(self):
        data=fixture();del data["variants"][0]["biosample"]
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_invalid_chr_name(self):
        data=fixture();data["variants"][0]["variant_id"]="chr9:100:A>G"
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_zero_coordinate(self):
        data=fixture();data["variants"][0]["variant_id"]="chrToy:0:A>G"
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_duplicate_variant(self):
        data=fixture();data["variants"][1]["variant_id"]=data["variants"][0]["variant_id"]
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_unchanged_allele(self):
        data=fixture();data["variants"][0]["variant_id"]="chrToy:100:A>A"
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_bad_biosample(self):
        data=fixture();data["variants"][0]["biosample"]="patient_neuron"
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_note_must_be_text(self):
        data=fixture();data["variants"][0]["teaching_note"]=123
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_note_capped(self):
        data=fixture();data["variants"][0]["teaching_note"]="a"*121
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_nonsensical_phyloP_range(self):
        data=fixture();data["variants"][0]["toy_phyloP_like"]=99.0
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_molecular_delta_range(self):
        data=fixture();data["variants"][0]["toy_splicing_delta"]=1.0001
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_molecular_inf_rejected(self):
        data=fixture();data["variants"][0]["toy_splicing_delta"]=math.inf
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_phyloP_bool_rejected(self):
        data=fixture();data["variants"][0]["toy_phyloP_like"]=False
        with self.assertRaises(ValueError):d.validate_fixture(data)
    def test_null_is_allowed(self):self.assertIsNone(d.validate_fixture(fixture())["variants"][4]["toy_phyloP_like"])
    def test_toy_no_clinical_identifiers(self):self.assertNotIn("patient_id",dumps_data(fixture()))


def dumps_data(data):return json.dumps(data,ensure_ascii=False)


class ClassificationTests(unittest.TestCase):
    def setUp(self):
        self.report=report()
        self.byid={v["variant_id"]:v for v in self.report["case_rows"]}
    def at(self, pos):return self.byid[next(k for k in self.byid if k.split(":")[1]==str(pos))]
    def test_total_count(self):self.assertEqual(self.report["variant_count"],7)
    def test_both(self):self.assertEqual(self.at(102)["evidence_class"],"BOTH_SUGGESTIVE")
    def test_conservation_only(self):self.assertEqual(self.at(100)["evidence_class"],"CONSERVATION_ONLY")
    def test_molecular_only(self):self.assertEqual(self.at(101)["evidence_class"],"MOLECULAR_ONLY")
    def test_neither(self):self.assertEqual(self.at(103)["evidence_class"],"NEITHER_TOY_SIGNAL")
    def test_acceleration(self):self.assertEqual(self.at(105)["evidence_class"],"ACCELERATION_FLAG")
    def test_missing_conservation(self):self.assertEqual(self.at(104)["evidence_class"],"INCOMPLETE_EVIDENCE")
    def test_missing_molecular(self):self.assertEqual(self.at(106)["evidence_class"],"INCOMPLETE_EVIDENCE")
    def test_missing_counts(self):self.assertEqual((self.report["missing_conservation_count"],self.report["missing_molecular_count"]),(1,1))
    def test_null_never_zero(self):self.assertIsNone(self.at(104)["toy_phyloP_like"])
    def test_negative_delta_preserved(self):self.assertEqual(self.at(101)["toy_splicing_delta"],-.42)
    def test_acceleration_not_functionless(self):self.assertIn("不能推定無功能",self.at(105)["next_research_question"])
    def test_both_not_causal(self):self.assertEqual(self.at(102)["manual_review"],"REQUIRED")
    def test_no_patient_data(self):self.assertFalse(self.report["clinical_inference"])
    def test_no_prediction(self):self.assertFalse(self.report["model_prediction"])
    def test_no_api(self):self.assertEqual(self.report["real_api_calls_this_run"],0)
    def test_class_counts(self):self.assertEqual(sum(self.report["classification_counts"].values()),7)
    def test_all_statuses(self):self.assertTrue(all(v["evidence_status"]=="SYNTHETIC_UNVALIDATED" for v in self.report["case_rows"]))
    def test_no_official_threshold(self):self.assertFalse(self.report["thresholds_are_official"])
    def test_exact_threshold_conservation(self):
        cfg=d.validate_fixture(fixture())["teaching_thresholds"];rec=copy.deepcopy(fixture()["variants"][3]);rec["toy_phyloP_like"]=2.0
        self.assertEqual(d.classify(rec,cfg)["evidence_class"],"CONSERVATION_ONLY")
    def test_exact_threshold_molecular(self):
        cfg=d.validate_fixture(fixture())["teaching_thresholds"];rec=copy.deepcopy(fixture()["variants"][3]);rec["toy_splicing_delta"]=-.25
        self.assertEqual(d.classify(rec,cfg)["evidence_class"],"MOLECULAR_ONLY")
    def test_exact_negative_acceleration(self):
        cfg=d.validate_fixture(fixture())["teaching_thresholds"];rec=copy.deepcopy(fixture()["variants"][3]);rec["toy_phyloP_like"]=-2.0
        self.assertEqual(d.classify(rec,cfg)["evidence_class"],"ACCELERATION_FLAG")
    def test_zero_is_not_null(self):
        cfg=d.validate_fixture(fixture())["teaching_thresholds"];rec=copy.deepcopy(fixture()["variants"][3]);rec["toy_splicing_delta"]=0.0
        self.assertEqual(d.classify(rec,cfg)["evidence_class"],"NEITHER_TOY_SIGNAL")
    def test_null_signal_boolean_undetermined(self):self.assertIsNone(self.at(104)["conservation_signal"])
    def test_missing_molecular_boolean_undetermined(self):self.assertIsNone(self.at(106)["molecular_signal"])
    def test_no_finite_probability(self):self.assertNotIn("pathogenic_probability",self.report)
    def test_source_fingerprint(self):self.assertEqual(self.report["fixture_sha256"],d.sha256_file(FIXTURE))


class RenderingTests(unittest.TestCase):
    def setUp(self):self.r=report();self.a=d.render_outputs(self.r)
    def test_artifact_count(self):self.assertEqual(len(self.a),5)
    def test_artifact_names(self):self.assertEqual(set(self.a),set(d.OUTPUT_FILES))
    def test_deterministic(self):self.assertEqual(self.a,d.render_outputs(self.r))
    def test_json_valid(self):self.assertEqual(json.loads(self.a[d.OUTPUT_FILES[0]])["variant_count"],7)
    def test_no_literal_nan(self):self.assertNotIn(b'NaN',self.a[d.OUTPUT_FILES[0]])
    def test_json_null(self):self.assertIn(b'"toy_phyloP_like": null',self.a[d.OUTPUT_FILES[0]])
    def test_csv_header(self):
        rows=list(csv.reader(io.StringIO(self.a[d.OUTPUT_FILES[1]].decode('utf-8'))))
        self.assertEqual(rows[0],list(d.CSV_COLUMNS))
    def test_csv_seven(self):
        rows=list(csv.reader(io.StringIO(self.a[d.OUTPUT_FILES[1]].decode('utf-8'))))
        self.assertEqual(len(rows)-1,7)
    def test_csv_empty_cell(self):
        rows=list(csv.DictReader(io.StringIO(self.a[d.OUTPUT_FILES[1]].decode('utf-8'))))
        self.assertEqual(rows[4]["toy_phyloP_like"],"")
    def test_svg_dimensions(self):self.assertIn(b'width="1920" height="1080"',self.a[d.OUTPUT_FILES[2]])
    def test_svg_background(self):self.assertIn(b'#E5DEDA',self.a[d.OUTPUT_FILES[2]])
    def test_svg_synthetic(self):self.assertIn(b'SYNTHETIC ONLY',self.a[d.OUTPUT_FILES[2]])
    def test_svg_seven_cards(self):self.assertEqual(self.a[d.OUTPUT_FILES[2]].count(b'phyloP-like  '),7)
    def test_svg_contains_ids(self):self.assertIn(b'chrToy:101:C&gt;T',self.a[d.OUTPUT_FILES[2]])
    def test_html_offline(self):self.assertIn(b'SYNTHETIC EDUCATIONAL DEMO',self.a[d.OUTPUT_FILES[3]])
    def test_html_interactivity(self):self.assertIn(b'button[data-filter]',self.a[d.OUTPUT_FILES[3]])
    def test_html_no_external_code(self):
        s=self.a[d.OUTPUT_FILES[3]].decode("utf-8")
        self.assertNotIn('src="http',s);self.assertNotIn('fetch(',s)
    def test_html_seven_table_rows(self):self.assertEqual(self.a[d.OUTPUT_FILES[3]].count(b'<tr data-kind='),7)
    def test_html_no_patient(self):self.assertNotIn(b'patient_id',self.a[d.OUTPUT_FILES[3]])
    def test_markdown_contains_all_cases(self):
        for row in self.r["case_rows"]:self.assertIn(row["variant_id"],self.a[d.OUTPUT_FILES[4]].decode("utf-8"))
    def test_markdown_no_real_api(self):self.assertIn("不呼叫外部 API",self.a[d.OUTPUT_FILES[4]].decode("utf-8"))
    def test_all_utf8(self):
        for blob in self.a.values(): self.assertIsInstance(blob.decode("utf-8"),str)
    def test_counterfactual_no_clinical(self):
        text=self.a[d.OUTPUT_FILES[4]].decode("utf-8")
        self.assertIn("臨床",text)


class EndToEndTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.output=self.root/"day03_output"
        self.fixture=self.root/FIXTURE.name
        self.fixture.write_bytes(FIXTURE.read_bytes())
    def make(self):return d.run(self.fixture,self.output)
    def test_run_pass(self):self.assertEqual(self.make()["status"],"PASS")
    def test_verify_pass(self):self.make();self.assertEqual(d.verify(self.fixture,self.output)["checksums"],"MATCH")
    def test_run_outputs_count(self):self.assertEqual(self.make()["artifact_count"],5)
    def test_run_never_network(self):self.assertEqual(self.make()["real_api_calls_this_run"],0)
    def test_five_outputs_plus_manifest(self):self.make();self.assertEqual(len(list(self.output.iterdir())),6)
    def test_source_sha_recorded(self):
        self.make();obj=json.loads((self.output/"day03_manifest.json").read_text())
        self.assertEqual(obj["source_sha256"],d.sha256_file(self.fixture))
    def test_unsigned_manifest(self):
        self.make();obj=json.loads((self.output/"day03_manifest.json").read_text())
        self.assertFalse(obj["manifest_is_signed"])
    def test_avoid_overwrite(self):
        self.make()
        with self.assertRaises(FileExistsError):self.make()
    def test_live_run_denied(self):
        with self.assertRaises(PermissionError):d.run(self.fixture,self.output,mode="live")
    def test_live_verify_denied(self):
        with self.assertRaises(PermissionError):d.verify(self.fixture,self.output,mode="live")
    def test_environment_key_does_not_enable_live(self):
        with mock.patch.dict(os.environ,{"ALPHA_GENOME_API_KEY":"SENSITIVE_DUMMY_NEVER_SHOW"}):
            status=self.make()
            self.assertNotIn("SENSITIVE_DUMMY",json.dumps(status))
    def test_no_sdk_imports(self):
        src=Path(d.__file__).read_text()
        self.assertNotIn("import alphagenome",src)
        self.assertNotIn("import requests",src)
    def test_modified_source_detected(self):
        self.make();self.fixture.write_text(self.fixture.read_text()+"\n")
        with self.assertRaisesRegex(ValueError,"Source"):d.verify(self.fixture,self.output)
    def test_different_fixture_filename_detected(self):
        self.make();renamed=self.root/"new.json";renamed.write_bytes(self.fixture.read_bytes())
        with self.assertRaisesRegex(ValueError,"Source"):d.verify(renamed,self.output)
    def test_missing_artifact_detected(self):
        self.make();(self.output/d.OUTPUT_FILES[1]).unlink()
        with self.assertRaises(ValueError):d.verify(self.fixture,self.output)
    def test_added_file_detected(self):
        self.make();(self.output/"extra.csv").write_text("x")
        with self.assertRaisesRegex(ValueError,"Unexpected"):d.verify(self.fixture,self.output)
    def test_altered_json_detected(self):
        self.make();p=self.output/d.OUTPUT_FILES[0];p.write_bytes(p.read_bytes()+b' ')
        with self.assertRaisesRegex(ValueError,"Artifact altered"):d.verify(self.fixture,self.output)
    def test_altered_csv_detected(self):
        self.make();p=self.output/d.OUTPUT_FILES[1];p.write_bytes(p.read_bytes()+b' ')
        with self.assertRaises(ValueError):d.verify(self.fixture,self.output)
    def test_altered_svg_detected(self):
        self.make();p=self.output/d.OUTPUT_FILES[2];p.write_bytes(p.read_bytes()+b' ')
        with self.assertRaises(ValueError):d.verify(self.fixture,self.output)
    def test_altered_html_detected(self):
        self.make();p=self.output/d.OUTPUT_FILES[3];p.write_bytes(p.read_bytes()+b' ')
        with self.assertRaises(ValueError):d.verify(self.fixture,self.output)
    def test_altered_markdown_detected(self):
        self.make();p=self.output/d.OUTPUT_FILES[4];p.write_bytes(p.read_bytes()+b' ')
        with self.assertRaises(ValueError):d.verify(self.fixture,self.output)
    def test_manifest_artifact_count(self):
        self.make();p=self.output/"day03_manifest.json";m=json.loads(p.read_text());m["artifact_count"]=42;p.write_text(json.dumps(m))
        with self.assertRaisesRegex(ValueError,"Artifact list"):d.verify(self.fixture,self.output)
    def test_manifest_signature_false_only(self):
        self.make();p=self.output/"day03_manifest.json";m=json.loads(p.read_text());m["manifest_is_signed"]=True;p.write_text(json.dumps(m))
        with self.assertRaisesRegex(ValueError,"Invalid unsigned"):d.verify(self.fixture,self.output)
    def test_manifest_checksum_forgery_does_not_pass(self):
        self.make();p=self.output/d.OUTPUT_FILES[0];p.write_bytes(p.read_bytes()+b' ')
        mpath=self.output/"day03_manifest.json";m=json.loads(mpath.read_text());m["artifacts_sha256"][d.OUTPUT_FILES[0]]=d.sha256_file(p)
        mpath.write_text(json.dumps(m))
        with self.assertRaisesRegex(ValueError,"diverges"):d.verify(self.fixture,self.output)
    def test_cli_run(self):
        p=subprocess.run([sys.executable,"-B",str(Path(d.__file__)),"run","--fixture",str(self.fixture),"--output-dir",str(self.output)],capture_output=True,text=True,check=True)
        self.assertEqual(json.loads(p.stdout)["variant_count"],7)
    def test_cli_verify(self):
        self.make();p=subprocess.run([sys.executable,"-B",str(Path(d.__file__)),"verify","--fixture",str(self.fixture),"--output-dir",str(self.output)],capture_output=True,text=True,check=True)
        self.assertEqual(json.loads(p.stdout)["files_verified"],5)
    def test_cli_live_fails(self):
        p=subprocess.run([sys.executable,str(Path(d.__file__)),"run","--mode","live","--fixture",str(self.fixture),"--output-dir",str(self.output)],capture_output=True,text=True)
        self.assertNotEqual(p.returncode,0);self.assertFalse(self.output.exists())
    def test_cli_help(self):
        p=subprocess.run([sys.executable,str(Path(d.__file__)),"--help"],capture_output=True,text=True,check=True)
        self.assertIn("offline only",p.stdout)
    def test_two_runs_byte_identical(self):
        self.make();another=self.root/"another_run";d.run(self.fixture,another)
        for name in (*d.OUTPUT_FILES,"day03_manifest.json"):
            self.assertEqual((self.output/name).read_bytes(),(another/name).read_bytes())

if __name__=="__main__":unittest.main()
