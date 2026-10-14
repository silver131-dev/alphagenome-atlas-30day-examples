"""Offline synthetic Day29 privacy gate tests; no patient data or network."""
from __future__ import annotations

import ast
import copy
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import day29_offline_privacy_gate_lab as lab

BASE = Path(__file__).resolve().parent
FIXTURE = BASE / "fixtures/day28_evidence_report.json"
MANIFEST = BASE / "fixtures/day28_report_manifest.json"
PROGRAM = BASE / "day29_offline_privacy_gate_lab.py"


def load():
    return json.loads(FIXTURE.read_text(encoding="utf-8"))


class Day29SchemaTests(unittest.TestCase):
    def setUp(self):
        self.report = load()

    def test_day28_actual_fixture_accepted(self):
        result = lab.validate_day28(self.report)
        self.assertEqual(result["candidate_count"], 8)
        self.assertEqual(result["selected_count"], 3)
        self.assertEqual(result["missing_score_count"], 1)

    def test_day28_unsigned_manifest_sha_matches(self):
        counts, prov = lab.load_and_validate(FIXTURE, MANIFEST)
        self.assertEqual(counts["day27_artifact_count_claimed_by_day28"], 23)
        self.assertEqual(prov["day28_report_sha256"],
                         hashlib.sha256(FIXTURE.read_bytes()).hexdigest())
        self.assertFalse(prov["day27_source_artifacts_independently_checked_this_run"])

    def test_missing_source(self):
        with self.assertRaises(ValueError):
            lab.load_and_validate(BASE / "missing.json", MANIFEST)

    def test_reject_extra_root_field(self):
        self.report["new_source_column"] = "synthetic"
        with self.assertRaises(ValueError): lab.validate_day28(self.report)

    def test_reject_patient_id_nested(self):
        self.report["annotations"][0]["patient_id"] = "SYNTHETIC_EXAMPLE"
        with self.assertRaisesRegex(ValueError, "DIRECT_IDENTIFIER_FIELD_NOT_ALLOWED"):
            lab.validate_day28(self.report)

    def test_reject_full_name_root(self):
        self.report["full_name"] = "Invented Person"
        with self.assertRaises(ValueError): lab.validate_day28(self.report)

    def test_reject_name_nested(self):
        self.report["source_provenance"]["name"] = "Invented"
        with self.assertRaises(ValueError): lab.validate_day28(self.report)

    def test_reject_non_synthetic(self):
        self.report["evidence_class"] = "REAL_PATIENT_DATA"
        with self.assertRaises(ValueError): lab.validate_day28(self.report)

    def test_reject_grch38_marker(self):
        self.report["assembly"] = "GRCh38"
        with self.assertRaises(ValueError): lab.validate_day28(self.report)

    def test_reject_live_mode(self):
        self.report["mode"] = "LIVE"
        with self.assertRaises(ValueError): lab.validate_day28(self.report)

    def test_reject_any_api_calls(self):
        self.report["real_api_calls_this_run"] = 1
        with self.assertRaises(ValueError): lab.validate_day28(self.report)

    def test_reject_boolean_api_call_fake_zero(self):
        self.report["real_api_calls_this_run"] = False
        with self.assertRaises(ValueError): lab.validate_day28(self.report)

    def test_reject_model_prediction(self):
        self.report["model_prediction"] = True
        with self.assertRaises(ValueError): lab.validate_day28(self.report)

    def test_reject_non_null_clinical_classification(self):
        self.report["clinical_classification"] = "pathogenic"
        with self.assertRaises(ValueError): lab.validate_day28(self.report)

    def test_reject_disease_probability(self):
        self.report["disease_probability"] = 0.9
        with self.assertRaises(ValueError): lab.validate_day28(self.report)

    def test_reject_signed_provenance_mismatch(self):
        self.report["source_provenance"]["signature_verified"] = True
        with self.assertRaises(ValueError): lab.validate_day28(self.report)

    def test_reject_wrong_day27_kind(self):
        self.report["source_provenance"]["source_kind"] = "REAL_API"
        with self.assertRaises(ValueError): lab.validate_day28(self.report)

    def test_reject_wrong_artifact_count(self):
        self.report["source_provenance"]["artifact_count_verified"] = 22
        with self.assertRaises(ValueError): lab.validate_day28(self.report)

    def test_reject_bad_sha_label(self):
        self.report["source_provenance"]["source_bundle_sha256"] = "bad"
        with self.assertRaises(ValueError): lab.validate_day28(self.report)

    def test_reject_mismatched_candidate_count(self):
        self.report["pipeline_counts"]["candidate_count"] = 9
        with self.assertRaises(ValueError): lab.validate_day28(self.report)

    def test_reject_missing_review_annotation(self):
        self.report["annotations"][1]["review_status"] = "APPROVED"
        with self.assertRaises(ValueError): lab.validate_day28(self.report)

    def test_reject_duplicate_variant(self):
        self.report["annotations"][1]["variant_id"] = self.report["annotations"][0]["variant_id"]
        with self.assertRaises(ValueError): lab.validate_day28(self.report)

    def test_reject_fake_missing_score(self):
        self.report["annotations"][2]["synthetic_score"] = 0.0
        with self.assertRaises(ValueError): lab.validate_day28(self.report)

    def test_zero_is_valid_non_missing_teaching_value(self):
        self.report["annotations"][2]["synthetic_score"] = 0.0
        self.report["annotations"][2]["score_status"] = "SYNTHETIC_SCORE"
        self.report["pipeline_counts"]["annotated_count"] = 3
        self.report["pipeline_counts"]["missing_score_count"] = 0
        self.report["missing_values"]["missing_score_variant_ids"] = []
        self.assertEqual(lab.validate_day28(self.report)["missing_score_count"], 0)

    def test_reject_nan_score(self):
        self.report["annotations"][0]["synthetic_score"] = float("nan")
        with self.assertRaises(ValueError): lab.validate_day28(self.report)

    def test_reject_bad_variant_identifiers(self):
        self.report["annotations"][0]["variant_id"] = "../../../patient.json"
        with self.assertRaises(ValueError): lab.validate_day28(self.report)

    def test_reject_combined_biological_claim(self):
        self.report["teaching_observations"]["combined_biological_score"] = 0.7
        with self.assertRaises(ValueError): lab.validate_day28(self.report)

    def test_reject_missing_scientific_limitation(self):
        self.report["model_limitations"] = [e for e in self.report["model_limitations"]
                                              if e["code"] != "PROVENANCE_NOT_SIGNED"]
        with self.assertRaises(ValueError): lab.validate_day28(self.report)


class Day29OutputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.folder = Path(self.temp.name)

    def tearDown(self): self.temp.cleanup()

    def create(self, name="out"):
        dst = self.folder/name
        result = lab.run(FIXTURE, MANIFEST, dst)
        return dst, result

    def test_run_creates_expected_files(self):
        dst, result = self.create()
        self.assertEqual(result["output_artifact_count"], 6)
        self.assertTrue(all((dst/fn).is_file() for fn in lab.ARTIFACT_FILES))
        self.assertTrue((dst/"day29_report_manifest.json").is_file())

    def test_privacy_gate_always_blocked(self):
        dst, result = self.create()
        audit = json.loads((dst/"day29_privacy_audit.json").read_text())
        self.assertEqual(result["release_gate"], "BLOCKED")
        self.assertEqual(audit["governance_review_status"], "MANUAL_REVIEW_REQUIRED")
        self.assertFalse(audit["demo_policy"]["external_transfer_allowed"])
        self.assertFalse(audit["demo_policy"]["publish_allowed"])
        self.assertFalse(audit["demo_policy"]["human_approval_completed"])
        self.assertFalse(audit["demo_policy"]["irb_or_ethics_review_verified"])

    def test_token_not_anonymization(self):
        dst, _ = self.create()
        audit = json.loads((dst/"day29_privacy_audit.json").read_text())
        self.assertRegex(audit["local_demo_token"], lab.DEMO_TOKEN_RE)
        self.assertFalse(audit["token_is_anonymization"])
        self.assertFalse(audit["reidentification_risk_eliminated"])

    def test_token_changes_each_run(self):
        a, _ = self.create("a")
        b, _ = self.create("b")
        x = json.loads((a/"day29_privacy_audit.json").read_text())
        y = json.loads((b/"day29_privacy_audit.json").read_text())
        self.assertNotEqual(x["local_demo_token"], y["local_demo_token"])

    def test_preview_excludes_variants_tokens_and_scores(self):
        dst, _ = self.create()
        for fn in ("day29_local_preview.json", "day29_privacy_dashboard.html",
                   "day29_privacy_flow.svg", "day29_release_decision.csv"):
            contents = (dst/fn).read_text(encoding="utf-8")
            for forbidden in ("chr1:", "DEMO_", "0.42", "-0.17", "patient_id", "full_name"):
                self.assertNotIn(forbidden, contents, fn)

    def test_html_preview_privacy_claim(self):
        dst, _ = self.create()
        contents=(dst/"day29_privacy_dashboard.html").read_text()
        self.assertIn("不是匿名化", contents)
        self.assertIn("外傳：禁止", contents)

    def test_svg_is_16_by_9_and_correct_color(self):
        dst, _ = self.create()
        from xml.etree import ElementTree as ET
        root=ET.fromstring((dst/"day29_privacy_flow.svg").read_text())
        self.assertEqual(root.attrib["width"], "1920")
        self.assertEqual(root.attrib["height"], "1080")
        self.assertTrue(any(e.attrib.get("fill")=="#E5DEDA" for e in root.iter()))

    def test_markdown_contains_limitations(self):
        dst, _ = self.create()
        contents=(dst/"day29_privacy_report.md").read_text()
        self.assertIn("不能證明資料來源", contents)
        self.assertIn("不保證基因資料無法被重新識別", contents)
        self.assertIn("本次並未逐一重驗", contents)

    def test_manifest_verify_succeeds(self):
        dst, _ = self.create()
        self.assertEqual(lab.verify(FIXTURE, MANIFEST, dst)["verified_artifacts"], 6)

    def test_verifier_catches_tamper(self):
        dst, _ = self.create()
        with (dst/"day29_local_preview.json").open("a") as f: f.write("TAMPER")
        with self.assertRaises(ValueError): lab.verify(FIXTURE, MANIFEST, dst)

    def test_verifier_catches_missing_file(self):
        dst, _ = self.create()
        (dst/"day29_privacy_report.md").unlink()
        with self.assertRaises(ValueError): lab.verify(FIXTURE, MANIFEST, dst)

    def test_verifier_catches_manifest_source_hash_change(self):
        dst, _ = self.create()
        path=dst/"day29_report_manifest.json"
        obj=json.loads(path.read_text())
        obj["source_report_sha256"]="0"*64
        path.write_text(json.dumps(obj))
        with self.assertRaises(ValueError): lab.verify(FIXTURE, MANIFEST, dst)

    def test_refuses_existing_output(self):
        dst, _=self.create()
        before=hashlib.sha256((dst/"day29_privacy_audit.json").read_bytes()).hexdigest()
        with self.assertRaises(ValueError): lab.run(FIXTURE, MANIFEST, dst)
        after=hashlib.sha256((dst/"day29_privacy_audit.json").read_bytes()).hexdigest()
        self.assertEqual(before, after)

    def test_preview_stable_across_runs_even_when_token_changes(self):
        a, _=self.create("a")
        b, _=self.create("b")
        for fn in ("day29_local_preview.json", "day29_privacy_dashboard.html",
                   "day29_privacy_flow.svg", "day29_release_decision.csv"):
            self.assertEqual((a/fn).read_bytes(), (b/fn).read_bytes(), fn)

    def test_source_not_modified(self):
        before=hashlib.sha256(FIXTURE.read_bytes()).hexdigest()
        self.create()
        self.assertEqual(before, hashlib.sha256(FIXTURE.read_bytes()).hexdigest())

    def test_rejects_sha_mismatched_report(self):
        dst=self.folder/"changed.json"
        obj=load()
        obj["pipeline_counts"]["candidate_count"] = 99
        dst.write_text(json.dumps(obj), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "DAY28_INPUT_HASH_MISMATCH"):
            lab.load_and_validate(dst, MANIFEST)

    def test_rejects_manifest_claimed_signed_true(self):
        bad=self.folder/"manifest.json"
        obj=json.loads(MANIFEST.read_text())
        obj["signed"]=True
        bad.write_text(json.dumps(obj), encoding="utf-8")
        with self.assertRaises(ValueError): lab.load_and_validate(FIXTURE,bad)

    def test_rejects_patient_data_even_if_rehash_manifest(self):
        report=self.folder/"bad_day28.json"
        manifest=self.folder/"bad_manifest.json"
        obj=load()
        obj["patient_id"] = "FICTIONAL-001"
        data=(json.dumps(obj,sort_keys=True)+"\n").encode()
        report.write_bytes(data)
        m=json.loads(MANIFEST.read_text())
        m["artifacts"]["day28_evidence_report.json"]["sha256"]=lab.sha256_bytes(data)
        m["artifacts"]["day28_evidence_report.json"]["bytes"]=len(data)
        manifest.write_text(json.dumps(m))
        with self.assertRaises(ValueError): lab.load_and_validate(report, manifest)

    def test_risk_decisions_do_not_depend_on_api_environment_variables(self):
        original=os.environ.get("ALPHAGENOME_API_KEY")
        try:
            os.environ["ALPHAGENOME_API_KEY"]="DUMMY_NOT_A_REAL_KEY"
            dst,result=self.create()
            self.assertFalse(result["external_transfer_allowed"])
            self.assertEqual(result["real_api_calls_this_run"], 0)
        finally:
            if original is None: os.environ.pop("ALPHAGENOME_API_KEY",None)
            else: os.environ["ALPHAGENOME_API_KEY"]=original


class Day29CliTests(unittest.TestCase):
    def test_cli_run_and_verify(self):
        with tempfile.TemporaryDirectory() as d:
            dst=Path(d)/"out"
            base=[sys.executable,str(PROGRAM)]
            args=["--day28-report",str(FIXTURE),"--day28-manifest",str(MANIFEST),"--output-dir",str(dst)]
            run=subprocess.run(base+["run","--mode","offline"]+args, capture_output=True, text=True)
            self.assertEqual(run.returncode, 0,run.stderr)
            self.assertEqual(json.loads(run.stdout)["release_gate"],"BLOCKED")
            verify=subprocess.run(base+["verify"]+args, capture_output=True, text=True)
            self.assertEqual(verify.returncode,0,verify.stderr)
            self.assertEqual(json.loads(verify.stdout)["verified_artifacts"],6)

    def test_cli_live_mode_refused(self):
        run=subprocess.run([sys.executable,str(PROGRAM),"run","--mode","live"],
                           capture_output=True,text=True)
        self.assertNotEqual(run.returncode,0)
        self.assertFalse(run.stdout)

    def test_cli_invalid_source_does_not_create_output(self):
        with tempfile.TemporaryDirectory() as d:
            dst=Path(d)/"out"
            run=subprocess.run([sys.executable,str(PROGRAM),"run","--day28-report",str(Path(d)/"missing.json"),
                                "--day28-manifest",str(MANIFEST),"--output-dir",str(dst)],
                               capture_output=True,text=True)
            self.assertNotEqual(run.returncode,0)
            self.assertFalse(dst.exists())

    def test_ast_no_network_or_sdk_imports(self):
        tree=ast.parse(PROGRAM.read_text(encoding="utf-8"))
        roots=set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import): roots.update(x.name.split(".")[0] for x in node.names)
            if isinstance(node, ast.ImportFrom) and node.module: roots.add(node.module.split(".")[0])
        self.assertFalse(roots.intersection({"alphagenome","requests","socket","http","urllib","aiohttp","boto3"}))


if __name__ == "__main__": unittest.main()
