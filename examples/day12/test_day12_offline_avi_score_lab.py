"""Unit and CLI tests for Day 12, no network or live API access."""
from __future__ import annotations

import ast
import copy
import csv
import hashlib
import json
import math
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import day12_offline_avi_score_lab as lesson

ROOT = Path(__file__).resolve().parent
TOY_FILE = ROOT / "day12_synthetic_avi_fixture.json"
COMPARISON_FILE = ROOT / "day12_synthetic_avi_comparison.json"


def sample_inputs():
    return copy.deepcopy(lesson.load_json(TOY_FILE)), copy.deepcopy(lesson.load_json(COMPARISON_FILE))


class RankMathTests(unittest.TestCase):
    def test_phred_10(self):
        self.assertAlmostEqual(lesson.phred_from_cdf(.9), 10)

    def test_phred_20(self):
        self.assertAlmostEqual(lesson.phred_from_cdf(.99), 20)

    def test_phred_30(self):
        self.assertAlmostEqual(lesson.phred_from_cdf(.999), 30)

    def test_top_percent_10(self):
        self.assertAlmostEqual(lesson.top_percent_from_cdf(.9), 10)

    def test_top_percent_20(self):
        self.assertAlmostEqual(lesson.top_percent_from_cdf(.99), 1)

    def test_top_percent_30(self):
        self.assertAlmostEqual(lesson.top_percent_from_cdf(.999), .1)

    def test_phred_zero_cdf(self):
        self.assertEqual(lesson.phred_from_cdf(0), 0)

    def test_top_percent_zero_cdf(self):
        self.assertEqual(lesson.top_percent_from_cdf(0), 100)

    def test_phred_cdf_one_rejected(self):
        with self.assertRaises(ValueError):lesson.phred_from_cdf(1)

    def test_negative_cdf_rejected(self):
        with self.assertRaises(ValueError):lesson.phred_from_cdf(-.1)

    def test_cdf_above_one_rejected(self):
        with self.assertRaises(ValueError):lesson.phred_from_cdf(1.2)

    def test_nan_rejected(self):
        with self.assertRaises(ValueError):lesson.phred_from_cdf(float("nan"))

    def test_infinity_rejected(self):
        with self.assertRaises(ValueError):lesson.phred_from_cdf(float("inf"))

    def test_true_not_a_float(self):
        with self.assertRaises(ValueError):lesson.phred_from_cdf(True)

    def test_non_numeric_rejected(self):
        with self.assertRaises(ValueError):lesson.phred_from_cdf("0.99")

    def test_higher_cdf_higher_phred(self):
        self.assertLess(lesson.phred_from_cdf(.95), lesson.phred_from_cdf(.995))

    def test_higher_cdf_lower_tail(self):
        self.assertGreater(lesson.top_percent_from_cdf(.95), lesson.top_percent_from_cdf(.995))

    def test_phred_and_top_fraction_inverse(self):
        cdf = .993690426555198
        p = lesson.phred_from_cdf(cdf)
        self.assertAlmostEqual(100 * 10 ** (-p/10), lesson.top_percent_from_cdf(cdf))


class SyntheticInputTests(unittest.TestCase):
    def setUp(self):
        self.toy, self.hist = sample_inputs()

    def invalid(self, key, value):
        self.toy[key] = value
        with self.assertRaises(ValueError):lesson.validate_synthetic(self.toy)

    def test_valid_fixture_has_four(self):
        self.assertEqual(len(lesson.validate_synthetic(self.toy)), 4)

    def test_not_synthetic_rejected(self):
        self.invalid("evidence_class", "HISTORICAL_LIVE_OUTPUT")

    def test_real_assembly_rejected(self):
        self.invalid("assembly", "GRCh38")

    def test_wrong_semantics_rejected(self):
        self.invalid("score_semantics", "CLINICAL_PROBABILITY")

    def test_clinical_true_rejected(self):
        self.invalid("clinical_interpretation", True)

    def test_candidates_not_list(self):
        self.invalid("candidates", "foo")

    def test_no_candidates_rejected(self):
        self.invalid("candidates", [])

    def test_many_candidates_rejected(self):
        self.invalid("candidates", self.toy["candidates"] * 9)

    def test_duplicate_candidate_rejected(self):
        self.toy["candidates"][1]["variant_id"] = "TOY_A"
        with self.assertRaises(ValueError):lesson.validate_synthetic(self.toy)

    def test_real_variant_name_rejected(self):
        self.toy["candidates"][0]["variant_id"] = "chr9:123:A>G"
        with self.assertRaises(ValueError):lesson.validate_synthetic(self.toy)

    def test_variant_newline_rejected(self):
        self.toy["candidates"][0]["variant_id"] = "TOY_A\n# injected"
        with self.assertRaises(ValueError):lesson.validate_synthetic(self.toy)

    def test_missing_with_attribute_rejected(self):
        self.toy["candidates"][3]["toy_attributions"] = [{"feature":"TEACHING_BAD", "value":1}]
        with self.assertRaises(ValueError):lesson.validate_synthetic(self.toy)

    def test_scored_without_attribute_rejected(self):
        self.toy["candidates"][0]["toy_attributions"] = []
        with self.assertRaises(ValueError):lesson.validate_synthetic(self.toy)

    def test_duplicate_attribute_rejected(self):
        self.toy["candidates"][0]["toy_attributions"][1]["feature"] = "TEACHING_PROTEIN_PROXY"
        with self.assertRaises(ValueError):lesson.validate_synthetic(self.toy)

    def test_attribute_name_must_be_teaching(self):
        self.toy["candidates"][0]["toy_attributions"][0]["feature"] = "MERGED_SPLICING"
        with self.assertRaises(ValueError):lesson.validate_synthetic(self.toy)

    def test_attribute_magnitude_bounded(self):
        self.toy["candidates"][0]["toy_attributions"][0]["value"] = 1.8
        with self.assertRaises(ValueError):lesson.validate_synthetic(self.toy)

    def test_attribute_nan_rejected(self):
        self.toy["candidates"][0]["toy_attributions"][0]["value"] = float("nan")
        with self.assertRaises(ValueError):lesson.validate_synthetic(self.toy)

    def test_attribute_string_rejected(self):
        self.toy["candidates"][0]["toy_attributions"][0]["value"] = "high"
        with self.assertRaises(ValueError):lesson.validate_synthetic(self.toy)

    def test_missing_score_kept_missing(self):
        r=lesson.build_report(self.toy, self.hist, "a"*64, "b"*64)
        self.assertIsNone(r["toy_candidates"][3]["phred"])
        self.assertEqual(r["toy_candidates"][3]["score_status"],"MISSING_SCORE")

    def test_same_phred_different_top_feature(self):
        r=lesson.build_report(self.toy, self.hist,"a"*64,"b"*64)
        a,b=r["toy_candidates"][:2]
        self.assertEqual(a["phred"],b["phred"])
        self.assertNotEqual(a["top_teaching_feature"],b["top_teaching_feature"])

    def test_cannot_convert_missing_to_zero(self):
        r=lesson.build_report(self.toy,self.hist,"a"*64,"b"*64)
        self.assertEqual(r["missing_score_count"],1)
        self.assertNotEqual(r["toy_candidates"][3]["phred"],0)


class SyntheticComparisonTests(unittest.TestCase):
    def setUp(self):
        self.toy,self.hist=sample_inputs()

    def test_comparison_valid(self):
        self.assertIs(lesson.validate_comparison(self.hist),self.hist)

    def test_comparison_phred(self):
        self.assertAlmostEqual(lesson.phred_from_cdf(self.hist["cdf_quantile"]),16.989700043360187)

    def test_comparison_tail_percent(self):
        self.assertAlmostEqual(lesson.top_percent_from_cdf(self.hist["cdf_quantile"]),2.0)

    def test_comparison_residual(self):
        self.assertAlmostEqual(self.hist["avi_raw"]-self.hist["attribution_sum"],0.10)

    def test_native_response_cannot_be_claimed(self):
        self.hist["native_response_retained"]=True
        with self.assertRaises(ValueError):lesson.validate_comparison(self.hist)

    def test_comparison_rejects_sdk_claim(self):
        self.hist["sdk_version"]="UNAUTHORIZED"
        with self.assertRaises(ValueError):lesson.validate_comparison(self.hist)

    def test_comparison_rejects_live_status_claim(self):
        self.hist["comparison_live_status"]="CURRENT_LIVE"
        with self.assertRaises(ValueError):lesson.validate_comparison(self.hist)

    def test_comparison_cdf_mismatch_rejected(self):
        self.hist["cdf_quantile"] = .9
        with self.assertRaises(ValueError):lesson.validate_comparison(self.hist)

    def test_comparison_residual_mismatch_rejected(self):
        self.hist["additivity_residual"] = 0
        with self.assertRaises(ValueError):lesson.validate_comparison(self.hist)

    def test_comparison_review_not_pass(self):
        self.hist["additivity_status"]="PASS"
        with self.assertRaises(ValueError):lesson.validate_comparison(self.hist)

    def test_comparison_wrong_variant(self):
        self.hist["variant_label"]="chr1:123:A>G"
        with self.assertRaises(ValueError):lesson.validate_comparison(self.hist)

    def test_comparison_wrong_assembly_claim(self):
        self.hist["assembly"]="GRCh38"
        with self.assertRaises(ValueError):lesson.validate_comparison(self.hist)

    def test_comparison_missing_required_role_rejected(self):
        self.hist.pop("dataset_role")
        with self.assertRaises(ValueError):lesson.validate_comparison(self.hist)

    def test_comparison_unrecognized_field_rejected(self):
        self.hist["native_api_claim"]="BAD"
        with self.assertRaises(ValueError):lesson.validate_comparison(self.hist)

    def test_comparison_feature_count(self):
        self.hist["feature_count"]=17
        with self.assertRaises(ValueError):lesson.validate_comparison(self.hist)

    def test_comparison_rejects_historical_output_claim(self):
        self.hist['evidence_class']='DERIVED_FROM_HISTORICAL_OUTPUT'
        with self.assertRaises(ValueError):lesson.validate_comparison(self.hist)

    def test_comparison_rejects_model_execution(self):
        self.hist['model_executed']=True
        with self.assertRaises(ValueError):lesson.validate_comparison(self.hist)

    def test_comparison_not_mixed_with_toys(self):
        r=lesson.build_report(self.toy,self.hist,"a"*64,"b"*64)
        self.assertNotIn("TOY_COMPARISON",[x["variant_id"] for x in r["toy_candidates"]])
        self.assertEqual(r["synthetic_comparison_reference"]["additivity_status"],"REVIEW")


class ArtifactIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base=Path(self.tmp.name)
        self.toy=self.base/"toy.json"
        self.hist=self.base/"comparison.json"
        self.toy.write_bytes(TOY_FILE.read_bytes())
        self.hist.write_bytes(COMPARISON_FILE.read_bytes())
        self.folder=self.base/"report"

    def execute(self):
        return lesson.run(self.toy,self.hist,self.folder)

    def test_run_pass(self):
        self.assertEqual(self.execute()["status"],"PASS")

    def test_output_artifacts_five(self):
        r=self.execute()
        self.assertEqual(r["output_artifact_count"],5)
        for name in lesson.ARTIFACTS:self.assertTrue((self.folder/name).is_file())

    def test_verify_match(self):
        self.execute()
        self.assertEqual(lesson.verify(self.toy,self.hist,self.folder)["checksums"],"MATCH")

    def test_manifest_is_not_signature(self):
        self.execute()
        self.assertTrue(lesson.load_json(self.folder/"day12_manifest.json")["manifest_is_not_a_signature"])

    def test_live_blocked(self):
        with self.assertRaises(ValueError):lesson.run(self.toy,self.hist,self.folder,"live")
        self.assertFalse(self.folder.exists())

    def test_existing_directory_not_overwritten(self):
        self.folder.mkdir()
        sentinel=self.folder/"sentinel.txt"
        sentinel.write_text("keep")
        with self.assertRaises(FileExistsError):self.execute()
        self.assertEqual(sentinel.read_text(),"keep")

    def test_toy_source_edit_detected(self):
        self.execute()
        self.toy.write_text(self.toy.read_text()+" ")
        with self.assertRaisesRegex(ValueError,"Source SHA"):
            lesson.verify(self.toy,self.hist,self.folder)

    def test_comparison_source_edit_detected(self):
        self.execute()
        self.hist.write_text(self.hist.read_text()+" ")
        with self.assertRaisesRegex(ValueError,"Source SHA"):
            lesson.verify(self.toy,self.hist,self.folder)

    def test_svg_edit_detected(self):
        self.execute()
        p=self.folder/lesson.ARTIFACTS[2]
        p.write_text(p.read_text()+"modified")
        with self.assertRaisesRegex(ValueError,"Output SHA"):
            lesson.verify(self.toy,self.hist,self.folder)

    def test_csv_edit_detected(self):
        self.execute()
        p=self.folder/lesson.ARTIFACTS[1]
        p.write_text(p.read_text().replace("TOY_A","TOY_X"))
        with self.assertRaisesRegex(ValueError,"Output SHA"):
            lesson.verify(self.toy,self.hist,self.folder)

    def test_json_edit_detected(self):
        self.execute()
        p=self.folder/lesson.ARTIFACTS[0]
        p.write_text(p.read_text().replace('"OFFLINE_ONLY"','"LIVE"'))
        with self.assertRaisesRegex(ValueError,"Output SHA"):
            lesson.verify(self.toy,self.hist,self.folder)

    def test_manifest_sha_tamper_rejected(self):
        self.execute()
        p=self.folder/"day12_manifest.json"
        obj=lesson.load_json(p)
        obj["artifact_sha256"][lesson.ARTIFACTS[0]]="0"*64
        lesson.save_json(p,obj)
        with self.assertRaisesRegex(ValueError,"Output SHA"):
            lesson.verify(self.toy,self.hist,self.folder)

    def test_manifest_gate_tamper_rejected(self):
        self.execute()
        p=self.folder/"day12_manifest.json"
        obj=lesson.load_json(p)
        obj["release_gate"]="APPROVED"
        lesson.save_json(p,obj)
        with self.assertRaisesRegex(ValueError,"Manifest status"):
            lesson.verify(self.toy,self.hist,self.folder)

    def test_csv_has_four_rows(self):
        self.execute()
        with (self.folder/lesson.ARTIFACTS[1]).open(encoding="utf-8",newline="") as f:
            self.assertEqual(len(list(csv.DictReader(f))),4)

    def test_csv_missing_blank_not_zero(self):
        self.execute()
        with (self.folder/lesson.ARTIFACTS[1]).open(encoding="utf-8",newline="") as f:
            rows=list(csv.DictReader(f))
        self.assertEqual(rows[-1]["phred"],"")
        self.assertEqual(rows[-1]["score_status"],"MISSING_SCORE")

    def test_svg_ratio(self):
        self.execute()
        self.assertIn('viewBox="0 0 1920 1080"',(self.folder/lesson.ARTIFACTS[2]).read_text())

    def test_svg_contains_toy_ids(self):
        self.execute()
        s=(self.folder/lesson.ARTIFACTS[2]).read_text()
        self.assertIn("TOY_A",s)
        self.assertIn("TOY_D",s)

    def test_html_contains_synthetic(self):
        self.execute()
        self.assertIn("SYNTHETIC EDUCATIONAL DEMO",(self.folder/lesson.ARTIFACTS[3]).read_text())

    def test_html_separates_history(self):
        self.execute()
        s=(self.folder/lesson.ARTIFACTS[3]).read_text()
        self.assertIn("合成對照範例",s)
        self.assertIn("additivity_status=REVIEW",s)

    def test_notes_prohibit_clinical(self):
        self.execute()
        notes=(self.folder/lesson.ARTIFACTS[4]).read_text()
        self.assertIn("非風險",notes)
        self.assertIn("release_gate=BLOCKED",notes)

    def test_report_no_model_api_calls(self):
        self.execute()
        r=lesson.load_json(self.folder/lesson.ARTIFACTS[0])
        self.assertFalse(r["model_inference_performed"])
        self.assertEqual(r["real_api_calls_this_run"],0)

    def test_report_comparison_independent(self):
        self.execute()
        r=lesson.load_json(self.folder/lesson.ARTIFACTS[0])
        self.assertEqual(r["synthetic_comparison_reference"]["feature_count"],18)
        self.assertEqual(r["toy_candidates"][1]["evidence_class"],"SYNTHETIC_FIXTURE")

    def test_second_run_deterministic(self):
        self.execute()
        other=self.base/"another"
        lesson.run(self.toy,self.hist,other)
        for artifact in (*lesson.ARTIFACTS,"day12_manifest.json"):
            self.assertEqual((self.folder/artifact).read_bytes(),(other/artifact).read_bytes())

    def test_invalid_toy_cannot_create_folder(self):
        self.toy.write_text('{"evidence_class":"REAL_PATIENT"}',encoding="utf-8")
        with self.assertRaises(ValueError):self.execute()
        self.assertFalse(self.folder.exists())

    def test_json_duplicate_key_rejected(self):
        self.toy.write_text('{"x":1,"x":2}',encoding="utf-8")
        with self.assertRaisesRegex(ValueError,"Duplicate JSON"):
            self.execute()

    def test_json_nan_rejected(self):
        self.toy.write_text('{"x":NaN}',encoding="utf-8")
        with self.assertRaisesRegex(ValueError,"Nonfinite JSON"):
            self.execute()

    def test_source_sha_retained(self):
        self.execute()
        report=lesson.load_json(self.folder/lesson.ARTIFACTS[0])
        self.assertEqual(report["input_provenance"]["toy_fixture_sha256"],hashlib.sha256(self.toy.read_bytes()).hexdigest())

    def test_no_external_network_import(self):
        code=ast.parse((ROOT/"day12_offline_avi_score_lab.py").read_text(encoding="utf-8"))
        imports=[n.module or "" for n in ast.walk(code) if isinstance(n,ast.ImportFrom)]
        imports += [alias.name for n in ast.walk(code) if isinstance(n,ast.Import) for alias in n.names]
        self.assertFalse(any(x.split(".")[0] in {"requests","urllib","http","socket","alphagenome","google","subprocess"} for x in imports))

    def test_cli_run_then_verify(self):
        args=[sys.executable,str(ROOT/"day12_offline_avi_score_lab.py")]
        common=["--fixture",str(self.toy),"--comparison-fixture",str(self.hist)]
        run=subprocess.run(args+["run"]+common+["--output-dir",str(self.folder)],capture_output=True,text=True)
        self.assertEqual(run.returncode,0,run.stderr+run.stdout)
        check=subprocess.run(args+["verify"]+common+["--report-dir",str(self.folder)],capture_output=True,text=True)
        self.assertEqual(check.returncode,0,check.stderr+check.stdout)
        self.assertEqual(json.loads(check.stdout)["checksums"],"MATCH")

    def test_cli_live_rejected_even_if_key_set(self):
        env=os.environ.copy()
        env["ALPHAGENOME_API_KEY"]="FAKE_DO_NOT_USE"
        done=subprocess.run([sys.executable,str(ROOT/"day12_offline_avi_score_lab.py"),"run","--mode","live"],capture_output=True,text=True,env=env,cwd=self.base)
        self.assertEqual(done.returncode,2)
        self.assertIn("ERROR",done.stdout)
        self.assertNotIn("FAKE_DO_NOT_USE",done.stdout)


if __name__ == "__main__":
    unittest.main()
