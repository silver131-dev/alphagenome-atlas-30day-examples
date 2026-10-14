"""Unit tests for the synthetic Day 14 splicing validation DESIGN lab."""
import ast
import copy
import csv
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import day14_offline_splicing_validation_lab as lab

HERE = Path(__file__).resolve().parent
FIXTURE = HERE / "day14_synthetic_splicing_design.json"


def fixture_copy():
    return json.loads(FIXTURE.read_text(encoding="utf-8"))


def plan_of(fixture=None):
    return lab.build_plan(fixture or fixture_copy(), "toy.json", "a" * 64)


class Day14FixtureTests(unittest.TestCase):
    def test_fixture_passes(self):
        self.assertEqual(lab.check_fixture(fixture_copy())["gene_symbol"], "DNM1")

    def test_planned_slot_count(self):
        self.assertEqual(plan_of()["planned_slot_count"], 24)

    def test_all_slots_are_not_performed(self):
        self.assertTrue(all(x["experiment_status"] == "NOT_PERFORMED" for x in plan_of()["slots"]))

    def test_all_measurements_are_none(self):
        self.assertTrue(all(x["rt_pcr_result"] is None and x["junction_sequence_result"] is None for x in plan_of()["slots"]))

    def test_no_observed_sample(self):
        self.assertEqual(plan_of()["observed_sample_count"], 0)

    def test_no_model_or_api(self):
        p = plan_of()
        self.assertEqual(p["live_api_calls"], 0)
        self.assertFalse(p["splicing_measured"])
        self.assertFalse(p["protein_function_verified"])
        self.assertFalse(p["clinical_interpretation"])

    def test_valid_constructs_and_roles(self):
        self.assertEqual({x["id"]: x["role"] for x in fixture_copy()["constructs"]}, lab.CONSTRUCTS)

    def test_valid_contexts(self):
        self.assertEqual({x["id"]: x["role"] for x in fixture_copy()["cell_contexts"]}, lab.CONTEXTS)

    def test_two_replicates(self):
        f = fixture_copy(); f["planned_biological_replicates"] = 2
        self.assertEqual(plan_of(f)["planned_slot_count"], 16)

    def test_six_replicates(self):
        f = fixture_copy(); f["planned_biological_replicates"] = 6
        self.assertEqual(plan_of(f)["planned_slot_count"], 48)

    def test_missing_field_rejected(self):
        f = fixture_copy(); del f["gene_symbol"]
        with self.assertRaises(ValueError): lab.check_fixture(f)

    def test_extra_patient_field_rejected(self):
        f = fixture_copy(); f["patient_id"] = "should-not-appear"
        with self.assertRaises(ValueError): lab.check_fixture(f)

    def test_wrong_data_kind_rejected(self):
        f = fixture_copy(); f["data_kind"] = "PATIENT_SAMPLE"
        with self.assertRaises(ValueError): lab.check_fixture(f)

    def test_wrong_scenario_rejected(self):
        f = fixture_copy(); f["scenario_id"] = "OTHER"
        with self.assertRaises(ValueError): lab.check_fixture(f)

    def test_wrong_gene_rejected(self):
        f = fixture_copy(); f["gene_symbol"] = "OTHER"
        with self.assertRaises(ValueError): lab.check_fixture(f)

    def test_wrong_reference_rejected(self):
        f = fixture_copy(); f["reference_allele"] = "C"
        with self.assertRaises(ValueError): lab.check_fixture(f)

    def test_wrong_alternative_rejected(self):
        f = fixture_copy(); f["alternative_allele"] = "T"
        with self.assertRaises(ValueError): lab.check_fixture(f)

    def test_real_coordinates_rejected(self):
        f = fixture_copy(); f["genomic_context"] = "GRCh38"
        with self.assertRaises(ValueError): lab.check_fixture(f)

    def test_wrong_assay_rejected(self):
        f = fixture_copy(); f["planned_assays"] = ["CLINICAL_DIAGNOSIS"]
        with self.assertRaises(ValueError): lab.check_fixture(f)

    def test_construct_missing_rejected(self):
        f = fixture_copy(); f["constructs"] = f["constructs"][1:]
        with self.assertRaises(ValueError): lab.check_fixture(f)

    def test_duplicate_construct_rejected(self):
        f = fixture_copy(); f["constructs"][-1] = copy.deepcopy(f["constructs"][0])
        with self.assertRaises(ValueError): lab.check_fixture(f)

    def test_wrong_construct_role_rejected(self):
        f = fixture_copy(); f["constructs"][0]["role"] = "UNKNOWN"
        with self.assertRaises(ValueError): lab.check_fixture(f)

    def test_extra_construct_field_rejected(self):
        f = fixture_copy(); f["constructs"][0]["secret"] = "bad"
        with self.assertRaises(ValueError): lab.check_fixture(f)

    def test_context_missing_rejected(self):
        f = fixture_copy(); f["cell_contexts"] = f["cell_contexts"][:1]
        with self.assertRaises(ValueError): lab.check_fixture(f)

    def test_duplicate_context_rejected(self):
        f = fixture_copy(); f["cell_contexts"][-1] = copy.deepcopy(f["cell_contexts"][0])
        with self.assertRaises(ValueError): lab.check_fixture(f)

    def test_wrong_context_role_rejected(self):
        f = fixture_copy(); f["cell_contexts"][0]["role"] = "PATIENT"
        with self.assertRaises(ValueError): lab.check_fixture(f)

    def test_bool_replicate_rejected(self):
        f = fixture_copy(); f["planned_biological_replicates"] = True
        with self.assertRaises(ValueError): lab.check_fixture(f)

    def test_low_replicate_rejected(self):
        f = fixture_copy(); f["planned_biological_replicates"] = 1
        with self.assertRaises(ValueError): lab.check_fixture(f)

    def test_large_replicate_rejected(self):
        f = fixture_copy(); f["planned_biological_replicates"] = 1000000
        with self.assertRaises(ValueError): lab.check_fixture(f)

    def test_nonint_replicate_rejected(self):
        f = fixture_copy(); f["planned_biological_replicates"] = "3"
        with self.assertRaises(ValueError): lab.check_fixture(f)

    def test_39_nt_conditional_codons(self):
        p = plan_of()
        self.assertTrue(p["in_frame_if_in_coding_region"])
        self.assertEqual(p["potential_added_codons_if_translated"], 13)
        self.assertFalse(p["protein_function_verified"])

    def test_40_nt_not_in_frame(self):
        f = fixture_copy(); f["putative_extension_nt"] = 40
        p = plan_of(f)
        self.assertFalse(p["in_frame_if_in_coding_region"])
        self.assertIsNone(p["potential_added_codons_if_translated"])

    def test_zero_nt_rejected(self):
        f = fixture_copy(); f["putative_extension_nt"] = 0
        with self.assertRaises(ValueError): lab.check_fixture(f)

    def test_bool_nt_rejected(self):
        f = fixture_copy(); f["putative_extension_nt"] = True
        with self.assertRaises(ValueError): lab.check_fixture(f)

    def test_overly_large_nt_rejected(self):
        f = fixture_copy(); f["putative_extension_nt"] = 10000
        with self.assertRaises(ValueError): lab.check_fixture(f)

    def test_list_not_mapping_rejected(self):
        with self.assertRaises(ValueError): lab.check_fixture([])


class Day14HypothesisTests(unittest.TestCase):
    def test_alt_relevant_model_is_hypothesis(self):
        self.assertEqual(lab.proposed_pattern("ALT", "NEURONAL_LIKE"), "HYPOTHESIS_ALT_MORE_CRYPTIC_JUNCTION_THAN_WT")

    def test_alt_comparator_unknown(self):
        self.assertEqual(lab.proposed_pattern("ALT", "COMPARATOR"), "CONTEXT_EFFECT_UNCERTAIN")

    def test_wt_baseline(self):
        self.assertEqual(lab.proposed_pattern("WT", "NEURONAL_LIKE"), "BASELINE_TO_MEASURE")

    def test_revert_rescue_is_not_asserted(self):
        self.assertEqual(lab.proposed_pattern("REVERT", "NEURONAL_LIKE"), "RESTORATION_TO_TEST")

    def test_empty_background(self):
        self.assertEqual(lab.proposed_pattern("EMPTY", "COMPARATOR"), "ASSAY_BACKGROUND_ONLY")

    def test_no_duplicate_planned_ids(self):
        xs = [r["slot_id"] for r in plan_of()["slots"]]
        self.assertEqual(len(xs), len(set(xs)))

    def test_two_contexts_each_has_twelve_slots(self):
        xs = plan_of()["slots"]
        self.assertEqual(sum(r["cell_context"] == "NEURONAL_LIKE" for r in xs), 12)
        self.assertEqual(sum(r["cell_context"] == "COMPARATOR" for r in xs), 12)

    def test_control_and_test_replication(self):
        xs = plan_of()["slots"]
        self.assertEqual(sum(r["construct"] == "WT" for r in xs), 6)
        self.assertEqual(sum(r["construct"] == "ALT" for r in xs), 6)

    def test_release_blocked(self):
        p = plan_of()
        self.assertIs(p["data_release_allowed"], False)
        self.assertIs(p["external_transfer_allowed"], False)
        self.assertEqual(p["governance_review_status"], "MANUAL_REVIEW_REQUIRED")

    def test_rejection_criteria_are_explicit(self):
        self.assertGreaterEqual(len(plan_of()["hypothesis_rejection_examples"]), 3)


class Day14RenderingTests(unittest.TestCase):
    def test_expected_artifact_count(self):
        self.assertEqual(len(lab.render_artifacts(plan_of())), 5)

    def test_csv_contains_24_rows(self):
        rows = list(csv.DictReader(io.StringIO(lab.csv_text(plan_of()))))
        self.assertEqual(len(rows), 24)

    def test_csv_missing_values_blank(self):
        rows = list(csv.DictReader(io.StringIO(lab.csv_text(plan_of()))))
        self.assertTrue(all(r["rt_pcr_result"] == "" for r in rows))

    def test_csv_not_performed(self):
        rows = list(csv.DictReader(io.StringIO(lab.csv_text(plan_of()))))
        self.assertTrue(all(r["experiment_status"] == "NOT_PERFORMED" for r in rows))

    def test_svg_dimensions(self):
        svg = lab.svg_text(plan_of())
        self.assertIn('width="1920" height="1080"', svg)

    def test_svg_requested_color(self):
        self.assertIn("#E5DEDA", lab.svg_text(plan_of()))

    def test_html_has_disclaimer(self):
        self.assertIn("SYNTHETIC EDUCATIONAL DESIGN", lab.html_text(plan_of(),lab.svg_text(plan_of())))

    def test_html_24_plan_rows(self):
        html = lab.html_text(plan_of(),lab.svg_text(plan_of()))
        self.assertEqual(html.count("<tr>"), 25)  # one header + 24 rows

    def test_markdown_says_not_measured(self):
        txt = lab.summary_markdown(plan_of())
        self.assertIn("NOT MEASURED", txt)

    def test_report_json_preserves_none(self):
        data = json.loads(lab.json_bytes(plan_of()))
        self.assertIsNone(data["slots"][0]["rt_pcr_result"])


class Day14IOTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.fixture = self.base / "day14_synthetic_splicing_design.json"
        self.fixture.write_bytes(FIXTURE.read_bytes())
        self.out = self.base / "out"

    def test_run_creates_files(self):
        summary = lab.run(self.fixture,self.out)
        self.assertEqual(summary["artifact_count"], 5)
        self.assertEqual(len(list(self.out.iterdir())), 6)  # manifest + 5 artifacts

    def test_verify_passes(self):
        lab.run(self.fixture,self.out)
        self.assertEqual(lab.verify(self.fixture,self.out)["checksums"],"MATCH")

    def test_refuse_overwrite(self):
        lab.run(self.fixture,self.out)
        with self.assertRaises(FileExistsError):lab.run(self.fixture,self.out)

    def test_manifest_source_hash(self):
        lab.run(self.fixture,self.out)
        manifest = json.loads((self.out / "day14_artifact_manifest.json").read_text())
        self.assertEqual(manifest["source_fixture"]["sha256"],lab.digest(self.fixture.read_bytes()))

    def test_modified_artifact_rejected(self):
        lab.run(self.fixture,self.out)
        (self.out / "day14_hypothesis_review.md").write_text("tampered",encoding="utf-8")
        with self.assertRaises(ValueError):lab.verify(self.fixture,self.out)

    def test_modified_manifest_rejected(self):
        lab.run(self.fixture,self.out)
        m=self.out / "day14_artifact_manifest.json"
        doc=json.loads(m.read_text()); doc["signature_present"]=True
        m.write_text(json.dumps(doc),encoding="utf-8")
        with self.assertRaises(ValueError):lab.verify(self.fixture,self.out)

    def test_modified_source_rejected(self):
        lab.run(self.fixture,self.out)
        fixture=fixture_copy(); fixture["putative_extension_nt"] = 40
        self.fixture.write_text(json.dumps(fixture),encoding="utf-8")
        with self.assertRaises(ValueError):lab.verify(self.fixture,self.out)

    def test_no_network_modules_required(self):
        tree=ast.parse((HERE / "day14_offline_splicing_validation_lab.py").read_text(encoding="utf-8"))
        imported = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.extend(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported.append(node.module.split(".")[0])
        self.assertFalse(set(imported) & {"requests", "socket", "urllib", "httpx", "alphagenome"})

    def test_default_command_in_subprocess(self):
        run=subprocess.run([sys.executable,str(HERE / 'day14_offline_splicing_validation_lab.py'),
                            'run','--mode','offline','--fixture',str(self.fixture),'--output-dir',str(self.out)],
                           capture_output=True,text=True,check=True)
        self.assertEqual(json.loads(run.stdout)["planned_slot_count"],24)

    def test_live_rejected(self):
        run=subprocess.run([sys.executable,str(HERE / 'day14_offline_splicing_validation_lab.py'),
                            'run','--mode','live','--fixture',str(self.fixture),'--output-dir',str(self.out)],
                           capture_output=True,text=True)
        self.assertNotEqual(run.returncode,0)
        self.assertFalse(self.out.exists())

    def test_key_presence_cannot_enable_live(self):
        with patch.dict(os.environ,{"ALPHAGENOME_API_KEY":"DUMMY_NOT_A_SECRET"}):
            run=subprocess.run([sys.executable,str(HERE / 'day14_offline_splicing_validation_lab.py'),
                                'run','--mode','live','--fixture',str(self.fixture),'--output-dir',str(self.out)],
                               capture_output=True,text=True)
        self.assertNotEqual(run.returncode,0)

    def test_verify_subprocess(self):
        lab.run(self.fixture,self.out)
        run=subprocess.run([sys.executable,str(HERE / 'day14_offline_splicing_validation_lab.py'),
                            'verify','--fixture',str(self.fixture),'--output-dir',str(self.out)],
                           capture_output=True,text=True,check=True)
        self.assertEqual(json.loads(run.stdout)["verified_artifact_count"],5)

    def test_input_file_not_modified(self):
        prior = self.fixture.read_bytes()
        lab.run(self.fixture,self.out)
        self.assertEqual(self.fixture.read_bytes(),prior)

    def test_file_hash_reproducible(self):
        lab.run(self.fixture,self.out)
        out2=self.base/'second'
        lab.run(self.fixture,out2)
        for fn in lab.OUTPUT_NAMES:
            self.assertEqual((self.out/fn).read_bytes(),(out2/fn).read_bytes())

    def test_run_no_silent_clinical_claim(self):
        lab.run(self.fixture,self.out)
        doc=json.loads((self.out/'day14_splicing_design.json').read_text())
        self.assertFalse(doc["clinical_interpretation"])
        self.assertFalse(doc["experiment_performed"])


if __name__ == "__main__":
    unittest.main()
