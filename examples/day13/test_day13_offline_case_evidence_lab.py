"""Python 3.12+ tests for the fully synthetic Day 13 evidence lesson."""
from __future__ import annotations

import ast
import copy
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import day13_offline_case_evidence_lab as lab

ROOT = Path(__file__).resolve().parent
FIXTURE = ROOT / "day13_synthetic_evidence_case.json"
SCRIPT = ROOT / "day13_offline_case_evidence_lab.py"


class ContractTests(unittest.TestCase):
    def setUp(self):
        self.record, self.hash = lab.read_fixture(FIXTURE)

    def test_load_synthetic(self):
        self.assertEqual(self.record["data_kind"], lab.KIND)

    def test_fixture_sha_is_hex(self):
        self.assertEqual(len(self.hash), 64)
        self.assertTrue(all(ch in "0123456789abcdef" for ch in self.hash))

    def test_six_distinct_stages(self):
        self.assertEqual([r["step_id"] for r in self.record["evidence_rows"]], list(lab.STAGES))

    def test_three_inconclusive(self):
        self.assertEqual(lab.packet(self.record, self.hash)["status_counts"]["INCONCLUSIVE"], 3)

    def test_no_clinical_output(self):
        p = lab.packet(self.record, self.hash)
        self.assertIsNone(p["clinical_classification"])
        self.assertFalse(p["experiment_performed"])
        self.assertFalse(p["atlas_query_performed"])

    def test_absent_model_score_not_zero(self):
        self.assertIsNone(lab.packet(self.record, self.hash)["model_score"])

    def test_null_blood_record(self):
        self.assertIsNone(self.record["evidence_rows"][3]["observed_value"])

    def test_outcomes_are_to_review_not_rank(self):
        p = lab.packet(self.record, self.hash)
        self.assertEqual(p["conclusion"], "RESEARCH_HYPOTHESIS_ONLY")
        self.assertEqual(p["review_status"], "REVIEW_REQUIRED")
        self.assertNotIn("rank", p)

    def test_wgs_is_explicitly_toy(self):
        self.assertEqual(self.record["evidence_rows"][2]["observed_value"], "TOY_INTRONIC_CANDIDATE")
        self.assertEqual(self.record["assembly"], "TOY_NOT_GRCH38")

    def test_output_names_unique(self):
        self.assertEqual(len(set(lab.OUTPUTS)), 5)

    def test_markdown_describes_inconclusive(self):
        self.assertIn("INCONCLUSIVE", lab.questions(lab.packet(self.record, self.hash)).decode())

    def test_csv_all_rows(self):
        table = lab.table_bytes(lab.packet(self.record, self.hash)).decode()
        rows = list(csv.DictReader(table.splitlines()))
        self.assertEqual(len(rows), 6)
        self.assertEqual(rows[4]["observed_value"], "")

    def test_svg_dimensions_and_disclaimers(self):
        svg = lab.svg_bytes(lab.packet(self.record, self.hash)).decode()
        self.assertIn('width="1920" height="1080"', svg)
        self.assertIn('#E5DEDA', svg)
        self.assertIn('API calls: 0', svg)

    def test_html_escape_malicious_note(self):
        r = copy.deepcopy(self.record)
        r["evidence_rows"][0]["review_question"] = "Could <script>alert(1)</script> change this?"
        lab.validate(r)
        report = lab.packet(r, self.hash)
        text = lab.html_bytes(report).decode()
        self.assertNotIn("<script>", text)
        self.assertIn("&lt;script&gt;", text)

    def test_only_standard_library_imports(self):
        tree = ast.parse(SCRIPT.read_text(encoding="utf-8"))
        imports = [n.module.split(".")[0] for n in tree.body if isinstance(n, ast.ImportFrom) and n.module]
        imports += [name.name.split(".")[0] for n in tree.body if isinstance(n, ast.Import) for name in n.names]
        self.assertTrue(set(imports).issubset({"__future__", "argparse", "csv", "hashlib", "html", "io", "json", "re", "pathlib", "typing"}))


# Each row is a distinct invalid input contract that should be rejected.
def mutator(name):
    def test(self):
        r, _ = lab.read_fixture(FIXTURE)
        r = copy.deepcopy(r)
        if name == "not_dict": r = []
        elif name == "extra_patient_id": r["patient_id"] = "private"
        elif name == "no_evidence": del r["evidence_rows"]
        elif name == "wrong_schema": r["schema_version"] = "OTHER"
        elif name == "wrong_kind": r["data_kind"] = "PATIENT"
        elif name == "wrong_scenario": r["scenario_id"] = "ANY"
        elif name == "wrong_gene": r["gene_symbol"] = "BRCA1"
        elif name == "wrong_build": r["assembly"] = "GRCh38"
        elif name == "wrong_variant": r["variant_id"] = "chr9:128225994:G>A"
        elif name == "wrong_phenotype": r["phenotype_label"] = "clinical seizure"
        elif name == "clinical_true": r["clinical_interpretation"] = True
        elif name == "api_one": r["real_api_calls"] = 1
        elif name == "api_bool": r["real_api_calls"] = False
        elif name == "experiment_true": r["experiment_performed"] = True
        elif name == "missing_row": r["evidence_rows"].pop()
        elif name == "extra_row": r["evidence_rows"].append(copy.deepcopy(r["evidence_rows"][0]))
        elif name == "row_reorder": r["evidence_rows"].reverse()
        elif name == "row_not_dict": r["evidence_rows"][0] = "PANEL"
        elif name == "row_extra": r["evidence_rows"][0]["patient"] = "x"
        elif name == "wrong_step": r["evidence_rows"][0]["step_id"] = "WGS"
        elif name == "wrong_type": r["evidence_rows"][0]["evidence_type"] = "SNP_ASSOCIATION"
        elif name == "status_overclaim": r["evidence_rows"][4]["status"] = "CONFIRMED"
        elif name == "status_bad_wes": r["evidence_rows"][1]["status"] = "NEGATIVE"
        elif name == "source_newline": r["evidence_rows"][0]["source_id"] = "TOY\nBREAK"
        elif name == "source_empty": r["evidence_rows"][0]["source_id"] = ""
        elif name == "source_lower": r["evidence_rows"][0]["source_id"] = "patient_file"
        elif name == "question_empty": r["evidence_rows"][0]["review_question"] = "?"
        elif name == "question_newline": r["evidence_rows"][0]["review_question"] = "Does it?\nOr not?"
        elif name == "question_not_string": r["evidence_rows"][0]["review_question"] = 100
        elif name == "wgs_not_toy": r["evidence_rows"][2]["observed_value"] = "real patient"
        elif name == "model_measured": r["evidence_rows"][4]["observed_value"] = 0.8
        elif name == "blood_measured": r["evidence_rows"][3]["observed_value"] = 0
        elif name == "validation_measured": r["evidence_rows"][5]["observed_value"] = 1.0
        else: self.fail("Unlisted test mutation")
        with self.assertRaises(ValueError):
            lab.validate(r)
    return test


class ValidationRejections(unittest.TestCase):
    pass


for kind in (
    "not_dict", "extra_patient_id", "no_evidence", "wrong_schema", "wrong_kind", "wrong_scenario",
    "wrong_gene", "wrong_build", "wrong_variant", "wrong_phenotype", "clinical_true", "api_one",
    "api_bool", "experiment_true", "missing_row", "extra_row", "row_reorder", "row_not_dict",
    "row_extra", "wrong_step", "wrong_type", "status_overclaim", "status_bad_wes", "source_newline",
    "source_empty", "source_lower", "question_empty", "question_newline", "question_not_string",
    "wgs_not_toy", "model_measured", "blood_measured", "validation_measured",
):
    setattr(ValidationRejections, "test_reject_" + kind, mutator(kind))


class ArtifactsTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.out = Path(self.tmp.name) / "generated"
        self.output = lab.run(str(FIXTURE), str(self.out))

    def test_run_reports_pass(self):
        self.assertEqual(self.output["status"], "PASS")

    def test_run_reports_real_api_zero(self):
        self.assertEqual(self.output["real_api_calls_this_run"], 0)

    def test_five_outputs_plus_manifest(self):
        self.assertEqual(len(list(self.out.iterdir())), 6)

    def test_packet_is_valid_json(self):
        r = json.loads((self.out / lab.OUTPUTS[0]).read_text())
        self.assertEqual(r["step_count"], 6)

    def test_verify_matches(self):
        self.assertEqual(lab.verify(str(FIXTURE), str(self.out))["checksums"], "MATCH")

    def test_repeat_refused(self):
        with self.assertRaises(ValueError): lab.run(str(FIXTURE), str(self.out))

    def test_modified_json_rejected(self):
        p = self.out / lab.OUTPUTS[0]
        p.write_bytes(p.read_bytes() + b" ")
        with self.assertRaises(ValueError): lab.verify(str(FIXTURE), str(self.out))

    def test_modified_svg_rejected(self):
        p = self.out / lab.OUTPUTS[3]
        p.write_bytes(p.read_bytes().replace(b"#E5DEDA", b"#FFFFFF"))
        with self.assertRaises(ValueError): lab.verify(str(FIXTURE), str(self.out))

    def test_missing_csv_rejected(self):
        (self.out / lab.OUTPUTS[1]).unlink()
        with self.assertRaises(ValueError): lab.verify(str(FIXTURE), str(self.out))

    def test_extra_file_rejected(self):
        (self.out / "unexpected.txt").write_text("x")
        with self.assertRaises(ValueError): lab.verify(str(FIXTURE), str(self.out))

    def test_manifest_wrong_sha_rejected(self):
        f = self.out / "day13_manifest.json"
        m = json.loads(f.read_text())
        m["fixture_sha256"] = "0"*64
        f.write_text(json.dumps(m))
        with self.assertRaises(ValueError): lab.verify(str(FIXTURE), str(self.out))

    def test_manifest_bad_json_rejected(self):
        (self.out / "day13_manifest.json").write_text("{")
        with self.assertRaises(ValueError): lab.verify(str(FIXTURE), str(self.out))

    def test_manifest_extra_key_rejected(self):
        f = self.out / "day13_manifest.json"
        m = json.loads(f.read_text()); m["extra"] = "bad"; f.write_text(json.dumps(m))
        with self.assertRaises(ValueError): lab.verify(str(FIXTURE), str(self.out))

    def test_source_mutation_rejected(self):
        other = Path(self.tmp.name) / "case.json"
        other.write_bytes(FIXTURE.read_bytes().replace(b"TOY_PANEL_NOTE", b"TOY_NEW_NOTE_1"))
        with self.assertRaises(ValueError): lab.verify(str(other), str(self.out))

    def test_manifest_absent_rejected(self):
        (self.out / "day13_manifest.json").unlink()
        with self.assertRaises(ValueError): lab.verify(str(FIXTURE), str(self.out))

    def test_no_patient_identifier_emitted(self):
        combined = "\n".join(p.read_text(encoding="utf-8") for p in self.out.iterdir())
        self.assertNotIn("patient_id", combined)
        self.assertNotIn("full_name", combined)

    def test_preview_data_marked_toy(self):
        webpage = (self.out / lab.OUTPUTS[4]).read_text()
        self.assertIn("TOY_NOT_GRCH38", webpage)
        self.assertIn("SYNTHETIC", webpage)

    def test_cli_run_stdout(self):
        out2 = Path(self.tmp.name) / "via_cli"
        r = subprocess.run([sys.executable, str(SCRIPT), "run", "--mode", "offline", "--fixture", str(FIXTURE), "--output-dir", str(out2)],
                           check=True, capture_output=True, text=True)
        self.assertEqual(json.loads(r.stdout)["step_count"], 6)

    def test_cli_live_refused_with_api_key(self):
        env = {**os.environ, "ALPHAGENOME_API_KEY": "SYNTHETIC_NOT_A_KEY"}
        r = subprocess.run([sys.executable, str(SCRIPT), "run", "--mode", "live", "--fixture", str(FIXTURE),
                            "--output-dir", str(Path(self.tmp.name)/"no_live")], env=env, capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0)
        self.assertFalse((Path(self.tmp.name) / "no_live").exists())
        self.assertIn("LIVE is not implemented", r.stderr)

    def test_cli_verify(self):
        r = subprocess.run([sys.executable, str(SCRIPT), "verify", "--fixture", str(FIXTURE), "--output-dir", str(self.out)],
                           capture_output=True, text=True, check=True)
        self.assertEqual(json.loads(r.stdout)["checksums"], "MATCH")

    def test_nonexistent_parent_rejected(self):
        with self.assertRaises(ValueError): lab.run(str(FIXTURE), str(Path(self.tmp.name)/"missing"/"new"))

    def test_symlink_fixture_rejected(self):
        f = Path(self.tmp.name)/"fixture_symlink.json"
        f.symlink_to(FIXTURE)
        with self.assertRaises(ValueError): lab.read_fixture(f)

    def test_non_json_fixture_rejected(self):
        f=Path(self.tmp.name)/"bad.json"; f.write_text("not json")
        with self.assertRaises(ValueError): lab.read_fixture(f)

    def test_missing_fixture_rejected(self):
        with self.assertRaises(ValueError): lab.read_fixture(Path(self.tmp.name)/"missing.json")

    def test_large_fixture_rejected(self):
        f=Path(self.tmp.name)/"large.json"; f.write_bytes(b" "*50001)
        with self.assertRaises(ValueError): lab.read_fixture(f)

    def test_repeat_bitwise_identical(self):
        other=Path(self.tmp.name)/"independent"
        lab.run(str(FIXTURE),str(other))
        for name in (*lab.OUTPUTS,"day13_manifest.json"):
            self.assertEqual((self.out/name).read_bytes(),(other/name).read_bytes())


if __name__ == "__main__":
    unittest.main(verbosity=2)
