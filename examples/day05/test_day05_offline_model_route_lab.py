"""Day 05 synthetic routing lab tests. No network/model/real patient data."""
from __future__ import annotations

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
from pathlib import Path
from unittest.mock import patch

import day05_offline_model_route_lab as lab

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "day05_synthetic_model_questions.json"


class FixtureMixin:
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.fixture = copy.deepcopy(json.loads(SOURCE.read_text(encoding="utf-8")))

    def save(self, data=None, name="fixture.json"):
        path = self.root / name
        path.write_text(json.dumps(self.fixture if data is None else data, ensure_ascii=False), encoding="utf-8")
        return path

    def invalid(self, data):
        with self.assertRaises(ValueError):
            lab.read_fixture(self.save(data))


class Day05FixtureTests(FixtureMixin, unittest.TestCase):
    def test_actual_synthetic_fixture_valid(self):
        self.assertEqual(len(lab.read_fixture(SOURCE)["cases"]), 7)

    def test_input_has_no_patient_fields(self):
        self.assertNotIn("patient_id", SOURCE.read_text(encoding="utf-8"))

    def test_refuse_missing_fixture(self):
        with self.assertRaises(ValueError):
            lab.read_fixture(self.root / "missing.json")

    def test_refuse_symlink_input(self):
        alias = self.root / "alias.json"
        alias.symlink_to(SOURCE)
        with self.assertRaises(ValueError):
            lab.read_fixture(alias)

    def test_refuse_broken_json(self):
        path = self.root / "broken.json"
        path.write_text("{bad", encoding="utf-8")
        with self.assertRaises(json.JSONDecodeError):
            lab.read_fixture(path)

    def test_refuse_fixture_larger_than_64k(self):
        p = self.root / "huge.json"
        p.write_text("x" * 65537, encoding="utf-8")
        with self.assertRaises(ValueError):
            lab.read_fixture(p)

    def test_refuse_top_level_array(self): self.invalid([])
    def test_refuse_missing_data_kind(self):
        d=copy.deepcopy(self.fixture); del d["data_kind"]; self.invalid(d)
    def test_refuse_unknown_top_field(self):
        d=copy.deepcopy(self.fixture); d["private_note"]="hello"; self.invalid(d)
    def test_refuse_wrong_data_kind(self):
        d=copy.deepcopy(self.fixture); d["data_kind"]="LIVE"; self.invalid(d)
    def test_refuse_wrong_assembly(self):
        d=copy.deepcopy(self.fixture); d["assembly"]="GRCh38"; self.invalid(d)
    def test_refuse_wrong_schema_version(self):
        d=copy.deepcopy(self.fixture); d["schema_version"]="2.0"; self.invalid(d)
    def test_refuse_missing_cases(self):
        d=copy.deepcopy(self.fixture); del d["cases"]; self.invalid(d)
    def test_refuse_cases_not_list(self):
        d=copy.deepcopy(self.fixture); d["cases"]={}; self.invalid(d)
    def test_refuse_no_cases(self):
        d=copy.deepcopy(self.fixture); d["cases"]=[]; self.invalid(d)
    def test_refuse_more_than_20_cases(self):
        d=copy.deepcopy(self.fixture); d["cases"]= [dict(self.fixture['cases'][0], case_id=f'TOY_{i:03d}') for i in range(21)];self.invalid(d)
    def test_refuse_non_dict_case(self):
        d=copy.deepcopy(self.fixture); d["cases"][0]=None; self.invalid(d)
    def test_refuse_extra_case_field(self):
        d=copy.deepcopy(self.fixture); d["cases"][0]["full_name"]="Alice"; self.invalid(d)
    def test_refuse_missing_case_field(self):
        d=copy.deepcopy(self.fixture); del d["cases"][0]["variant"]; self.invalid(d)
    def test_refuse_duplicate_ids(self):
        d=copy.deepcopy(self.fixture); d["cases"][1]["case_id"]=d["cases"][0]["case_id"]; self.invalid(d)
    def test_refuse_lowercase_case_id(self):
        d=copy.deepcopy(self.fixture); d["cases"][0]["case_id"]="toy_one"; self.invalid(d)
    def test_refuse_numeric_case_id(self):
        d=copy.deepcopy(self.fixture); d["cases"][0]["case_id"]=10; self.invalid(d)
    def test_refuse_empty_case_id(self):
        d=copy.deepcopy(self.fixture); d["cases"][0]["case_id"]=""; self.invalid(d)
    def test_refuse_variant_not_str(self):
        d=copy.deepcopy(self.fixture); d["cases"][0]["variant"]=12; self.invalid(d)
    def test_refuse_real_chromosome(self):
        d=copy.deepcopy(self.fixture); d["cases"][0]["variant"]="chr1:101:A>G"; self.invalid(d)
    def test_refuse_zero_coordinate(self):
        d=copy.deepcopy(self.fixture); d["cases"][0]["variant"]="chrToy:0:A>G"; self.invalid(d)
    def test_refuse_ref_same_as_alt(self):
        d=copy.deepcopy(self.fixture); d["cases"][0]["variant"]="chrToy:101:A>A"; self.invalid(d)
    def test_refuse_indel(self):
        d=copy.deepcopy(self.fixture); d["cases"][0]["variant"]="chrToy:101:A>AG"; self.invalid(d)
    def test_refuse_broken_notation(self):
        d=copy.deepcopy(self.fixture); d["cases"][0]["variant"]="chrToy:101:A:G"; self.invalid(d)
    def test_refuse_unrecognized_consequence(self):
        d=copy.deepcopy(self.fixture); d["cases"][0]["consequence"]="PATHOGENIC"; self.invalid(d)
    def test_refuse_consequence_list(self):
        d=copy.deepcopy(self.fixture); d["cases"][0]["consequence"]=["MISSENSE"]; self.invalid(d)
    def test_refuse_unrecognized_question(self):
        d=copy.deepcopy(self.fixture); d["cases"][0]["research_question"]="IS_CAUSAL"; self.invalid(d)
    def test_refuse_question_list(self):
        d=copy.deepcopy(self.fixture); d["cases"][0]["research_question"]=["RNA_SPLICING"]; self.invalid(d)
    def test_accept_synthetic_other_position(self):
        d=copy.deepcopy(self.fixture); d["cases"][0]["variant"]="chrToy:9:T>G";self.assertEqual(lab.read_fixture(self.save(d))["cases"][0]["variant"],"chrToy:9:T>G")


class Day05RoutingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = lab.read_fixture(SOURCE)
        cls.report = lab.build_report(cls.fixture, lab.sha256_file(SOURCE))
        cls.byid = {r["case_id"]: r for r in cls.report["rows"]}

    def test_schema_mode(self): self.assertEqual(self.report["mode"],"OFFLINE_ONLY")
    def test_assembly_toy_only(self): self.assertEqual(self.report["assembly"],"TOY_NOT_GRCH38")
    def test_count_seven(self): self.assertEqual(self.report["case_count"],7)
    def test_structure_count(self): self.assertEqual(self.report["structure_route_count"],2)
    def test_regulation_count(self): self.assertEqual(self.report["regulatory_route_count"],3)
    def test_review_count(self): self.assertEqual(self.report["question_review_count"],2)
    def test_counts_add_up(self):
        self.assertEqual(sum(self.report[k] for k in ("structure_route_count","regulatory_route_count","question_review_count")),7)
    def test_missense_structure(self): self.assertEqual(self.byid["TOY_MISSENSE_SHAPE"]["model_family_to_study"],"ALPHAFOLD_FAMILY")
    def test_missense_complex(self): self.assertEqual(self.byid["TOY_MISSENSE_COMPLEX"]["next_evidence"],"EXPERIMENTAL_COMPLEX_AND_FUNCTION")
    def test_promoter_expression(self): self.assertEqual(self.byid["TOY_PROMOTER_EXPR"]["candidate_route"],"REGULATORY_EVIDENCE_START")
    def test_intron_splice(self): self.assertEqual(self.byid["TOY_INTRON_SPLICE"]["next_evidence"],"TISSUE_RELEVANT_RNA_OR_SPLICE_ASSAY")
    def test_missense_can_still_require_splice_question(self): self.assertEqual(self.byid["TOY_MISSENSE_SPLICE"]["model_family_to_study"],"ALPHAGENOME_MODEL")
    def test_promoter_not_automatically_protein_structure(self): self.assertEqual(self.byid["TOY_PROMOTER_STRUCTURE"]["candidate_route"],"QUESTION_REVIEW_REQUIRED")
    def test_unspecified_is_manual_review(self): self.assertEqual(self.byid["TOY_UNKNOWN"]["candidate_route"],"QUESTION_REVIEW_REQUIRED")
    def test_no_model_called(self): self.assertTrue(all(r["model_called"] is False for r in self.report["rows"]))
    def test_no_predicted_effect(self): self.assertTrue(all(r["predicted_effect"] is None for r in self.report["rows"]))
    def test_no_clinical_risk(self): self.assertTrue(all(r["clinical_risk"] is None for r in self.report["rows"]))
    def test_all_manual_review(self): self.assertTrue(all(r["manual_review"] == "REQUIRED" for r in self.report["rows"]))
    def test_all_planned_only(self): self.assertTrue(all(r["analysis_status"] == "PLANNED_ONLY" for r in self.report["rows"]))
    def test_all_synthetic_provenance(self): self.assertTrue(all(r["evidence_class"] == "SYNTHETIC_FIXTURE" for r in self.report["rows"]))
    def test_source_sha_is_actual(self): self.assertEqual(self.report["input_sha256"],lab.sha256_file(SOURCE))
    def test_mention_alphafold3_complexes(self): self.assertTrue(any("AlphaFold 3" in n for n in self.report["notices"]))
    def test_mention_regulation_not_timecourse(self): self.assertTrue(any("time-course" in n for n in self.report["notices"]))
    def test_no_api_counter(self): self.assertEqual(self.report["real_api_calls_this_run"],0)
    def test_gene_expression_missense_route_is_question_driven(self):
        r=lab.route_case(dict(self.fixture['cases'][0],research_question="GENE_EXPRESSION"))
        self.assertEqual(r["candidate_route"],"REGULATORY_EVIDENCE_START")
    def test_chromatin_question_routes_regulatory(self):
        r=lab.route_case(dict(self.fixture['cases'][0],research_question="CHROMATIN_ACCESSIBILITY"))
        self.assertEqual(r["next_evidence"],"TISSUE_RELEVANT_CHROMATIN_ASSAY")
    def test_unknown_structure_needs_review(self):
        r=lab.route_case(dict(self.fixture['cases'][0],consequence="UNKNOWN"))
        self.assertEqual(r["model_family_to_study"],"NONE_SELECTED")


class Day05ArtifactTests(FixtureMixin, unittest.TestCase):
    def run_fixture(self):
        src=self.save();out=self.root/"out"; summary=lab.run_lab(src,out)
        return src,out,summary

    def test_run_five_artifacts(self):
        _,out,s=self.run_fixture();self.assertEqual(s["artifact_count"],5);self.assertTrue(all((out/n).exists() for n in lab.ARTIFACTS))
    def test_verify_pass(self):
        src,out,_=self.run_fixture();self.assertEqual(lab.verify_lab(src,out)["checksums"],"MATCH")
    def test_manifest_schema(self):
        _,out,_=self.run_fixture();d=json.loads((out/"day05_sha256_manifest.json").read_text());self.assertEqual(set(d['files']),set(lab.ARTIFACTS))
    def test_csv_rows(self):
        _,out,_=self.run_fixture();r=list(csv.DictReader(io.StringIO((out/"day05_routing_table.csv").read_text(encoding="utf-8"))));self.assertEqual(len(r),7)
    def test_csv_structure_data(self):
        _,out,_=self.run_fixture();t=(out/"day05_routing_table.csv").read_text();self.assertIn("ALPHAFOLD_FAMILY",t)
    def test_json_null(self):
        _,out,_=self.run_fixture();data=json.loads((out/"day05_routing_report.json").read_text());self.assertIsNone(data['rows'][0]['predicted_effect'])
    def test_svg_dimensions(self):
        _,out,_=self.run_fixture();svg=(out/"day05_model_map.svg").read_text();self.assertIn('width="1920" height="1080"',svg)
    def test_svg_background(self):
        _,out,_=self.run_fixture();s=(out/"day05_model_map.svg").read_text();self.assertIn('#E5DEDA',s)
    def test_svg_no_external_assets(self):
        _,out,_=self.run_fixture();s=(out/"day05_model_map.svg").read_text();self.assertNotIn('http://www.w3.org/1999/xlink',s)
    def test_svg_valid_xml(self):
        import xml.etree.ElementTree as ET
        _,out,_=self.run_fixture();ET.fromstring((out/"day05_model_map.svg").read_text())
    def test_html_doctype(self):
        _,out,_=self.run_fixture();self.assertTrue((out/"day05_routing_dashboard.html").read_text().startswith("<!doctype html>"))
    def test_html_seven_rows(self):
        _,out,_=self.run_fixture();self.assertEqual((out/"day05_routing_dashboard.html").read_text().count('<tr>'),8)
    def test_html_synthetic_warning(self):
        _,out,_=self.run_fixture();self.assertIn("No AlphaFold or AlphaGenome model was executed",(out/"day05_routing_dashboard.html").read_text())
    def test_html_xss_escaping(self):
        report=lab.build_report(lab.read_fixture(SOURCE),lab.sha256_file(SOURCE));report['rows'][0]['case_id']='<img onerror=alert(1)>'
        h=lab.render_html(report);self.assertNotIn('<img onerror',h);self.assertIn('&lt;img',h)
    def test_markdown_cards(self):
        _,out,_=self.run_fixture();s=(out/"day05_learning_card.md").read_text();self.assertIn('TOY_MISSENSE_SHAPE',s)
    def test_no_overwrite_output(self):
        src,out,_=self.run_fixture();
        with self.assertRaises(FileExistsError): lab.run_lab(src,out)
    def test_no_overwrite_empty_existing_dir(self):
        src=self.save();out=self.root/"existing";out.mkdir()
        with self.assertRaises(FileExistsError):lab.run_lab(src,out)
    def test_refuse_output_symlink(self):
        src=self.save(); out=self.root/"alias";out.symlink_to(self.root,target_is_directory=True)
        with self.assertRaises(FileExistsError):lab.run_lab(src,out)
    def test_reject_tampered_json(self):
        src,out,_=self.run_fixture();(out/lab.ARTIFACTS[0]).write_text('[]')
        with self.assertRaises(ValueError):lab.verify_lab(src,out)
    def test_reject_tampered_csv(self):
        src,out,_=self.run_fixture();(out/lab.ARTIFACTS[1]).write_text('wrong')
        with self.assertRaises(ValueError):lab.verify_lab(src,out)
    def test_reject_tampered_svg(self):
        src,out,_=self.run_fixture();(out/lab.ARTIFACTS[2]).write_text('wrong')
        with self.assertRaises(ValueError):lab.verify_lab(src,out)
    def test_reject_tampered_html(self):
        src,out,_=self.run_fixture();(out/lab.ARTIFACTS[3]).write_text('wrong')
        with self.assertRaises(ValueError):lab.verify_lab(src,out)
    def test_reject_tampered_md(self):
        src,out,_=self.run_fixture();(out/lab.ARTIFACTS[4]).write_text('wrong')
        with self.assertRaises(ValueError):lab.verify_lab(src,out)
    def test_refuse_modified_input(self):
        src,out,_=self.run_fixture();d=json.loads(src.read_text());d['cases'][0]['variant']='chrToy:9:A>T';src.write_text(json.dumps(d))
        with self.assertRaises(ValueError):lab.verify_lab(src,out)
    def test_refuse_unexpected_output(self):
        src,out,_=self.run_fixture();(out/"EXTRA.txt").write_text('untracked')
        with self.assertRaises(ValueError):lab.verify_lab(src,out)
    def test_refuse_removed_artifact(self):
        src,out,_=self.run_fixture();(out/lab.ARTIFACTS[0]).unlink()
        with self.assertRaises(ValueError):lab.verify_lab(src,out)
    def test_refuse_swapped_manifest_hash(self):
        src,out,_=self.run_fixture();f=out/'day05_sha256_manifest.json';m=json.loads(f.read_text());m['files'][lab.ARTIFACTS[0]]='0'*64;f.write_text(json.dumps(m))
        with self.assertRaises(ValueError):lab.verify_lab(src,out)
    def test_refuse_manifest_with_extra_file(self):
        src,out,_=self.run_fixture();f=out/'day05_sha256_manifest.json';m=json.loads(f.read_text());m['files']['../../secret']='0'*64;f.write_text(json.dumps(m))
        with self.assertRaises(ValueError):lab.verify_lab(src,out)
    def test_refuse_modified_manifest_source_name(self):
        src,out,_=self.run_fixture();f=out/'day05_sha256_manifest.json';m=json.loads(f.read_text());m['source_fixture_name']='other.json';f.write_text(json.dumps(m))
        with self.assertRaises(ValueError):lab.verify_lab(src,out)
    def test_refuse_modified_manifest_api_count(self):
        src,out,_=self.run_fixture();f=out/'day05_sha256_manifest.json';m=json.loads(f.read_text());m['real_api_calls_this_run']=1;f.write_text(json.dumps(m))
        with self.assertRaises(ValueError):lab.verify_lab(src,out)
    def test_refuse_linked_output_file(self):
        src,out,_=self.run_fixture();p=out/lab.ARTIFACTS[0];orig=p.read_bytes();p.unlink();other=self.root/'elsewhere';other.write_bytes(orig);p.symlink_to(other)
        with self.assertRaises(ValueError):lab.verify_lab(src,out)
    def test_byte_for_byte_repeat(self):
        src,out,_=self.run_fixture();p=self.root/'again';lab.run_lab(src,p)
        self.assertTrue(all((out/n).read_bytes()==(p/n).read_bytes() for n in lab.ARTIFACTS))
    def test_manifest_byte_repeat(self):
        src,out,_=self.run_fixture();p=self.root/'again';lab.run_lab(src,p)
        self.assertEqual((out/'day05_sha256_manifest.json').read_bytes(),(p/'day05_sha256_manifest.json').read_bytes())


class Day05CliTests(FixtureMixin, unittest.TestCase):
    def cmd(self, *args, environ=None):
        env=os.environ.copy();env.update(environ or {})
        return subprocess.run([sys.executable,"-B",str(HERE/'day05_offline_model_route_lab.py'),*args],cwd=self.root,env=env,capture_output=True,text=True,timeout=15)

    def test_cli_run(self):
        r=self.cmd('run','--mode','offline','--fixture',str(SOURCE),'--output-dir','out');self.assertEqual(r.returncode,0,r.stderr);self.assertEqual(json.loads(r.stdout)['case_count'],7)
    def test_cli_verify(self):
        r=self.cmd('run','--mode','offline','--fixture',str(SOURCE),'--output-dir','out');self.assertEqual(r.returncode,0)
        v=self.cmd('verify','--mode','offline','--fixture',str(SOURCE),'--output-dir','out');self.assertEqual(v.returncode,0,v.stderr);self.assertEqual(json.loads(v.stdout)['files_verified'],5)
    def test_cli_live_rejected(self):
        r=self.cmd('run','--mode','live','--fixture',str(SOURCE),'--output-dir','out');self.assertNotEqual(r.returncode,0);self.assertFalse((self.root/'out').exists())
    def test_cli_live_with_key_rejected(self):
        r=self.cmd('run','--mode','live','--fixture',str(SOURCE),'--output-dir','out',environ={'ALPHAGENOME_API_KEY':'fake-placeholder'});self.assertNotEqual(r.returncode,0);self.assertFalse((self.root/'out').exists())
    def test_cli_missing_fixture_rejected(self):
        r=self.cmd('run','--mode','offline','--fixture','nonexistent','--output-dir','out');self.assertNotEqual(r.returncode,0)
    def test_cli_missing_mode_rejected(self):
        r=self.cmd('run','--fixture',str(SOURCE),'--output-dir','out');self.assertNotEqual(r.returncode,0)
    def test_cli_missing_output_rejected(self):
        r=self.cmd('run','--mode','offline','--fixture',str(SOURCE));self.assertNotEqual(r.returncode,0)
    def test_cli_secret_never_printed(self):
        fake_secret='TEST-NEVER-PRINT-ALPHAGENOME-KEY-1234567';r=self.cmd('run','--mode','offline','--fixture',str(SOURCE),'--output-dir','out',environ={'ALPHAGENOME_API_KEY':fake_secret});self.assertEqual(r.returncode,0);self.assertNotIn(fake_secret,r.stdout+r.stderr)
    def test_script_does_not_import_network_libraries(self):
        content=(HERE/'day05_offline_model_route_lab.py').read_text(encoding='utf-8')
        import ast
        tree=ast.parse(content)
        imports=[]
        for node in ast.walk(tree):
            if isinstance(node,ast.Import): imports.extend(a.name.split('.')[0] for a in node.names)
            if isinstance(node,ast.ImportFrom) and node.module: imports.append(node.module.split('.')[0])
        self.assertTrue({'requests','urllib','http','socket','alphagenome','google'}.isdisjoint(imports))
    def test_cli_no_model_claim(self):
        r=self.cmd('run','--mode','offline','--fixture',str(SOURCE),'--output-dir','out');self.assertEqual(json.loads(r.stdout)['model_predictions_made'],0)


if __name__ == '__main__':
    unittest.main(verbosity=2)
