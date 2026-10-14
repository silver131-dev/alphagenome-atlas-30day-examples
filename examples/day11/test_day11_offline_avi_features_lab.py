"""Day11 offline feature catalogue, evidence boundaries, CLI and reproducibility tests."""
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
from pathlib import Path
from unittest.mock import patch

import day11_offline_avi_features_lab as lab

ROOT=Path(__file__).resolve().parent
FIXTURE_PATH=ROOT/'day11_synthetic_avi_features.json'
COMPARISON_PATH=ROOT/'day11_synthetic_avi_comparison.json'
SCRIPT=ROOT/'day11_offline_avi_features_lab.py'


def fixture():
    return json.loads(FIXTURE_PATH.read_text(encoding='utf-8'))


def comparison():
    return json.loads(COMPARISON_PATH.read_text(encoding='utf-8'))


def call_cli(*args, cwd=None, env=None):
    return subprocess.run([sys.executable, '-B', str(SCRIPT), *map(str, args)],cwd=cwd,env=env,capture_output=True,text=True,check=False)


class CatalogueTests(unittest.TestCase):
    def test_canonical_name_count(self):
        self.assertEqual(len(lab.FEATURE_NAMES),18)

    def test_unique_feature_names(self):
        self.assertEqual(len(set(lab.FEATURE_NAMES)),18)

    def test_canonical_group_count(self):
        self.assertEqual(len(lab.GROUP_ORDER),5)

    def test_groups_match_research_counts(self):
        self.assertEqual(lab.GROUP_COUNT,{'AlphaGenome':10,'AlphaMissense':1,'VEP consequences':3,'Conservation':2,'Indel flags':2})

    def test_group_count_sums_to_18(self):
        self.assertEqual(sum(lab.GROUP_COUNT.values()),18)

    def test_indel_flags_are_both_present(self):
        self.assertIn('IS_INSERTION',lab.FEATURE_NAMES)
        self.assertIn('IS_DELETION',lab.FEATURE_NAMES)

    def test_splicing_is_merged(self):
        self.assertIn('MERGED_SPLICING',lab.FEATURE_NAMES)
        self.assertEqual(sum(n=='MERGED_SPLICING' for n in lab.FEATURE_NAMES),1)

    def test_conservation_is_two_inputs(self):
        self.assertEqual([n for n,g,_ in lab.FEATURES if g=='Conservation'],['CACTUS_241_WAY','PHASTCONS_470_WAY'])

    def test_catalog_no_empty_descriptions(self):
        self.assertTrue(all(d.strip() for _,_,d in lab.FEATURES))

    def test_binary_feature_count(self):
        self.assertEqual(len(lab.BINARY_FEATURES),5)

    def test_no_external_package_import_in_lab(self):
        tree=ast.parse(SCRIPT.read_text(encoding='utf-8'))
        roots={n.names[0].name.split('.')[0] for n in ast.walk(tree) if isinstance(n,ast.Import)}
        roots|={n.module.split('.')[0] for n in ast.walk(tree) if isinstance(n,ast.ImportFrom) and n.module}
        self.assertFalse(roots & {'alphagenome','requests','httpx','urllib3','aiohttp','google','dotenv'})


class FixtureTests(unittest.TestCase):
    def test_fixture_three_cases(self):
        self.assertEqual(len(lab.validate_fixture(fixture())),3)

    def test_toy_assembly(self):
        self.assertEqual(fixture()['assembly'],'TOY_NOT_GRCH38')

    def test_non_synthetic_rejected(self):
        x=fixture();x['data_kind']='REAL_PATIENT'
        with self.assertRaises(ValueError):lab.validate_fixture(x)

    def test_unrecognized_top_level_field_rejected(self):
        x=fixture();x['patient_id']='123'
        with self.assertRaises(ValueError):lab.validate_fixture(x)

    def test_dataset_role_mismatch_rejected(self):
        x=fixture();x['dataset_role']='REAL_API_RESPONSE'
        with self.assertRaises(ValueError):lab.validate_fixture(x)

    def test_empty_cases_rejected(self):
        x=fixture();x['cases']=[]
        with self.assertRaises(ValueError):lab.validate_fixture(x)

    def test_too_many_cases_rejected(self):
        x=fixture();x['cases']=[copy.deepcopy(x['cases'][0]) for _ in range(21)]
        with self.assertRaises(ValueError):lab.validate_fixture(x)

    def test_duplicate_ids_rejected(self):
        x=fixture();x['cases'][1]['case_id']=x['cases'][0]['case_id']
        with self.assertRaises(ValueError):lab.validate_fixture(x)

    def test_bad_variant_name_rejected(self):
        x=fixture();x['cases'][0]['variant_label']='chr9:128226027:G>A'
        with self.assertRaises(ValueError):lab.validate_fixture(x)

    def test_newline_case_id_rejected(self):
        x=fixture();x['cases'][0]['case_id']='TOY_A\n# injected'
        with self.assertRaises(ValueError):lab.validate_fixture(x)

    def test_boolean_score_rejected(self):
        x=fixture();x['cases'][0]['raw_score']=True
        with self.assertRaises(ValueError):lab.validate_fixture(x)

    def test_nan_score_rejected(self):
        x=fixture();x['cases'][0]['raw_score']=float('nan')
        with self.assertRaises(ValueError):lab.validate_fixture(x)

    def test_infinite_input_rejected(self):
        x=fixture();x['cases'][0]['feature_values']['MAX_ABS_ATAC']=float('inf')
        with self.assertRaises(ValueError):lab.validate_fixture(x)

    def test_missing_feature_name_rejected(self):
        x=fixture();del x['cases'][0]['feature_values']['MAX_ABS_ATAC']
        with self.assertRaises(ValueError):lab.validate_fixture(x)

    def test_extra_feature_name_rejected(self):
        x=fixture();x['cases'][0]['feature_values']['PATIENT_ID']=1
        with self.assertRaises(ValueError):lab.validate_fixture(x)

    def test_negative_max_abs_rejected(self):
        x=fixture();x['cases'][0]['feature_values']['MAX_ABS_ATAC']=-0.5
        with self.assertRaises(ValueError):lab.validate_fixture(x)

    def test_negative_cactus_is_allowed(self):
        x=fixture();x['cases'][0]['feature_values']['CACTUS_241_WAY']=-0.5
        self.assertEqual(lab.validate_fixture(x)[0]['feature_values']['CACTUS_241_WAY'],-0.5)

    def test_phastcons_outside_probability_range_rejected(self):
        x=fixture();x['cases'][0]['feature_values']['PHASTCONS_470_WAY']=1.5
        with self.assertRaises(ValueError):lab.validate_fixture(x)

    def test_indel_flag_outside_binary_rejected(self):
        x=fixture();x['cases'][0]['feature_values']['IS_INSERTION']=0.4
        with self.assertRaises(ValueError):lab.validate_fixture(x)

    def test_two_indel_flags_cannot_both_be_true(self):
        x=fixture();x['cases'][0]['feature_values']['IS_INSERTION']=1;x['cases'][0]['feature_values']['IS_DELETION']=1
        with self.assertRaises(ValueError):lab.validate_fixture(x)

    def test_attribution_for_indicator_is_not_binary(self):
        x=fixture();x['cases'][0]['attributions']['START_LOST']=0.2
        lab.validate_fixture(x)

    def test_attribution_boolean_rejected(self):
        x=fixture();x['cases'][0]['attributions']['MERGED_SPLICING']=False
        with self.assertRaises(ValueError):lab.validate_fixture(x)

    def test_attribution_omitted_with_raw_score_rejected(self):
        x=fixture();x['cases'][0]['attributions']=None
        with self.assertRaises(ValueError):lab.validate_fixture(x)

    def test_attribution_present_with_null_score_rejected(self):
        x=fixture();x['cases'][2]['attributions']=x['cases'][0]['attributions']
        with self.assertRaises(ValueError):lab.validate_fixture(x)

    def test_null_inside_complete_attributions_rejected(self):
        x=fixture();x['cases'][0]['attributions']['MERGED_SPLICING']=None
        with self.assertRaises(ValueError):lab.validate_fixture(x)

    def test_missing_case_keeps_nulls(self):
        out=lab.validate_fixture(fixture())[-1]
        self.assertIsNone(out['raw_score'])
        self.assertIsNone(out['attributions'])
        self.assertTrue(all(v is None for v in out['feature_values'].values()))

    def test_file_has_no_real_genome_coordinates(self):
        text=FIXTURE_PATH.read_text(encoding='utf-8')
        self.assertNotIn('chr22:36201698',text)
        self.assertNotIn('chr9:128225',text)


class SyntheticComparisonTests(unittest.TestCase):
    def test_comparison_summary_class(self):
        x=lab.validate_comparison(comparison())
        self.assertEqual(x['evidence_class'],'SYNTHETIC_EDUCATIONAL_DEMO')

    def test_comparison_is_not_native(self):
        self.assertFalse(lab.validate_comparison(comparison())['native_response_retained'])

    def test_comparison_residual(self):
        x=lab.validate_comparison(comparison())
        self.assertAlmostEqual(x['additivity_residual'],0.10)

    def test_comparison_additivity_review(self):
        self.assertEqual(lab.validate_comparison(comparison())['additivity_status'],'REVIEW')

    def test_comparison_claiming_sdk_is_rejected(self):
        x=comparison();x['sdk_version']='UNAUTHORIZED_DEMO'
        with self.assertRaises(ValueError):lab.validate_comparison(x)

    def test_comparison_native_false_required(self):
        x=comparison();x['native_response_retained']=True
        with self.assertRaises(ValueError):lab.validate_comparison(x)

    def test_comparison_clinical_false_required(self):
        x=comparison();x['clinical_interpretation']=True
        with self.assertRaises(ValueError):lab.validate_comparison(x)

    def test_comparison_wrong_variant_rejected(self):
        x=comparison();x['variant_label']='chr1:123:A>G'
        with self.assertRaises(ValueError):lab.validate_comparison(x)

    def test_comparison_wrong_feature_count_rejected(self):
        x=comparison();x['feature_count']=17
        with self.assertRaises(ValueError):lab.validate_comparison(x)

    def test_comparison_residual_corrupted_rejected(self):
        x=comparison();x['additivity_residual']=1.0
        with self.assertRaises(ValueError):lab.validate_comparison(x)

    def test_comparison_missing_score_rejected(self):
        x=comparison();x['avi_raw']=None
        with self.assertRaises(ValueError):lab.validate_comparison(x)

    def test_comparison_optional(self):
        self.assertIsNone(lab.validate_comparison(None))

    def test_comparison_rejects_historical_output_claim(self):
        x=comparison();x['evidence_class']='DERIVED_FROM_HISTORICAL_OUTPUT'
        with self.assertRaises(ValueError):lab.validate_comparison(x)

    def test_comparison_rejects_live_model_claim(self):
        x=comparison();x['model_executed']=True
        with self.assertRaises(ValueError):lab.validate_comparison(x)


class ReportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report=lab.build_report(fixture(),comparison())
        cls.cases=lab.validate_fixture(fixture())

    def test_offline_mode(self):
        self.assertEqual(self.report['mode'],'OFFLINE_ONLY')

    def test_evidence_level(self):
        self.assertEqual(self.report['evidence_class'],'SYNTHETIC_EDUCATIONAL_DEMO')

    def test_case_count(self):
        self.assertEqual(self.report['case_count'],3)

    def test_feature_count_18(self):
        self.assertEqual(self.report['feature_count'],18)

    def test_group_count_5(self):
        self.assertEqual(len(self.report['feature_group_counts']),5)

    def test_toy_splicing_additivity_pass(self):
        self.assertEqual(self.report['cases'][0]['additivity_status'],'PASS')
        self.assertAlmostEqual(self.report['cases'][0]['additivity_sum'] if 'additivity_sum' in self.report['cases'][0] else self.report['cases'][0]['attribution_sum'],1.07)

    def test_toy_conservation_additivity_review(self):
        self.assertEqual(self.report['cases'][1]['additivity_status'],'REVIEW')
        self.assertAlmostEqual(self.report['cases'][1]['additivity_residual'],0.05)

    def test_missing_case_not_applicable(self):
        self.assertEqual(self.report['cases'][2]['additivity_status'],'NOT_APPLICABLE_MISSING')

    def test_missing_raw_not_zero(self):
        self.assertIsNone(self.report['cases'][2]['raw_score'])

    def test_missing_count(self):
        self.assertEqual(self.report['missing_case_count'],1)

    def test_review_count(self):
        self.assertEqual(self.report['review_count'],1)

    def test_toy_top_splicing(self):
        self.assertEqual(self.report['cases'][0]['top_attribution']['feature'],'MERGED_SPLICING')

    def test_toy_top_conservation(self):
        self.assertEqual(self.report['cases'][1]['top_attribution']['feature'],'CACTUS_241_WAY')

    def test_comparison_kept_separate_from_toy(self):
        self.assertNotIn(self.report['comparison_summary']['variant_label'],json.dumps(self.report['cases']))

    def test_comparison_review_kept(self):
        self.assertEqual(self.report['comparison_additivity_status'],'REVIEW')

    def test_manual_review(self):
        self.assertEqual(self.report['manual_review'],'REQUIRED')

    def test_release_not_granted(self):
        self.assertEqual(self.report['publication_permission'],'NOT_GRANTED')

    def test_no_real_api(self):
        self.assertEqual(self.report['real_api_calls_this_run'],0)
        self.assertFalse(self.report['official_avi_model_executed'])

    def test_no_key_read(self):
        self.assertFalse(self.report['api_key_read'])

    def test_no_clinical_claim(self):
        self.assertFalse(self.report['clinical_interpretation'])

    def test_abs_attribution_sorting(self):
        x=self.report['cases'][0]['top5_attributions']
        self.assertTrue(all(abs(x[i]['contribution']) >= abs(x[i+1]['contribution']) for i in range(len(x)-1)))

    def test_csv_value_rows_54(self):
        data=lab.make_csv(self.cases,attribution=False)
        self.assertEqual(len(list(csv.DictReader(io.StringIO(data)))),54)

    def test_csv_attribution_rows_54(self):
        data=lab.make_csv(self.cases,attribution=True)
        self.assertEqual(len(list(csv.DictReader(io.StringIO(data)))),54)

    def test_csv_null_is_blank(self):
        rows=list(csv.DictReader(io.StringIO(lab.make_csv(self.cases,attribution=True))))
        self.assertEqual([x['value'] for x in rows if x['case_id']=='TOY_MISSING'],['']*18)

    def test_csv_no_comparison_variant(self):
        self.assertNotIn('TOY_COMPARISON',lab.make_csv(self.cases,attribution=True))

    def test_svg_dimensions_1920x1080(self):
        svg=lab.make_svg(self.report,self.cases)
        self.assertIn('width="1920"',svg)
        self.assertIn('height="1080"',svg)

    def test_svg_background(self):
        self.assertIn('#E5DEDA',lab.make_svg(self.report,self.cases))

    def test_svg_comparison_review(self):
        svg=lab.make_svg(self.report,self.cases)
        self.assertIn('SYNTHETIC COMPARISON',svg)
        self.assertIn('REVIEW',svg)

    def test_svg_xml_parses(self):
        import xml.etree.ElementTree as ET
        ET.fromstring(lab.make_svg(self.report,self.cases))

    def test_html_contains_attribution_table(self):
        page=lab.make_html(self.report,self.cases)
        self.assertIn('MERGED_SPLICING',page)
        self.assertIn('CACTUS_241_WAY',page)

    def test_html_no_remote_js(self):
        page=lab.make_html(self.report,self.cases)
        self.assertNotIn('<script',page)
        self.assertNotIn('http:',page)
        self.assertNotIn('https:',page)

    def test_html_has_toy_disclaimer(self):
        self.assertIn('TOY_NOT_GRCH38',lab.make_html(self.report,self.cases))

    def test_notes_no_clinical_inference(self):
        self.assertIn('不是疾病風險',lab.make_notes(self.report))

    def test_comparison_optional_status(self):
        r=lab.build_report(fixture())
        self.assertEqual(r['comparison_additivity_status'],'NOT_PROVIDED')


class ArtifactTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.out=Path(self.temp.name)/'output'

    def generate(self):
        return lab.write_files(self.out,FIXTURE_PATH,COMPARISON_PATH)

    def test_run_creates_7_files(self):
        self.generate()
        self.assertEqual(len(list(self.out.iterdir())),7)

    def test_run_json_has_evidence(self):
        self.generate();r=json.loads((self.out/'day11_feature_report.json').read_text())
        self.assertEqual(r['comparison_additivity_status'],'REVIEW')

    def test_run_manifest_has_inputs(self):
        self.generate();man=json.loads((self.out/'day11_manifest.json').read_text())
        self.assertEqual(man['inputs']['fixture_sha256'],lab.sha256_file(FIXTURE_PATH))

    def test_run_manifest_unsigned(self):
        self.generate();man=json.loads((self.out/'day11_manifest.json').read_text())
        self.assertIs(man['signed'],False)

    def test_run_manifest_release_blocked(self):
        self.generate();man=json.loads((self.out/'day11_manifest.json').read_text())
        self.assertEqual(man['release_gate'],'BLOCKED')

    def test_verify_passes(self):
        self.generate();self.assertEqual(lab.verify_files(self.out,FIXTURE_PATH,COMPARISON_PATH)['status'],'PASS')

    def test_refuses_existing_output(self):
        self.generate()
        with self.assertRaises(FileExistsError):self.generate()

    def test_verify_detects_tampered_csv(self):
        self.generate();p=self.out/'day11_feature_values.csv';p.write_text('TAMPERED',encoding='utf-8')
        with self.assertRaises(ValueError):lab.verify_files(self.out,FIXTURE_PATH,COMPARISON_PATH)

    def test_verify_detects_tampered_svg(self):
        self.generate();p=self.out/'day11_feature_map.svg';p.write_text('TAMPERED',encoding='utf-8')
        with self.assertRaises(ValueError):lab.verify_files(self.out,FIXTURE_PATH,COMPARISON_PATH)

    def test_verify_detects_tampered_html(self):
        self.generate();p=self.out/'day11_feature_dashboard.html';p.write_text('TAMPERED',encoding='utf-8')
        with self.assertRaises(ValueError):lab.verify_files(self.out,FIXTURE_PATH,COMPARISON_PATH)

    def test_verify_detects_tampered_json(self):
        self.generate();p=self.out/'day11_feature_report.json';p.write_text('TAMPERED',encoding='utf-8')
        with self.assertRaises(ValueError):lab.verify_files(self.out,FIXTURE_PATH,COMPARISON_PATH)

    def test_verify_detects_tampered_manifest(self):
        self.generate();p=self.out/'day11_manifest.json';m=json.loads(p.read_text());m['files']['day11_feature_report.json']='0'*64;p.write_text(json.dumps(m))
        with self.assertRaises(ValueError):lab.verify_files(self.out,FIXTURE_PATH,COMPARISON_PATH)

    def test_verify_detects_missing_output(self):
        self.generate();(self.out/'day11_feature_map.svg').unlink()
        with self.assertRaises(ValueError):lab.verify_files(self.out,FIXTURE_PATH,COMPARISON_PATH)

    def test_verify_detects_added_output(self):
        self.generate();(self.out/'surprise.txt').write_text('extra')
        with self.assertRaises(ValueError):lab.verify_files(self.out,FIXTURE_PATH,COMPARISON_PATH)

    def test_verify_detects_modified_fixture(self):
        self.generate();p=Path(self.temp.name)/'other.json';x=fixture();x['cases'][0]['raw_score']=0.2;p.write_text(json.dumps(x))
        with self.assertRaises(ValueError):lab.verify_files(self.out,p,COMPARISON_PATH)

    def test_run_is_bit_reproducible(self):
        self.generate();one={p.name:p.read_bytes() for p in self.out.iterdir()}
        other=Path(self.temp.name)/'second'
        lab.write_files(other,FIXTURE_PATH,COMPARISON_PATH)
        self.assertEqual(one,{p.name:p.read_bytes() for p in other.iterdir()})

    def test_abs_paths_not_embedded_in_manifest(self):
        self.generate();self.assertNotIn(self.temp.name,(self.out/'day11_manifest.json').read_text())

    def test_verify_missing_manifest_raises(self):
        self.generate();(self.out/'day11_manifest.json').unlink()
        with self.assertRaises(FileNotFoundError):lab.verify_files(self.out,FIXTURE_PATH,COMPARISON_PATH)

    def test_sha_method_matches_hashlib(self):
        self.assertEqual(lab.sha256_file(FIXTURE_PATH),hashlib.sha256(FIXTURE_PATH.read_bytes()).hexdigest())


class CliTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.out=Path(self.temp.name)/'day11_output'

    def run_cli(self,*custom):
        return call_cli('run','--mode','offline','--fixture',FIXTURE_PATH,'--comparison-fixture',COMPARISON_PATH,'--output-dir',self.out,*custom,cwd=self.temp.name)

    def verify_cli(self,*custom):
        return call_cli('verify','--mode','offline','--fixture',FIXTURE_PATH,'--comparison-fixture',COMPARISON_PATH,'--output-dir',self.out,*custom,cwd=self.temp.name)

    def test_cli_run_pass(self):
        c=self.run_cli();self.assertEqual(c.returncode,0,c.stderr)
        self.assertEqual(json.loads(c.stdout)['status'],'PASS')

    def test_cli_verify_pass(self):
        self.assertEqual(self.run_cli().returncode,0)
        c=self.verify_cli();self.assertEqual(c.returncode,0,c.stderr)
        self.assertEqual(json.loads(c.stdout)['checksums'],'MATCH')

    def test_cli_real_api_calls_zero(self):
        c=self.run_cli();self.assertEqual(json.loads(c.stdout)['real_api_calls_this_run'],0)

    def test_cli_live_mode_refused(self):
        c=call_cli('run','--mode','live','--fixture',FIXTURE_PATH,'--comparison-fixture',COMPARISON_PATH,'--output-dir',self.out)
        self.assertNotEqual(c.returncode,0)
        self.assertIn('LIVE_DISABLED',c.stderr)
        self.assertFalse(self.out.exists())

    def test_cli_live_refused_even_with_key(self):
        env=os.environ.copy();env['ALPHAGENOME_API_KEY']='INTENTIONALLY_FAKE_NEVER_USABLE'
        c=call_cli('run','--mode','live','--fixture',FIXTURE_PATH,'--comparison-fixture',COMPARISON_PATH,'--output-dir',self.out,env=env)
        self.assertNotEqual(c.returncode,0)
        self.assertNotIn('INTENTIONALLY_FAKE',c.stderr)
        self.assertFalse(self.out.exists())

    def test_cli_invalid_fixture_rejected(self):
        bad=Path(self.temp.name)/'bad.json';bad.write_text('{"data_kind":"PATIENT"}')
        c=call_cli('run','--fixture',bad,'--comparison-fixture',COMPARISON_PATH,'--output-dir',self.out)
        self.assertNotEqual(c.returncode,0)
        self.assertFalse(self.out.exists())

    def test_cli_existing_output_rejected(self):
        self.assertEqual(self.run_cli().returncode,0)
        self.assertNotEqual(self.run_cli().returncode,0)

    def test_cli_verify_missing_file_fails(self):
        self.assertEqual(self.run_cli().returncode,0)
        (self.out/'day11_feature_values.csv').unlink()
        self.assertNotEqual(self.verify_cli().returncode,0)

    def test_cli_from_unrelated_cwd(self):
        c=call_cli('run','--fixture',FIXTURE_PATH,'--comparison-fixture',COMPARISON_PATH,'--output-dir',self.out,cwd='/tmp')
        self.assertEqual(c.returncode,0,c.stderr)

    def test_cli_no_api_key_required(self):
        env=os.environ.copy();env.pop('ALPHAGENOME_API_KEY',None);env.pop('ALPHA_GENOME_API_KEY',None)
        c=call_cli('run','--fixture',FIXTURE_PATH,'--comparison-fixture',COMPARISON_PATH,'--output-dir',self.out,env=env)
        self.assertEqual(c.returncode,0,c.stderr)

    def test_cli_verify_without_run_fails(self):
        self.assertNotEqual(self.verify_cli().returncode,0)

    def test_cli_help(self):
        c=call_cli('--help');self.assertEqual(c.returncode,0);self.assertIn('OFFLINE synthetic AVI',c.stdout)


# Each of the 18 official-style feature keys must survive the fixture -> report -> CSV contract.
def _make_feature_contract_test(name):
    def test(self):
        x=fixture();r=lab.build_report(x,comparison())
        self.assertEqual(sum(z['name']==name for z in r['feature_catalog']),1)
        self.assertIn(name,lab.make_csv(lab.validate_fixture(x),attribution=False))
        self.assertIn(name,lab.make_csv(lab.validate_fixture(x),attribution=True))
    return test

for _name in lab.FEATURE_NAMES:
    setattr(CatalogueTests,'test_feature_contract_'+_name.lower(),_make_feature_contract_test(_name))


if __name__=='__main__':
    unittest.main(verbosity=2)
