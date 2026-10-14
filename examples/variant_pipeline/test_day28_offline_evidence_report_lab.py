"""Day28 standalone tests; synthetic fixtures only; stdlib unittest."""
from __future__ import annotations

import copy
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

import day28_offline_evidence_report_lab as mod

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'example_output'/'day27_output'
SCRIPT=ROOT/'day28_offline_evidence_report_lab.py'


class Day28SchemaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report=mod.report_from_day27(SOURCE)
        cls.rows=copy.deepcopy(mod.d27.load_json(SOURCE/'day23'/'day23_annotation_report.json')['annotations'])

    def test_schema(self):self.assertEqual(self.report['schema_version'],mod.SCHEMA)
    def test_mode_offline_only(self):self.assertEqual(self.report['mode'],'OFFLINE_ONLY')
    def test_assembly_is_explicit_toy(self):self.assertEqual(self.report['assembly'],'TOY_NOT_GRCH38')
    def test_annotation_counts(self):self.assertEqual((len(self.report['annotations']),self.report['pipeline_counts']['annotated_count']), (3,2))
    def test_day21_counts(self):
        c=self.report['pipeline_counts'];self.assertEqual((c['candidate_count'],c['qc_selected'],c['qc_review_required'],c['qc_rejected']),(8,3,2,3))
    def test_original_source_artifacts_23(self):self.assertEqual(self.report['source_provenance']['artifact_count_verified'],23)
    def test_day27_manifest_hash(self):
        self.assertEqual(self.report['source_provenance']['day27_run_manifest_sha256'],mod.sha256(SOURCE/'day27_run_manifest.json'))
    def test_day27_summary_hash(self):
        self.assertEqual(self.report['source_provenance']['day27_summary_sha256'],mod.sha256(SOURCE/'day27_pipeline_summary.json'))
    def test_day27_source_bundle_sha_64(self):self.assertRegex(self.report['source_provenance']['source_bundle_sha256'],r'^[0-9a-f]{64}$')
    def test_source_manifest_unsigned(self):self.assertIs(self.report['source_provenance']['signature_verified'],False)
    def test_third_missing_not_zero(self):
        a=self.report['annotations'];self.assertEqual(a[-1]['synthetic_score'],None);self.assertEqual(a[-1]['score_status'],'MISSING_SCORE')
    def test_missing_variant_id(self):self.assertEqual(self.report['missing_values']['missing_score_variant_ids'],['chr1:110:A>C'])
    def test_all_review_required(self):self.assertTrue(all(a['review_status']=='REVIEW_REQUIRED' for a in self.report['annotations']))
    def test_score_zero_not_missing(self):
        rows=copy.deepcopy(self.rows);rows[0]['synthetic_score']=0.0
        result,missing=mod.annotations_with_missing(rows)
        self.assertEqual(result[0]['synthetic_score'],0.0);self.assertEqual(len(missing),1)
    def test_missing_score_error(self):
        rows=copy.deepcopy(self.rows);rows[2]['score_status']='SYNTHETIC_SCORE'
        with self.assertRaises(ValueError):mod.annotations_with_missing(rows)
    def test_double_missing_error(self):
        rows=copy.deepcopy(self.rows);rows[0]['synthetic_score']=None;rows[0]['score_status']='MISSING_SCORE'
        with self.assertRaises(ValueError):mod.annotations_with_missing(rows)
    def test_invalid_boolean_score(self):
        rows=copy.deepcopy(self.rows);rows[0]['synthetic_score']=True
        with self.assertRaises(ValueError):mod.annotations_with_missing(rows)
    def test_invalid_inf_score(self):
        rows=copy.deepcopy(self.rows);rows[0]['synthetic_score']=float('inf')
        with self.assertRaises(ValueError):mod.annotations_with_missing(rows)
    def test_invalid_nan_score(self):
        rows=copy.deepcopy(self.rows);rows[0]['synthetic_score']=float('nan')
        with self.assertRaises(ValueError):mod.annotations_with_missing(rows)
    def test_duplicate_variant_rejected(self):
        rows=copy.deepcopy(self.rows);rows[1]['variant_id']=rows[0]['variant_id']
        with self.assertRaises(ValueError):mod.annotations_with_missing(rows)
    def test_non_toy_variant_rejected(self):
        rows=copy.deepcopy(self.rows);rows[1]['variant_id']='chrX:123:A>G'
        with self.assertRaises(ValueError):mod.annotations_with_missing(rows)
    def test_source_status_required(self):
        rows=copy.deepcopy(self.rows);rows[0]['score_source']='LIVE_MODEL'
        with self.assertRaises(ValueError):mod.annotations_with_missing(rows)
    def test_clinical_flags_explicit_false(self):
        self.assertIsNone(self.report['clinical_classification']);self.assertFalse(self.report['confidence_calibrated']);self.assertFalse(self.report['clinical_interpretation'])
    def test_score_scales_not_combined(self):self.assertIsNone(self.report['teaching_observations']['combined_biological_score'])
    def test_motif_two_deltas(self):self.assertEqual(self.report['teaching_observations']['motif_ppm_like']['MOTIF_B'],.12)
    def test_limitations_nonempty(self):self.assertGreaterEqual(len(self.report['model_limitations']),6)
    def test_no_real_api(self):self.assertEqual(self.report['real_api_calls_this_run'],0)


class Day28RenderingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.report=mod.report_from_day27(SOURCE)
    def test_markdown_contains_null(self):self.assertIn('`null`',mod.markdown_report(self.report))
    def test_markdown_contains_sha(self):self.assertIn('SHA-256',mod.markdown_report(self.report))
    def test_markdown_contains_uncertainty(self):self.assertIn('不能相加',mod.markdown_report(self.report))
    def test_html_escapes_source(self):
        r=copy.deepcopy(self.report);r['annotations'][0]['variant_id']='<b>&'
        h=mod.html_report(r)
        self.assertIn('&lt;b&gt;&amp;',h);self.assertNotIn('<b>&</code>',h)
    def test_html_no_external_assets(self):
        html=mod.html_report(self.report)
        self.assertNotIn('<script',html);self.assertNotIn('https://',html)
    def test_html_links(self):
        h=mod.html_report(self.report)
        for f in ('day28_evidence_report.json','day28_evidence_map.svg','day28_evidence_report.md'):
            self.assertIn(f,h)
    def test_csv_three_rows(self):self.assertEqual(len(list(csv.DictReader(io.StringIO(mod.report_csv(self.report))))),3)
    def test_csv_missing_blank(self):
        rows=list(csv.DictReader(io.StringIO(mod.report_csv(self.report))))
        self.assertEqual(rows[-1]['synthetic_score'],'');self.assertEqual(rows[-1]['score_status'],'MISSING_SCORE')
    def test_svg_xml_parse(self):self.assertTrue(ET.fromstring(mod.draw_evidence_svg(self.report)).tag.endswith('svg'))
    def test_svg_16_9(self):
        svg=ET.fromstring(mod.draw_evidence_svg(self.report));self.assertEqual((svg.get('width'),svg.get('height')) ,('1920','1080'))
    def test_svg_background(self):self.assertIn('#E5DEDA',mod.draw_evidence_svg(self.report))
    def test_svg_manual_review(self):self.assertIn('REQUIRED',mod.draw_evidence_svg(self.report))


class Day28IntegrationTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix='day28-test-')
        self.root=Path(self.tmp.name)
    def tearDown(self):self.tmp.cleanup()
    def test_emit_five_artifacts_plus_manifest(self):
        target=self.root/'report';mod.write_outputs(SOURCE,target)
        self.assertEqual(set(p.name for p in target.iterdir()), set(mod.ARTIFACTS)|{'day28_report_manifest.json'})
    def test_output_manifest_sha_matches(self):
        target=self.root/'report';mod.write_outputs(SOURCE,target)
        manifest=mod.d27.load_json(target/'day28_report_manifest.json')
        for name,v in manifest['artifacts'].items():self.assertEqual(mod.sha256(target/name),v['sha256'])
    def test_verify_day28_against_day27(self):
        target=self.root/'report';mod.write_outputs(SOURCE,target)
        self.assertEqual(mod.verify_report(SOURCE,target)['source_alignment'],'MATCH')
    def test_output_never_overwritten(self):
        target=self.root/'report';mod.write_outputs(SOURCE,target)
        with self.assertRaises(ValueError):mod.write_outputs(SOURCE,target)
    def test_invalid_source_no_output(self):
        target=self.root/'report'
        with self.assertRaises(ValueError):mod.write_outputs(self.root/'missing',target)
        self.assertFalse(target.exists())
    def test_report_tampering_rejected(self):
        target=self.root/'report';mod.write_outputs(SOURCE,target)
        with (target/'day28_evidence_report.md').open('a') as f:f.write('forged')
        with self.assertRaises(ValueError):mod.verify_report(SOURCE,target)
    def test_report_json_tampering_rejected(self):
        target=self.root/'report';mod.write_outputs(SOURCE,target)
        with (target/'day28_evidence_report.json').open('a') as f:f.write(' ')
        with self.assertRaises(ValueError):mod.verify_report(SOURCE,target)
    def test_report_missing_artifact_rejected(self):
        target=self.root/'report';mod.write_outputs(SOURCE,target)
        (target/'day28_evidence_table.csv').unlink()
        with self.assertRaises(ValueError):mod.verify_report(SOURCE,target)
    def test_manifest_bad_schema_rejected(self):
        target=self.root/'report';mod.write_outputs(SOURCE,target)
        p=target/'day28_report_manifest.json'; a=json.loads(p.read_text());a['schema_version']='broken';p.write_text(json.dumps(a))
        with self.assertRaises(ValueError):mod.verify_report(SOURCE,target)
    def test_source_tampering_rejected(self):
        run=self.root/'src';shutil.copytree(SOURCE,run)
        p=run/'day23/day23_annotation_report.json'
        with p.open('a') as f:f.write('tampered')
        with self.assertRaises(ValueError):mod.report_from_day27(run)
    def test_source_missing_artifact_rejected(self):
        run=self.root/'src';shutil.copytree(SOURCE,run)
        (run/'day22/day22_qc_report.json').unlink()
        with self.assertRaises(ValueError):mod.report_from_day27(run)
    def test_source_path_traversal_rejected(self):
        run=self.root/'src';shutil.copytree(SOURCE,run)
        p=run/'day27_run_manifest.json';a=json.loads(p.read_text());a['artifacts']['../escape']={'sha256':'0'*64,'bytes':0};a['artifact_count']+=1;p.write_text(json.dumps(a))
        with self.assertRaises(ValueError):mod.report_from_day27(run)
    def test_source_symlink_rejected(self):
        link=self.root/'sym';link.symlink_to(SOURCE,target_is_directory=True)
        with self.assertRaises(ValueError):mod.report_from_day27(link)
    def test_report_symlink_rejected(self):
        target=self.root/'report';mod.write_outputs(SOURCE,target)
        link=self.root/'sym';link.symlink_to(target,target_is_directory=True)
        with self.assertRaises(ValueError):mod.verify_report(SOURCE,link)
    def test_clean_two_runs_same_bytes(self):
        a=self.root/'a';b=self.root/'b';mod.write_outputs(SOURCE,a);mod.write_outputs(SOURCE,b)
        for name in (*mod.ARTIFACTS,'day28_report_manifest.json'):
            self.assertEqual((a/name).read_bytes(),(b/name).read_bytes(),name)
    def test_source_files_unchanged(self):
        files=(SOURCE/'day27_pipeline_summary.json',SOURCE/'day27_run_manifest.json')
        before=[mod.sha256(p) for p in files];mod.write_outputs(SOURCE,self.root/'report')
        self.assertEqual(before,[mod.sha256(p) for p in files])


class Day28CliTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix='day28-cli-'); self.root=Path(self.tmp.name)
    def tearDown(self):self.tmp.cleanup()
    def invoke(self,*args,cwd=None,env=None):
        return subprocess.run([sys.executable,'-B',str(SCRIPT),*args],cwd=cwd or ROOT,capture_output=True,text=True,env=env)
    def test_main_command(self):
        r=self.invoke('run','--mode','offline','--run-dir',str(SOURCE),'--output-dir',str(self.root/'report'))
        self.assertEqual(r.returncode,0,r.stderr)
        self.assertEqual(json.loads(r.stdout)['missing_score_count'],1)
    def test_verify_command(self):
        mod.write_outputs(SOURCE,self.root/'report')
        r=self.invoke('verify','--run-dir',str(SOURCE),'--report-dir',str(self.root/'report'))
        self.assertEqual(r.returncode,0,r.stderr)
        self.assertEqual(json.loads(r.stdout)['checksums'],'MATCH')
    def test_live_rejected(self):
        r=self.invoke('run','--mode','live','--run-dir',str(SOURCE),'--output-dir',str(self.root/'report'))
        self.assertNotEqual(r.returncode,0);self.assertFalse((self.root/'report').exists())
    def test_no_mode_rejected(self):
        self.assertNotEqual(self.invoke('run','--run-dir',str(SOURCE)).returncode,0)
    def test_key_in_environment_does_not_enable_live(self):
        env=os.environ.copy();env['ALPHAGENOME_API_KEY']='DUMMY_NOT_REAL';env['ALPHA_GENOME_API_KEY']='DUMMY_NOT_REAL'
        r=self.invoke('run','--mode','live','--output-dir',str(self.root/'report'),env=env)
        self.assertNotEqual(r.returncode,0);self.assertFalse((self.root/'report').exists())
    def test_working_directory_independent(self):
        target=self.root/'report'
        r=self.invoke('run','--mode','offline','--run-dir',str(SOURCE),'--output-dir',str(target),cwd=self.root)
        self.assertEqual(r.returncode,0,r.stderr)
    def test_no_api_key_needed(self):
        env=os.environ.copy();env.pop('ALPHAGENOME_API_KEY',None);env.pop('ALPHA_GENOME_API_KEY',None)
        r=self.invoke('run','--mode','offline','--run-dir',str(SOURCE),'--output-dir',str(self.root/'report'),env=env)
        self.assertEqual(r.returncode,0,r.stderr)


if __name__=='__main__':unittest.main(verbosity=2)
