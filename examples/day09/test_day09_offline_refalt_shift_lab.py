"""Behavioral tests for standalone synthetic Day09 teaching lab; no network."""
import copy
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

import day09_offline_refalt_shift_lab as lab

HERE=Path(__file__).resolve().parent
SOURCE=HERE/'day09_synthetic_refalt_shift.json'
MODULE=HERE/'day09_offline_refalt_shift_lab.py'


def fresh():
    return json.loads(SOURCE.read_text(encoding='utf-8'))


class Day09FixtureTests(unittest.TestCase):
    def test_fixture_loads(self): self.assertEqual(lab.load_fixture(SOURCE)['assembly'],lab.ASSEMBLY)
    def test_source_label(self): self.assertIn('shifted',fresh()['label'])
    def test_valid_fixture(self): self.assertIs(lab.validated_fixture(fresh()) is not None,True)
    def test_reference_at_104(self):
        f=fresh();self.assertEqual(f['reference_sequence'][f['positions_1based'].index(104)],'A')
    def test_artifact_names_unique(self): self.assertEqual(len(set(lab.ARTIFACT_NAMES)),6)
    def test_digest_matches_hashlib(self): self.assertEqual(lab.digest(b'test'),hashlib.sha256(b'test').hexdigest())
    def test_canonical_json_deterministic(self):
        self.assertEqual(lab.canonical_json({'a':1,'b':2}),lab.canonical_json({'b':2,'a':1}))
    def test_canonical_rejects_nan(self):
        with self.assertRaises(ValueError): lab.canonical_json({'a':float('nan')})


class Day09ValidationTests(unittest.TestCase):
    def expect_bad(self,fn):
        f=fresh();fn(f)
        with self.assertRaises((ValueError,TypeError)):
            lab.validated_fixture(f)
    def test_not_synthetic(self):self.expect_bad(lambda d:d.update(data_kind='REAL_DATA'))
    def test_wrong_assembly(self):self.expect_bad(lambda d:d.update(assembly='GRCh38'))
    def test_variant_wrong_chrom(self):self.expect_bad(lambda d:d.update(variant='chr9:104:A>G'))
    def test_variant_zero_position(self):self.expect_bad(lambda d:d.update(variant='chrToy:0:A>G'))
    def test_variant_ref_equals_alt(self):self.expect_bad(lambda d:d.update(variant='chrToy:104:A>A'))
    def test_variant_non_base_alt(self):self.expect_bad(lambda d:d.update(variant='chrToy:104:A>AT'))
    def test_variant_ref_mismatch(self):self.expect_bad(lambda d:d.update(variant='chrToy:104:C>G'))
    def test_variant_outside_coords(self):self.expect_bad(lambda d:d.update(variant='chrToy:999:A>G'))
    def test_variant_missing(self):self.expect_bad(lambda d:d.pop('variant'))
    def test_positions_duplicate(self):self.expect_bad(lambda d:d['positions_1based'].__setitem__(3,102))
    def test_positions_wrong_order(self):self.expect_bad(lambda d:d['positions_1based'].reverse())
    def test_positions_gap(self):self.expect_bad(lambda d:d['positions_1based'].__setitem__(5,110))
    def test_positions_negative(self):self.expect_bad(lambda d:d['positions_1based'].__setitem__(0,-1))
    def test_positions_bool(self):self.expect_bad(lambda d:d['positions_1based'].__setitem__(0,True))
    def test_positions_short(self):self.expect_bad(lambda d:d.update(positions_1based=[100]))
    def test_sequence_wrong_length(self):self.expect_bad(lambda d:d.update(reference_sequence='ACGT'))
    def test_sequence_invalid_base(self):self.expect_bad(lambda d:d.update(reference_sequence='TTGAANCGT'))
    def test_ref_short_signal(self):self.expect_bad(lambda d:d['reference_track'].pop())
    def test_alt_short_signal(self):self.expect_bad(lambda d:d['alternate_track'].pop())
    def test_ref_negative_signal(self):self.expect_bad(lambda d:d['reference_track'].__setitem__(1,-0.1))
    def test_alt_signal_gt_one(self):self.expect_bad(lambda d:d['alternate_track'].__setitem__(1,1.1))
    def test_alt_signal_nan(self):self.expect_bad(lambda d:d['alternate_track'].__setitem__(1,float('nan')))
    def test_ref_signal_inf(self):self.expect_bad(lambda d:d['reference_track'].__setitem__(1,float('inf')))
    def test_ref_signal_bool(self):self.expect_bad(lambda d:d['reference_track'].__setitem__(1,False))
    def test_ref_signal_string(self):self.expect_bad(lambda d:d['reference_track'].__setitem__(1,'0.2'))
    def test_empty_label(self):self.expect_bad(lambda d:d.update(label='   '))
    def test_shift_missing(self):self.expect_bad(lambda d:d.pop('shift_only_windows'))
    def test_shift_ref_wrong_start(self):self.expect_bad(lambda d:d['shift_only_windows']['ref_positions_1based'].__setitem__(0,101))
    def test_shift_alt_wrong_start(self):self.expect_bad(lambda d:d['shift_only_windows']['alt_positions_1based'].__setitem__(0,100))
    def test_shift_alt_corrupt_value(self):self.expect_bad(lambda d:d['shift_only_windows']['alt_values'].__setitem__(3,0.1))
    def test_shift_ref_corrupt_value(self):self.expect_bad(lambda d:d['shift_only_windows']['ref_values'].__setitem__(2,0.8))
    def test_shift_ref_short(self):self.expect_bad(lambda d:d['shift_only_windows']['ref_values'].pop())
    def test_shift_alt_short(self):self.expect_bad(lambda d:d['shift_only_windows']['alt_values'].pop())
    def test_shift_has_to_be_dict(self):self.expect_bad(lambda d:d.update(shift_only_windows='wrong'))
    def test_shift_alt_nan(self):self.expect_bad(lambda d:d['shift_only_windows']['alt_values'].__setitem__(0,float('nan')))


class Day09MathTests(unittest.TestCase):
    def setUp(self):self.report=lab.build_report(fresh(),lab.digest(SOURCE.read_bytes()))
    def test_evidence_class(self):self.assertEqual(self.report['evidence_class'],lab.KIND)
    def test_genome_build(self):self.assertEqual(self.report['assembly'],lab.ASSEMBLY)
    def test_variant(self):self.assertEqual(self.report['variant'],'chrToy:104:A>G')
    def test_effect_point_count(self):self.assertEqual(self.report['effect_comparison']['point_count'],9)
    def test_ref_alt_max_position(self):self.assertEqual(self.report['effect_comparison']['largest_absolute_delta_position_1based'],104)
    def test_ref_alt_max_delta(self):self.assertEqual(self.report['effect_comparison']['largest_absolute_delta'],-0.3)
    def test_delta_positive_on_105(self):self.assertEqual(self.report['effect_comparison']['rows'][5]['delta'],0.15)
    def test_delta_zero_on_100(self):self.assertEqual(self.report['effect_comparison']['rows'][0]['delta'],0.0)
    def test_delta_zero_on_108(self):self.assertEqual(self.report['effect_comparison']['rows'][-1]['delta'],0.0)
    def test_naive_shift_max(self):self.assertEqual(self.report['shift_only_probe']['naive_max_absolute_delta'],0.36)
    def test_shift_no_biological_change(self):self.assertIs(self.report['shift_only_probe']['biological_signal_changed'],False)
    def test_coordinate_aligned_max_zero(self):self.assertEqual(self.report['shift_only_probe']['aligned_max_absolute_delta'],0)
    def test_coordinate_aligned_shared_count(self):self.assertEqual(self.report['shift_only_probe']['coordinate_aligned_shared_position_count'],7)
    def test_no_diff_for_shared_positions(self):self.assertEqual({r['aligned_delta'] for r in self.report['shift_only_probe']['aligned_rows']},{0.0})
    def test_shift_unmatched_left(self):self.assertEqual(self.report['shift_only_probe']['unmatched_ref_positions'],[100])
    def test_shift_unmatched_right(self):self.assertEqual(self.report['shift_only_probe']['unmatched_alt_positions'],[108])
    def test_naive_compares_mismatched_coordinates(self):
        x=self.report['shift_only_probe']['naive_rows'][0];self.assertNotEqual(x['ref_position_1based'],x['alt_position_1based'])
    def test_no_model_prediction(self):self.assertIs(self.report['model_prediction'],False)
    def test_no_reference_verification(self):self.assertIs(self.report['genome_reference_verified'],False)
    def test_no_clinical_interpretation(self):self.assertIs(self.report['clinical_interpretation'],False)
    def test_api_calls_zero(self):self.assertEqual(self.report['real_api_calls_this_run'],0)
    def test_fixture_sha(self):self.assertEqual(self.report['input_provenance']['fixture_sha256'],lab.digest(SOURCE.read_bytes()))
    def test_limits_include_augmentation_caveat(self):
        self.assertIn('Shift Augmentation',' '.join(self.report['scientific_limits']))


class Day09PresentationTests(unittest.TestCase):
    def setUp(self):
        self.r=lab.build_report(fresh(),lab.digest(SOURCE.read_bytes()))
        self.svg=lab.svg_text(self.r)
    def test_svg_xml_parses(self):self.assertTrue(ET.fromstring(self.svg).tag.endswith('svg'))
    def test_svg_16_9(self):
        root=ET.fromstring(self.svg);self.assertEqual((root.get('width'),root.get('height')),('1920','1080'))
    def test_svg_palette(self):self.assertIn('#E5DEDA',self.svg)
    def test_svg_ref_alt(self):
        self.assertIn('data-track="REF"',self.svg);self.assertIn('data-track="ALT"',self.svg)
    def test_svg_synthetic_disclaimer(self):self.assertIn('SYNTHETIC EDUCATIONAL DEMO',self.svg)
    def test_svg_shift_caveat(self):self.assertIn('Shift Augmentation',self.svg)
    def test_html_contains_table(self):self.assertIn('<table>',lab.html_text(self.r,self.svg))
    def test_html_disclaimer(self):self.assertIn('NO API',lab.html_text(self.r,self.svg))
    def test_html_source_sha(self):self.assertIn(self.r['input_provenance']['fixture_sha256'],lab.html_text(self.r,self.svg))
    def test_notes_disclaimer(self):self.assertIn('not a model window-augmentation test',lab.note_text(self.r))
    def test_csv_effect_header(self):
        csvs=lab.csv_text(('position_1based','ref_proxy','alt_proxy','delta'),self.r['effect_comparison']['rows']);self.assertEqual(len(list(csv.DictReader(io.StringIO(csvs)))),9)
    def test_csv_shift_count(self):
        csvs=lab.csv_text(('ref_position_1based','alt_position_1based','naive_delta'),self.r['shift_only_probe']['naive_rows']);self.assertEqual(len(list(csv.DictReader(io.StringIO(csvs)))),8)


class Day09RunAndVerifyTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.path=Path(self.tmp.name)/'out';self.fix=Path(self.tmp.name)/SOURCE.name
        self.fix.write_bytes(SOURCE.read_bytes())
    def tearDown(self):self.tmp.cleanup()
    def run_once(self):return lab.run(self.fix,self.path)
    def test_run_creates_manifest(self):
        self.run_once();self.assertTrue((self.path/'day09_sha256_manifest.json').is_file())
    def test_run_creates_six_artifacts(self):
        self.run_once();self.assertEqual(sum((self.path/n).is_file() for n in lab.ARTIFACT_NAMES),6)
    def test_run_stdout_contract(self):
        result=self.run_once();self.assertEqual((result['status'],result['mode'],result['real_api_calls_this_run']),('PASS','OFFLINE_ONLY',0))
    def test_verify_ok(self):
        self.run_once();self.assertEqual(lab.verify(self.fix,self.path)['checksums'],'MATCH')
    def test_reject_existing_output(self):
        self.run_once()
        with self.assertRaises(FileExistsError):self.run_once()
    def test_reject_file_tamper(self):
        self.run_once();p=self.path/'day09_effect_delta.csv';p.write_text('malicious!',encoding='utf8')
        with self.assertRaises(ValueError):lab.verify(self.fix,self.path)
    def test_reject_manifest_tamper(self):
        self.run_once();p=self.path/'day09_sha256_manifest.json';obj=json.loads(p.read_text());obj['files'].pop('day09_ref_alt_report.json');p.write_text(json.dumps(obj))
        with self.assertRaises(ValueError):lab.verify(self.fix,self.path)
    def test_reject_source_change(self):
        self.run_once();self.fix.write_bytes(self.fix.read_bytes()+b' ')
        with self.assertRaises(ValueError):lab.verify(self.fix,self.path)
    def test_reject_source_renaming(self):
        self.run_once();new=self.fix.with_name('renamed.json');self.fix.rename(new)
        with self.assertRaises(ValueError):lab.verify(new,self.path)
    def test_reject_report_contract_even_if_hash_forged(self):
        self.run_once();p=self.path/'day09_ref_alt_report.json';obj=json.loads(p.read_text());obj['effect_comparison']['largest_absolute_delta']=0.99;p.write_bytes(lab.canonical_json(obj))
        manifest_path=self.path/'day09_sha256_manifest.json';m=json.loads(manifest_path.read_text());m['files'][p.name]=lab.digest(p.read_bytes());manifest_path.write_bytes(lab.canonical_json(m))
        with self.assertRaises(ValueError):lab.verify(self.fix,self.path)
    def test_deterministic_repeated_run(self):
        self.run_once()
        b={n:(self.path/n).read_bytes() for n in [*lab.ARTIFACT_NAMES,'day09_sha256_manifest.json']}
        second=Path(self.tmp.name)/'other';lab.run(self.fix,second)
        self.assertEqual(b,{n:(second/n).read_bytes() for n in b})
    def test_json_report_parses(self):
        self.run_once();r=json.loads((self.path/lab.ARTIFACT_NAMES[0]).read_text());self.assertEqual(r['variant'],'chrToy:104:A>G')
    def test_manifest_fingerprint_format(self):
        self.run_once();m=json.loads((self.path/'day09_sha256_manifest.json').read_text());self.assertEqual(len(m['source_fixture_sha256']),64)
    def test_cli_rejects_live(self):
        cp=subprocess.run([sys.executable,str(MODULE),'run','--mode','live','--fixture',str(self.fix),'--output-dir',str(self.path)],capture_output=True,text=True)
        self.assertNotEqual(cp.returncode,0);self.assertFalse(self.path.exists())
    def test_cli_run_and_verify(self):
        for cmd in ('run','verify'):
            p=subprocess.run([sys.executable,str(MODULE),cmd,'--fixture',str(self.fix),'--output-dir',str(self.path)],capture_output=True,text=True)
            self.assertEqual(p.returncode,0,p.stderr)
            self.assertEqual(json.loads(p.stdout)['real_api_calls_this_run'],0)
    def test_api_key_env_does_not_enable_live(self):
        env=os.environ.copy();env['ALPHAGENOME_API_KEY']='FAKE_DO_NOT_USE'
        p=subprocess.run([sys.executable,str(MODULE),'run','--mode','live','--fixture',str(self.fix),'--output-dir',str(self.path)],env=env,capture_output=True,text=True)
        self.assertNotEqual(p.returncode,0)
        self.assertFalse(self.path.exists())
    def test_only_stdlib_imports(self):
        import ast
        tree=ast.parse(MODULE.read_text(encoding='utf-8'))
        imports={alias.name.split('.')[0] for x in ast.walk(tree) if isinstance(x,ast.Import) for alias in x.names}
        imports|={x.module.split('.')[0] for x in ast.walk(tree) if isinstance(x,ast.ImportFrom) and x.module}
        self.assertFalse(imports&{'requests','httpx','socket','urllib','alphagenome','google','pandas','numpy'})

if __name__=='__main__':unittest.main()
