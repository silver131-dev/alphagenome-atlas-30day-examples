"""Day 04 isolated, no-network tests. Python 3.12+ standard library only."""
import copy
import json
import os
import socket
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock
from xml.etree import ElementTree as ET

import day04_offline_mvp_bridge_lab as lab

ROOT = Path(__file__).resolve().parent
FIXTURE = ROOT / 'day04_synthetic_mvp.json'
PYTHON = ROOT / 'day04_offline_mvp_bridge_lab.py'


class FixtureBase(unittest.TestCase):
    def setUp(self):
        self.data = lab.read_spec(FIXTURE)

    def changed(self, update):
        data=copy.deepcopy(self.data)
        update(data)
        return data

    def bad(self, update):
        with self.assertRaises(ValueError):
            lab.validate_spec(self.changed(update))


class VariantCoordinateTests(FixtureBase):
    def test_variant_valid(self): self.assertEqual(lab.parse_variant('chrToy:110:A>G'),(110,'A','G'))
    def test_variant_pos_one(self): self.assertEqual(lab.parse_variant('chrToy:1:C>T'),(1,'C','T'))
    def test_variant_wrong_chromosome(self):
        with self.assertRaises(ValueError):lab.parse_variant('chr1:110:A>G')
    def test_variant_same_alleles(self):
        with self.assertRaises(ValueError):lab.parse_variant('chrToy:110:A>A')
    def test_variant_multibase(self):
        with self.assertRaises(ValueError):lab.parse_variant('chrToy:110:AC>G')
    def test_variant_lowercase(self):
        with self.assertRaises(ValueError):lab.parse_variant('chrToy:110:a>g')
    def test_variant_zero(self):
        with self.assertRaises(ValueError):lab.parse_variant('chrToy:0:A>G')
    def test_variant_negative(self):
        with self.assertRaises(ValueError):lab.parse_variant('chrToy:-1:A>G')
    def test_variant_spaces(self):
        with self.assertRaises(ValueError):lab.parse_variant(' chrToy:110:A>G ')
    def test_variant_wrong_type(self):
        with self.assertRaises(ValueError):lab.parse_variant(110)
    def test_region_valid(self): self.assertEqual(lab.parse_region('chrToy:100-120'),(100,120))
    def test_region_single_base(self):self.assertEqual(lab.parse_region('chrToy:110-110'),(110,110))
    def test_region_reversed(self):
        with self.assertRaises(ValueError):lab.parse_region('chrToy:120-100')
    def test_region_too_large(self):
        with self.assertRaises(ValueError):lab.parse_region('chrToy:1-100002')
    def test_region_zero(self):
        with self.assertRaises(ValueError):lab.parse_region('chrToy:0-120')
    def test_region_nonstring(self):
        with self.assertRaises(ValueError):lab.parse_region(None)
    def test_region_incorrect_chromosome(self):
        with self.assertRaises(ValueError):lab.parse_region('chr9:100-120')
    def test_region_missing_end(self):
        with self.assertRaises(ValueError):lab.parse_region('chrToy:100-')


class SchemaTests(FixtureBase):
    def test_valid_source(self):lab.validate_spec(self.data)
    def test_not_mapping(self):
        with self.assertRaises(ValueError):lab.validate_spec([])
    def test_unexpected_top_level(self):self.bad(lambda d:d.update({'patient_name':'x'}))
    def test_missing_top_level(self):self.bad(lambda d:d.pop('assembly'))
    def test_wrong_evidence(self):self.bad(lambda d:d.update({'evidence_class':'REAL_MODEL'}))
    def test_wrong_assembly(self):self.bad(lambda d:d.update({'assembly':'GRCh38'}))
    def test_wrong_schema(self):self.bad(lambda d:d.update({'schema_version':'V2'}))
    def test_actual_sample_label_forbidden(self):self.bad(lambda d:d.update({'biosample_label':'patient'}))
    def test_fake_avi_rejected(self):self.bad(lambda d:d['variant_view'].update({'avi_phred':32.2}))
    def test_fake_score_status_rejected(self):self.bad(lambda d:d['variant_view'].update({'score_status':'OFFICIAL'}))
    def test_fake_importance_rejected(self):self.bad(lambda d:d['variant_view'].update({'feature_importance_available':True}))
    def test_extra_view_field(self):self.bad(lambda d:d['variant_view'].update({'clinical_risk':0.9}))
    def test_variant_outside_region(self):self.bad(lambda d:d.update({'variant_query':'chrToy:200:A>G'}))
    def test_variant_included_in_positions(self):self.bad(lambda d:d['locus_view'].update({'positions_1based':[100,101,102,103,104]}))
    def test_duplicate_positions(self):self.bad(lambda d:d['locus_view'].update({'positions_1based':[100,105,105,115,120]}))
    def test_unsorted_positions(self):self.bad(lambda d:d['locus_view'].update({'positions_1based':[100,110,105,115,120]}))
    def test_boolean_position(self):self.bad(lambda d:d['locus_view'].update({'positions_1based':[True,105,110,115,120]}))
    def test_position_outside(self):self.bad(lambda d:d['locus_view'].update({'positions_1based':[99,105,110,115,120]}))
    def test_missing_track(self):self.bad(lambda d:d['locus_view'].update({'tracks':d['locus_view']['tracks'][:1]}))
    def test_track_wrong_order(self):self.bad(lambda d:d['locus_view']['tracks'].reverse())
    def test_track_wrong_id(self):self.bad(lambda d:d['locus_view']['tracks'][0].update({'track_id':'AVI_SCORE'}))
    def test_track_bad_length(self):self.bad(lambda d:d['locus_view']['tracks'][0].update({'values':[0.5]}))
    def test_track_nan(self):self.bad(lambda d:d['locus_view']['tracks'][0]['values'].__setitem__(0,float('nan')))
    def test_track_inf(self):self.bad(lambda d:d['locus_view']['tracks'][0]['values'].__setitem__(0,float('inf')))
    def test_track_negative(self):self.bad(lambda d:d['locus_view']['tracks'][0]['values'].__setitem__(0,-0.001))
    def test_track_over_one(self):self.bad(lambda d:d['locus_view']['tracks'][0]['values'].__setitem__(0,1.1))
    def test_track_boolean(self):self.bad(lambda d:d['locus_view']['tracks'][0]['values'].__setitem__(0,False))
    def test_missing_value_allowed(self):self.assertIsNone(self.data['locus_view']['tracks'][1]['values'][3]);lab.validate_spec(self.data)
    def test_zero_value_allowed(self):self.assertEqual(self.data['locus_view']['tracks'][1]['values'][2],0.0);lab.validate_spec(self.data)
    def test_region_edge_check(self):
        d=self.changed(lambda d:d.update({'region_query':'chrToy:110-110'}))
        with self.assertRaises(ValueError):lab.validate_spec(d)


class ArtifactTests(FixtureBase):
    def test_csv_retains_zero(self):
        self.assertIn('110,0.42,0\n',lab.locus_csv(self.data))
    def test_csv_retains_missing(self):
        self.assertIn('115,0.04,\n',lab.locus_csv(self.data))
    def test_contract_variant_null(self):
        v=lab.build_contract(self.data,'a'*64)['variant_view']
        self.assertIsNone(v['avi_phred']);self.assertEqual(v['score_status'],'NOT_QUERIED')
    def test_contract_region_conversion(self):
        r=lab.build_contract(self.data,'a'*64)['locus_view']['region_0based_half_open']
        self.assertEqual(r,{'start':99,'end':120})
    def test_contract_flags(self):
        c=lab.build_contract(self.data,'a'*64)
        self.assertIs(c['real_model_predictions'],False)
        self.assertIs(c['external_transfer_allowed'],False)
        self.assertEqual(c['real_api_calls_this_run'],0)
    def test_svg_is_valid_xml(self):
        root=ET.fromstring(lab.build_svg(self.data));self.assertEqual(root.attrib.get('width'),'1920')
    def test_svg_is_169(self):
        root=ET.fromstring(lab.build_svg(self.data));self.assertEqual(root.attrib.get('viewBox'),'0 0 1920 1080')
    def test_svg_has_both_views(self):
        svg=lab.build_svg(self.data);self.assertIn('Variant View',svg);self.assertIn('Locus View',svg)
    def test_svg_contains_synthetic_warning(self):self.assertIn('SYNTHETIC EDUCATIONAL DEMO',lab.build_svg(self.data))
    def test_html_switches_views(self):
        h=lab.build_html(self.data,{});self.assertIn('aria-controls="view-locus"',h);self.assertIn('addEventListener',h)
    def test_html_disallows_external_network(self):
        h=lab.build_html(self.data,{})
        self.assertIn("default-src 'none'",h)
        self.assertNotIn('fetch(',h);self.assertNotIn('XMLHttpRequest',h)
    def test_html_not_fake_ai_studio(self):self.assertIn('不是 Google AI Studio 實際匯出結果',lab.build_html(self.data,{}))
    def test_mapping_mentions_node_vs_python(self):
        m=lab.reverse_mapping_md(self.data);self.assertIn('React',m);self.assertIn('Python',m);self.assertIn('Node.js',m)
    def test_mapping_no_api(self):self.assertIn('real_api_calls_this_run=0',lab.reverse_mapping_md(self.data))


class RunVerifyTests(FixtureBase):
    def setUp(self):
        super().setUp()
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.source=self.root/'day04_synthetic_mvp.json'
        self.source.write_bytes(FIXTURE.read_bytes())
        self.out=self.root/'out'

    def run_first(self):return lab.run_lab(self.source,self.out)

    def test_run_status(self):self.assertEqual(self.run_first()['status'],'PASS')
    def test_run_artifact_count(self):self.assertEqual(self.run_first()['artifact_count'],5)
    def test_run_has_six_files_including_manifest(self):
        self.run_first();self.assertEqual(len(list(self.out.iterdir())),6)
    def test_verify(self):
        self.run_first();self.assertEqual(lab.verify_lab(self.source,self.out)['checksums'],'MATCH')
    def test_source_sha(self):
        self.run_first();m=json.loads((self.out/'day04_manifest.json').read_text());self.assertEqual(m['source_fixture_sha256'],lab.sha256_file(self.source))
    def test_no_overwrite(self):
        self.run_first()
        with self.assertRaises(FileExistsError):self.run_first()
    def test_fixture_changed(self):
        self.run_first();self.source.write_bytes(self.source.read_bytes()+b' ')
        with self.assertRaises(ValueError):lab.verify_lab(self.source,self.out)
    def test_artifact_tamper(self):
        self.run_first();(self.out/'day04_mvp_dashboard.html').write_text('tampered')
        with self.assertRaises(ValueError):lab.verify_lab(self.source,self.out)
    def test_csv_tamper(self):
        self.run_first();(self.out/'day04_locus_tracks.csv').write_text('fake')
        with self.assertRaises(ValueError):lab.verify_lab(self.source,self.out)
    def test_manifest_fake_file(self):
        self.run_first();p=self.out/'day04_manifest.json';m=json.loads(p.read_text());m['files']['../../leak']='x';p.write_text(json.dumps(m))
        with self.assertRaises(ValueError):lab.verify_lab(self.source,self.out)
    def test_manifest_wrong_mode(self):
        self.run_first();p=self.out/'day04_manifest.json';m=json.loads(p.read_text());m['mode']='LIVE';p.write_text(json.dumps(m))
        with self.assertRaises(ValueError):lab.verify_lab(self.source,self.out)
    def test_manifest_deleted_artifact(self):
        self.run_first();(self.out/'day04_mvp_two_views.svg').unlink()
        with self.assertRaises(FileNotFoundError):lab.verify_lab(self.source,self.out)
    def test_contract_tamper_with_matching_manifest_hash(self):
        self.run_first();p=self.out/'day04_mvp_contract.json';d=json.loads(p.read_text());d['real_api_calls_this_run']=5;p.write_text(json.dumps(d))
        mfile=self.out/'day04_manifest.json';m=json.loads(mfile.read_text());m['files']['day04_mvp_contract.json']=lab.sha256_file(p);mfile.write_text(json.dumps(m))
        with self.assertRaises(ValueError):lab.verify_lab(self.source,self.out)
    def test_csv_tamper_with_matching_manifest_hash(self):
        self.run_first();p=self.out/'day04_locus_tracks.csv';p.write_text('position_1based,splicing_proxy,accessibility_proxy\n110,0.42,0.9\n')
        mfile=self.out/'day04_manifest.json';m=json.loads(mfile.read_text());m['files']['day04_locus_tracks.csv']=lab.sha256_file(p);mfile.write_text(json.dumps(m))
        with self.assertRaises(ValueError):lab.verify_lab(self.source,self.out)
    def test_outputs_reproducible(self):
        self.run_first();other=self.root/'out2';lab.run_lab(self.source,other)
        for name in [*lab.ARTIFACTS,'day04_manifest.json']:
            with self.subTest(name=name):self.assertEqual((self.out/name).read_bytes(),(other/name).read_bytes())
    def test_no_network_io_in_lab(self):
        with mock.patch.object(socket,'create_connection',side_effect=AssertionError('Network call')):
            self.run_first();lab.verify_lab(self.source,self.out)
    def test_cli_run(self):
        p=subprocess.run([sys.executable,str(PYTHON),'run','--mode','offline','--fixture',str(self.source),'--output-dir',str(self.out)],capture_output=True,text=True)
        self.assertEqual(p.returncode,0,p.stderr)
        self.assertEqual(json.loads(p.stdout)['status'],'PASS')
    def test_cli_verify(self):
        self.run_first();p=subprocess.run([sys.executable,str(PYTHON),'verify','--mode','offline','--fixture',str(self.source),'--output-dir',str(self.out)],capture_output=True,text=True)
        self.assertEqual(p.returncode,0,p.stderr)
        self.assertEqual(json.loads(p.stdout)['files_verified'],5)
    def test_live_mode_rejected(self):
        p=subprocess.run([sys.executable,str(PYTHON),'run','--mode','live','--fixture',str(self.source),'--output-dir',str(self.out)],capture_output=True,text=True)
        self.assertNotEqual(p.returncode,0)
        self.assertFalse(self.out.exists())
    def test_real_key_presence_not_used(self):
        with mock.patch.dict(os.environ,{'ALPHA_GENOME_API_KEY':'FAKE_DO_NOT_USE','ALPHAGENOME_API_KEY':'FAKE_DO_NOT_USE'}):
            res=self.run_first();self.assertEqual(res['real_api_calls_this_run'],0)
            for f in self.out.iterdir():
                self.assertNotIn('FAKE_DO_NOT_USE',f.read_text(encoding='utf-8'))


if __name__=='__main__':unittest.main()
