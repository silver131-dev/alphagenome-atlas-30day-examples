"""Offline tests for Day 10's tiny, explicitly synthetic index."""
import ast
import copy
import csv
import gzip
import hashlib
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest import mock
from xml.etree import ElementTree as ET

import day10_offline_atlas_index_lab as m

HERE=Path(__file__).resolve().parent
FIXTURE=HERE/'day10_synthetic_atlas_catalog.json'
SCRIPT=HERE/'day10_offline_atlas_index_lab.py'


class CoordinateTests(unittest.TestCase):
    def test_region_valid(self):
        self.assertEqual(m.parse_region('chrToy:100-120'),('chrToy',100,120))
    def test_region_single_base(self):
        self.assertEqual(m.parse_region('chrDemo:1-1'),('chrDemo',1,1))
    def test_region_bad_chrom(self):
        with self.assertRaises(ValueError):m.parse_region('chr1:100-120')
    def test_region_bad_syntax(self):
        for value in ['chrToy:100', 'chrToy:2>3', 'chrToy:010-20', '', 'chrToy:5.5-7']:
            with self.subTest(value=value), self.assertRaises(ValueError):m.parse_region(value)
    def test_region_zero(self):
        with self.assertRaises(ValueError):m.parse_region('chrToy:0-2')
    def test_region_reverse(self):
        with self.assertRaises(ValueError):m.parse_region('chrToy:120-100')
    def test_region_limit(self):
        with self.assertRaises(ValueError):m.parse_region('chrToy:1-201')
    def test_variant_valid(self):
        self.assertEqual(m.parse_variant('chrToy:110:C>G'),('chrToy',110,'C','G'))
    def test_variant_bad_allele(self):
        with self.assertRaises(ValueError):m.parse_variant('chrToy:110:CA>G')
    def test_variant_bad_pos(self):
        with self.assertRaises(ValueError):m.parse_variant('chrToy:0:C>G')
    def test_variant_no_change(self):
        with self.assertRaises(ValueError):m.parse_variant('chrToy:110:C>C')
    def test_bucket_first_last(self):
        self.assertEqual([m.bucket_id(p) for p in [1,10,11,20,21,100,101,110,111,120]],
                         [0,0,1,1,2,9,10,10,11,11])
    def test_bucket_invalid(self):
        for p in [0, -1, True]:
            with self.subTest(p=p),self.assertRaises(ValueError):m.bucket_id(p)
    def test_bucket_width(self):
        self.assertEqual(m.bucket_id(21,20),1)


class FixtureTests(unittest.TestCase):
    def setUp(self):
        self.doc=json.loads(FIXTURE.read_text())
    def test_fixture_count(self):
        self.assertEqual(len(m.validate_catalog(self.doc)),16)
    def test_explicit_synthetic_and_assembly(self):
        self.assertEqual(self.doc['data_kind'],m.DATA_KIND)
        self.assertEqual(self.doc['assembly'],m.ASSEMBLY)
    def test_fixture_unordered_sorts(self):
        x=copy.deepcopy(self.doc)
        x['rows'].reverse()
        self.assertEqual(m.validate_catalog(x),m.validate_catalog(self.doc))
    def test_duplicate_variant(self):
        x=copy.deepcopy(self.doc);x['rows'].append(x['rows'][0])
        with self.assertRaises(ValueError):m.validate_catalog(x)
    def test_conflicting_ref_same_position(self):
        x=copy.deepcopy(self.doc);x['rows'][0]['ref']='A'
        with self.assertRaises(ValueError):m.validate_catalog(x)
    def test_forbid_real_assembly(self):
        x=copy.deepcopy(self.doc);x['assembly']='GRCh38'
        with self.assertRaises(ValueError):m.validate_catalog(x)
    def test_forbid_mislabeled_synthetic(self):
        x=copy.deepcopy(self.doc);x['data_kind']='REAL_PATIENT'
        with self.assertRaises(ValueError):m.validate_catalog(x)
    def test_forbid_extra_patient_column(self):
        x=copy.deepcopy(self.doc);x['rows'][0]['patient_id']='1234'
        with self.assertRaises(ValueError):m.validate_catalog(x)
    def test_forbid_unknown_chrom(self):
        x=copy.deepcopy(self.doc);x['rows'][0]['chrom']='chr9'
        with self.assertRaises(ValueError):m.validate_catalog(x)
    def test_forbid_boolean_position(self):
        x=copy.deepcopy(self.doc);x['rows'][0]['pos_1based']=True
        with self.assertRaises(ValueError):m.validate_catalog(x)
    def test_forbid_negative_position(self):
        x=copy.deepcopy(self.doc);x['rows'][0]['pos_1based']=-1
        with self.assertRaises(ValueError):m.validate_catalog(x)
    def test_forbid_same_ref_alt(self):
        x=copy.deepcopy(self.doc);x['rows'][0]['alt']=x['rows'][0]['ref']
        with self.assertRaises(ValueError):m.validate_catalog(x)
    def test_forbid_nan_score(self):
        x=copy.deepcopy(self.doc);x['rows'][0]['toy_proxy_score']=float('nan')
        with self.assertRaises(ValueError):m.validate_catalog(x)
    def test_forbid_infinite_score(self):
        x=copy.deepcopy(self.doc);x['rows'][0]['toy_proxy_score']=float('inf')
        with self.assertRaises(ValueError):m.validate_catalog(x)
    def test_forbid_negative_score(self):
        x=copy.deepcopy(self.doc);x['rows'][0]['toy_proxy_score']=-0.1
        with self.assertRaises(ValueError):m.validate_catalog(x)
    def test_forbid_score_above_one(self):
        x=copy.deepcopy(self.doc);x['rows'][0]['toy_proxy_score']=1.1
        with self.assertRaises(ValueError):m.validate_catalog(x)
    def test_zero_is_not_missing(self):
        x=copy.deepcopy(self.doc);x['rows'][0]['toy_proxy_score']=0.0
        self.assertEqual(next(r['toy_proxy_score'] for r in m.validate_catalog(x) if r['chrom']=='chrToy' and r['pos_1based']==95),0.0)
    def test_null_score_is_allowed(self):
        self.assertEqual(sum(r['toy_proxy_score'] is None for r in m.validate_catalog(self.doc)),1)
    def test_no_rows_rejected(self):
        x=copy.deepcopy(self.doc);x['rows']=[]
        with self.assertRaises(ValueError):m.validate_catalog(x)
    def test_load_hash(self):
        rows,digest=m.load_fixture(FIXTURE)
        self.assertEqual(len(rows),16)
        self.assertEqual(digest,hashlib.sha256(FIXTURE.read_bytes()).hexdigest())
    def test_missing_description(self):
        x=copy.deepcopy(self.doc);x['description']=''
        with self.assertRaises(ValueError):m.validate_catalog(x)


class BuiltIndexTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp=tempfile.TemporaryDirectory()
        cls.folder=Path(cls.tmp.name)/'output'
        cls.result=m.run(FIXTURE,'chrToy:100-120','chrToy:110:C>G',cls.folder)
        cls.index=json.loads((cls.folder/m.INDEX_FILE).read_text())
        cls.report=json.loads((cls.folder/m.REPORT_FILE).read_text())
    @classmethod
    def tearDownClass(cls):cls.tmp.cleanup()
    def test_cli_summary(self):self.assertEqual(self.result['status'],'PASS')
    def test_fixture_total_16(self):self.assertEqual(self.index['total_records'],16)
    def test_total_buckets_8(self):self.assertEqual(self.index['bucket_count'],8)
    def test_selected_buckets_3(self):self.assertEqual(self.report['buckets_read'],3)
    def test_selected_buckets_correct_ids(self):self.assertEqual(self.report['bucket_ids_read'],[9,10,11])
    def test_read_rows_10(self):self.assertEqual(self.report['rows_read_from_buckets'],10)
    def test_region_hits_9(self):self.assertEqual(self.report['region_hit_count'],9)
    def test_variant_exact_hit(self):self.assertEqual(self.report['variant_hits'],['chrToy:110:C>G'])
    def test_window_conversion(self):self.assertEqual(self.report['query_region_0based_half_open'],{'chrom':'chrToy','start':99,'end':120})
    def test_hit_positions_included(self):self.assertEqual([r['pos_1based'] for r in self.report['rows']],[100,101,108,109,110,110,111,119,120])
    def test_null_not_zero(self):
        rows=self.report['rows']
        self.assertEqual(sum(r['score_status']=='MISSING_SCORE' for r in rows),1)
        self.assertIsNone(next(r['toy_proxy_score'] for r in rows if r['pos_1based']==109))
    def test_no_api_and_no_clinical(self):
        self.assertFalse(self.report['atlas_query_performed'])
        self.assertFalse(self.report['clinical_interpretation'])
        self.assertEqual(self.report['real_api_calls_this_run'],0)
    def test_two_variants_same_position(self):
        rows=[r for r in self.report['rows'] if r['pos_1based']==110]
        self.assertEqual({r['alt'] for r in rows},{'G','T'})
    def test_next_single_base(self):
        rpt=m.query_index(self.folder,self.index,'chrToy:121-121',None)
        self.assertEqual(rpt['region_hit_count'],1)
        self.assertEqual(rpt['rows'][0]['pos_1based'],121)
    def test_edge_120_in_and_121_out(self):
        self.assertIn(120,[r['pos_1based'] for r in self.report['rows']])
        self.assertNotIn(121,[r['pos_1based'] for r in self.report['rows']])
    def test_other_chrom_only(self):
        rpt=m.query_index(self.folder,self.index,'chrDemo:100-110',None)
        self.assertEqual(rpt['region_hit_count'],2)
        self.assertTrue(all(r['chrom']=='chrDemo' for r in rpt['rows']))
    def test_missing_chrom_window(self):
        rpt=m.query_index(self.folder,self.index,'chrDemo:111-112',None)
        self.assertEqual(rpt['region_hit_count'],0)
        self.assertEqual(rpt['buckets_read'],0)
    def test_reject_variant_outside_region(self):
        with self.assertRaises(ValueError):m.query_index(self.folder,self.index,'chrToy:100-120','chrToy:121:A>T')
    def test_reject_variant_other_chrom(self):
        with self.assertRaises(ValueError):m.query_index(self.folder,self.index,'chrToy:100-120','chrDemo:110:G>T')
    def test_exact_variant_not_found(self):
        rpt=m.query_index(self.folder,self.index,'chrToy:100-120','chrToy:102:A>G')
        self.assertEqual(rpt['variant_hit_count'],0)
    def test_query_does_not_read_unrelated_chunk(self):
        rows=m.query_index(self.folder,self.index,'chrToy:100-100',None)
        self.assertEqual(rows['buckets_read'],1)
        self.assertEqual(rows['region_hit_count'],1)
    def test_chunk_hashes(self):
        for meta in self.index['buckets']:
            with self.subTest(meta=meta['path']):
                self.assertEqual(m.sha_file(self.folder/meta['path']),meta['sha256'])
    def test_chunk_roundtrip(self):
        recovered=[r for b in self.index['buckets'] for r in m.read_bucket(self.folder,b)]
        self.assertEqual(sorted(recovered,key=lambda r:(r['chrom'],r['pos_1based'],r['ref'],r['alt'])),m.load_fixture(FIXTURE)[0])
    def test_blocks_are_gzip(self):
        meta=self.index['buckets'][0]
        self.assertTrue((self.folder/meta['path']).read_bytes().startswith(b'\x1f\x8b'))
    def test_unsafe_chunk_path_rejected(self):
        meta=dict(self.index['buckets'][0],path='../external.json.gz')
        with self.assertRaises(ValueError):m.read_bucket(self.folder,meta)
    def test_csv_structure(self):
        with (self.folder/m.CSV_FILE).open(newline='',encoding='utf-8') as f:rows=list(csv.DictReader(f))
        self.assertEqual(len(rows),9)
        self.assertEqual(rows[0]['variant_id'],'chrToy:100:A>C')
    def test_csv_missing_is_empty(self):
        with (self.folder/m.CSV_FILE).open(newline='',encoding='utf-8') as f:rows=list(csv.DictReader(f))
        row=next(r for r in rows if r['variant_id']=='chrToy:109:A>G')
        self.assertEqual(row['toy_proxy_score'],'')
        self.assertEqual(row['score_status'],'MISSING_SCORE')
    def test_html_has_disclaimer(self):
        page=(self.folder/m.HTML_FILE).read_text()
        self.assertIn('NOT THE GOOGLE INDEX',page)
        self.assertIn('not',page.lower())
    def test_html_has_query_and_9_hits(self):
        page=(self.folder/m.HTML_FILE).read_text()
        self.assertIn('chrToy:100-120',page)
        self.assertIn('9 hits',page)
    def test_svg_dimensions_palette(self):
        root=ET.fromstring((self.folder/m.SVG_FILE).read_text())
        self.assertEqual((root.attrib['width'],root.attrib['height']),('1920','1080'))
        self.assertIn('#E5DEDA',(self.folder/m.SVG_FILE).read_text())
    def test_svg_declares_not_google(self):
        self.assertIn('NOT GOOGLE INTERNALS',(self.folder/m.SVG_FILE).read_text())
    def test_manifest_has_files(self):
        manifest=json.loads((self.folder/m.MANIFEST_FILE).read_text())
        self.assertEqual(manifest['file_count'],13)
    def test_verify_accepts(self):
        self.assertEqual(m.verify(FIXTURE,self.folder)['checksums'],'MATCH')


class FailureModeTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.folder=self.root/'generated';m.run(FIXTURE,'chrToy:100-120','chrToy:110:C>G',self.folder)
    def tearDown(self):self.tmp.cleanup()
    def test_refuse_overwrite(self):
        with self.assertRaises(ValueError):m.run(FIXTURE,'chrToy:100-120',None,self.folder)
    def test_bad_region_has_no_output(self):
        new=self.root/'new'
        with self.assertRaises(ValueError):m.run(FIXTURE,'chrToy:120-100',None,new)
        self.assertFalse(new.exists())
    def test_bad_variant_has_no_output(self):
        new=self.root/'new'
        with self.assertRaises(ValueError):m.run(FIXTURE,'chrToy:100-120','chrToy:110:C>C',new)
        self.assertFalse(new.exists())
    def test_tamper_chunk_detected(self):
        index=json.loads((self.folder/m.INDEX_FILE).read_text())
        f=self.folder/index['buckets'][0]['path'];f.write_bytes(f.read_bytes()+b'bad')
        with self.assertRaises(ValueError):m.verify(FIXTURE,self.folder)
    def test_tamper_report_detected(self):
        path=self.folder/m.REPORT_FILE;report=json.loads(path.read_text());report['region_hit_count']=99
        path.write_bytes(m.canonical_json(report))
        with self.assertRaises(ValueError):m.verify(FIXTURE,self.folder)
    def test_tamper_csv_detected(self):
        with (self.folder/m.CSV_FILE).open('a') as f:f.write('forged')
        with self.assertRaises(ValueError):m.verify(FIXTURE,self.folder)
    def test_deleted_artifact_detected(self):
        (self.folder/m.SVG_FILE).unlink()
        with self.assertRaises(ValueError):m.verify(FIXTURE,self.folder)
    def test_extra_artifact_detected(self):
        (self.folder/'extra.txt').write_text('test')
        with self.assertRaises(ValueError):m.verify(FIXTURE,self.folder)
    def test_changed_fixture_detected(self):
        bad=self.root/'catalog.json';bad.write_text(FIXTURE.read_text()+'\n')
        with self.assertRaises(ValueError):m.verify(bad,self.folder)
    def test_manifest_path_traversal_detected(self):
        f=self.folder/m.MANIFEST_FILE
        doc=json.loads(f.read_text());doc['files'][0]['path']='../outside.txt';f.write_bytes(m.canonical_json(doc))
        with self.assertRaises(ValueError):m.verify(FIXTURE,self.folder)
    def test_directory_missing(self):
        with self.assertRaises(ValueError):m.verify(FIXTURE,self.root/'none')
    def test_test_fixture_has_no_patient_fields(self):
        catalog=json.loads(FIXTURE.read_text())
        self.assertTrue(all(not {'patient_id','name','age','phenotype'} & set(r) for r in catalog['rows']))


class CommandTests(unittest.TestCase):
    def setUp(self):self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
    def tearDown(self):self.tmp.cleanup()
    def command(self,*args,env=None):
        return subprocess.run([sys.executable,'-B',str(SCRIPT),*map(str,args)],capture_output=True,text=True,env=env)
    def test_run_exit_and_stdout(self):
        r=self.command('run','--fixture',FIXTURE,'--output-dir',self.root/'run')
        self.assertEqual(r.returncode,0,r.stderr)
        self.assertEqual(json.loads(r.stdout)['region_hit_count'],9)
    def test_verify_exit_and_stdout(self):
        self.command('run','--fixture',FIXTURE,'--output-dir',self.root/'run')
        r=self.command('verify','--fixture',FIXTURE,'--output-dir',self.root/'run')
        self.assertEqual(r.returncode,0,r.stderr)
        self.assertEqual(json.loads(r.stdout)['checksums'],'MATCH')
    def test_live_rejected_even_when_key_env_set(self):
        env=dict(os.environ,ALPHAGENOME_API_KEY='DUMMY_NOT_REAL')
        r=self.command('run','--mode','live','--output-dir',self.root/'not_created',env=env)
        self.assertEqual(r.returncode,2)
        self.assertFalse((self.root/'not_created').exists())
        self.assertIn('LIVE_API_NOT_AVAILABLE',r.stderr)
    def test_no_network_client_imports(self):
        tree=ast.parse(SCRIPT.read_text(encoding='utf-8'))
        names=[]
        for n in ast.walk(tree):
            if isinstance(n,ast.Import): names.extend(x.name.split('.')[0] for x in n.names)
            if isinstance(n,ast.ImportFrom) and n.module:names.append(n.module.split('.')[0])
        self.assertFalse({'socket','urllib','requests','http','alphagenome','subprocess'} & set(names))
    def test_allow_no_api_key(self):
        env=dict(os.environ)
        env.pop('ALPHAGENOME_API_KEY',None);env.pop('ALPHA_GENOME_API_KEY',None)
        r=self.command('run','--fixture',FIXTURE,'--output-dir',self.root/'run',env=env)
        self.assertEqual(r.returncode,0,r.stderr)
    def test_run_from_other_workdir(self):
        r=subprocess.run([sys.executable,'-B',str(SCRIPT),'run','--fixture',str(FIXTURE),'--output-dir',str(self.root/'output')],capture_output=True,text=True,cwd=self.root)
        self.assertEqual(r.returncode,0,r.stderr)
    def test_script_error_is_machine_readable(self):
        r=self.command('run','--fixture',self.root/'missing','--output-dir',self.root/'new')
        self.assertEqual(r.returncode,2)
        self.assertEqual(json.loads(r.stderr)['status'],'FAIL')
    def test_two_clean_runs_same_sha(self):
        self.command('run','--fixture',FIXTURE,'--output-dir',self.root/'a')
        self.command('run','--fixture',FIXTURE,'--output-dir',self.root/'b')
        docs=[json.loads((self.root/part/m.MANIFEST_FILE).read_text()) for part in ('a','b')]
        self.assertEqual(docs[0],docs[1])
        self.assertEqual(m.sha_file(self.root/'a'/m.MANIFEST_FILE),m.sha_file(self.root/'b'/m.MANIFEST_FILE))


if __name__=='__main__':unittest.main()
