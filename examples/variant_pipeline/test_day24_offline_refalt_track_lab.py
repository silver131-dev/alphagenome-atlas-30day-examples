"""Day24 offline fixture/shape/delta/provenance/HTML/CLI regression tests."""
from __future__ import annotations
import ast
import copy
import csv
import io
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
import xml.etree.ElementTree as ET

import day24_offline_refalt_track_lab as lab

ROOT = Path(__file__).resolve().parent
TRACK = ROOT / "day24_synthetic_tracks.json"
PACKET = ROOT / "day23_annotation_report.json"


class Day24TrackContractTests(unittest.TestCase):
    def setUp(self):
        self.track = lab.load_json(TRACK)
        self.packet = lab.load_json(PACKET)

    def test_default_track_contract(self):
        self.assertIs(lab.check_track_fixture(self.track), self.track)

    def test_baseline_point_count_and_variant(self):
        r = lab.compare_tracks(self.track)
        self.assertEqual(r['point_count'], 5)
        self.assertEqual(r['variant_id'], 'chr1:110:A>G')

    def test_exact_delta_sequence(self):
        r = lab.compare_tracks(self.track)
        self.assertEqual([row['delta_alt_minus_ref'] for row in r['point_rows']], [0.0, 0.25, -0.3, 0.1, 0.0])

    def test_max_absolute_delta_and_coordinate(self):
        r = lab.compare_tracks(self.track)
        self.assertEqual(r['largest_absolute_delta_position_1based'], 110)
        self.assertEqual(r['largest_absolute_delta'], -0.3)

    def test_variant_here_marker(self):
        r = lab.compare_tracks(self.track)
        self.assertEqual([x['position_1based'] for x in r['point_rows'] if x['variant_here']], [110])

    def test_scientific_boundary_flags(self):
        r = lab.compare_tracks(self.track)
        self.assertIs(r['model_prediction'], False)
        self.assertIs(r['official_atlas_data'], False)
        self.assertIs(r['clinical_interpretation'], False)
        self.assertIs(r['reference_allele_verified'], False)
        self.assertEqual(r['real_api_calls_this_run'], 0)
        self.assertEqual(r['mechanism_status'], 'HYPOTHESIS_ONLY_NOT_VALIDATED')

    def test_ref_equals_alt_gives_zero_delta(self):
        d = copy.deepcopy(self.track); d['alt'] = d['ref'][:]
        r = lab.compare_tracks(d)
        self.assertEqual([a['delta_alt_minus_ref'] for a in r['point_rows']], [0] * 5)
        self.assertEqual(r['largest_absolute_delta'], 0)

    def test_no_mutation_of_fixture(self):
        old = copy.deepcopy(self.track)
        lab.compare_tracks(self.track, day23_packet=self.packet)
        self.assertEqual(self.track, old)

    def test_wrong_length_rejected(self):
        d=copy.deepcopy(self.track); d['alt'].pop()
        with self.assertRaises(ValueError): lab.check_track_fixture(d)

    def test_too_few_points_rejected(self):
        d=copy.deepcopy(self.track)
        for name in ('ref', 'alt', 'positions_1based'): d[name] = d[name][:1]
        with self.assertRaises(ValueError): lab.check_track_fixture(d)

    def test_too_many_points_rejected(self):
        d=copy.deepcopy(self.track)
        d['positions_1based']=list(range(1,1003));d['ref']=[0.5]*1002;d['alt']=[0.5]*1002
        with self.assertRaises(ValueError): lab.check_track_fixture(d)

    def test_nonconsecutive_positions_rejected(self):
        d=copy.deepcopy(self.track);d['positions_1based'][3]=114
        with self.assertRaises(ValueError): lab.check_track_fixture(d)

    def test_noninteger_position_rejected(self):
        d=copy.deepcopy(self.track);d['positions_1based'][2]=110.0
        with self.assertRaises(ValueError): lab.check_track_fixture(d)

    def test_negative_position_rejected(self):
        d=copy.deepcopy(self.track);d['positions_1based']=[-1,0,1,2,3]
        with self.assertRaises(ValueError): lab.check_track_fixture(d)

    def test_variant_outside_window_rejected(self):
        d=copy.deepcopy(self.track);d['variant_id']='chr1:200:A>G'
        with self.assertRaises(ValueError): lab.check_track_fixture(d)

    def test_same_ref_alt_rejected(self):
        d=copy.deepcopy(self.track);d['variant_id']='chr1:110:A>A'
        with self.assertRaises(ValueError): lab.check_track_fixture(d)

    def test_bool_numeric_rejected(self):
        d=copy.deepcopy(self.track);d['ref'][1]=True
        with self.assertRaises(ValueError): lab.check_track_fixture(d)

    def test_out_of_range_rejected(self):
        d=copy.deepcopy(self.track);d['ref'][0]=1.1
        with self.assertRaises(ValueError): lab.check_track_fixture(d)

    def test_nan_inf_rejected(self):
        for val in (math.nan, math.inf, -math.inf):
            d=copy.deepcopy(self.track);d['alt'][2]=val
            with self.subTest(val=val), self.assertRaises(ValueError):lab.check_track_fixture(d)

    def test_synthetic_only_assembly(self):
        d=copy.deepcopy(self.track);d['assembly']='GRCh38'
        with self.assertRaises(ValueError): lab.check_track_fixture(d)

    def test_biosample_must_be_toy(self):
        d=copy.deepcopy(self.track);d['biosample']='human_neuron'
        with self.assertRaises(ValueError): lab.check_track_fixture(d)

    def test_model_prediction_true_rejected(self):
        d=copy.deepcopy(self.track);d['model_prediction']=True
        with self.assertRaises(ValueError): lab.check_track_fixture(d)

    def test_fake_atlas_true_rejected(self):
        d=copy.deepcopy(self.track);d['official_atlas_data']=True
        with self.assertRaises(ValueError): lab.check_track_fixture(d)

    def test_real_api_count_cannot_be_nonzero(self):
        d=copy.deepcopy(self.track);d['real_api_calls_this_run']=1
        with self.assertRaises(ValueError): lab.check_track_fixture(d)

    def test_wrong_schema_and_type_rejected(self):
        d=copy.deepcopy(self.track);d['schema_version']='not-approved'
        with self.assertRaises(ValueError): lab.check_track_fixture(d)
        with self.assertRaises(ValueError): lab.check_track_fixture([])

    def test_day23_link_correct(self):
        r=lab.compare_tracks(self.track,day23_packet=self.packet)
        self.assertEqual(r['day23_annotation_link']['synthetic_score'],-0.17)
        self.assertEqual(r['day23_annotation_link']['review_status'],'REVIEW_REQUIRED')
        self.assertEqual(r['day23_annotation_link']['score_status'],'SYNTHETIC_SCORE')
        self.assertIn('DO NOT combine',r['day23_annotation_link']['note'])

    def test_day23_absent_is_optional(self):
        r=lab.compare_tracks(self.track)
        self.assertIsNone(r['day23_annotation_link'])
        self.assertIn('No Day23',lab.render_html(r,lab.render_svg(r)))

    def test_day23_wrong_variant_rejected(self):
        d=copy.deepcopy(self.packet)
        d['annotations']=[row for row in d['annotations'] if row['variant_id'] != self.track['variant_id']]
        d['selected_count']=len(d['annotations'])
        with self.assertRaises(ValueError): lab.check_day23_packet(d,self.track['variant_id'])

    def test_day23_live_claim_rejected(self):
        d=copy.deepcopy(self.packet);d['mode']='LIVE'
        with self.assertRaises(ValueError): lab.check_day23_packet(d,self.track['variant_id'])

    def test_day23_false_atlas_claim_rejected(self):
        d=copy.deepcopy(self.packet);d['official_avi']=True
        with self.assertRaises(ValueError): lab.check_day23_packet(d,self.track['variant_id'])

    def test_day23_unreviewed_selected_rejected(self):
        d=copy.deepcopy(self.packet);d['annotations'][0]['review_status']='APPROVED'
        with self.assertRaises(ValueError): lab.check_day23_packet(d,self.track['variant_id'])

    def test_day23_duplicate_rejected(self):
        d=copy.deepcopy(self.packet);d['annotations'].append(copy.deepcopy(d['annotations'][0]));d['selected_count']+=1
        with self.assertRaises(ValueError): lab.check_day23_packet(d,self.track['variant_id'])

    def test_day23_missing_score_contract(self):
        d=copy.deepcopy(self.packet);d['annotations'][2]['synthetic_score']=0.0
        with self.assertRaises(ValueError): lab.check_day23_packet(d,self.track['variant_id'])

    def test_day23_count_mismatch_rejected(self):
        d=copy.deepcopy(self.packet);d['selected_count']=4
        with self.assertRaises(ValueError): lab.check_day23_packet(d,self.track['variant_id'])

    def test_day23_invalid_score_rejected(self):
        d=copy.deepcopy(self.packet);d['annotations'][1]['synthetic_score']=float('nan')
        with self.assertRaises(ValueError): lab.check_day23_packet(d,self.track['variant_id'])


class Day24OutputsTests(unittest.TestCase):
    def setUp(self):
        self.track=lab.load_json(TRACK)
        self.packet=lab.load_json(PACKET)
        self.report=lab.compare_tracks(self.track,day23_packet=self.packet)

    def test_svg_xml_parse_and_two_traces(self):
        tree=ET.fromstring(lab.render_svg(self.report))
        ns={'s':lab.SVG_NS}
        self.assertEqual(len(tree.findall('.//s:polyline',ns)),2)
        tracks=[p.attrib.get('data-track') for p in tree.findall('.//s:polyline',ns)]
        self.assertEqual(tracks,['REF','ALT'])

    def test_svg_16_by_9(self):
        tree=ET.fromstring(lab.render_svg(self.report))
        self.assertEqual(tree.attrib.get('viewBox'),'0 0 960 540')

    def test_svg_has_synthetic_label(self):
        self.assertIn('SYNTHETIC EDUCATIONAL DATA',lab.render_svg(self.report))

    def test_svg_renders_delta_lines(self):
        tree=ET.fromstring(lab.render_svg(self.report));ns={'s':lab.SVG_NS}
        lines=tree.findall('.//s:line',ns)
        self.assertEqual(len([x for x in lines if 'data-delta-pos' in x.attrib]),5)

    def test_html_independently_contains_svg(self):
        html=lab.render_html(self.report,lab.render_svg(self.report))
        self.assertIn('<svg',html)
        self.assertIn('TOY_NOT_GRCH38',html)
        self.assertIn('Day23 SELECTED',html)
        self.assertIn('HYPOTHESIS',html.upper())

    def test_html_escapes_extra_question(self):
        r=copy.deepcopy(self.report)
        r['mechanism_questions']=['<script>alert("x")</script>']
        html=lab.render_html(r,lab.render_svg(r))
        self.assertIn('&lt;script&gt;',html)
        self.assertNotIn('<script>',html)

    def test_content_hash_matches_file(self):
        import hashlib
        x=lab.source_record(TRACK)
        self.assertEqual(x['sha256'],hashlib.sha256(TRACK.read_bytes()).hexdigest())

    def test_missing_input_fails_without_writing(self):
        with tempfile.TemporaryDirectory() as tmp:
            output=Path(tmp)/'out'
            with self.assertRaises(ValueError):lab.emit_artifacts(Path(tmp)/'missing.json',PACKET,output)
            self.assertFalse(output.exists())

    def test_invalid_json_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'bad.json';p.write_text('{bad',encoding='utf-8')
            with self.assertRaises(ValueError):lab.load_json(p)

    def test_input_limit(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'big.json';p.write_text(' '*100)
            with self.assertRaises(ValueError):lab.load_json(p,max_bytes=25)

    def test_output_files_and_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary=lab.emit_artifacts(TRACK,PACKET,Path(tmp)/'out')
            self.assertEqual(summary['status'],'PASS')
            self.assertTrue(summary['day23_linked'])
            for key in ('json_file','csv_file','svg_file','html_file'):
                self.assertTrue(Path(summary[key]).is_file())
            report=json.loads(Path(summary['json_file']).read_text())
            self.assertEqual(report['input_provenance']['track_fixture']['sha256'],lab.source_record(TRACK)['sha256'])
            self.assertEqual(report['input_provenance']['day23_packet']['sha256'],lab.source_record(PACKET)['sha256'])

    def test_output_csv_rows_and_delta(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary=lab.emit_artifacts(TRACK,PACKET,Path(tmp)/'out')
            with open(summary['csv_file'],newline='',encoding='utf-8') as stream:
                rows=list(csv.DictReader(stream))
            self.assertEqual(len(rows),5)
            self.assertEqual(rows[2]['delta_alt_minus_ref'],'-0.3')

    def test_cli_stdout_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            buff=io.StringIO()
            with redirect_stdout(buff):
                out=lab.main(['--mode','offline','--track-fixture',str(TRACK),
                              '--day23-report',str(PACKET),'--output-dir',str(Path(tmp)/'out')])
            self.assertEqual(out,json.loads(buff.getvalue()))
            self.assertEqual(out['real_api_calls_this_run'],0)

    def test_live_mode_rejected(self):
        with self.assertRaises(SystemExit) as ctx:
            with redirect_stdout(io.StringIO()):
                lab.main(['--mode','live'])
        self.assertEqual(ctx.exception.code,2)

    def test_no_sdk_network_imports(self):
        import ast
        path=ROOT/'day24_offline_refalt_track_lab.py'
        tree=ast.parse(path.read_text(encoding='utf-8'))
        allowed={'__future__','argparse','csv','hashlib','html','json','math','pathlib','re','xml'}
        for node in ast.walk(tree):
            if isinstance(node,ast.Import):
                for a in node.names:self.assertIn(a.name.split('.')[0],allowed)
            elif isinstance(node,ast.ImportFrom):
                self.assertIn((node.module or '').split('.')[0],allowed)


if __name__=='__main__':
    unittest.main(verbosity=2)
