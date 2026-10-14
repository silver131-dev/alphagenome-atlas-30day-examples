"""Independent unit tests for Day15 educational-only saturation design."""
from __future__ import annotations

import copy
import csv
import io
import json
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

import day15_offline_saturation_design_lab as d

ROOT=Path(__file__).resolve().parent
FIXTURE=ROOT/'day15_synthetic_design.json'


def spec():
    return json.loads(FIXTURE.read_text(encoding='utf-8'))


class SpecValidationTests(unittest.TestCase):
    def test_valid_fixture(self):
        self.assertEqual(d.validate_spec(spec())['reference_sequence'],'ACGT')

    def test_not_a_mapping(self):
        with self.assertRaises(ValueError): d.validate_spec([])

    def test_missing_field(self):
        s=spec(); del s['assembly']
        with self.assertRaises(ValueError): d.validate_spec(s)

    def test_rejects_patient_field(self):
        s=spec(); s['patient_id']='NOT_ALLOWED'
        with self.assertRaises(ValueError): d.validate_spec(s)

    def test_rejects_non_synthetic_marker(self):
        s=spec(); s['data_kind']='REAL_PATIENT'
        with self.assertRaises(ValueError): d.validate_spec(s)

    def test_rejects_grch38(self):
        s=spec(); s['assembly']='GRCh38'
        with self.assertRaises(ValueError): d.validate_spec(s)

    def test_rejects_real_chr(self):
        s=spec(); s['chromosome']='chr1'
        with self.assertRaises(ValueError): d.validate_spec(s)

    def test_rejects_true_experiment_flag(self):
        s=spec(); s['experiment_performed']=True
        with self.assertRaises(ValueError): d.validate_spec(s)

    def test_rejects_clinical_flag(self):
        s=spec(); s['clinical_interpretation']=True
        with self.assertRaises(ValueError): d.validate_spec(s)

    def test_rejects_bad_schema(self):
        s=spec(); s['schema_version']='DAY15_LIVE_V1'
        with self.assertRaises(ValueError): d.validate_spec(s)

    def test_rejects_non_acgt_sequence(self):
        s=spec(); s['reference_sequence']='ACNT'
        with self.assertRaises(ValueError): d.validate_spec(s)

    def test_rejects_lowercase_sequence(self):
        s=spec(); s['reference_sequence']='acgt'
        with self.assertRaises(ValueError): d.validate_spec(s)

    def test_rejects_empty_sequence(self):
        s=spec(); s['reference_sequence']=''
        with self.assertRaises(ValueError): d.validate_spec(s)

    def test_rejects_overlong_sequence(self):
        s=spec(); s['reference_sequence']='A'*13
        with self.assertRaises(ValueError): d.validate_spec(s)

    def test_rejects_bad_position(self):
        for p in (0, -1, 1000001, 3.5, True, '101'):
            with self.subTest(p=p), self.assertRaises(ValueError):
                s=spec(); s['start_1based']=p; d.validate_spec(s)

    def test_rejects_bad_repeats(self):
        for repeats in (0,7,-1,1.5,True):
            with self.subTest(repeats=repeats), self.assertRaises(ValueError):
                s=spec(); s['planned_repeats']=repeats; d.validate_spec(s)

    def test_rejects_bad_seed(self):
        for seed in (-1,100000000,'0',False):
            with self.subTest(seed=seed), self.assertRaises(ValueError):
                s=spec(); s['randomization_seed']=seed; d.validate_spec(s)

    def test_read_spec_hashes_bytes(self):
        s, digest=d.read_spec(FIXTURE)
        self.assertEqual(digest,d.sha256_bytes(FIXTURE.read_bytes()))
        self.assertEqual(s['assembly'],'TOY_NOT_GRCH38')

    def test_malformed_json(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'x.json';p.write_text('{',encoding='utf8')
            with self.assertRaises(ValueError): d.read_spec(p)

    def test_missing_file(self):
        with self.assertRaises(ValueError): d.read_spec(ROOT/'not_here_123.json')


class CandidateTests(unittest.TestCase):
    def setUp(self):
        self.s=spec(); self.c=d.enumerate_snvs(self.s)

    def test_four_positions_12_snvs(self):
        self.assertEqual(len(self.c),12)

    def test_three_alternatives_per_base(self):
        from collections import Counter
        self.assertEqual(set(Counter(c['position_1based'] for c in self.c).values()),{3})

    def test_never_ref_equals_alt(self):
        self.assertTrue(all(c['ref']!=c['alt'] for c in self.c))

    def test_variants_are_unique(self):
        self.assertEqual(len({c['variant_id'] for c in self.c}),12)

    def test_all_positions_are_toy(self):
        self.assertTrue(all(c['variant_id'].startswith('chrToy:') for c in self.c))

    def test_all_measurements_null(self):
        self.assertTrue(all(c['observed_effect'] is None and c['measurement_status']=='NOT_MEASURED' for c in self.c))

    def test_explicit_model_not_run(self):
        self.assertTrue(all(c['model_prediction']=='NOT_RUN' for c in self.c))

    def test_first_position_variants(self):
        ids=[x['variant_id'] for x in self.c if x['position_1based']==101]
        self.assertEqual(ids,['chrToy:101:A>C','chrToy:101:A>G','chrToy:101:A>T'])

    def test_single_base_always_three_snvs(self):
        s=spec();s['reference_sequence']='A'
        self.assertEqual(len(d.enumerate_snvs(s)),3)


class AllocationTests(unittest.TestCase):
    def setUp(self):
        self.s=spec();self.c=d.enumerate_snvs(self.s);self.a=d.allocate_slots(self.s,self.c)

    def test_45_planned_slots(self):
        self.assertEqual(len(self.a),45)

    def test_15_per_run(self):
        for i in range(1,4):
            self.assertEqual(sum(x['planned_run']==f'R{i}' for x in self.a),15)

    def test_every_candidate_once_per_run(self):
        wanted={c['candidate_id'] for c in self.c}|set(d.CONTROLS)
        for r in ('R1','R2','R3'):
            self.assertEqual({a['candidate_id'] for a in self.a if a['planned_run']==r},wanted)

    def test_control_count_is_nine(self):
        self.assertEqual(sum(x['item_type']=='CONTROL_PLACEHOLDER' for x in self.a),9)

    def test_no_unmeasured_values_filled(self):
        self.assertTrue(all(x['observed_readout'] is None and x['measurement_status']=='NOT_MEASURED' for x in self.a))

    def test_not_biological_replicates(self):
        self.assertTrue(all(x['biological_independence_verified'] is False for x in self.a))

    def test_unique_slots(self):
        self.assertEqual(len({x['slot_id'] for x in self.a}),45)

    def test_deterministic_allocation(self):
        self.assertEqual(self.a,d.allocate_slots(self.s,self.c))

    def test_changed_seed_changes_order(self):
        s=spec();s['randomization_seed']+=1
        self.assertNotEqual(self.a,d.allocate_slots(s,d.enumerate_snvs(s)))

    def test_edited_candidate_rejected(self):
        fake=copy.deepcopy(self.c);fake[0]['variant_id']='chrToy:999:A>C'
        with self.assertRaises(ValueError): d.allocate_slots(self.s,fake)

    def test_duplicate_candidate_rejected(self):
        fake=copy.deepcopy(self.c);fake[0]=copy.deepcopy(fake[1])
        with self.assertRaises(ValueError): d.allocate_slots(self.s,fake)

    def test_two_repeats(self):
        s=spec();s['planned_repeats']=2
        self.assertEqual(len(d.allocate_slots(s,d.enumerate_snvs(s))),30)


class OutputContracts(unittest.TestCase):
    def setUp(self):
        self.s=spec();self.p=d.build_plan(self.s,'a'*64)

    def test_coverage(self):
        self.assertEqual((self.p['expected_possible_snvs'],self.p['candidate_count'],self.p['toy_snv_coverage_fraction']),(12,12,1.0))

    def test_model_not_performed(self):
        self.assertFalse(self.p['experiment_performed'])
        self.assertFalse(self.p['saturation_genome_editing_performed'])
        self.assertFalse(self.p['minigene_assay_performed'])
        self.assertFalse(self.p['reporter_assay_performed'])

    def test_no_api_calls(self):
        self.assertEqual(self.p['real_api_calls_this_run'],0)

    def test_no_clinical_interpretation(self):
        self.assertFalse(self.p['clinical_interpretation'])
        self.assertEqual(self.p['research_ethics_approval'],'NOT_ASSESSED')

    def test_invalid_digest_rejected(self):
        for s in ('', 'xyz', 'g'*64, 'A'*64):
            with self.subTest(s=s), self.assertRaises(ValueError): d.build_plan(self.s,s)

    def test_all_artifact_types(self):
        self.assertEqual(set(d.artifact_bytes(self.p)),set(d.ARTIFACT_NAMES))

    def test_candidate_csv_12_rows(self):
        data=d.artifact_bytes(self.p)['day15_snv_candidates.csv'].decode('utf8')
        self.assertEqual(len(list(csv.DictReader(io.StringIO(data)))),12)

    def test_slot_csv_45_rows(self):
        data=d.artifact_bytes(self.p)['day15_planned_slots.csv'].decode('utf8')
        self.assertEqual(len(list(csv.DictReader(io.StringIO(data)))),45)

    def test_svg_is_parseable_and_16_9(self):
        svg=d.render_svg(self.p);root=ET.fromstring(svg)
        self.assertEqual((root.attrib['width'],root.attrib['height']),('1920','1080'))
        self.assertIn('#E5DEDA',svg)

    def test_svg_not_wetlab_evidence(self):
        self.assertIn('SYNTHETIC / IN SILICO ONLY',d.render_svg(self.p))

    def test_html_has_explicit_disclaimer(self):
        t=d.render_html(self.p,d.render_svg(self.p))
        self.assertIn('NOT_MEASURED / NOT_PERFORMED',t)
        self.assertIn('chrToy:101:A&gt;C',t)

    def test_notes_warn_no_biological_replicates(self):
        self.assertIn('not verified biological replicates',d.render_notes(self.p))

    def test_json_null(self):
        self.assertIn('"observed_effect": null',d.to_json(self.p))


class CliRunAndVerifyTests(unittest.TestCase):
    def test_run_and_verify(self):
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)/'out';r=d.run(FIXTURE,out)
            self.assertEqual(r['status'],'PASS')
            self.assertEqual(r['artifact_count'],7)
            self.assertEqual(d.verify(FIXTURE,out)['checksums'],'MATCH')

    def test_run_uses_real_fixture_sha(self):
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)/'out';d.run(FIXTURE,out)
            p=json.loads((out/'day15_design_plan.json').read_text())
            self.assertEqual(p['fixture_provenance']['sha256'],d.sha256_bytes(FIXTURE.read_bytes()))

    def test_reject_existing_output(self):
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)/'out';out.mkdir()
            with self.assertRaises(FileExistsError):d.run(FIXTURE,out)

    def test_missing_parent_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            with self.assertRaises(ValueError):d.run(FIXTURE,Path(td)/'missing'/'out')

    def test_tampered_svg_fails_verify(self):
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)/'out';d.run(FIXTURE,out)
            with (out/'day15_saturation_map.svg').open('a',encoding='utf8') as f:f.write('TAMPER')
            with self.assertRaises(ValueError):d.verify(FIXTURE,out)

    def test_missing_artifact_fails(self):
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)/'out';d.run(FIXTURE,out);(out/'day15_design_notes.md').unlink()
            with self.assertRaises(ValueError):d.verify(FIXTURE,out)

    def test_modified_fixture_fails_verify(self):
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)/'out';d.run(FIXTURE,out)
            source=Path(td)/'modified.json'
            s=spec();s['randomization_seed']+=1;source.write_text(json.dumps(s))
            with self.assertRaises(ValueError):d.verify(source,out)

    def test_rerun_produces_bitwise_same_artifacts(self):
        with tempfile.TemporaryDirectory() as td:
            a,b=Path(td)/'a',Path(td)/'b';d.run(FIXTURE,a);d.run(FIXTURE,b)
            for filename in (*d.ARTIFACT_NAMES,'day15_manifest.json'):
                self.assertEqual((a/filename).read_bytes(),(b/filename).read_bytes(),filename)

    def test_cli_rejects_live(self):
        with tempfile.TemporaryDirectory() as td:
            p=subprocess.run([sys.executable,str(ROOT/'day15_offline_saturation_design_lab.py'),'run','--mode','live','--output-dir',str(Path(td)/'out')],capture_output=True,text=True)
            self.assertNotEqual(p.returncode,0)
            self.assertFalse((Path(td)/'out').exists())

    def test_cli_run_stdout(self):
        with tempfile.TemporaryDirectory() as td:
            p=subprocess.run([sys.executable,str(ROOT/'day15_offline_saturation_design_lab.py'),'run','--fixture',str(FIXTURE),'--output-dir',str(Path(td)/'out')],capture_output=True,text=True,check=True)
            res=json.loads(p.stdout)
            self.assertEqual((res['candidate_count'],res['planned_slot_count'],res['real_api_calls_this_run']),(12,45,0))

    def test_cli_verify_stdout(self):
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)/'out';d.run(FIXTURE,out)
            p=subprocess.run([sys.executable,str(ROOT/'day15_offline_saturation_design_lab.py'),'verify','--fixture',str(FIXTURE),'--output-dir',str(out)],capture_output=True,text=True,check=True)
            self.assertEqual(json.loads(p.stdout)['checksums'],'MATCH')


if __name__=='__main__':
    unittest.main()
