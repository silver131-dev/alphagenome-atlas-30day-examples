"""Day22 offline-only unit and contract tests; Python 3.12+."""
from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from day21_offline_vcf_review_lab import (
    ASSEMBLY, QC_ORDER, check_qc, in_bed, parse_synthetic_vcf,
    load_toy_reference, parse_toy_bed,
)
from day22_offline_qc_explorer_lab import (
    FAIL_REASON, EVIDENCE_CLASS, execute_offline, extra_probe_rows,
    make_probe, qc_report, render_html, trace_qc, validate_min_qual,
)

HERE = Path(__file__).resolve().parent
VCF_PATH = HERE / "day21_synthetic.vcf"
REF_PATH = HERE / "day21_toy_reference.json"
BED_PATH = HERE / "day21_synthetic_regions.bed"
REF = load_toy_reference(REF_PATH.read_text(encoding="utf-8"))
BED = parse_toy_bed(BED_PATH.read_text(encoding="utf-8"))
VCF = VCF_PATH.read_text(encoding="utf-8")

EXPECTED_SHA = {
    "day21_synthetic.vcf": "adfc480a659e581a4d69f6b7b01890016ddc6ad0d6de6e709521fc6c2804800d",
    "day21_toy_reference.json": "43e8cda00040afb123073f628b1249826f37d4fb4732aed578bf64ab28f9cd5a",
    "day21_synthetic_regions.bed": "219e3264fc7b05d0d3c5ee8af57e07d6cb80df4e7618bcf11be4af2b7640a910",
}


class VcfAndFixtureContinuity(unittest.TestCase):
    def test_day21_vcf_sha_exact(self):
        self.assertEqual(hashlib.sha256(VCF_PATH.read_bytes()).hexdigest(), EXPECTED_SHA[VCF_PATH.name])

    def test_day21_reference_sha_exact(self):
        self.assertEqual(hashlib.sha256(REF_PATH.read_bytes()).hexdigest(), EXPECTED_SHA[REF_PATH.name])

    def test_day21_bed_sha_exact(self):
        self.assertEqual(hashlib.sha256(BED_PATH.read_bytes()).hexdigest(), EXPECTED_SHA[BED_PATH.name])

    def test_seven_vcf_rows_expand_to_eight_alts(self):
        self.assertEqual(len([l for l in VCF.splitlines() if l and not l.startswith('#')]), 7)
        self.assertEqual(len(parse_synthetic_vcf(VCF)), 8)

    def test_multiallelic_candidate_identifiers(self):
        vs = [r['variant_id'] for r in parse_synthetic_vcf(VCF)]
        self.assertIn('chr1:110:A>G', vs)
        self.assertIn('chr1:110:A>C', vs)
        self.assertEqual(len(vs), len(set(vs)))

    def test_synthetic_only_fixture(self):
        with self.assertRaises(ValueError):
            parse_synthetic_vcf(VCF.replace('##source=SYNTHETIC_EDUCATIONAL_DEMO', '##source=clinical'))

    def test_toy_not_real_reference(self):
        self.assertEqual(ASSEMBLY, 'TOY_NOT_GRCH38')
        self.assertEqual(len(REF['chr1']), 200)
        self.assertEqual(set(REF['chr1']), {'A'})

    def test_bed_exact_0based_intervals(self):
        self.assertEqual(BED, [('chr1',90,115),('chr1',145,165)])

    def test_no_real_scores_needed_by_day22(self):
        self.assertNotIn('scores', qc_report(VCF, REF, BED))


class DecisionPathTests(unittest.TestCase):
    def test_qc_order_matches_day21(self):
        self.assertEqual(QC_ORDER, ['FILTER','BED','REF','INDEL','MISSING_QUAL','LOW_QUAL'])

    def test_all_known_fail_reasons_present(self):
        self.assertEqual(set(FAIL_REASON), set(QC_ORDER))

    def test_filter_stops_first(self):
        row=make_probe('A',130,'G','T',qual=10,filt='q10')
        r=trace_qc(row,REF,BED)
        self.assertEqual((r['qc_status'],r['reason']),('REJECTED','FILTER_NOT_PASS'))
        self.assertEqual([s['check'] for s in r['checked_steps']], ['FILTER'])
        self.assertEqual(r['not_checked_steps'][0],'BED')

    def test_out_of_bed_low_qual_is_masked(self):
        r=trace_qc(make_probe('B',120,'A','T',qual=10),REF,BED)
        self.assertEqual(r['reason'],'OUTSIDE_TEACHING_BED')
        self.assertIn('LOW_QUAL',r['not_checked_steps'])

    def test_out_of_bed_ref_mismatch_is_masked(self):
        r=trace_qc(make_probe('C',140,'G','T'),REF,BED)
        self.assertEqual(r['reason'],'OUTSIDE_TEACHING_BED')
        self.assertIn('REF',r['not_checked_steps'])

    def test_in_bed_low_qual_is_reached(self):
        r=trace_qc(make_probe('D',105,'A','T',qual=10),REF,BED)
        self.assertEqual(r['reason'],'LOW_QUAL')
        self.assertEqual(r['checked_steps'][-1]['check'],'LOW_QUAL')

    def test_in_bed_ref_mismatch_is_reached(self):
        r=trace_qc(make_probe('E',106,'G','T'),REF,BED)
        self.assertEqual(r['reason'],'REF_MISMATCH')
        self.assertEqual(r['checked_steps'][-1]['check'],'REF')
        self.assertNotIn('LOW_QUAL',[s['check'] for s in r['checked_steps']])

    def test_indel_sent_for_review(self):
        r=trace_qc(make_probe('F',150,'A','AT'),REF,BED)
        self.assertEqual((r['qc_status'],r['reason']),('REVIEW_REQUIRED','INDEL_NOT_NORMALIZED'))
        self.assertEqual(r['checked_steps'][-1]['check'],'INDEL')

    def test_missing_qual_sent_for_review(self):
        r=trace_qc(make_probe('G',160,'A','T',qual=None),REF,BED)
        self.assertEqual((r['qc_status'],r['reason']),('REVIEW_REQUIRED','MISSING_QUAL'))
        self.assertEqual(r['checked_steps'][-1]['check'],'MISSING_QUAL')

    def test_bed_coordinate_boundary_below_start(self):
        self.assertFalse(in_bed('chr1',90,BED))
        self.assertEqual(trace_qc(make_probe('H',90,'A','C'),REF,BED)['reason'],'OUTSIDE_TEACHING_BED')

    def test_bed_first_included_one_based_91(self):
        self.assertTrue(in_bed('chr1',91,BED))
        r=trace_qc(make_probe('I',91,'A','C'),REF,BED)
        self.assertEqual((r['qc_status'],r['bed_anchor_0based']),('SELECTED',90))

    def test_bed_last_included_one_based_115(self):
        self.assertEqual(trace_qc(make_probe('J',115,'A','C'),REF,BED)['qc_status'],'SELECTED')

    def test_bed_first_excluded_one_based_116(self):
        self.assertEqual(trace_qc(make_probe('K',116,'A','C'),REF,BED)['reason'],'OUTSIDE_TEACHING_BED')

    def test_second_bed_first_included_one_based_146(self):
        self.assertEqual(trace_qc(make_probe('L',146,'A','C'),REF,BED)['qc_status'],'SELECTED')

    def test_chromosome_separation(self):
        self.assertEqual(trace_qc(make_probe('M',105,'A','C',chrom='chr2'),REF,BED)['reason'],'OUTSIDE_TEACHING_BED')

    def test_all_selected_reach_final_qual_step(self):
        rows=parse_synthetic_vcf(VCF)
        sel=[trace_qc(r,REF,BED) for r in rows if check_qc(r,REF,BED)[0]=='SELECTED']
        self.assertEqual(len(sel),3)
        self.assertTrue(all(len(x['checked_steps'])==6 for x in sel))

    def test_validate_threshold_nan_inf_negative(self):
        for value in [float('nan'),float('inf'),float('-inf'),-1,True,'20',None]:
            with self.subTest(value=repr(value)),self.assertRaises(ValueError):
                validate_min_qual(value)

    def test_changed_threshold_respects_decisions(self):
        r=qc_report(VCF,REF,BED,min_qual=55)
        self.assertEqual((r['selected_count'],r['rejected_count'],r['review_required_count']), (1,5,2))

    def test_source_day21_qc_parity_for_every_probe(self):
        for row in extra_probe_rows():
            with self.subTest(label=row['label']):
                new=trace_qc(row,REF,BED)
                self.assertEqual((new['qc_status'],new['reason']),check_qc(row,REF,BED))


class ReportAndArtifactTests(unittest.TestCase):
    def test_report_counts_are_day21_unchanged(self):
        r=qc_report(VCF,REF,BED)
        self.assertEqual((r['candidate_count'],r['selected_count'],r['review_required_count'],r['rejected_count']), (8,3,2,3))

    def test_reason_distribution(self):
        r=qc_report(VCF,REF,BED)
        self.assertEqual(r['reason_counts'],{
            'TEACHING_QC_PASS':3,'OUTSIDE_TEACHING_BED':2,'FILTER_NOT_PASS':1,
            'INDEL_NOT_NORMALIZED':1,'MISSING_QUAL':1})

    def test_supplementary_probe_independent_of_default_counts(self):
        r=qc_report(VCF,REF,BED)
        self.assertEqual(len(r['supplementary_probe_cases']),10)
        self.assertEqual(r['candidate_count'],8)
        self.assertEqual(len(r['qc_rows']),8)

    def test_evidence_labels_and_live_disabled(self):
        r=qc_report(VCF,REF,BED)
        self.assertEqual(r['data_kind'],EVIDENCE_CLASS)
        self.assertEqual(r['mode'],'OFFLINE_ONLY')
        self.assertFalse(r['real_genome_used'])
        self.assertEqual(r['real_api_calls_this_run'],0)
        self.assertFalse(r['atlas_query_performed'])
        self.assertFalse(r['clinical_interpretation'])

    def test_report_can_roundtrip_json(self):
        r=qc_report(VCF,REF,BED)
        self.assertEqual(json.loads(json.dumps(r))['qc_rows'][0]['reason'],'TEACHING_QC_PASS')

    def test_render_html_mentions_synthetic_and_qc(self):
        h=render_html(qc_report(VCF,REF,BED))
        self.assertIn('SYNTHETIC EDUCATIONAL DEMO',h)
        self.assertIn('TOY_NOT_GRCH38',h)
        self.assertIn('REF_MISMATCH_INSIDE_BED',h)
        self.assertIn('atlas_query_performed=false',h)
        self.assertIn('<details>',h)

    def test_html_escapes_report_user_strings(self):
        r=qc_report(VCF,REF,BED)
        r['qc_rows'][0]['variant_id']='<script>alert(1)</script>'
        html=render_html(r)
        self.assertNotIn('<script>',html)
        self.assertIn('&lt;script&gt;',html)

    def test_file_artifacts_created_with_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            s=execute_offline(VCF_PATH,REF_PATH,BED_PATH,Path(tmp))
            self.assertEqual(s['status'],'PASS')
            self.assertEqual(s['real_api_calls_this_run'],0)
            obj=json.loads((Path(tmp)/'day22_qc_report.json').read_text())
            self.assertEqual(obj['candidate_count'],8)
            self.assertEqual(obj['input_provenance']['vcf']['sha256'],EXPECTED_SHA[VCF_PATH.name])
            self.assertTrue((Path(tmp)/'day22_qc_dashboard.html').is_file())

    def test_cli_outputs_match_contract(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=subprocess.run([sys.executable,str(HERE/'day22_offline_qc_explorer_lab.py'),
                              '--mode','offline','--vcf',str(VCF_PATH),
                              '--reference',str(REF_PATH),'--regions',str(BED_PATH),
                              '--output-dir',str(Path(tmp)/'out')],
                             cwd=HERE,capture_output=True,text=True,check=True)
            s=json.loads(p.stdout)
            self.assertEqual((s['candidate_count'],s['selected_count'],s['review_required_count'],s['rejected_count']), (8,3,2,3))
            self.assertEqual(s['low_qual_probe'],'LOW_QUAL')
            self.assertEqual(s['ref_mismatch_probe'],'REF_MISMATCH')
            self.assertTrue(Path(s['json_file']).is_file())

    def test_explicit_live_rejected(self):
        p=subprocess.run([sys.executable,str(HERE/'day22_offline_qc_explorer_lab.py'),'--mode','live'],
                         cwd=HERE,capture_output=True,text=True)
        self.assertNotEqual(p.returncode,0)
        self.assertIn('invalid choice',p.stderr)

    def test_invalid_threshold_rejected(self):
        with self.assertRaises(ValueError):
            qc_report(VCF,REF,BED,min_qual=-0.5)

    def test_python_312_ast_parse(self):
        for name in ('day21_offline_vcf_review_lab.py','day22_offline_qc_explorer_lab.py',
                     'test_day22_offline_qc_explorer_lab.py'):
            with self.subTest(name=name):
                ast.parse((HERE/name).read_text(encoding='utf-8'),feature_version=(3,12))

    def test_day22_has_no_network_or_alphagenome_imports(self):
        tree=ast.parse((HERE/'day22_offline_qc_explorer_lab.py').read_text(encoding='utf-8'))
        imports=[i.module for i in ast.walk(tree) if isinstance(i,ast.ImportFrom)]
        imports.extend(a.name for i in ast.walk(tree) if isinstance(i,ast.Import) for a in i.names)
        self.assertNotIn('alphagenome',imports)
        self.assertNotIn('requests',imports)
        self.assertNotIn('httpx',imports)
        self.assertNotIn('urllib.request',imports)

    def test_byte_deterministic_same_inputs(self):
        with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
            execute_offline(VCF_PATH,REF_PATH,BED_PATH,a)
            execute_offline(VCF_PATH,REF_PATH,BED_PATH,b)
            for name in ('day22_qc_report.json','day22_qc_dashboard.html'):
                self.assertEqual((Path(a)/name).read_bytes(),(Path(b)/name).read_bytes())


if __name__=='__main__':
    unittest.main()
