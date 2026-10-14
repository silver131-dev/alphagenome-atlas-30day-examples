"""Contract and scientific-boundary tests for the Day 07 synthetic-only lab."""
from __future__ import annotations

import ast
import copy
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

import day07_offline_cnn_transformer_lab as lab

HERE = Path(__file__).resolve().parent
FIXTURE_PATH = HERE / "day07_synthetic_multiscale.json"
SCRIPT_PATH = HERE / "day07_offline_cnn_transformer_lab.py"


class FixtureContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture, cls.sha = lab.load_fixture(FIXTURE_PATH)

    def test_fixture_loads(self):
        self.assertEqual(self.fixture["chrom"], "chrToy")

    def test_input_sha256_matches_bytes(self):
        self.assertEqual(self.sha, hashlib.sha256(FIXTURE_PATH.read_bytes()).hexdigest())

    def test_requires_synthetic_marker(self):
        bad = copy.deepcopy(self.fixture)
        bad["evidence_class"] = "HISTORICAL_LIVE_OUTPUT"
        self.assert_bad(bad)

    def test_rejects_real_assembly(self):
        bad = copy.deepcopy(self.fixture)
        bad["assembly"] = "GRCh38"
        self.assert_bad(bad)

    def test_rejects_extra_patient_id(self):
        bad = copy.deepcopy(self.fixture)
        bad["patient_id"] = "DO_NOT_IMPORT"
        self.assert_bad(bad)

    def test_rejects_missing_source(self):
        bad = copy.deepcopy(self.fixture)
        bad.pop("sequence")
        self.assert_bad(bad)

    def test_rejects_non_object(self):
        self.assert_bad([self.fixture])

    def test_rejects_non_toy_chrom(self):
        bad = copy.deepcopy(self.fixture)
        bad["chrom"] = "chr1"
        self.assert_bad(bad)

    def test_rejects_nonpositive_start(self):
        bad = copy.deepcopy(self.fixture)
        bad["sequence_start_1based"] = 0
        self.assert_bad(bad)

    def test_rejects_boolean_start(self):
        bad = copy.deepcopy(self.fixture)
        bad["sequence_start_1based"] = True
        self.assert_bad(bad)

    def test_rejects_short_dna(self):
        bad = copy.deepcopy(self.fixture)
        bad["sequence"] = bad["sequence"][:-1]
        self.assert_bad(bad)

    def test_rejects_ambiguous_base(self):
        bad = copy.deepcopy(self.fixture)
        bad["sequence"] = bad["sequence"][:-1]+"N"
        self.assert_bad(bad)

    def test_rejects_lowercase_dna(self):
        bad = copy.deepcopy(self.fixture)
        bad["sequence"] = bad["sequence"].lower()
        self.assert_bad(bad)

    def test_rejects_bin_size_16(self):
        bad = copy.deepcopy(self.fixture)
        bad["bin_size"] = 16
        self.assert_bad(bad)

    def test_rejects_noninteger_bin_size(self):
        bad = copy.deepcopy(self.fixture)
        bad["bin_size"] = 8.0
        self.assert_bad(bad)

    def test_rejects_longer_local_motif(self):
        bad = copy.deepcopy(self.fixture)
        bad["local_motif"] = "ACGT"
        self.assert_bad(bad)

    def test_rejects_motif_mismatch(self):
        bad = copy.deepcopy(self.fixture)
        bad["local_motif_start_1based"] += 1
        self.assert_bad(bad)

    def test_rejects_wrong_remote_distance(self):
        bad = copy.deepcopy(self.fixture)
        bad["distal_motif_start_1based"] += 1
        self.assert_bad(bad)

    def test_rejects_extra_variant_field(self):
        bad = copy.deepcopy(self.fixture)
        bad["variants"][0]["patient_id"] = "X"
        self.assert_bad(bad)

    def test_rejects_missing_variant(self):
        bad = copy.deepcopy(self.fixture)
        bad["variants"].pop()
        self.assert_bad(bad)

    def test_rejects_swapped_variant_order(self):
        bad = copy.deepcopy(self.fixture)
        bad["variants"].reverse()
        self.assert_bad(bad)

    def test_rejects_real_variant_source(self):
        bad = copy.deepcopy(self.fixture)
        bad["variants"][0]["source"] = "PATIENT_VCF"
        self.assert_bad(bad)

    def test_rejects_ref_mismatch(self):
        bad = copy.deepcopy(self.fixture)
        bad["variants"][0]["ref"] = "A"
        self.assert_bad(bad)

    def test_rejects_equal_ref_alt(self):
        bad = copy.deepcopy(self.fixture)
        bad["variants"][0]["alt"] = "T"
        self.assert_bad(bad)

    def test_rejects_out_of_window_variant(self):
        bad = copy.deepcopy(self.fixture)
        bad["variants"][0]["position_1based"] = 3000
        self.assert_bad(bad)

    def test_rejects_invalid_variant_id(self):
        bad = copy.deepcopy(self.fixture)
        bad["variants"][0]["variant_id"] = "chr22:36201698:A>C"
        self.assert_bad(bad)

    def test_rejects_non_snv(self):
        bad = copy.deepcopy(self.fixture)
        bad["variants"][0]["alt"] = "CG"
        self.assert_bad(bad)

    def test_rejects_duplicate_local_hit(self):
        bad = copy.deepcopy(self.fixture)
        seq = list(bad["sequence"])
        seq[25:28] = list("ACG")
        bad["sequence"] = "".join(seq)
        self.assert_bad(bad)

    def test_rejects_variant_that_preserves_motif(self):
        bad = copy.deepcopy(self.fixture)
        bad["variants"][0]["position_1based"] = 1052
        bad["variants"][0]["ref"] = "G"
        bad["variants"][0]["alt"] = "A"
        bad["variants"][0]["variant_id"] = "chrToy:1052:G>A"
        self.assert_bad(bad)

    def assert_bad(self, data):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "fixture.json"
            p.write_text(json.dumps(data), encoding="utf-8")
            with self.assertRaises((ValueError, TypeError, KeyError)):
                lab.load_fixture(p)


class LocalKernelAndCoordinateTests(unittest.TestCase):
    def test_find_exact_hits(self):
        self.assertEqual(lab.find_exact_hits("AACGA", "ACG"), [1])

    def test_find_overlapping_hits(self):
        self.assertEqual(lab.find_exact_hits("AAAA", "AA"), [0,1,2])

    def test_no_motif_hit(self):
        self.assertEqual(lab.find_exact_hits("AAAA", "CGT"), [])

    def test_invalid_motif(self):
        with self.assertRaises(ValueError): lab.find_exact_hits("AAAA", "")

    def test_invalid_dna(self):
        with self.assertRaises(ValueError): lab.find_exact_hits("AACN", "AC")

    def test_fixed_kernel_full_match(self):
        self.assertEqual(lab.fixed_kernel_scan("ACG", "ACG"), [1.0])

    def test_fixed_kernel_partial_match(self):
        self.assertAlmostEqual(lab.fixed_kernel_scan("ATG", "ACG")[0], 2/3)

    def test_fixed_kernel_no_match(self):
        self.assertEqual(lab.fixed_kernel_scan("TTT", "ACG"), [0.0])

    def test_fixed_kernel_invalid(self):
        with self.assertRaises(ValueError): lab.fixed_kernel_scan("ATX", "ACG")

    def test_pool_eight_bins(self):
        f,_=lab.load_fixture(FIXTURE_PATH)
        p,h=lab.pool_full_hits(f["sequence"], "ACG", 8)
        self.assertEqual(p, [0,1,0,0,0,0,0,0])
        self.assertEqual(h, [10])

    def test_pool_distal(self):
        f,_=lab.load_fixture(FIXTURE_PATH)
        p,h=lab.pool_full_hits(f["sequence"], "TGA", 8)
        self.assertEqual(p, [0,0,0,0,0,0,1,0]);self.assertEqual(h,[50])

    def test_pool_invalid_size(self):
        with self.assertRaises(ValueError): lab.pool_full_hits("AAAA", "AA", 3)

    def test_snv_changes_one_base(self):
        self.assertEqual(lab.apply_snv("ACG", 1012, "C", "T", 1011), "ATG")

    def test_snv_does_not_modify_reference(self):
        seq="ACG";lab.apply_snv(seq,1012,"C","T",1011);self.assertEqual(seq,"ACG")

    def test_snv_wrong_ref_rejected(self):
        with self.assertRaises(ValueError): lab.apply_snv("ACG",1012,"G","T",1011)

    def test_snv_outside_rejected(self):
        with self.assertRaises(ValueError): lab.apply_snv("ACG",1010,"C","T",1011)

    def test_snv_equal_rejected(self):
        with self.assertRaises(ValueError): lab.apply_snv("ACG",1012,"C","C",1011)


class AttentionAndMultiscaleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture,cls.sha=lab.load_fixture(FIXTURE_PATH)
        cls.report=lab.build_report(cls.fixture,cls.sha,FIXTURE_PATH.name)

    def test_softmax_uniform(self):
        self.assertEqual(lab.softmax([0,0,0,0]),[.25]*4)

    def test_softmax_one_high(self):
        w=lab.softmax([0,0,3,0]);self.assertGreater(w[2],0.8);self.assertAlmostEqual(sum(w),1.0)

    def test_softmax_rejects_empty(self):
        with self.assertRaises(ValueError): lab.softmax([])

    def test_softmax_rejects_nan(self):
        with self.assertRaises(ValueError): lab.softmax([0,float('nan')])

    def test_softmax_rejects_inf(self):
        with self.assertRaises(ValueError): lab.softmax([0,float('inf')])

    def test_attention_selects_bin_six(self):
        a=[0,1,0,0,0,0,0,0];b=[0,0,0,0,0,0,1,0]
        weights,context=lab.scalar_attention(a,b,1)
        self.assertEqual(max(range(8),key=lambda x:weights[x]),6)
        self.assertAlmostEqual(weights[6],context)

    def test_attention_weights_sum_one(self):
        a=[0,1,0,0,0,0,0,0];b=[0,0,0,0,0,0,1,0]
        weights,_=lab.scalar_attention(a,b,1)
        self.assertAlmostEqual(sum(weights),1)

    def test_attention_no_distal_yields_zero_context(self):
        a=[0,1,0,0,0,0,0,0];b=[0]*8
        _,context=lab.scalar_attention(a,b,1)
        self.assertEqual(context,0)

    def test_attention_invalid_bin(self):
        with self.assertRaises(ValueError):lab.scalar_attention([1],[1],1)

    def test_attention_invalid_shape(self):
        with self.assertRaises(ValueError):lab.scalar_attention([1],[1,0],0)

    def test_attention_rejects_non_binary(self):
        with self.assertRaises(ValueError):lab.scalar_attention([.5,0],[1,0],0)

    def test_reference_exact_local_hit(self):
        r=self.report["cases"][0]["local_context"]
        self.assertEqual(r["local_kernel_full_hits_1based"],[1011])

    def test_reference_exact_distal_hit(self):
        r=self.report["cases"][0]["local_context"]
        self.assertEqual(r["distal_kernel_full_hits_1based"],[1051])

    def test_reference_fused_score(self):
        r=self.report["cases"][0]["local_context"]
        self.assertAlmostEqual(r["toy_fused_output"],0.87078,places=5)

    def test_distal_change_preserves_local_hit(self):
        self.assertTrue(self.report["cases"][1]["local_cnn_unchanged"])
        self.assertEqual(self.report["cases"][1]["local_context"]["local_kernel_full_hits_1based"],[1011])

    def test_distal_change_erases_remote_hit(self):
        self.assertEqual(self.report["cases"][1]["local_context"]["distal_kernel_full_hits_1based"],[])

    def test_distal_change_proxy_is_half(self):
        self.assertEqual(self.report["cases"][1]["local_context"]["toy_fused_output"],0.5)

    def test_local_change_erases_local_hit(self):
        self.assertEqual(self.report["cases"][2]["local_context"]["local_kernel_full_hits_1based"],[])

    def test_local_change_proxy_zero(self):
        self.assertEqual(self.report["cases"][2]["local_context"]["toy_fused_output"],0.0)

    def test_local_change_preserves_remote_hit(self):
        self.assertEqual(self.report["cases"][2]["local_context"]["distal_kernel_full_hits_1based"],[1051])

    def test_reference_scores_are_not_probabilities(self):
        self.assertFalse(self.report["model_prediction"])
        self.assertFalse(self.report["clinical_interpretation"])

    def test_reference_coarse_bin_size(self):
        self.assertEqual(self.report["bin_count"],8)
        self.assertEqual(self.report["bin_size_bp"],8)

    def test_decoder_one_bp_position(self):
        self.assertEqual(self.report["cases"][0]["local_context"]["decoder_proxy_1bp_nonzero"][0]["position_1based"],1011)

    def test_distance_and_case_count(self):
        self.assertEqual(self.report["linear_distance_bp"],40)
        self.assertEqual(len(self.report["cases"]),3)


class OutputAndSafetyTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.folder=Path(self.tmp.name)
        self.fixture=FIXTURE_PATH
        self.out=self.folder/'output'
        lab.run(self.fixture,self.out)

    def test_run_outputs_five_files(self):
        self.assertEqual({p.name for p in self.out.iterdir()},set(lab.OUTPUT_FILES)|{lab.MANIFEST_NAME})

    def test_run_refuses_existing_dir(self):
        with self.assertRaises(FileExistsError):lab.run(self.fixture,self.out)

    def test_verify_succeeds(self):
        self.assertEqual(lab.verify(self.fixture,self.out)["checksums"],"MATCH")

    def test_json_valid(self):
        self.assertEqual(json.loads((self.out/lab.OUTPUT_FILES[0]).read_text())["bin_count"],8)

    def test_csv_has_24_rows(self):
        with (self.out/lab.OUTPUT_FILES[1]).open(newline='',encoding='utf-8') as f:
            a=list(csv.DictReader(f))
        self.assertEqual(len(a),24)

    def test_csv_three_distinct_cases(self):
        with (self.out/lab.OUTPUT_FILES[1]).open(newline='',encoding='utf-8') as f:
            a=list(csv.DictReader(f))
        self.assertEqual({r["case"] for r in a},{"REFERENCE","DISTAL_CHANGE","LOCAL_CHANGE"})

    def test_csv_attention_from_local_bin(self):
        with (self.out/lab.OUTPUT_FILES[1]).open(newline='',encoding='utf-8') as f:
            rows=list(csv.DictReader(f))
        ref=[r for r in rows if r['case']=='REFERENCE']
        self.assertEqual(max(ref,key=lambda x:float(x['attention_from_local_bin']))['bin_index_0based'],'6')

    def test_csv_coordinate_mapping(self):
        with (self.out/lab.OUTPUT_FILES[1]).open(newline='',encoding='utf-8') as f:
            rows=list(csv.DictReader(f))
        self.assertEqual(rows[0]["bin_start_1based"],"1001")
        self.assertEqual(rows[0]["bin_end_1based_inclusive"],"1008")

    def test_svg_xml_valid(self):
        root=ET.fromstring((self.out/lab.OUTPUT_FILES[2]).read_text())
        self.assertTrue(root.tag.endswith('svg'))

    def test_svg_dimensions(self):
        root=ET.parse(self.out/lab.OUTPUT_FILES[2]).getroot()
        self.assertEqual(root.attrib['viewBox'],'0 0 1920 1080')

    def test_svg_uses_requested_background(self):
        self.assertIn('#E5DEDA',(self.out/lab.OUTPUT_FILES[2]).read_text())

    def test_html_is_local_standalone(self):
        doc=(self.out/lab.OUTPUT_FILES[3]).read_text(encoding='utf-8')
        self.assertIn('<svg ',doc);self.assertNotIn('<script',doc)

    def test_html_synthetic_disclaimer(self):
        doc=(self.out/lab.OUTPUT_FILES[3]).read_text(encoding='utf-8')
        self.assertIn('SYNTHETIC',doc);self.assertIn('OFFLINE ONLY',doc)

    def test_manifest_contains_source_sha(self):
        m=json.loads((self.out/lab.MANIFEST_NAME).read_text())
        self.assertEqual(m['source']['sha256'],hashlib.sha256(self.fixture.read_bytes()).hexdigest())

    def test_manifest_not_signed(self):
        m=json.loads((self.out/lab.MANIFEST_NAME).read_text())
        self.assertFalse(m['signed_provenance'])

    def test_manifests_has_exactly_four_outputs(self):
        m=json.loads((self.out/lab.MANIFEST_NAME).read_text())
        self.assertEqual(set(m['output_sha256']),set(lab.OUTPUT_FILES))

    def test_tamper_json_detected(self):
        p=self.out/lab.OUTPUT_FILES[0];p.write_text(p.read_text()+" ",encoding='utf-8')
        with self.assertRaises(ValueError):lab.verify(self.fixture,self.out)

    def test_tamper_csv_detected(self):
        p=self.out/lab.OUTPUT_FILES[1];p.write_text(p.read_text().replace('REFERENCE','OTHER',1),encoding='utf-8')
        with self.assertRaises(ValueError):lab.verify(self.fixture,self.out)

    def test_tamper_svg_detected(self):
        p=self.out/lab.OUTPUT_FILES[2];p.write_text(p.read_text().replace('#E5DEDA','#FFFFFF'),encoding='utf-8')
        with self.assertRaises(ValueError):lab.verify(self.fixture,self.out)

    def test_tamper_html_detected(self):
        p=self.out/lab.OUTPUT_FILES[3];p.write_text(p.read_text()+'fake',encoding='utf-8')
        with self.assertRaises(ValueError):lab.verify(self.fixture,self.out)

    def test_tamper_manifest_detected(self):
        p=self.out/lab.MANIFEST_NAME;m=json.loads(p.read_text());m['signed_provenance']=True;p.write_text(json.dumps(m))
        with self.assertRaises(ValueError):lab.verify(self.fixture,self.out)

    def test_tamper_fixture_detected(self):
        p=self.folder/'fixture.json';d=json.loads(self.fixture.read_text());d['sequence_start_1based']+=1
        p.write_text(json.dumps(d))
        with self.assertRaises(ValueError):lab.verify(p,self.out)

    def test_missing_file_detected(self):
        (self.out/lab.OUTPUT_FILES[1]).unlink()
        with self.assertRaises(FileNotFoundError):lab.verify(self.fixture,self.out)

    def test_byte_stable_run(self):
        second=self.folder/'new';lab.run(self.fixture,second)
        self.assertEqual({p.name:p.read_bytes() for p in self.out.iterdir()},
                         {p.name:p.read_bytes() for p in second.iterdir()})

    def test_imports_no_network_or_sdk(self):
        tree=ast.parse(SCRIPT_PATH.read_text(encoding='utf-8'))
        mods=set()
        for node in ast.walk(tree):
            if isinstance(node,ast.Import):mods.update(alias.name.split('.')[0] for alias in node.names)
            if isinstance(node,ast.ImportFrom):mods.add((node.module or '').split('.')[0])
        self.assertFalse(mods & {'socket','requests','urllib','httpx','alphagenome','google','subprocess'})

    def test_cli_offline_in_api_key_environment(self):
        env=os.environ.copy();env['ALPHAGENOME_API_KEY']='SYNTHETIC_DO_NOT_USE';env['ALPHA_GENOME_API_KEY']='SYNTHETIC_DO_NOT_USE'
        cmd=[sys.executable, str(SCRIPT_PATH),'verify','--mode','offline',
             '--fixture',str(self.fixture),'--output-dir',str(self.out)]
        p=subprocess.run(cmd,env=env,capture_output=True,text=True,check=True)
        self.assertEqual(json.loads(p.stdout)['real_api_calls_this_run'],0)

    def test_cli_refuses_live(self):
        cmd=[sys.executable,str(SCRIPT_PATH),'run','--mode','live',
             '--fixture',str(self.fixture),'--output-dir',str(self.folder/'live')]
        p=subprocess.run(cmd,capture_output=True,text=True)
        self.assertNotEqual(p.returncode,0)
        self.assertFalse((self.folder/'live').exists())

    def test_cli_help(self):
        p=subprocess.run([sys.executable,str(SCRIPT_PATH),'--help'],capture_output=True,text=True)
        self.assertEqual(p.returncode,0);self.assertIn('--mode',p.stdout)


if __name__ == '__main__':
    unittest.main(verbosity=2)
