"""Day 27 offline-only orchestration acceptance tests (Python 3.12+)."""
from __future__ import annotations

import csv
import hashlib
from importlib import import_module
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

import day27_alpha_variant_cli as cli


class Day27PipelineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory(prefix="day27-unit-tests-")
        cls.base=Path(cls.temp.name)
        cls.output=cls.base/"run"
        cls.result=cli.run_pipeline(cls.output)
        cls.summary=cli.load_json(cls.output/"day27_pipeline_summary.json")
        cls.manifest=cli.load_json(cls.output/"day27_run_manifest.json")

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_pipeline_complete(self):
        self.assertEqual(self.result["status"],"PASS")
        self.assertEqual(self.result["pipeline_stage_count"],6)
        self.assertEqual(self.result["real_api_calls_this_run"],0)

    def test_exact_synthetic_provenance(self):
        self.assertEqual(self.summary["evidence_class"],cli.EVIDENCE)
        self.assertEqual(self.summary["assembly"],"TOY_NOT_GRCH38")
        self.assertEqual(self.summary["mode"],"OFFLINE_ONLY")
        self.assertFalse(self.summary["official_atlas_data"])
        self.assertFalse(self.summary["model_prediction"])
        self.assertFalse(self.summary["clinical_interpretation"])

    def test_stage_order_and_source(self):
        self.assertEqual([s["day"] for s in self.summary["stages"]],list(cli.STAGES))
        self.assertEqual([s["module"] for s in self.summary["stages"]],
                         [cli.SOURCE_MODULES[d]+".py" for d in cli.STAGES])

    def test_qc_counts(self):
        self.assertEqual((self.summary["candidate_count"],self.summary["selected_count"],
                          self.summary["qc_review_count"],self.summary["qc_rejected_count"]),
                         (8,3,2,3))

    def test_day22_true_shadowed_branches(self):
        report=cli.load_json(self.output/cli.RUN_FILES[22])
        probes={r["label"]:r["reason"] for r in report["supplementary_probe_cases"]}
        self.assertEqual(probes["LOW_QUAL_INSIDE_BED"],"LOW_QUAL")
        self.assertEqual(probes["REF_MISMATCH_INSIDE_BED"],"REF_MISMATCH")
        default={r["variant_id"]:r["reason"] for r in report["qc_rows"]}
        self.assertEqual(default["chr1:120:A>T"],"OUTSIDE_TEACHING_BED")
        self.assertEqual(default["chr1:140:G>T"],"OUTSIDE_TEACHING_BED")

    def test_day23_annotations_and_null(self):
        report=cli.load_json(self.output/cli.RUN_FILES[23])
        self.assertEqual((report["annotated_count"],report["missing_score_count"]),(2,1))
        self.assertIsNone(report["annotations"][-1]["synthetic_score"])
        self.assertEqual(report["annotations"][-1]["score_status"],"MISSING_SCORE")

    def test_day23_source_qc_rows_identical(self):
        d21=cli.load_json(self.output/cli.RUN_FILES[21])
        d23=cli.load_json(self.output/cli.RUN_FILES[23])
        self.assertEqual(d21["qc_rows"],d23["qc_rows"])

    def test_review_flags_all_required(self):
        d23=cli.load_json(self.output/cli.RUN_FILES[23])
        self.assertTrue(all(r["review_status"]=="REVIEW_REQUIRED" for r in d23["annotations"]))

    def test_day24_peak_and_no_score_merge(self):
        d24=cli.load_json(self.output/cli.RUN_FILES[24])
        self.assertEqual(d24["largest_absolute_delta_position_1based"],110)
        self.assertEqual(d24["largest_absolute_delta"],-0.3)
        self.assertEqual([r["delta_alt_minus_ref"] for r in d24["point_rows"]],
                         [0.0,0.25,-0.3,0.1,0.0])
        self.assertEqual(d24["day23_annotation_link"]["synthetic_score"],-0.17)
        self.assertFalse(self.summary["score_combination_performed"])

    def test_day25_links_day24_and_two_traces(self):
        d25=cli.load_json(self.output/cli.RUN_FILES[25])
        self.assertEqual(d25["point_count"],5)
        self.assertFalse(d25["day23_score_combined_with_delta"])
        tree=ET.parse(self.output/"day25/day25_locus_tracks.svg")
        traces=[el.get("data-track") for el in tree.iter() if el.tag.endswith("polyline")]
        self.assertEqual(traces,["REF","ALT"])

    def test_day26_ppm_opposite_deltas(self):
        self.assertEqual(self.summary["ppm_anchor_deltas"],{"MOTIF_A":-0.12,"MOTIF_B":0.12})
        d26=cli.load_json(self.output/cli.RUN_FILES[26])
        self.assertTrue(d26["day25_linked"])
        self.assertFalse(d26["day25_delta_combined_with_ppm"])
        self.assertFalse(d26["ppm_definition"]["log_odds_pwm_performed"])
        self.assertEqual(d26["window_rows_count"],26)

    def test_every_stage_json_is_present(self):
        self.assertTrue(all((self.output/cli.RUN_FILES[d]).is_file() for d in cli.STAGES))

    def test_all_artifact_hashes_verified(self):
        result=cli.verify_run(self.output)
        self.assertEqual(result["status"],"PASS")
        self.assertEqual(result["checksums"],"MATCH")
        self.assertGreaterEqual(result["artifact_count"],20)

    def test_stage_json_sha_matches_manifest(self):
        for entry in self.summary["stages"]:
            self.assertEqual(entry["json_sha256"],self.manifest["artifacts"][entry["json_file"]]["sha256"])

    def test_csv_six_stages(self):
        with (self.output/"day27_stages.csv").open(encoding="utf-8",newline="") as f:
            rows=list(csv.DictReader(f))
        self.assertEqual(len(rows),6)
        self.assertTrue(all(r["status"]=="PASS" for r in rows))

    def test_synthetic_visualization_is_valid_svg(self):
        tree=ET.parse(self.output/"day27_pipeline_flow.svg")
        self.assertEqual(tree.getroot().attrib["viewBox"],"0 0 1920 1080")
        self.assertIn("SYNTHETIC",(self.output/"day27_review_dashboard.html").read_text(encoding="utf-8"))

    def test_dashboard_uses_local_relative_links(self):
        h=(self.output/"day27_review_dashboard.html").read_text(encoding="utf-8")
        for day in cli.STAGES:
            self.assertIn(cli.RUN_FILES[day],h)
        self.assertNotIn("<script",h)
        self.assertNotIn("https://",h)

    def test_source_bundle_sha_matches_files(self):
        verified=cli.verify_bundle(Path(cli.__file__).parent)
        self.assertEqual(verified["files_checked"],12)
        self.assertEqual(self.summary["source_manifest_sha256"],verified["manifest_sha256"])

    def test_no_network_modules_in_cli_or_reused_source(self):
        import ast
        root=Path(cli.__file__).parent
        for mod in ["day27_alpha_variant_cli.py",*(m+".py" for m in cli.SOURCE_MODULES.values())]:
            tree=ast.parse((root/mod).read_text(encoding="utf-8"))
            modules=[]
            for n in ast.walk(tree):
                if isinstance(n,ast.Import): modules.extend(alias.name for alias in n.names)
                elif isinstance(n,ast.ImportFrom) and n.module: modules.append(n.module)
            self.assertTrue(all(x.split(".")[0] not in ("requests","httpx","urllib","socket","alphagenome","google")
                                for x in modules),mod)

    def test_original_day21_function_still_used(self):
        module=import_module(cli.SOURCE_MODULES[21])
        self.assertEqual(module.build_review_packet.__name__,"build_review_packet")
        self.assertEqual(module.QC_ORDER,["FILTER","BED","REF","INDEL","MISSING_QUAL","LOW_QUAL"])

    def test_existing_output_does_not_get_overwritten(self):
        before=(self.output/"day27_pipeline_summary.json").read_bytes()
        with self.assertRaisesRegex(ValueError,"already exists"):
            cli.run_pipeline(self.output)
        self.assertEqual((self.output/"day27_pipeline_summary.json").read_bytes(),before)

    def test_mutated_artifact_rejected(self):
        dest=self.base/"corrupted"
        shutil.copytree(self.output,dest)
        path=dest/"day26/day26_motif_scan.csv"
        path.write_bytes(path.read_bytes()+b"#corrupted")
        with self.assertRaisesRegex(ValueError,"checksum mismatch"):
            cli.verify_run(dest)

    def test_mutated_summary_rejected(self):
        dest=self.base/"corrupted_summary"
        shutil.copytree(self.output,dest)
        path=dest/"day27_pipeline_summary.json"
        obj=cli.load_json(path)
        obj["selected_count"]=999
        path.write_text(json.dumps(obj),encoding="utf-8")
        with self.assertRaisesRegex(ValueError,"checksum mismatch"):
            cli.verify_run(dest)

    def test_manifest_traversal_rejected(self):
        dest=self.base/"bad_manifest"
        shutil.copytree(self.output,dest)
        mpath=dest/"day27_run_manifest.json"
        manifest=cli.load_json(mpath)
        manifest["artifacts"]["../illegal.json"]={"sha256":"0"*64,"bytes":100}
        manifest["artifact_count"]+=1
        mpath.write_text(json.dumps(manifest),encoding="utf-8")
        with self.assertRaisesRegex(ValueError,"Manifest path/checksum"):
            cli.verify_run(dest)

    def test_bundle_drift_is_rejected(self):
        root=self.base/"drift_bundle"
        root.mkdir()
        shutil.copy2(Path(cli.__file__).parent/"DAY27_SOURCE_MANIFEST.json",root/"DAY27_SOURCE_MANIFEST.json")
        for rel in cli.load_json(root/"DAY27_SOURCE_MANIFEST.json")["files"]:
            src=Path(cli.__file__).parent/rel
            dst=root/rel
            dst.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(src,dst)
        (root/"fixtures/day21_synthetic.vcf").write_text("tampered\n",encoding="utf-8")
        with self.assertRaisesRegex(ValueError,"Bundle drift"):
            cli.verify_bundle(root)

    def test_output_missing_when_bundle_invalid(self):
        root=self.base/"drift_bundle"
        out=self.base/"blocked_new_run"
        with self.assertRaises(ValueError):
            cli.run_pipeline(out,bundle_root=root)
        self.assertFalse(out.exists())

    def test_failure_rolls_back_partial_stage_dir(self):
        root=self.base/"bad_fixture_bundle"
        shutil.copytree(Path(cli.__file__).parent,root,ignore=shutil.ignore_patterns("__pycache__","day27_output"))
        f=root/"fixtures/day21_synthetic.vcf"
        f.write_text("BAD DATA",encoding="utf-8")
        man=root/"DAY27_SOURCE_MANIFEST.json"
        obj=cli.load_json(man)
        obj["files"]["fixtures/day21_synthetic.vcf"]={"sha256":cli.file_hash(f),"bytes":f.stat().st_size}
        man.write_text(json.dumps(obj),encoding="utf-8")
        out=self.base/"bad_pipeline_output"
        with self.assertRaises(ValueError):
            cli.run_pipeline(out,bundle_root=root)
        self.assertFalse(out.exists())
        self.assertEqual(list(self.base.glob(".day27-working-*")),[])

    def test_live_not_supported_even_if_env_var_present(self):
        import os
        e=dict(os.environ)
        e["ALPHAGENOME_API_KEY"]="test-placeholder-not-a-real-key"
        proc=subprocess.run([sys.executable,str(Path(cli.__file__).resolve()),"run","--mode","live","--output-dir",str(self.base/"never")],
                            text=True,capture_output=True,cwd=self.base,env=e)
        self.assertNotEqual(proc.returncode,0)
        self.assertIn("invalid choice",proc.stderr)
        self.assertFalse((self.base/"never").exists())

    def test_no_live_subcommand_and_help(self):
        proc=subprocess.run([sys.executable,str(Path(cli.__file__).resolve()),"--help"],
                            capture_output=True,text=True,cwd=self.base,check=True)
        self.assertIn("run",proc.stdout)
        self.assertIn("verify",proc.stdout)
        self.assertNotIn("--authorize-live",proc.stdout)

    def test_cli_verify_subcommand(self):
        proc=subprocess.run([sys.executable,str(Path(cli.__file__).resolve()),"verify","--run-dir",str(self.output)],
                            text=True,capture_output=True,cwd=self.base,check=True)
        result=json.loads(proc.stdout)
        self.assertEqual(result["checksums"],"MATCH")
        self.assertEqual(result["real_api_calls_this_run"],0)

    def test_cli_runs_when_cwd_elsewhere(self):
        out=self.base/"another_run"
        proc=subprocess.run([sys.executable,str(Path(cli.__file__).resolve()),"run","--mode","offline","--output-dir",str(out)],
                            text=True,capture_output=True,cwd=self.base,check=True)
        result=json.loads(proc.stdout)
        self.assertEqual(result["status"],"PASS")
        self.assertTrue((out/"day27_pipeline_summary.json").is_file())

    def test_run_csv_json_remain_synthetic(self):
        for day in cli.STAGES:
            doc=cli.load_json(self.output/cli.RUN_FILES[day])
            self.assertNotEqual(doc.get("assembly"),"GRCh38")
            self.assertEqual(doc["real_api_calls_this_run"],0)

    def test_duplicate_json_key_rejected(self):
        p=self.base/"duplicate.json"
        p.write_text('{"a":1,"a":2}',encoding="utf-8")
        with self.assertRaisesRegex(ValueError,"Duplicate"):
            cli.load_json(p)


if __name__=="__main__":
    unittest.main()
