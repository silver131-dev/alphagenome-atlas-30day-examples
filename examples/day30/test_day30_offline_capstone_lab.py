"""Day30 offline evidence chain tests: standard library, no API or SDK."""
from __future__ import annotations

import copy
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import day30_offline_capstone_lab as d

ROOT=Path(__file__).resolve().parent
D27=ROOT/"example_output/day27_output"
D28=ROOT/"example_output/day28_output"
D29=ROOT/"example_output/day29_output"
TITLES=ROOT/"author_approved_day01_day28_titles.csv"


class EvidenceSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=d.verify_source_chain(D27,D28,D29)

    def test_source_manifest_counts(self):
        self.assertEqual(self.r["source_artifact_counts"],{"day27":23,"day28":5,"day29":6,"total":34})

    def test_six_pipeline_stages(self):
        self.assertEqual(self.r["day27_pipeline_stages"],6)

    def test_qc_counts(self):
        self.assertEqual([self.r[x] for x in ("candidate_count","qc_selected","qc_review_required","qc_rejected")],[8,3,2,3])

    def test_annotations_and_missing(self):
        self.assertEqual([self.r["annotated_count"],self.r["missing_score_count"]],[2,1])

    def test_false_live(self):
        self.assertFalse(self.r["models_called_in_this_capstone"])
        self.assertEqual(self.r["real_api_calls_this_run"],0)

    def test_historical_live_not_imported(self):
        self.assertFalse(self.r["historical_live_evidence_imported"])

    def test_reidentification_risk_not_eliminated(self):
        self.assertFalse(self.r["reidentification_risk_eliminated"])

    def test_signature_not_verified(self):
        self.assertFalse(self.r["source_manifest_digital_signatures_verified"])

    def test_assembly_is_toy(self):
        self.assertEqual(self.r["genome_assembly"],"TOY_NOT_GRCH38")

    def test_publication_blocked(self):
        self.assertEqual(self.r["privacy_release_gate"],"BLOCKED")
        self.assertFalse(self.r["external_transfer_allowed"])
        self.assertFalse(self.r["human_approval_completed"])

    def test_human_review_pending(self):
        self.assertEqual(self.r["manual_review"],"REQUIRED")

    def test_report_does_not_contain_snp_or_private_token(self):
        s=json.dumps(self.r,ensure_ascii=False)
        for val in ("chr1:110", "chr22:362", "DEMO_", "local_demo_token", "patient_id"):
            self.assertNotIn(val,s)

    def test_hash_links_day27_day28(self):
        src=d.read_json(D28/"day28_evidence_report.json")["source_provenance"]
        self.assertEqual(src["day27_run_manifest_sha256"],d.digest(D27/"day27_run_manifest.json"))

    def test_hash_links_day28_day29(self):
        src=d.read_json(D29/"day29_privacy_audit.json")["source"]
        self.assertEqual(src["day28_manifest_sha256"],d.digest(D28/"day28_report_manifest.json"))

    def test_hash_digest_matches_sha256lib(self):
        content=(D29/"day29_privacy_audit.json").read_bytes()
        self.assertEqual(d.digest(D29/"day29_privacy_audit.json"),hashlib.sha256(content).hexdigest())

    def test_day27_source_manifest_23(self):
        m,h=d.verify_manifest(D27,"day27_run_manifest.json",23)
        self.assertEqual(len(m["artifacts"]),23)
        self.assertEqual(len(h),64)

    def test_day28_source_manifest_5(self):
        m,_=d.verify_manifest(D28,"day28_report_manifest.json",5)
        self.assertIs(m["signed"],False)

    def test_day29_source_manifest_6(self):
        m,_=d.verify_manifest(D29,"day29_report_manifest.json",6)
        self.assertIs(m["signed"],False)

    def test_day29_preview_is_aggregate_only(self):
        preview=d.read_json(D29/"day29_local_preview.json")
        self.assertIs(preview["per_variant_data_included"],False)
        self.assertIs(preview["random_token_included"],False)

    def test_day28_score_null_untouched(self):
        ann=d.read_json(D28/"day28_evidence_report.json")["annotations"]
        self.assertEqual(sum(item["synthetic_score"] is None for item in ann),1)

    def test_minimal_python_import_surface(self):
        source=(ROOT/"day30_offline_capstone_lab.py").read_text(encoding="utf-8")
        for dangerous in ("import alphagenome","import requests","import urllib", "import socket", "http.client", "get_variant_annotation("):
            self.assertNotIn(dangerous,source)


class TitleAndVisualizationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=d.verify_source_chain(D27,D28,D29)
        cls.titles=d.read_approved_titles(TITLES)

    def test_28_author_titles(self):
        self.assertEqual(len(self.titles),28)

    def test_title_day01_matches_author(self):
        self.assertTrue(self.titles[0]["approved_title"].startswith("破解 98% 基因暗物質"))

    def test_title_day28_matches_author(self):
        self.assertIn("學會說清楚我們不知道什麼",self.titles[-1]["approved_title"])

    def test_30_day_index_complete(self):
        m=d.make_map(self.r,self.titles)
        self.assertEqual(sum(line.startswith("| "+f"{i:02d}+"[:-1]) for i in range(1,31) for line in m.splitlines()),30)

    def test_map_has_day29_and_day30(self):
        m=d.make_map(self.r,self.titles)
        self.assertIn("| 29 | 前瞻展望",m)
        self.assertIn("| 30 | 30 天實戰總結",m)

    def test_map_clearly_separates_live(self):
        m=d.make_map(self.r,self.titles)
        self.assertIn("HISTORICAL LIVE",m)
        self.assertIn("CURRENT LIVE",m)
        self.assertIn("0.10.0",m)

    def test_svg_16_by_9(self):
        svg=d.make_svg(self.r)
        self.assertIn('width="1920" height="1080"',svg)
        self.assertIn('fill="#E5DEDA"',svg)

    def test_svg_contains_source_totals(self):
        svg=d.make_svg(self.r)
        self.assertIn("共 34 個來源檔",svg)
        self.assertIn("PRIVATE",svg)

    def test_html_standalone_no_external_url(self):
        html=d.make_html(self.r,d.make_svg(self.r))
        self.assertIn("<svg",html)
        self.assertNotIn('src="https://',html)
        self.assertNotIn('<script',html)

    def test_checklist_release_unchecked(self):
        ch=d.make_checklist(self.r)
        self.assertIn("- [ ] GitHub 公開",ch)
        self.assertIn("- [x]",ch)


class LifecycleTests(unittest.TestCase):
    def test_run_verify_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/"new"
            result=d.create_report(out,D27,D28,D29,TITLES)
            self.assertEqual(result["status"],"PASS")
            self.assertEqual(d.verify_report(out,D27,D28,D29,TITLES)["checksums"],"MATCH")

    def test_all_output_files_expected(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/"out"
            d.create_report(out,D27,D28,D29,TITLES)
            self.assertEqual(len(list(out.glob("*"))),6)
            self.assertTrue((out/"day30_capstone_dashboard.html").is_file())

    def test_existing_output_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(d.EvidenceError):
                d.create_report(Path(tmp),D27,D28,D29,TITLES)

    def test_deterministic_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            a,b=Path(tmp)/"a",Path(tmp)/"b"
            d.create_report(a,D27,D28,D29,TITLES)
            d.create_report(b,D27,D28,D29,TITLES)
            self.assertEqual({x.name:d.digest(x) for x in a.iterdir()},
                             {x.name:d.digest(x) for x in b.iterdir()})

    def test_manifest_refuses_modified_dashboard(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/"out"
            d.create_report(out,D27,D28,D29,TITLES)
            with (out/"day30_capstone_dashboard.html").open("a") as f:f.write("TAMPER")
            with self.assertRaises(d.EvidenceError):
                d.verify_report(out,D27,D28,D29,TITLES)

    def test_title_mutation_invalidates_verification(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/"out"
            d.create_report(out,D27,D28,D29,TITLES)
            original=TITLES.read_text(encoding="utf-8-sig")
            mut=Path(tmp)/"titles.csv";mut.write_text(original.replace("破解 98%","破解 97%"),encoding="utf-8")
            with self.assertRaises(d.EvidenceError):
                d.verify_report(out,D27,D28,D29,mut)

    def test_malformed_approved_title_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"bad.csv"
            p.write_text("day,approved_title,title_status\n01,bad,AUTHOR_APPROVED\n")
            with self.assertRaises(d.EvidenceError):d.read_approved_titles(p)

    def test_run_cli_offline(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/"cli"
            command=[sys.executable,str(ROOT/"day30_offline_capstone_lab.py"),"run","--mode","offline",
              "--day27-dir",str(D27),"--day28-dir",str(D28),"--day29-dir",str(D29),
              "--titles-file",str(TITLES),"--output-dir",str(out)]
            r=subprocess.run(command,text=True,capture_output=True)
            self.assertEqual(r.returncode,0,r.stderr)
            self.assertEqual(json.loads(r.stdout)["status"],"PASS")
            v=subprocess.run([sys.executable,str(ROOT/"day30_offline_capstone_lab.py"),"verify",
              "--day27-dir",str(D27),"--day28-dir",str(D28),"--day29-dir",str(D29),
              "--titles-file",str(TITLES),"--report-dir",str(out)],text=True,capture_output=True)
            self.assertEqual(v.returncode,0,v.stderr)
            self.assertEqual(json.loads(v.stdout)["checksums"],"MATCH")

    def test_live_mode_is_rejected_with_env_key_set(self):
        with tempfile.TemporaryDirectory() as tmp:
            env=os.environ.copy();env["ALPHAGENOME_API_KEY"]="THIS_IS_NOT_A_KEY"
            r=subprocess.run([sys.executable,str(ROOT/"day30_offline_capstone_lab.py"),
               "run","--mode","live","--output-dir",str(Path(tmp)/"out")],text=True,capture_output=True,env=env)
            self.assertNotEqual(r.returncode,0)
            self.assertFalse((Path(tmp)/"out").exists())
            self.assertNotIn("THIS_IS_NOT_A_KEY",r.stderr+r.stdout)

    def test_nonexistent_manifest_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(d.EvidenceError):
                d.verify_source_chain(Path(tmp),D28,D29)


class TamperingTests(unittest.TestCase):
    def duplicate_source(self):
        tmp=tempfile.TemporaryDirectory()
        base=Path(tmp.name)
        for name,src in (("d27",D27),("d28",D28),("d29",D29)):
            shutil.copytree(src,base/name)
        return tmp,base

    def test_day27_data_modification_refused(self):
        ctx,base=self.duplicate_source()
        with ctx:
            f=base/"d27/day22/day22_qc_report.json"
            f.write_text(f.read_text()+"\n",encoding="utf-8")
            with self.assertRaises(d.EvidenceError):d.verify_source_chain(base/"d27",base/"d28",base/"d29")

    def test_day28_data_modification_refused(self):
        ctx,base=self.duplicate_source()
        with ctx:
            f=base/"d28/day28_evidence_report.json"
            f.write_text(f.read_text()+"\n",encoding="utf-8")
            with self.assertRaises(d.EvidenceError):d.verify_source_chain(base/"d27",base/"d28",base/"d29")

    def test_day29_data_modification_refused(self):
        ctx,base=self.duplicate_source()
        with ctx:
            f=base/"d29/day29_privacy_audit.json"
            f.write_text(f.read_text()+"\n",encoding="utf-8")
            with self.assertRaises(d.EvidenceError):d.verify_source_chain(base/"d27",base/"d28",base/"d29")

    def test_day27_file_deleted_refused(self):
        ctx,base=self.duplicate_source()
        with ctx:
            (base/"d27/day25/day25_locus_tracks.svg").unlink()
            with self.assertRaises(d.EvidenceError):d.verify_source_chain(base/"d27",base/"d28",base/"d29")

    def test_day29_bad_claim_even_with_manifest_updated(self):
        ctx,base=self.duplicate_source()
        with ctx:
            p=base/"d29/day29_privacy_audit.json";x=json.loads(p.read_text())
            x["release_gate"]="APPROVED"
            p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
            m=base/"d29/day29_report_manifest.json";man=json.loads(m.read_text())
            man["artifacts"][p.name]["sha256"]=d.digest(p)
            man["artifacts"][p.name]["bytes"]=p.stat().st_size
            m.write_text(json.dumps(man,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
            with self.assertRaisesRegex(d.EvidenceError,"release_gate"):
                d.verify_source_chain(base/"d27",base/"d28",base/"d29")

    def test_manifest_path_traversal_refused(self):
        ctx,base=self.duplicate_source()
        with ctx:
            mp=base/"d29/day29_report_manifest.json";m=json.loads(mp.read_text())
            detail=m["artifacts"].pop("day29_release_decision.csv")
            m["artifacts"]["../day29_release_decision.csv"]=detail
            mp.write_text(json.dumps(m,ensure_ascii=False)+"\n",encoding="utf-8")
            with self.assertRaises(d.EvidenceError):d.verify_source_chain(base/"d27",base/"d28",base/"d29")

    def test_manifest_symlink_artifact_refused(self):
        ctx,base=self.duplicate_source()
        with ctx:
            p=base/"d29/day29_privacy_dashboard.html"
            copy=p.read_bytes();p.unlink()
            q=base/"d29/independent_copy.html";q.write_bytes(copy)
            p.symlink_to(q)
            with self.assertRaises(d.EvidenceError):d.verify_source_chain(base/"d27",base/"d28",base/"d29")

    def test_day29_private_preview_tampered_and_rehashed_refused(self):
        ctx,base=self.duplicate_source()
        with ctx:
            p=base/"d29/day29_local_preview.json";v=json.loads(p.read_text())
            v["random_token_included"]=True
            p.write_text(json.dumps(v,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
            m=base/"d29/day29_report_manifest.json";man=json.loads(m.read_text())
            man["artifacts"][p.name]={"bytes":p.stat().st_size,"sha256":d.digest(p)}
            m.write_text(json.dumps(man,ensure_ascii=False)+"\n",encoding="utf-8")
            with self.assertRaisesRegex(d.EvidenceError,"random_token_included"):
                d.verify_source_chain(base/"d27",base/"d28",base/"d29")


if __name__=="__main__":unittest.main()
