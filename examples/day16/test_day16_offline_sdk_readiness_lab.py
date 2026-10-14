"""Day 16 offline-only tests: no network, no API keys, no installed SDK needed."""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import types
import unittest
from unittest import mock

import day16_offline_sdk_readiness_lab as lab


def report(**updates):
    data = {
        "installed_sdk": "0.10.0", "key_names_present": set(),
        "python_version": (3, 12, 3), "in_virtualenv": True,
    }
    data.update(updates)
    return lab.build_report(**data)


class Day16OfflineContractTests(unittest.TestCase):
    def test_expected_sdk_constant(self):
        self.assertEqual(lab.EXPECTED_SDK, "0.10.0")

    def test_supported_python(self):
        self.assertTrue(report()["python_supported"])

    def test_old_unsupported_python(self):
        r = report(python_version=(3, 9, 19))
        self.assertEqual(r["status"], "REVIEW_REQUIRED")
        self.assertFalse(r["python_supported"])

    def test_310_supported(self):
        self.assertTrue(report(python_version=(3, 10, 0))["python_supported"])

    def test_python_version_string(self):
        self.assertEqual(report()["python_version"], "3.12.3")

    def test_reject_bad_version(self):
        with self.assertRaises(ValueError):
            report(python_version=(3, -1, 0))

    def test_reject_bad_version_length(self):
        with self.assertRaises(ValueError):
            report(python_version=(3, 12))

    def test_virtualenv_true(self):
        self.assertTrue(report()["virtualenv_active"])

    def test_virtualenv_false_not_failure(self):
        r = report(in_virtualenv=False)
        self.assertFalse(r["virtualenv_active"])
        self.assertEqual(r["status"], "PASS")

    def test_sdk_match(self):
        self.assertEqual(report()["sdk"]["status"], "VERSION_MATCH")

    def test_sdk_different(self):
        self.assertEqual(report(installed_sdk="0.9.0")["sdk"]["status"], "VERSION_DIFFERS")

    def test_sdk_missing(self):
        with mock.patch.object(lab, "sdk_distribution_version", return_value=None):
            r = report(installed_sdk=None)
        self.assertEqual(r["sdk"]["status"], "NOT_INSTALLED")
        self.assertIsNone(r["sdk"]["installed_version"])

    def test_no_key(self):
        s = report()["api_key_presence_only"]
        self.assertEqual(set(s), set(lab.KEY_VARIABLES))
        self.assertTrue(all(v == "NOT_SET_IN_PROCESS" for v in s.values()))

    def test_alpha_genome_key_presence_only(self):
        r = report(key_names_present={"ALPHA_GENOME_API_KEY"})
        self.assertEqual(r["api_key_presence_only"]["ALPHA_GENOME_API_KEY"], "PRESENT_NOT_VALIDATED")
        self.assertEqual(r["api_key_presence_only"]["ALPHAGENOME_API_KEY"], "NOT_SET_IN_PROCESS")

    def test_alphagenome_key_presence_only(self):
        r = report(key_names_present={"ALPHAGENOME_API_KEY"})
        self.assertEqual(r["api_key_presence_only"]["ALPHAGENOME_API_KEY"], "PRESENT_NOT_VALIDATED")

    def test_both_keys_presence_only(self):
        r = report(key_names_present=set(lab.KEY_VARIABLES))
        self.assertEqual(list(r["api_key_presence_only"].values()), ["PRESENT_NOT_VALIDATED"] * 2)

    def test_no_key_value_in_serialized_report(self):
        secret = "TOPSECRET-DO-NOT-PUBLISH-ABC123"
        with mock.patch.dict(os.environ, {"ALPHAGENOME_API_KEY": secret}):
            with mock.patch.object(lab, "sdk_distribution_version", return_value=None):
                text = json.dumps(lab.build_report())
        self.assertNotIn(secret, text)

    def test_sdk_presence_not_live_auth(self):
        self.assertFalse(report()["api_key_validated"])
        self.assertEqual(report()["current_live_access"], "NOT_TESTED")

    def test_no_client_created(self):
        self.assertIs(report()["api_client_created"], False)

    def test_no_live_calls(self):
        self.assertEqual(report()["real_api_calls_this_run"], 0)

    def test_no_model_results(self):
        self.assertFalse(report()["model_or_atlas_response_obtained"])

    def test_no_git_visibility_claim(self):
        self.assertEqual(report()["git_visibility_check"], "NOT_PERFORMED")

    def test_offline_mode(self):
        self.assertEqual(report()["mode"], "OFFLINE_ONLY")

    def test_offline_never_imports_sdk(self):
        with mock.patch.object(lab, "api_interface_snapshot", side_effect=AssertionError("SDK inspection called")):
            r = report(mode="offline")
        self.assertEqual(r["real_api_calls_this_run"], 0)

    def test_reject_live_mode(self):
        with self.assertRaises(ValueError):
            report(mode="live")

    def test_reject_empty_expected_version(self):
        with self.assertRaises(ValueError):
            report(expected_sdk="")

    def test_inspect_only_with_injected_snapshot(self):
        x = {"inspection": "LOCAL_IMPORT_ONLY", "all_checked_exports_available": True}
        r = report(mode="inspect-sdk", inspected_api=x)
        self.assertEqual(r["sdk_interface_check"], x)
        self.assertEqual(r["real_api_calls_this_run"], 0)

    def test_inspect_missing_sdk_requires_review(self):
        with mock.patch.object(lab, "sdk_distribution_version", return_value=None):
            r = report(mode="inspect-sdk", installed_sdk=None,
                       inspected_api={"all_checked_exports_available": True})
        self.assertEqual(r["status"], "REVIEW_REQUIRED")

    def test_inspect_missing_interface_requires_review(self):
        r = report(mode="inspect-sdk", inspected_api={"all_checked_exports_available": False})
        self.assertEqual(r["status"], "REVIEW_REQUIRED")

    def test_inspect_different_sdk_requires_review(self):
        r = report(mode="inspect-sdk", installed_sdk="0.9.0",
                   inspected_api={"all_checked_exports_available": True})
        self.assertEqual(r["status"], "REVIEW_REQUIRED")


class Day16SDKStaticInspectionTests(unittest.TestCase):
    def test_import_inspects_without_invoking_factories(self):
        create_atlas = mock.Mock()
        create_model = mock.Mock()
        a = types.SimpleNamespace(create=create_atlas)
        m = types.SimpleNamespace(create=create_model, OutputType=types.SimpleNamespace(RNA_SEQ="RNA_SEQ"))
        v = types.SimpleNamespace(Variant=types.SimpleNamespace(from_str=lambda x: x))
        mappable = {"alphagenome.atlas.atlas": a, "alphagenome.models.dna_client": m, "alphagenome.data.genome": v}
        with mock.patch.object(lab.importlib, "import_module", side_effect=lambda name: mappable[name]):
            status = lab.api_interface_snapshot()
        self.assertTrue(status["all_checked_exports_available"])
        create_atlas.assert_not_called()
        create_model.assert_not_called()

    def test_missing_one_interface_reported(self):
        a = types.SimpleNamespace(create=lambda x: x)
        m = types.SimpleNamespace(create=lambda x: x, OutputType=types.SimpleNamespace())
        v = types.SimpleNamespace(Variant=types.SimpleNamespace(from_str=lambda x: x))
        mods = [a, m, v]
        with mock.patch.object(lab.importlib, "import_module", side_effect=mods):
            status = lab.api_interface_snapshot()
        self.assertFalse(status["model_rna_seq_enum_available"])
        self.assertFalse(status["all_checked_exports_available"])

    def test_import_error_opaque_no_secrets(self):
        with mock.patch.object(lab.importlib, "import_module", side_effect=ImportError("secret!")):
            status = lab.api_interface_snapshot()
        self.assertEqual(status, {"inspection": "UNAVAILABLE", "error_class": "ImportError"})
        self.assertNotIn("secret!", json.dumps(status))

    def test_import_os_error_opaque(self):
        with mock.patch.object(lab.importlib, "import_module", side_effect=OSError("a private path")):
            status = lab.api_interface_snapshot()
        self.assertEqual(status["error_class"], "OSError")


class Day16ArtifactTests(unittest.TestCase):
    def test_html_is_standalone(self):
        h = lab.render_html(report())
        self.assertIn("<!doctype html>", h)
        self.assertIn("#E5DEDA", h)
        self.assertIn("LOCAL ONLY", h)

    def test_html_escape_untrusted_version(self):
        h = lab.render_html(report(installed_sdk="<script>alert('x')</script>"))
        self.assertNotIn("<script>", h)
        self.assertIn("&lt;script&gt;", h)

    def test_html_no_key_value(self):
        secret = "VERY_PRIVATE_KEY_123456"
        with mock.patch.dict(os.environ, {"ALPHAGENOME_API_KEY": secret}):
            with mock.patch.object(lab, "sdk_distribution_version", return_value=None):
                h = lab.render_html(lab.build_report())
        self.assertNotIn(secret, h)
        self.assertIn("PRESENT_NOT_VALIDATED", h)

    def test_json_html_written(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "new_output"
            files = lab.write_outputs(out, report())
            self.assertTrue(Path(files["json_file"]).is_file())
            self.assertTrue(Path(files["html_file"]).is_file())
            self.assertEqual(json.loads(Path(files["json_file"]).read_text())["status"], "PASS")

    def test_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "existing"
            p.mkdir()
            (p / "important.txt").write_text("KEEP")
            with self.assertRaises(FileExistsError):
                lab.write_outputs(p, report())
            self.assertEqual((p / "important.txt").read_text(), "KEEP")

    def test_cli_outputs_expected_files(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "out"
            with mock.patch.object(lab, "sdk_distribution_version", return_value=None):
                res = lab.main(["--mode", "offline", "--output-dir", str(out)])
            self.assertEqual(res, 0)
            self.assertTrue((out / lab.REPORT_NAME).exists())
            self.assertTrue((out / lab.HTML_NAME).exists())

    def test_cli_rejects_live(self):
        with self.assertRaises(SystemExit) as exc:
            lab.main(["--mode", "live"])
        self.assertEqual(exc.exception.code, 2)

    def test_cli_rejects_second_run_same_folder(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "out"
            with mock.patch.object(lab, "sdk_distribution_version", return_value=None):
                lab.main(["--output-dir", str(out)])
                with self.assertRaises(SystemExit) as cm:
                    lab.main(["--output-dir", str(out)])
            self.assertEqual(cm.exception.code, 2)

    def test_cli_subprocess_no_network_and_no_key_leak(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "out"
            env = {**os.environ, "ALPHAGENOME_API_KEY": "SECRET_DAY16_SENTINEL"}
            res = subprocess.run(
                [sys.executable, "-B", str(Path(lab.__file__)), "--mode", "offline", "--output-dir", str(out)],
                check=False, capture_output=True, text=True, env=env, timeout=10,
            )
            self.assertEqual(res.returncode, 0, res.stderr)
            self.assertNotIn("SECRET_DAY16_SENTINEL", res.stdout)
            self.assertNotIn("SECRET_DAY16_SENTINEL", (out / lab.REPORT_NAME).read_text())
            self.assertFalse(json.loads((out / lab.REPORT_NAME).read_text())["api_key_validated"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
