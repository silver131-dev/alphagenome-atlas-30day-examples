"""Offline synthetic fixtures: behavior, provenance, and fail-closed tampering."""
from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from fixture_helpers import FixtureError, load_verified_fixture  # noqa: E402
from day17_historical_avi_reader import build_report as avi_report  # noqa: E402
from day18_historical_splicing_reader import build_report as splice_report  # noqa: E402


class SyntheticReadersTest(unittest.TestCase):
    def test_day17_fictional_avi(self) -> None:
        report = avi_report()
        self.assertEqual(report['status'], 'SYNTHETIC_EDUCATIONAL_DEMO')
        self.assertEqual(report['avi_interpretation']['fictional_raw_score'], 0.85)
        self.assertEqual(report['avi_interpretation']['fictional_quantile'], 0.99)
        self.assertEqual(report['avi_interpretation']['fictional_phred_like'], 20.0)
        self.assertEqual(report['avi_interpretation']['percentile_percent'], 99)
        self.assertEqual(report['avi_interpretation']['illustrative_feature_sum'], 0.75)
        self.assertEqual(report['avi_interpretation']['illustrative_score_minus_feature_sum'], 0.1)
        self.assertEqual(report['avi_interpretation']['toy_additivity_review'], 'REVIEW')

    def test_day18_splicing_delta_direction(self) -> None:
        report = splice_report()
        sites = report['splice_site_ref_alt_comparison']
        self.assertEqual(sites[0]['delta_alt_minus_ref'], -0.7)
        self.assertEqual(sites[0]['toy_direction'], 'decrease')
        self.assertEqual(sites[1]['delta_alt_minus_ref'], 0.25)
        self.assertEqual(sites[1]['toy_direction'], 'increase')
        self.assertEqual(report['rna_like_ref_alt_comparison'][0]['delta_alt_minus_ref'], -3.0)

    def test_both_days_share_same_fictional_variant(self) -> None:
        self.assertEqual(avi_report()['variant'], splice_report()['variant'])
        self.assertEqual(avi_report()['variant']['chromosome'], 'chrToy')

    def test_source_provenance(self) -> None:
        for r in (avi_report(), splice_report()):
            self.assertEqual(r['provenance']['source_class'], 'SYNTHETIC_FIXTURE')
            self.assertEqual(r['provenance']['real_api_calls_this_run'], 0)
            self.assertEqual(r['provenance']['integrity_status'], 'FIXTURE_SHA256_VERIFIED')
            self.assertFalse(r['provenance']['real_model_predictions'])
            self.assertEqual(r['provenance']['reference_data_authenticity'], 'NOT_CLAIMED')

    def _copy(self, tmp: str) -> Path:
        folder = Path(tmp)
        for path in HERE.glob('*.json'):
            shutil.copy2(path, folder / path.name)
        return folder

    def _rewrite_manifest_hash(self, folder: Path, file_name: str) -> None:
        man = folder / 'SYNTHETIC_FIXTURE_MANIFEST.json'
        data = json.loads(man.read_text(encoding='utf-8'))
        data['files'][file_name] = hashlib.sha256((folder / file_name).read_bytes()).hexdigest()
        man.write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')

    def test_tampering_fails_closed_before_read(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            folder = self._copy(tmp)
            path = folder / 'day17_synthetic_avi_fixture.json'
            path.write_bytes(path.read_bytes() + b' ')
            with self.assertRaisesRegex(FixtureError, 'FIXTURE_SHA256_MISMATCH'):
                load_verified_fixture(path.name, 17, folder)

    def test_non_synthetic_source_fails_even_with_new_hash(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            folder = self._copy(tmp)
            path = folder / 'day18_synthetic_splicing_fixture.json'
            data = json.loads(path.read_text(encoding='utf-8'))
            data['source']['real_model_output'] = True
            path.write_text(json.dumps(data), encoding='utf-8')
            self._rewrite_manifest_hash(folder, path.name)
            with self.assertRaisesRegex(FixtureError, 'Source declaration invalid'):
                splice_report(folder)

    def test_invalid_avi_quantile_fails_even_with_new_hash(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            folder = self._copy(tmp)
            path = folder / 'day17_synthetic_avi_fixture.json'
            data = json.loads(path.read_text(encoding='utf-8'))
            data['avi']['quantile'] = 1.0
            path.write_text(json.dumps(data), encoding='utf-8')
            self._rewrite_manifest_hash(folder, path.name)
            with self.assertRaisesRegex(FixtureError, 'quantile must be strictly'):
                avi_report(folder)

    def test_mismatched_ref_alt_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            folder = self._copy(tmp)
            path = folder / 'day18_synthetic_splicing_fixture.json'
            data = json.loads(path.read_text(encoding='utf-8'))
            data['variant']['alt'] = 'T'
            path.write_text(json.dumps(data), encoding='utf-8')
            self._rewrite_manifest_hash(folder, path.name)
            with self.assertRaisesRegex(FixtureError, 'identities differ'):
                splice_report(folder)

    def test_bad_position_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            folder = self._copy(tmp)
            path = folder / 'day18_synthetic_splicing_fixture.json'
            data = json.loads(path.read_text(encoding='utf-8'))
            data['splice_sites'][0]['position_1based'] = -1
            path.write_text(json.dumps(data), encoding='utf-8')
            self._rewrite_manifest_hash(folder, path.name)
            with self.assertRaisesRegex(FixtureError, 'Invalid 1-based'):
                splice_report(folder)

    def test_no_network_or_model_imports(self) -> None:
        prohibited = ('socket', 'requests', 'httpx', 'alphagenome', 'google', 'urllib', 'aiohttp', 'subprocess')
        for name in ['fixture_helpers.py', 'day17_historical_avi_reader.py', 'day18_historical_splicing_reader.py']:
            tree = ast.parse((HERE / name).read_text(encoding='utf-8'))
            imports = []
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    imports.extend(alias.name for alias in node.names)
                if isinstance(node, ast.ImportFrom) and node.module:
                    imports.append(node.module)
            self.assertFalse([x for x in imports if x.split('.')[0] in prohibited], name)

    def test_cli_runs_without_working_directory_dependency(self) -> None:
        for day in ('17', '18'):
            name = 'day17_historical_avi_reader.py' if day == '17' else 'day18_historical_splicing_reader.py'
            with tempfile.TemporaryDirectory() as work:
                result = subprocess.run([sys.executable, str(HERE / name), '--verify-only'],
                                        cwd=work, capture_output=True, text=True, timeout=10)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn(f'DAY{day}_FIXTURE_QA=PASS', result.stdout)

    def test_report_export(self) -> None:
        for name, stem in [('day17_historical_avi_reader.py', 'day17_synthetic_avi_report'),
                           ('day18_historical_splicing_reader.py', 'day18_synthetic_splicing_report')]:
            with tempfile.TemporaryDirectory() as tmp:
                result = subprocess.run([sys.executable, str(HERE / name), '--output-dir', tmp],
                                        capture_output=True, text=True, timeout=10)
                self.assertEqual(result.returncode, 0, result.stderr)
                payload = json.loads((Path(tmp) / f'{stem}.json').read_text())
                self.assertEqual(payload['provenance']['real_api_calls_this_run'], 0)
                self.assertIn('SYNTHETIC EDUCATIONAL FIXTURE', (Path(tmp) / f'{stem}.md').read_text())

    def test_no_old_historical_values_or_coordinates(self) -> None:
        # Never embed original model-derived decimals/digests in public tests either.
        import re
        corpus = "\n".join(p.read_text(encoding="utf-8") for p in HERE.glob("*.json"))
        for forbidden in ("DERIVED_FROM_HISTORICAL_OUTPUT", "GRCh38", "chr22:",
                          "UBERON:", "source_native_api_response", "model_executed_true"):
            self.assertNotIn(forbidden, corpus)
        self.assertFalse(re.search(r"\b[0-9]+\.[0-9]{9,}\b", corpus),
                         "Unexpected overly precise value in hand-authored toy fixture")

    def test_unreviewed_extra_field_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            folder = self._copy(tmp)
            p = folder / "day18_synthetic_splicing_fixture.json"
            obj = json.loads(p.read_text(encoding="utf-8"))
            obj["source"]["real_ontology"] = "not_allowed"
            p.write_text(json.dumps(obj), encoding="utf-8")
            self._rewrite_manifest_hash(folder, p.name)
            with self.assertRaisesRegex(FixtureError, "Unexpected synthetic source metadata"):
                splice_report(folder)

if __name__ == '__main__':
    unittest.main()
