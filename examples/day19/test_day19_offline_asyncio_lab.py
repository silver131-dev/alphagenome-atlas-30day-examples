"""All tests run OFFLINE using synthetic request/response scripts."""
import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

from day19_offline_asyncio_lab import (
    FakeAPI, delay_for, main, query_with_retry, run_batch,
)
import asyncio


async def no_wait(_seconds):
    """Tests record delays without spending time sleeping between attempts."""
    await asyncio.sleep(0)


class Day19AsyncTests(unittest.IsolatedAsyncioTestCase):
    async def test_default_summary(self):
        result = await run_batch(sleeper=no_wait)
        self.assertEqual((result["scenario_count"], result["success_count"],
                          result["failure_count"], result["total_attempts"]),
                         (5, 3, 2, 11))
        self.assertEqual(result["real_api_calls_this_run"], 0)
        self.assertEqual(result["evidence_class"], "SYNTHETIC_FIXTURE")
        self.assertAlmostEqual(result["sum_scheduled_retry_wait_s"], 0.10)

    async def test_429_respects_simulated_retry_after(self):
        result = await run_batch(sleeper=no_wait)
        row = result["results"][0]
        self.assertEqual(row["history"], [429, 200])
        self.assertEqual(row["retry_wait_s"], [0.03])
        self.assertEqual(row["status"], "OK")

    async def test_503_exponential_backoff(self):
        row = (await run_batch(sleeper=no_wait))["results"][1]
        self.assertEqual(row["history"], [503, 503, 200])
        self.assertEqual(row["retry_wait_s"], [0.01, 0.02])

    async def test_bad_request_not_retried(self):
        row = (await run_batch(sleeper=no_wait))["results"][2]
        self.assertEqual((row["status"], row["attempts"]), ("NOT_RETRYABLE", 1))
        self.assertEqual(row["retry_wait_s"], [])

    async def test_retry_limit(self):
        row = (await run_batch(sleeper=no_wait))["results"][3]
        self.assertEqual((row["status"], row["attempts"]), ("RETRY_EXHAUSTED", 3))
        self.assertEqual(row["history"], [503, 503, 503])

    async def test_timeout_retry(self):
        row = (await run_batch(sleeper=no_wait))["results"][4]
        self.assertEqual(row["history"], ["TIMEOUT", 200])
        self.assertEqual(row["synthetic_score"], 0.07)

    async def test_semaphore_limits_concurrent_mock_calls(self):
        result = await run_batch(max_concurrency=2, sleeper=no_wait)
        self.assertTrue(1 <= result["observed_max_mock_in_flight"] <= 2)
        result = await run_batch(max_concurrency=1, sleeper=no_wait)
        self.assertEqual(result["observed_max_mock_in_flight"], 1)

    async def test_401_not_retried(self):
        api = FakeAPI({"chrToy:999:A>C": [{"code": 401}, {"code": 200, "synthetic_score": 0.2}]})
        row = await query_with_retry("chrToy:999:A>C", api, asyncio.Semaphore(1), sleeper=no_wait)
        self.assertEqual(row["status"], "NOT_RETRYABLE")
        self.assertEqual(row["attempts"], 1)

    async def test_invalid_config(self):
        with self.assertRaises(ValueError):
            await run_batch(max_concurrency=0, sleeper=no_wait)
        with self.assertRaises(ValueError):
            await run_batch(max_attempts=0, sleeper=no_wait)

    def test_backoff_cap(self):
        self.assertEqual(delay_for(503, 5), 0.08)
        self.assertEqual(delay_for(429, 1, 0.07), 0.07)


class Day19CliTests(unittest.TestCase):
    def test_cli_writes_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            stream = io.StringIO()
            with contextlib.redirect_stdout(stream):
                main(["--output-dir", tmp])
            summary = json.loads(stream.getvalue())
            payload = json.loads((Path(tmp) / "day19_mock_batch.json").read_text())
            self.assertEqual(summary["real_api_calls_this_run"], 0)
            self.assertEqual(payload["success_count"], 3)

    def test_live_mode_not_available(self):
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as exc:
                main(["--mode", "live"])
        self.assertEqual(exc.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
