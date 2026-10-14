"""Day 19: deterministic OFFLINE-only asyncio retry/limit laboratory.

Research education only. No SDK, HTTP client, credentials, or live API access.
HTTP-like 429/503/400 statuses and TIMEOUT are entirely simulated.
"""
from __future__ import annotations

import argparse
import asyncio
import json
from pathlib import Path

# These variants and numbers are invented. They are NOT GRCh38 or API outputs.
SCENARIOS = {
    "chrToy:101:A>G": [
        {"code": 429, "retry_after_s": 0.03},
        {"code": 200, "synthetic_score": 0.12},
    ],
    "chrToy:102:C>T": [
        {"code": 503}, {"code": 503},
        {"code": 200, "synthetic_score": -0.21},
    ],
    "chrToy:103:G>A": [{"code": 400}],
    "chrToy:104:T>C": [{"code": 503}, {"code": 503}, {"code": 503}],
    "chrToy:105:A>C": [
        {"code": "TIMEOUT"}, {"code": 200, "synthetic_score": 0.07},
    ],
}
RETRYABLE = {429, 503, "TIMEOUT"}


class FakeAPI:
    """One preprogrammed local response per attempt; never makes network calls."""

    def __init__(self, scenarios=None):
        self.scenarios = SCENARIOS if scenarios is None else scenarios
        self.in_flight = 0
        self.max_in_flight = 0

    async def query(self, variant, attempt):
        self.in_flight += 1
        self.max_in_flight = max(self.max_in_flight, self.in_flight)
        try:
            # A tiny local delay permits asyncio tasks to overlap for teaching.
            await asyncio.sleep(0.005)
            outcomes = self.scenarios[variant]
            if attempt > len(outcomes):
                raise ValueError("The teaching scenario has no such attempt")
            return outcomes[attempt - 1]
        finally:
            self.in_flight -= 1


def delay_for(code, attempt, retry_after_s=None):
    """Teaching backoff: 0.01, 0.02,... seconds, max 0.08 seconds."""
    delay = min(0.08, 0.01 * (2 ** (attempt - 1)))
    if code == 429 and isinstance(retry_after_s, (int, float)) and retry_after_s >= 0:
        delay = max(delay, retry_after_s)
    return round(delay, 3)


async def query_with_retry(variant, fake_api, semaphore, max_attempts=3, sleeper=asyncio.sleep):
    """Retry only retryable mock errors; release semaphore while waiting."""
    history = []
    waits = []
    for attempt in range(1, max_attempts + 1):
        # This semaphore limits simultaneous requests, NOT requests per second.
        async with semaphore:
            response = await fake_api.query(variant, attempt)
        code = response["code"]
        history.append(code)

        if code == 200:
            return {"variant": variant, "status": "OK", "attempts": attempt,
                    "history": history, "retry_wait_s": waits,
                    "synthetic_score": response["synthetic_score"]}
        if code not in RETRYABLE:
            return {"variant": variant, "status": "NOT_RETRYABLE", "attempts": attempt,
                    "history": history, "retry_wait_s": waits,
                    "synthetic_score": None}
        if attempt == max_attempts:
            return {"variant": variant, "status": "RETRY_EXHAUSTED", "attempts": attempt,
                    "history": history, "retry_wait_s": waits,
                    "synthetic_score": None}

        wait = delay_for(code, attempt, response.get("retry_after_s"))
        waits.append(wait)
        await sleeper(wait)

    raise AssertionError("Unreachable")


async def run_batch(max_concurrency=2, max_attempts=3, scenarios=None, sleeper=asyncio.sleep):
    """Simulate a batch without keys, SDK, internet, or patient records."""
    if not 1 <= max_concurrency <= 8:
        raise ValueError("max_concurrency must be 1..8")
    if not 1 <= max_attempts <= 5:
        raise ValueError("max_attempts must be 1..5")
    api = FakeAPI(scenarios)
    semaphore = asyncio.Semaphore(max_concurrency)
    # gather() returns results in input order even if tasks complete differently.
    rows = await asyncio.gather(*(
        query_with_retry(v, api, semaphore, max_attempts, sleeper)
        for v in api.scenarios
    ))
    return {
        "evidence_class": "SYNTHETIC_FIXTURE",
        "mode": "OFFLINE_ONLY",
        "sdk_used": False,
        "real_api_calls_this_run": 0,
        "clinical_interpretation": False,
        "scenario_count": len(rows),
        "success_count": sum(r["status"] == "OK" for r in rows),
        "failure_count": sum(r["status"] != "OK" for r in rows),
        "total_attempts": sum(r["attempts"] for r in rows),
        "max_concurrency_limit": max_concurrency,
        "observed_max_mock_in_flight": api.max_in_flight,
        "sum_scheduled_retry_wait_s": round(
            sum(sum(r["retry_wait_s"]) for r in rows), 3),
        "results": rows,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description="Day19 OFFLINE asyncio mock retry lab")
    parser.add_argument("--mode", choices=["offline"], default="offline")
    parser.add_argument("--max-concurrency", type=int, default=2)
    parser.add_argument("--max-attempts", type=int, default=3)
    parser.add_argument("--output-dir", default="day19_output")
    args = parser.parse_args(argv)
    try:
        payload = asyncio.run(run_batch(args.max_concurrency, args.max_attempts))
    except ValueError as exc:
        parser.error(str(exc))
    folder = Path(args.output_dir)
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "day19_mock_batch.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    summary = {"mode": payload["mode"], "scenario_count": payload["scenario_count"],
               "success_count": payload["success_count"],
               "failure_count": payload["failure_count"],
               "total_attempts": payload["total_attempts"],
               "real_api_calls_this_run": payload["real_api_calls_this_run"],
               "output_file": str(folder / "day19_mock_batch.json")}
    print(json.dumps(summary, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
