# AlphaGenome Atlas｜Day 17–18 Synthetic-only 整合驗證報告

**Date:** 2026-10-14  ·  **Author:** ARiN  ·  **PUBLIC_RELEASE=HOLD (historical QA snapshot before publishing)**

> This is a prepublication evidence record. Its HOLD status reflects the time of the original audit, not an assertion about the repository's present GitHub visibility. See current GitHub repository status and rights notice.

## Input baseline
- User supplied: `AlphaGenome_Atlas_30Day_SYNTHETIC_DAY11_12_20261014(1).zip`
- Input SHA-256: `f9015b9f517c5dfd9fc78f792f6eced6eb668e044b5b5152c240ce99f31b5f60`
- Input content: 165 files, preconverted Day11/12 synthetic-only, validated by ZIP CRC.

## Day 17–18 changes
- Replace hard-coded historical AVI scores/source fingerprints and historical splicing values/ontology/genomic positions with **hand-authored fictional chrToy** fixtures.
- Day17: raw toy AVI 0.85; quantile 0.99; PHRED-like 20.0 computed by `-10 log10(1-q)`; toy feature-contribution sum 0.75 and residual 0.10 recorded as `REVIEW` (illustrative only).
- Day18: fictional donor REF 0.8 / ALT 0.1 (Δ -0.7), acceptor REF 0.25 / ALT 0.5 (Δ +0.25), RNA-like arbitrary-unit REF/ALT pairs; direction does not prove biological effects.
- SHA-256 manifest verifies **bundled synthetic input bytes**, not authentic Google DeepMind model provenance. No historical API response or patient data was used.
- The old `day17_historical_avi_reader.py` and `day18_historical_splicing_reader.py` filenames remain for compatibility; their behavior is completely synthetic-only.
- Preserve Day11/12 source bytes and all 30 article URLs; update only reader files, tests, documentation, and new synthetic fixtures/outputs.

## Observed tests
| Verification | Observed result |
|---|---|
| Focused historical_readers | **14/14 PASS** |
| Full offline `python scripts/run_tests.py` | **21 directories, 1842 tests PASS, 0 failures** |
| `scripts/check_repo_safety.py` | PASS — 0 common credential-regex hits; limited scan only |
| Day17 / Day18 `--verify-only` | PASS / PASS, `API_CALLS=0` |
| Fixture SHA-256 match | 2/2 PASS |
| Tampering, source labels, invalid score, REF/ALT identity | Negative cases PASS |
| Original Day11/12 files | Byte-identical to baseline |
| README 30 iThome links | Byte-identical `article_links.csv`, 30/30 links in README |
| Real model API calls during tests | 0 |
| Legacy Day17/18 historical SHA + quantitative extracts in new reader tree | Not detected in static scan |

## File boundaries
- Unchanged original files: **159 / 165**.
- Modified files: **6**.
- Added synthetic-only files: **9** (before this report).
- Deleted original baseline files: 0.
- No original AlphaGenome native API artifacts are bundled.

### Modified paths
- `NOTICE.md`
- `README.md`
- `docs/PROVENANCE.md`
- `examples/historical_readers/day17_historical_avi_reader.py`
- `examples/historical_readers/day18_historical_splicing_reader.py`
- `examples/historical_readers/test_historical_readers.py`

### Added paths
- `examples/historical_readers/README.md`
- `examples/historical_readers/SYNTHETIC_FIXTURE_MANIFEST.json`
- `examples/historical_readers/day17_synthetic_avi_fixture.json`
- `examples/historical_readers/day18_synthetic_splicing_fixture.json`
- `examples/historical_readers/example_output/day17_synthetic_avi_report.json`
- `examples/historical_readers/example_output/day17_synthetic_avi_report.md`
- `examples/historical_readers/example_output/day18_synthetic_splicing_report.json`
- `examples/historical_readers/example_output/day18_synthetic_splicing_report.md`
- `examples/historical_readers/fixture_helpers.py`

## Publication gate
- `DAY17_SYNTHETIC_ONLY=PASS` (this packaged candidate)
- `DAY18_SYNTHETIC_ONLY=PASS` (this packaged candidate)
- `FIXTURE_SOURCE_QA=PASS` (integrity + explicit synthetic-label verification, **not actual model authenticity**)
- `NO_LIVE_API=PASS` (offline scripts and tests)
- `REPOSITORY_TESTS=PASS` (this environment, Python 3.13; WSL 3.12/CI should rerun)
- `PUBLIC_RELEASE=HOLD`: no GitHub push/visibility change, author must approve publishing; perform independent source/rights review.
- WSL original staging checkout is **not automatically modified** by this downloaded ZIP.
- Keep historical AlphaGenome reference artifacts in private, separate research storage only.
