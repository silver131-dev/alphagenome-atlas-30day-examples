# AlphaGenome Atlas | Day 17–18 — SYNTHETIC ONLY (offline teaching example)

> **NOT REAL ALPHAGENOME OUTPUT.** All scores, positions, tracks and REF/ALT
> differences are fictional teaching values, including `chrToy:110:A>G`.
> No AlphaGenome API calls or patient data. Not for clinical interpretation.

These readers are integrated into a **separately reviewed 30-day archive candidate**.
They replace the old historical-source readers without changing the Day01–30 article index.
This conversion is not authorization to publish and must be retested in any other checkout.

## Learning outcomes

* Day 17: toy AVI raw score, quantile, percentile, PHRED-like transformation,
  feature contributions, and the distinction between ranking and clinical
  classification. The toy score is **not** a pathogenicity probability.
* Day 18: fictive splice donor/acceptor scores, strand/1-based positions,
  RNA-like arbitrary-unit signals, and paired REF/ALT `delta = ALT - REF`.
  A negative toy donor delta does **not** prove an actual splicing defect.
* Both: bundled `SYNTHETIC_FIXTURE_MANIFEST.json` with SHA-256 hashes, verified
  against bundled fixture bytes. Integrity is **not** evidence that Google
  or DeepMind generated a score. The origin is `hand_authored`.

The two entrypoint names retain `historical_` only to support the staged public
repository's previous file layout; scripts and reports explicitly reject any
historical-output claim. Adapt all existing callers and tests in the original
repository before publishing.

## Run (Python 3.12+ standard library)

From anywhere:

```bash
python3.12 examples/historical_readers/day17_historical_avi_reader.py --verify-only
python3.12 examples/historical_readers/day18_historical_splicing_reader.py --verify-only
python3.12 examples/historical_readers/day17_historical_avi_reader.py --output-dir /tmp/day17-demo
python3.12 examples/historical_readers/day18_historical_splicing_reader.py --output-dir /tmp/day18-demo
python3.12 -m unittest discover -s examples/historical_readers -p 'test_*.py' -v
```

The programs do not accept API keys, issue network calls, import AlphaGenome SDK,
or require any credentials. Failed SHA-256 or source-label checks exit nonzero.

## Provenance limitations

- The manifest demonstrates that *these exact files* have not changed since its
  checksums were authored. It does not validate any historical API response.
- All values are hand-authored fiction, not observations, not model estimates,
  and not medically meaningful. No real source values are inherited.
- Authors must review Day 17–18 articles, README links, `docs/PROVENANCE.md`,
  source lists and CI references for references to obsolete historical results.
- Earlier standalone QA recorded `PUBLIC_RELEASE=HOLD` before full repository tests.
  That historical state is not a current GitHub visibility assertion. Full
  repository QA, rights review and maintainer-controlled publication are distinct.
