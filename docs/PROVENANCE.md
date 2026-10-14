# Provenance and scientific boundaries

This repository contains original offline teaching examples from the revised Day 01–30 course mapping. It is not the hosted React/Node Research Agent or an official Google DeepMind distribution.

- Original 2026-10-13 candidate SHA-256: `fbc0ec7e724b1a1e11df4dbe7915d8d65df29c6b55cbf29a4f5180f997d670c8`
- 2026-10-14 merged candidate SHA-256: `858c8afc4617a750fea1cc8405cadde1e5dc6b461d8c018e2b9c203da11b608d`
- Original Day11–12 historical AVI context files were removed; Day11–12 code and tests were updated to accept only synthetic comparison data. Day17–18 source and tests were also replaced with synthetic-only fixture readers. All other existing teaching source and fixture files were retained unchanged relative to the input Day11–12 synthetic baseline (per byte-level diff audit).
- **Required archived example_output fixtures retained: 61 files.** Day 28 tests read the Day 27 example run directory. Day 30 also includes saved reference outputs. These are part of the reproducible teaching example and may be required for offline tests; their existence does not mean real API requests were made.
- Internal release files excluded: 2.
- No live AlphaGenome, Gemini, PubMed or clinical API inference is performed by this build.
- DEMO / SYNTHETIC / FIXTURE data must not be presented as live model outputs, validated ACMG classifications or clinical diagnoses.
- Days 17–18 were converted to **hand-authored synthetic-only** teaching readers. Their fixtures use the invented `chrToy` genome, toy AVI raw/quantile/PHRED-like values, artificial feature contributions, and illustrative RNA/splicing REF/ALT tracks. `examples/historical_readers/SYNTHETIC_FIXTURE_MANIFEST.json` pins SHA-256 of the two synthetic inputs; hashes establish internal byte integrity only, not model authenticity or external provenance. No real historical numerical extracts, coordinates, ontologies, or source-file fingerprints remain in these two reader scripts or their bundled fixtures.
- Day11–12 synthetic-only comparison inputs remove that specific historical AVI redistribution concern. The Day17/18 refactor removes analogous source-code extracts. **A public release decision requires maintainer authorization and separate checks of software QA, licenses, scientific claims and private information.** Offline fixture tests do not independently establish legal or clinical clearance. Consult the current GitHub repository visibility and release record to determine publication status.
- MIT applies only to the author's original code and rights-cleared materials. Third-party services, brand marks and hosted model outputs are not relicensed. See NOTICE.md.

Local tests / SHA-256 checks are not scientific, privacy, licensing or clinical certification.
