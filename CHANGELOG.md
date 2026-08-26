# Changelog

All notable changes to this project are documented in this file.

## [1.0.0] - 2026-01-01

Initial public release.

- Engine Block 1: classic IRFinder-principle ratio (`IRratio_classic`), corrected
  (`best_perf`) and traditional (`classic`) mask policies, recount3 adapter.
- Engine Block 2: coverage-shape layer for `long`-tier introns (>=1000bp) —
  uniformity, cliff detection, GC/mappability bias flags, IPA candidate ranking.
- Engine Block 2b/2c: IPA-driven down-weight removed from `IRratio_corrected`
  (no validated basis); `retention_confidence` decoupled from `shape_bias_flag`.
- Engine Block 3: two-point donor/acceptor coverage-shape layer for
  `awkward_mid`-tier introns (100-999bp).
- `irshape` CLI (`build-reference`, `run`), external reference-bundle resolution
  (`--ref-dir` / `$IRSHAPE_REF`), Docker and Apptainer packaging.
- `bam` input adapter interface (stubbed; raises `NotImplementedError`).
