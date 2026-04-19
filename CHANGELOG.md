# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `src/indices.py` module with pure biodiversity index functions:
  Shannon, Simpson, Gini–Simpson, Pielou evenness, species richness,
  and Margalef richness — all with Google-style docstrings and
  explicit edge-case handling.
- `BiodiversityScoreEngine.compute_indices` and `compute_per_plot`
  methods for pooled and per-plot diversity metrics.
- Pielou's evenness index as the featured new metric, with the
  documented `J' = 1.0` convention for the degenerate single-species
  case.
- `tests/` suite (56 tests) covering known-value checks for every
  index, engine integration, and edge cases (empty input, all-zero
  counts, single species, NaN values, missing columns).
- `demo/sample_data.csv` with 20 realistic rows across five plots
  and three sites, using tropical and temperate forest species
  (e.g. *Shorea robusta*, *Quercus robur*, *Dipterocarpus turbinatus*).
- Public-API re-exports in `src/__init__.py`.
- Domain reference section and quick-start example in `README.md`.

### Changed
- `BiodiversityScoreEngine.analyze` now returns biodiversity indices
  (previously only descriptive statistics).
- `validate` enforces the presence of required `species` and `count`
  columns with an actionable error message.
- `preprocess` no longer mutates its input DataFrame.
- `load_data` raises `FileNotFoundError` eagerly for missing files.

## [0.1.0] - Initial release

### Added
- `BiodiversityScoreEngine` skeleton with CSV/Excel ingestion and
  descriptive-statistic analysis.
- `src/data_generator.py` for synthetic sample data.
- `examples/basic_usage.py` demo script.
