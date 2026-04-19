# Biodiversity Score Engine

A lightweight Python toolkit for computing standard biodiversity indices from
species-count survey data produced by reforestation and nature-based-solution
(NbS) monitoring programmes.

## What it does

Given a tidy table of observations (one row per species per plot), the engine
aggregates records and reports:

- **Shannon diversity (H')** — information-theoretic diversity, natural log
- **Simpson's index (D)** — probability two individuals share a species
- **Gini–Simpson (1 - D)** — complementary diversity form
- **Pielou's evenness (J')** — how equally individuals are distributed
- **Species richness (S)** — count of distinct species observed
- **Margalef richness** — richness normalised by sample size

Metrics are produced for the pooled sample and, when a `plot_id` column is
present, broken out per plot.

## Installation

```bash
pip install -r requirements.txt
```

Python 3.9+ is supported.

## Quick start

```python
from src.main import BiodiversityScoreEngine

engine = BiodiversityScoreEngine()
df = engine.load_data("demo/sample_data.csv")

# One-shot pipeline
result = engine.run("demo/sample_data.csv")
print(result["indices"])
# {'shannon': 2.14..., 'simpson': 0.13..., 'gini_simpson': 0.86...,
#  'pielou_evenness': 0.89..., 'species_richness': 11,
#  'margalef_richness': 2.69..., 'total_individuals': 201.0}

# Per-plot breakdown
per_plot = engine.compute_per_plot(df)
print(per_plot)
```

The pure index functions are also importable directly:

```python
from src.indices import shannon_index, pielou_evenness

shannon_index([10, 10, 10, 10])   # → ln(4) ≈ 1.3863
pielou_evenness([97, 1, 1, 1])    # → ~0.129 (heavy dominance)
```

## Data format

Expected CSV/Excel columns:

| column          | required | notes                                   |
|-----------------|----------|-----------------------------------------|
| `species`       | yes      | Scientific or common name               |
| `count`         | yes      | Non-negative individual count           |
| `plot_id`       | optional | Enables per-plot index breakdown        |
| `date_surveyed` | optional | ISO date string                         |
| `site_name`     | optional | Free text                               |
| `observer`      | optional | Free text                               |

Column names are case- and whitespace-insensitive (`"Plot ID"` → `plot_id`).

See [`demo/sample_data.csv`](demo/sample_data.csv) for a small working
example drawn from tropical and temperate forest species.

## Features

- Data ingestion from CSV and Excel (`.xlsx` / `.xls`)
- Input validation with explicit, actionable `ValueError` messages
- Pooled and per-plot biodiversity index suite
- Pure, dependency-light index functions (`src/indices.py`)
- Immutable preprocessing (input DataFrames are never mutated)
- Edge-case handling: empty input, all-zero counts, single species, NaN values
- Sample-data generator for synthetic testing (`src/data_generator.py`)
- Comprehensive pytest suite covering known-value and edge cases

## Running the tests

```bash
pip install pytest
pytest tests/
```

## Project structure

```
biodiversity-score-engine/
├── src/
│   ├── __init__.py         # Public API re-exports
│   ├── main.py             # BiodiversityScoreEngine pipeline
│   ├── indices.py          # Pure index functions
│   └── data_generator.py   # Synthetic data generator
├── tests/
│   ├── test_indices.py     # Index function unit tests
│   └── test_engine.py      # Engine integration tests
├── demo/
│   └── sample_data.csv     # Small worked example
├── data/                   # User data (gitignored)
├── examples/
│   └── basic_usage.py
├── requirements.txt
├── CHANGELOG.md
└── README.md
```

## Domain reference

The metrics implemented here are textbook ecological indices:

- Shannon, C. E. (1948). *A Mathematical Theory of Communication.*
- Simpson, E. H. (1949). Measurement of Diversity. *Nature*, 163, 688.
- Pielou, E. C. (1966). The measurement of diversity in different types of
  biological collections. *J. Theor. Biol.*, 13, 131–144.
- Margalef, R. (1958). Information theory in ecology. *General Systems*, 3.
- Magurran, A. E. (2004). *Measuring Biological Diversity.* Blackwell.

This engine makes no novel claims — it is a clean, tested implementation of
standard biodiversity metrics suitable for reforestation / NbS monitoring.

## License

MIT License — free to use, modify, and distribute.
