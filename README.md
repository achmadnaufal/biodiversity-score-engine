# Biodiversity Score Engine

Biodiversity index calculator for reforestation and nature-based solution projects

## Features
- Data ingestion from CSV/Excel input files
- Automated analysis and KPI calculation
- Summary statistics and trend reporting
- Sample data generator for testing and development

## Installation

```bash
pip install -r requirements.txt
```

## Quick Start

```python
from src.main import BiodiversityScoreEngine

analyzer = BiodiversityScoreEngine()
df = analyzer.load_data("data/sample.csv")
result = analyzer.analyze(df)
print(result)
```

## Data Format

Expected CSV columns: `plot_id, species, count, period, survey_date`

## Project Structure

```
biodiversity-score-engine/
├── src/
│   ├── main.py          # Core analysis logic
│   └── data_generator.py # Sample data generator
├── data/                # Data directory (gitignored for real data)
├── examples/            # Usage examples
├── requirements.txt
└── README.md
```

## License

MIT License — free to use, modify, and distribute.
