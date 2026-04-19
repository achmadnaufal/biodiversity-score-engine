"""Integration tests for :class:`BiodiversityScoreEngine`."""
import math
from pathlib import Path

import pandas as pd
import pytest

from src.main import BiodiversityScoreEngine


DEMO_CSV = Path(__file__).resolve().parent.parent / "demo" / "sample_data.csv"


@pytest.fixture
def tiny_df():
    """A minimal two-plot survey with known expected indices."""
    return pd.DataFrame(
        {
            "plot_id": ["P1", "P1", "P1", "P1", "P2", "P2"],
            "species": ["A", "B", "C", "D", "A", "B"],
            "count": [10, 10, 10, 10, 50, 50],
        }
    )


class TestValidation:
    def test_empty_dataframe_raises(self):
        engine = BiodiversityScoreEngine()
        with pytest.raises(ValueError, match="empty"):
            engine.validate(pd.DataFrame())

    def test_missing_species_column(self):
        engine = BiodiversityScoreEngine()
        df = pd.DataFrame({"count": [1, 2, 3]})
        with pytest.raises(ValueError, match="species"):
            engine.validate(df)

    def test_missing_count_column(self):
        engine = BiodiversityScoreEngine()
        df = pd.DataFrame({"species": ["A", "B"]})
        with pytest.raises(ValueError, match="count"):
            engine.validate(df)

    def test_accepts_uppercase_columns(self):
        engine = BiodiversityScoreEngine()
        df = pd.DataFrame({"Species": ["A", "B"], "Count": [1, 2]})
        # preprocess lowers the column names; validate should then pass.
        assert engine.validate(engine.preprocess(df)) is True


class TestPreprocess:
    def test_does_not_mutate_input(self):
        engine = BiodiversityScoreEngine()
        df = pd.DataFrame({"Species ": ["A"], " Count": [1]})
        original = df.copy(deep=True)
        _ = engine.preprocess(df)
        pd.testing.assert_frame_equal(df, original)

    def test_normalises_column_names(self):
        engine = BiodiversityScoreEngine()
        df = pd.DataFrame({"Plot ID": ["P1"], "Species": ["A"], "Count": [1]})
        out = engine.preprocess(df)
        assert list(out.columns) == ["plot_id", "species", "count"]


class TestComputeIndices:
    def test_uniform_four_species(self, tiny_df):
        engine = BiodiversityScoreEngine()
        p1 = tiny_df[tiny_df["plot_id"] == "P1"]
        result = engine.compute_indices(p1)
        assert result["shannon"] == pytest.approx(math.log(4))
        assert result["simpson"] == pytest.approx(0.25)
        assert result["gini_simpson"] == pytest.approx(0.75)
        assert result["pielou_evenness"] == pytest.approx(1.0)
        assert result["species_richness"] == 4
        assert result["total_individuals"] == 40.0

    def test_aggregates_duplicate_species(self):
        # Two rows for the same species should sum.
        engine = BiodiversityScoreEngine()
        df = pd.DataFrame(
            {
                "species": ["A", "A", "B", "B"],
                "count": [3, 7, 4, 6],
            }
        )
        result = engine.compute_indices(df)
        # Aggregated: A=10, B=10 → uniform 2 species → H = ln(2).
        assert result["shannon"] == pytest.approx(math.log(2))
        assert result["species_richness"] == 2
        assert result["total_individuals"] == 20.0

    def test_all_zero_counts_raises(self):
        engine = BiodiversityScoreEngine()
        df = pd.DataFrame({"species": ["A", "B"], "count": [0, 0]})
        with pytest.raises(ValueError):
            engine.compute_indices(df)

    def test_single_species(self):
        engine = BiodiversityScoreEngine()
        df = pd.DataFrame({"species": ["Shorea robusta"], "count": [25]})
        result = engine.compute_indices(df)
        assert result["shannon"] == pytest.approx(0.0)
        assert result["simpson"] == pytest.approx(1.0)
        assert result["species_richness"] == 1


class TestPerPlot:
    def test_per_plot_returns_row_per_plot(self, tiny_df):
        engine = BiodiversityScoreEngine()
        per_plot = engine.compute_per_plot(tiny_df)
        assert set(per_plot.index) == {"P1", "P2"}
        assert per_plot.loc["P1", "shannon"] == pytest.approx(math.log(4))
        assert per_plot.loc["P2", "shannon"] == pytest.approx(math.log(2))

    def test_per_plot_requires_plot_id(self):
        engine = BiodiversityScoreEngine()
        df = pd.DataFrame({"species": ["A"], "count": [1]})
        with pytest.raises(ValueError, match="plot_id"):
            engine.compute_per_plot(df)


class TestAnalyze:
    def test_analyze_includes_indices(self, tiny_df):
        engine = BiodiversityScoreEngine()
        result = engine.analyze(tiny_df)
        assert "indices" in result
        assert "per_plot" in result
        assert result["indices"]["species_richness"] == 4  # A, B, C, D pooled

    def test_analyze_empty_raises(self):
        engine = BiodiversityScoreEngine()
        with pytest.raises(ValueError):
            engine.analyze(pd.DataFrame())

    def test_to_dataframe_flattens(self, tiny_df):
        engine = BiodiversityScoreEngine()
        result = engine.analyze(tiny_df)
        flat = engine.to_dataframe(result)
        assert set(flat.columns) == {"metric", "value"}
        assert any(m.startswith("indices.shannon") for m in flat["metric"])


class TestDemoData:
    def test_demo_csv_exists(self):
        assert DEMO_CSV.exists(), f"expected demo data at {DEMO_CSV}"

    def test_engine_runs_on_demo_csv(self):
        engine = BiodiversityScoreEngine()
        result = engine.run(str(DEMO_CSV))
        assert result["indices"]["species_richness"] >= 5
        assert result["indices"]["shannon"] > 0
        assert 0 < result["indices"]["simpson"] <= 1
        assert "per_plot" in result and len(result["per_plot"]) >= 2

    def test_load_data_missing_file_raises(self, tmp_path):
        engine = BiodiversityScoreEngine()
        with pytest.raises(FileNotFoundError):
            engine.load_data(str(tmp_path / "does_not_exist.csv"))


class TestNaNHandling:
    def test_nan_counts_are_ignored(self):
        engine = BiodiversityScoreEngine()
        df = pd.DataFrame(
            {
                "species": ["A", "B", "C"],
                "count": [10.0, 10.0, float("nan")],
            }
        )
        # NaN is dropped by groupby(min_count=1); A and B remain uniform.
        result = engine.compute_indices(df)
        assert result["species_richness"] == 2
        assert result["shannon"] == pytest.approx(math.log(2))
