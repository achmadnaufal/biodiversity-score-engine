"""
Biodiversity index calculator for reforestation and nature-based solution projects.

This module exposes :class:`BiodiversityScoreEngine`, a high-level
orchestrator that loads tabular survey data, validates it, and computes
a suite of standard biodiversity metrics (Shannon, Simpson, Gini-Simpson,
Pielou evenness, species richness, Margalef richness) per plot and for
the pooled sample.

Author: github.com/achmadnaufal
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd

from src.indices import (
    gini_simpson_index,
    margalef_richness,
    pielou_evenness,
    shannon_index,
    simpson_index,
    species_richness,
)


REQUIRED_COLUMNS = ("species", "count")


class BiodiversityScoreEngine:
    """Multi-index biodiversity calculator.

    Loads species-count survey data and computes standard ecological
    diversity metrics. Supports per-plot aggregation when a ``plot_id``
    column is present.

    Attributes:
        config: Optional configuration dictionary. Currently unused;
            reserved for future extension (custom log base, rarefaction
            depth, etc.).

    Example:
        >>> engine = BiodiversityScoreEngine()
        >>> df = engine.load_data("demo/sample_data.csv")
        >>> indices = engine.compute_indices(df)
        >>> indices["shannon"]  # doctest: +SKIP
        2.14
    """

    def __init__(self, config: Optional[Dict[str, Any]] = None) -> None:
        self.config = dict(config) if config else {}

    def load_data(self, filepath: str) -> pd.DataFrame:
        """Load survey data from a CSV or Excel file.

        Args:
            filepath: Path to the input file. ``.xlsx`` / ``.xls`` are
                parsed with :func:`pandas.read_excel`; all other suffixes
                fall through to :func:`pandas.read_csv`.

        Returns:
            The loaded DataFrame.

        Raises:
            FileNotFoundError: If ``filepath`` does not exist.
            ValueError: If the file cannot be parsed.
        """
        p = Path(filepath)
        if not p.exists():
            raise FileNotFoundError(f"Input file not found: {filepath}")
        if p.suffix.lower() in (".xlsx", ".xls"):
            return pd.read_excel(filepath)
        return pd.read_csv(filepath)

    def validate(self, df: pd.DataFrame) -> bool:
        """Validate that the DataFrame is fit for index calculation.

        Args:
            df: Input DataFrame.

        Returns:
            ``True`` if validation passes.

        Raises:
            ValueError: If ``df`` is empty or is missing required
                columns (``species``, ``count``).
        """
        if df is None or df.empty:
            raise ValueError("Input DataFrame is empty")
        cols = {c.lower().strip() for c in df.columns}
        missing = [c for c in REQUIRED_COLUMNS if c not in cols]
        if missing:
            raise ValueError(
                f"Input is missing required column(s): {missing}. "
                f"Required columns are: {list(REQUIRED_COLUMNS)}"
            )
        return True

    def preprocess(self, df: pd.DataFrame) -> pd.DataFrame:
        """Normalise column names and drop fully-empty rows.

        Args:
            df: Raw input DataFrame.

        Returns:
            A new DataFrame (the input is not mutated) with lower-cased,
            whitespace-normalised column names and empty rows dropped.
        """
        if df is None:
            raise ValueError("df must not be None")
        out = df.copy()
        out.dropna(how="all", inplace=True)
        out.columns = [str(c).lower().strip().replace(" ", "_") for c in out.columns]
        return out

    def species_abundance(self, df: pd.DataFrame) -> pd.Series:
        """Aggregate raw records into a species -> total count series.

        Args:
            df: DataFrame containing ``species`` and ``count`` columns.

        Returns:
            A pandas Series indexed by species name with integer/float
            total counts, sorted in descending abundance order.

        Raises:
            ValueError: If required columns are missing or the
                resulting abundance vector is empty.
        """
        df = self.preprocess(df)
        self.validate(df)
        grouped = (
            df.groupby("species", dropna=True)["count"]
            .sum(min_count=1)
            .dropna()
            .sort_values(ascending=False)
        )
        if grouped.empty:
            raise ValueError("No species counts available after aggregation")
        return grouped

    def compute_indices(self, df: pd.DataFrame) -> Dict[str, float]:
        """Compute the full suite of diversity indices for a dataset.

        Args:
            df: DataFrame containing ``species`` and ``count`` columns.

        Returns:
            A dictionary with keys:
                - ``shannon``: Shannon diversity (natural log).
                - ``simpson``: Simpson's D.
                - ``gini_simpson``: 1 - Simpson's D.
                - ``pielou_evenness``: Pielou's J'.
                - ``species_richness``: Number of distinct species.
                - ``margalef_richness``: Margalef's richness index.
                - ``total_individuals``: Sum of observed counts.

        Raises:
            ValueError: If the dataset is empty or all counts are zero.
        """
        abundance = self.species_abundance(df)
        counts = abundance.to_numpy(dtype=float)
        return {
            "shannon": shannon_index(counts),
            "simpson": simpson_index(counts),
            "gini_simpson": gini_simpson_index(counts),
            "pielou_evenness": pielou_evenness(counts),
            "species_richness": species_richness(counts),
            "margalef_richness": margalef_richness(counts),
            "total_individuals": float(counts.sum()),
        }

    def compute_per_plot(self, df: pd.DataFrame) -> pd.DataFrame:
        """Compute diversity indices per ``plot_id`` group.

        Args:
            df: DataFrame with at minimum ``plot_id``, ``species`` and
                ``count`` columns.

        Returns:
            A DataFrame indexed by ``plot_id`` with one column per
            index. Plots with zero total abundance are skipped.

        Raises:
            ValueError: If required columns (including ``plot_id``) are
                missing.
        """
        df = self.preprocess(df)
        if "plot_id" not in df.columns:
            raise ValueError("compute_per_plot requires a 'plot_id' column")
        self.validate(df)

        rows = []
        for plot_id, sub in df.groupby("plot_id", dropna=True):
            try:
                idx = self.compute_indices(sub)
            except ValueError:
                # Skip plots that yield no valid counts.
                continue
            idx["plot_id"] = plot_id
            rows.append(idx)

        if not rows:
            raise ValueError("No plots produced valid diversity indices")

        return pd.DataFrame(rows).set_index("plot_id").sort_index()

    def analyze(self, df: pd.DataFrame) -> Dict[str, Any]:
        """Run the full analysis pipeline on a DataFrame.

        Produces both descriptive table stats and the biodiversity index
        suite in a single result dictionary.

        Args:
            df: Input DataFrame.

        Returns:
            A dictionary with descriptive metadata (``total_records``,
            ``columns``, ``missing_pct``) and a nested ``indices`` key
            holding the diversity-metric results.

        Raises:
            ValueError: If the DataFrame cannot be validated or indices
                cannot be computed.
        """
        df = self.preprocess(df)
        self.validate(df)

        result: Dict[str, Any] = {
            "total_records": int(len(df)),
            "columns": list(df.columns),
            "missing_pct": (df.isnull().sum() / len(df) * 100).round(1).to_dict(),
            "indices": self.compute_indices(df),
        }

        numeric_df = df.select_dtypes(include="number")
        if not numeric_df.empty:
            result["summary_stats"] = numeric_df.describe().round(3).to_dict()
            result["totals"] = numeric_df.sum().round(2).to_dict()
            result["means"] = numeric_df.mean().round(3).to_dict()

        if "plot_id" in df.columns:
            try:
                result["per_plot"] = self.compute_per_plot(df).round(4).to_dict(orient="index")
            except ValueError:
                result["per_plot"] = {}

        return result

    def run(self, filepath: str) -> Dict[str, Any]:
        """Full pipeline: ``load -> validate -> analyze``.

        Args:
            filepath: Path to the input survey file.

        Returns:
            The analysis result dictionary (see :meth:`analyze`).

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If the file cannot be parsed or validated.
        """
        df = self.load_data(filepath)
        self.validate(df)
        return self.analyze(df)

    def to_dataframe(self, result: Dict[str, Any]) -> pd.DataFrame:
        """Flatten a nested analysis result into a tidy DataFrame.

        Args:
            result: The dictionary returned by :meth:`analyze`.

        Returns:
            A long-format DataFrame with ``metric`` and ``value``
            columns. Nested dictionaries are flattened with
            dot-separated keys.
        """
        rows = []
        for k, v in result.items():
            if isinstance(v, dict):
                for kk, vv in v.items():
                    if isinstance(vv, dict):
                        for kkk, vvv in vv.items():
                            rows.append({"metric": f"{k}.{kk}.{kkk}", "value": vvv})
                    else:
                        rows.append({"metric": f"{k}.{kk}", "value": vv})
            else:
                rows.append({"metric": k, "value": v})
        return pd.DataFrame(rows)
