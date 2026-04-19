"""Package: biodiversity-score-engine.

Public API:
    BiodiversityScoreEngine: high-level pipeline orchestrator.
    shannon_index, simpson_index, gini_simpson_index, pielou_evenness,
    species_richness, margalef_richness: pure index functions.
"""
from src.indices import (
    gini_simpson_index,
    margalef_richness,
    pielou_evenness,
    shannon_index,
    simpson_index,
    species_richness,
)
from src.main import BiodiversityScoreEngine

__all__ = [
    "BiodiversityScoreEngine",
    "gini_simpson_index",
    "margalef_richness",
    "pielou_evenness",
    "shannon_index",
    "simpson_index",
    "species_richness",
]
