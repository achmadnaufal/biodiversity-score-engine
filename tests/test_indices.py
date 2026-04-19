"""Unit tests for the pure biodiversity index functions."""
import math

import numpy as np
import pytest

from src.indices import (
    gini_simpson_index,
    margalef_richness,
    pielou_evenness,
    shannon_index,
    simpson_index,
    species_richness,
)


# ---------------------------------------------------------------------------
# Shannon index
# ---------------------------------------------------------------------------

class TestShannonIndex:
    def test_uniform_distribution_equals_ln_n(self):
        # For n equally abundant species, H' = ln(n).
        assert shannon_index([10, 10, 10, 10]) == pytest.approx(math.log(4))
        assert shannon_index([1, 1, 1, 1, 1]) == pytest.approx(math.log(5))

    def test_single_species_is_zero(self):
        assert shannon_index([42]) == pytest.approx(0.0)

    def test_base_two_gives_bits(self):
        # Uniform 4 species: log2(4) = 2.
        assert shannon_index([1, 1, 1, 1], base=2) == pytest.approx(2.0)

    def test_base_ten(self):
        assert shannon_index([1, 1, 1, 1], base=10) == pytest.approx(math.log10(4))

    def test_known_value(self):
        # Worked example: [5, 3, 2] with N=10. Proportions: 0.5, 0.3, 0.2.
        expected = -(0.5 * math.log(0.5) + 0.3 * math.log(0.3) + 0.2 * math.log(0.2))
        assert shannon_index([5, 3, 2]) == pytest.approx(expected)

    def test_zero_counts_are_ignored(self):
        # Zeros should not contribute (0 * log(0) is defined as 0).
        assert shannon_index([10, 10, 10, 10, 0, 0]) == pytest.approx(math.log(4))

    def test_accepts_mapping(self):
        counts = {"A. mangium": 10, "S. robusta": 10}
        assert shannon_index(counts) == pytest.approx(math.log(2))

    def test_accepts_numpy_array(self):
        assert shannon_index(np.array([1.0, 1.0, 1.0, 1.0])) == pytest.approx(math.log(4))

    # -- error paths -------------------------------------------------------

    def test_empty_raises(self):
        with pytest.raises(ValueError, match="empty"):
            shannon_index([])

    def test_all_zero_raises(self):
        with pytest.raises(ValueError, match="positive"):
            shannon_index([0, 0, 0])

    def test_negative_raises(self):
        with pytest.raises(ValueError, match="non-negative"):
            shannon_index([1, -2, 3])

    def test_all_nan_raises(self):
        with pytest.raises(ValueError, match="NaN"):
            shannon_index([float("nan"), float("nan")])

    def test_nan_mixed_with_values(self):
        # NaNs are dropped; remaining values still produce a valid result.
        assert shannon_index([10, 10, float("nan")]) == pytest.approx(math.log(2))

    def test_non_numeric_raises(self):
        with pytest.raises(ValueError, match="numeric"):
            shannon_index(["a", "b", "c"])

    def test_invalid_base(self):
        with pytest.raises(ValueError, match="base"):
            shannon_index([1, 1], base=1)
        with pytest.raises(ValueError, match="base"):
            shannon_index([1, 1], base=0)
        with pytest.raises(ValueError, match="base"):
            shannon_index([1, 1], base=-2)

    def test_none_raises(self):
        with pytest.raises(ValueError):
            shannon_index(None)


# ---------------------------------------------------------------------------
# Simpson / Gini-Simpson
# ---------------------------------------------------------------------------

class TestSimpsonIndex:
    def test_uniform_distribution(self):
        # D = sum(p^2) = n * (1/n)^2 = 1/n.
        assert simpson_index([1, 1, 1, 1]) == pytest.approx(0.25)
        assert simpson_index([5, 5]) == pytest.approx(0.5)

    def test_single_species_is_one(self):
        assert simpson_index([99]) == pytest.approx(1.0)

    def test_known_value(self):
        # [5, 3, 2], N=10 → 0.25 + 0.09 + 0.04 = 0.38
        assert simpson_index([5, 3, 2]) == pytest.approx(0.38)

    def test_gini_simpson_is_complement(self):
        assert gini_simpson_index([1, 1, 1, 1]) == pytest.approx(0.75)
        assert gini_simpson_index([5, 3, 2]) == pytest.approx(0.62)

    def test_empty_raises(self):
        with pytest.raises(ValueError):
            simpson_index([])


# ---------------------------------------------------------------------------
# Pielou evenness (the newly-added feature)
# ---------------------------------------------------------------------------

class TestPielouEvenness:
    def test_uniform_is_one(self):
        assert pielou_evenness([10, 10, 10, 10]) == pytest.approx(1.0)

    def test_single_species_is_one_by_convention(self):
        # S=1 makes ln(S)=0; we return 1.0 as the documented convention.
        assert pielou_evenness([42]) == pytest.approx(1.0)

    def test_dominance_reduces_evenness(self):
        even = pielou_evenness([10, 10, 10, 10])
        uneven = pielou_evenness([97, 1, 1, 1])
        assert uneven < even
        assert 0.0 <= uneven <= 1.0

    def test_known_value(self):
        # Worked example: [97, 1, 1, 1] with ln(4) denominator.
        n = 100
        p = np.array([97, 1, 1, 1]) / n
        h = -np.sum(p * np.log(p))
        expected = h / math.log(4)
        assert pielou_evenness([97, 1, 1, 1]) == pytest.approx(expected)

    def test_bounded_zero_to_one(self):
        rng = np.random.default_rng(0)
        for _ in range(20):
            counts = rng.integers(1, 100, size=rng.integers(2, 10)).tolist()
            j = pielou_evenness(counts)
            assert 0.0 <= j <= 1.0 + 1e-9

    def test_empty_raises(self):
        with pytest.raises(ValueError):
            pielou_evenness([])

    def test_all_zero_raises(self):
        with pytest.raises(ValueError):
            pielou_evenness([0, 0, 0])

    def test_negative_raises(self):
        with pytest.raises(ValueError):
            pielou_evenness([5, -1, 2])

    def test_all_nan_raises(self):
        with pytest.raises(ValueError):
            pielou_evenness([float("nan")])


# ---------------------------------------------------------------------------
# Richness metrics
# ---------------------------------------------------------------------------

class TestRichness:
    def test_species_richness_counts_nonzero(self):
        assert species_richness([5, 0, 3, 2]) == 3
        assert species_richness([1, 1, 1, 1]) == 4
        assert species_richness([9]) == 1

    def test_species_richness_empty_raises(self):
        with pytest.raises(ValueError):
            species_richness([])

    def test_margalef_known_value(self):
        # [10, 10, 10, 10]: S=4, N=40 → (4-1)/ln(40)
        expected = 3 / math.log(40)
        assert margalef_richness([10, 10, 10, 10]) == pytest.approx(expected)

    def test_margalef_single_observation_degenerate(self):
        # N=1 is degenerate; we document it returns 0.0.
        assert margalef_richness([1]) == pytest.approx(0.0)

    def test_margalef_single_species_many_individuals(self):
        # S=1, N>1: (1-1)/ln(N) = 0.
        assert margalef_richness([50]) == pytest.approx(0.0)

    def test_margalef_empty_raises(self):
        with pytest.raises(ValueError):
            margalef_richness([])

    def test_margalef_negative_raises(self):
        with pytest.raises(ValueError):
            margalef_richness([5, -3])
