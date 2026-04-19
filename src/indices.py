"""
Biodiversity indices used by the Biodiversity Score Engine.

All functions in this module are pure: they accept a sequence of species
counts (or a mapping from species name to count) and return a scalar
diversity metric. They raise ``ValueError`` on invalid input (empty,
non-numeric, negative values, or all-zero counts) so callers can fail
fast at system boundaries.

References:
    - Shannon, C. E. (1948). A Mathematical Theory of Communication.
    - Simpson, E. H. (1949). Measurement of Diversity. Nature.
    - Pielou, E. C. (1966). The measurement of diversity in different
      types of biological collections. J. Theor. Biol.
    - Margalef, R. (1958). Information theory in ecology.
"""
from __future__ import annotations

import math
from typing import Iterable, Mapping, Sequence, Union

import numpy as np


CountsLike = Union[Sequence[float], Iterable[float], Mapping[str, float], np.ndarray]


def _to_counts_array(counts: CountsLike) -> np.ndarray:
    """Coerce user input into a validated 1-D numpy array of counts.

    Args:
        counts: A sequence, iterable, mapping, or numpy array of species
            counts. Mappings use their values only.

    Returns:
        A 1-D float numpy array with NaNs removed, containing only the
        non-zero species counts.

    Raises:
        ValueError: If ``counts`` is empty, contains negative values or
            NaNs that cannot be dropped to leave any counts, or sums to
            zero after filtering.
        TypeError: If ``counts`` is not an iterable/mapping/array.
    """
    if counts is None:
        raise ValueError("counts must not be None")

    if isinstance(counts, Mapping):
        values = list(counts.values())
    elif isinstance(counts, np.ndarray):
        values = counts.tolist()
    elif hasattr(counts, "__iter__"):
        values = list(counts)
    else:
        raise TypeError(
            f"counts must be a sequence, mapping, or ndarray, got {type(counts).__name__}"
        )

    if len(values) == 0:
        raise ValueError("counts is empty; cannot compute a biodiversity index")

    try:
        arr = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"counts must be numeric; got: {exc}") from exc

    # Drop NaNs defensively — these are common in ecological field data.
    arr = arr[~np.isnan(arr)]
    if arr.size == 0:
        raise ValueError("counts contains only NaN values")

    if np.any(arr < 0):
        raise ValueError("counts must be non-negative")

    if arr.sum() <= 0:
        raise ValueError("counts must contain at least one positive observation")

    # Zero counts represent species not observed in this sample and
    # should be excluded from entropy/richness calculations.
    return arr[arr > 0]


def shannon_index(counts: CountsLike, base: float = math.e) -> float:
    """Compute the Shannon diversity index (H').

    ``H' = -sum(p_i * log(p_i))`` where ``p_i`` is the proportional
    abundance of species ``i``. Higher values indicate greater diversity.

    Args:
        counts: Species counts (see :func:`_to_counts_array`).
        base: Logarithm base. Defaults to natural log (``e``). Use
            ``2`` for bits, ``10`` for decimal digits.

    Returns:
        The Shannon index as a non-negative float. Returns ``0.0`` for a
        single-species community.

    Raises:
        ValueError: If ``counts`` is empty, all-zero, contains negatives
            or only NaNs, or if ``base`` is not positive / equals 1.

    Example:
        >>> shannon_index([10, 10, 10, 10])  # uniform → ln(4)
        1.3862943611198906
    """
    if base <= 0 or base == 1:
        raise ValueError("base must be positive and not equal to 1")

    arr = _to_counts_array(counts)
    proportions = arr / arr.sum()
    return float(-np.sum(proportions * np.log(proportions) / math.log(base)))


def simpson_index(counts: CountsLike) -> float:
    """Compute Simpson's diversity index ``D = sum(p_i**2)``.

    This is the probability that two individuals drawn at random belong
    to the same species. Lower values indicate greater diversity.

    Args:
        counts: Species counts.

    Returns:
        Simpson's ``D`` in the range ``(0, 1]``.

    Raises:
        ValueError: For empty, negative, all-NaN or all-zero input.

    Example:
        >>> round(simpson_index([10, 10, 10, 10]), 4)
        0.25
    """
    arr = _to_counts_array(counts)
    proportions = arr / arr.sum()
    return float(np.sum(proportions ** 2))


def gini_simpson_index(counts: CountsLike) -> float:
    """Compute the Gini-Simpson index ``1 - D``.

    Args:
        counts: Species counts.

    Returns:
        A value in ``[0, 1)`` where higher means more diverse.

    Raises:
        ValueError: For empty, negative, all-NaN or all-zero input.

    Example:
        >>> round(gini_simpson_index([10, 10, 10, 10]), 4)
        0.75
    """
    return 1.0 - simpson_index(counts)


def pielou_evenness(counts: CountsLike) -> float:
    """Compute Pielou's evenness index ``J' = H' / ln(S)``.

    Evenness measures how equally individuals are distributed among the
    species present. ``J' = 1`` indicates perfectly even abundance;
    lower values indicate dominance by a few species.

    For a single-species community evenness is mathematically undefined
    (``ln(1) = 0``). By ecological convention this function returns
    ``1.0`` in that degenerate case, since a single species is trivially
    "evenly" distributed with itself.

    Args:
        counts: Species counts.

    Returns:
        Pielou's ``J'`` in the range ``[0, 1]``.

    Raises:
        ValueError: For empty, negative, all-NaN or all-zero input.

    Example:
        >>> round(pielou_evenness([10, 10, 10, 10]), 4)
        1.0
        >>> round(pielou_evenness([97, 1, 1, 1]), 3)
        0.129
    """
    arr = _to_counts_array(counts)
    species_richness = arr.size
    if species_richness == 1:
        return 1.0
    h = shannon_index(arr)
    return float(h / math.log(species_richness))


def species_richness(counts: CountsLike) -> int:
    """Return the number of distinct species with at least one observation.

    Args:
        counts: Species counts.

    Returns:
        The count of species with non-zero abundance.

    Raises:
        ValueError: For empty, negative, all-NaN or all-zero input.

    Example:
        >>> species_richness([5, 0, 3, 2])
        3
    """
    return int(_to_counts_array(counts).size)


def margalef_richness(counts: CountsLike) -> float:
    """Compute Margalef's richness index ``D_Mg = (S - 1) / ln(N)``.

    Normalises species richness (``S``) by sample size (``N``) using a
    logarithmic correction. Undefined for ``N = 1``; this function
    returns ``0.0`` in that degenerate case.

    Args:
        counts: Species counts.

    Returns:
        Margalef's richness as a non-negative float.

    Raises:
        ValueError: For empty, negative, all-NaN or all-zero input.

    Example:
        >>> round(margalef_richness([10, 10, 10, 10]), 4)
        0.8123
    """
    arr = _to_counts_array(counts)
    n_total = float(arr.sum())
    s = arr.size
    if n_total <= 1:
        return 0.0
    return float((s - 1) / math.log(n_total))
