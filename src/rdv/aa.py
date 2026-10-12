"""Amino-acid helpers: code conversion and two substitution scores that need no data download.

Grantham (1974) and BLOSUM62 are the chemistry-only baselines of the plan (feature group C and
the control signals of RQ1): neither has been calibrated for ACMG PP3/BP4, and neither saw ClinVar.
"""

from __future__ import annotations

import math

import polars as pl

ONE = "ARNDCQEGHILKMFPSTWYV"  # BLOSUM row/column order
THREE = ["Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
         "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"]
THREE_TO_ONE = dict(zip(THREE, ONE))

# Grantham composition c, polarity p, molecular volume v.
_CPV = {
    "A": (0.00, 8.1, 31.0), "R": (0.65, 10.5, 124.0), "N": (1.33, 11.6, 56.0),
    "D": (1.38, 13.0, 54.0), "C": (2.75, 5.5, 55.0), "Q": (0.89, 10.5, 85.0),
    "E": (0.92, 12.3, 83.0), "G": (0.74, 9.0, 3.0), "H": (0.58, 10.4, 96.0),
    "I": (0.00, 5.2, 111.0), "L": (0.00, 4.9, 111.0), "K": (0.33, 11.3, 119.0),
    "M": (0.00, 5.7, 105.0), "F": (0.00, 5.2, 132.0), "P": (0.39, 8.0, 32.5),
    "S": (1.42, 9.2, 32.0), "T": (0.71, 8.6, 61.0), "W": (0.13, 5.4, 170.0),
    "Y": (0.20, 6.2, 136.0), "V": (0.00, 5.9, 84.0),
}
_ALPHA, _BETA, _GAMMA, _SCALE = 1.833, 0.1018, 0.000399, 50.723


def grantham(a: str, b: str) -> float:
    """Grantham distance between two one-letter amino acids (0 = identical, max about 215)."""
    (c1, p1, v1), (c2, p2, v2) = _CPV[a], _CPV[b]
    return _SCALE * math.sqrt(_ALPHA * (c1 - c2) ** 2 + _BETA * (p1 - p2) ** 2 + _GAMMA * (v1 - v2) ** 2)


# Standard BLOSUM62 in ONE order (identical to Biopython's substitution_matrices.load("BLOSUM62")).
_BLOSUM62 = (
    (4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0),
    (-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3),
    (-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3),
    (-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3),
    (0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1),
    (-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2),
    (-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2),
    (0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3),
    (-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3),
    (-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3),
    (-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1),
    (-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2),
    (-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1),
    (-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1),
    (-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2),
    (1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2),
    (0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0),
    (-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3),
    (-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1),
    (0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4),
)


def blosum62(a: str, b: str) -> int:
    return _BLOSUM62[ONE.index(a)][ONE.index(b)]


def substitution_scores(ref: pl.Series, alt: pl.Series) -> pl.DataFrame:
    """Grantham and BLOSUM62 for paired one-letter columns, via a 400-entry lookup."""
    table = pl.DataFrame(
        [(a, b, grantham(a, b), blosum62(a, b)) for a in ONE for b in ONE],
        schema={"ref": pl.Utf8, "alt": pl.Utf8, "grantham": pl.Float64, "blosum62": pl.Int8},
        orient="row",
    )
    return (
        pl.DataFrame({"ref": ref, "alt": alt})
        .join(table, on=["ref", "alt"], how="left")
        .select("grantham", "blosum62")
    )
