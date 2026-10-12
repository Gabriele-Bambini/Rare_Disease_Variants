import polars as pl
import pytest

from rdv.aa import ONE, THREE, THREE_TO_ONE, blosum62, grantham, substitution_scores


@pytest.mark.parametrize(
    "a,b,published",
    [  # Grantham (1974) table, rounded
        ("I", "L", 5), ("R", "K", 26), ("D", "E", 45), ("F", "Y", 22), ("L", "M", 15),
        ("A", "G", 60), ("S", "R", 110), ("C", "S", 112), ("C", "R", 180), ("G", "W", 184), ("C", "W", 215),
    ],
)
def test_grantham_matches_published_table(a, b, published):
    assert round(grantham(a, b)) == pytest.approx(published, abs=1)


def test_grantham_is_a_distance():
    assert all(grantham(a, a) == 0 for a in ONE)
    assert all(grantham(a, b) == pytest.approx(grantham(b, a)) for a in ONE for b in ONE)
    assert max(grantham(a, b) for a in ONE for b in ONE) == pytest.approx(215, abs=1)


def test_blosum62_known_entries_and_symmetry():
    assert (blosum62("W", "W"), blosum62("C", "C"), blosum62("A", "R"), blosum62("I", "V")) == (11, 9, -1, 3)
    assert all(blosum62(a, b) == blosum62(b, a) for a in ONE for b in ONE)


def test_blosum62_matches_biopython():
    sm = pytest.importorskip("Bio.Align.substitution_matrices")
    m = sm.load("BLOSUM62")
    assert all(blosum62(a, b) == m[a][b] for a in ONE for b in ONE)


def test_code_tables_are_consistent():
    assert len(ONE) == len(THREE) == len(set(THREE)) == 20
    assert THREE_TO_ONE["Arg"] == "R" and THREE_TO_ONE["Trp"] == "W"


def test_substitution_scores_vectorised():
    out = substitution_scores(pl.Series(["I", "C"]), pl.Series(["L", "W"]))
    assert out["grantham"].to_list() == pytest.approx([5, 215], abs=1)
    assert out["blosum62"].to_list() == [2, -2]
