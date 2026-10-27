import datetime as dt

import polars as pl
import pytest

from rdv.scv_text import (
    FLAGS, PATTERNS, aggregate, composition, decompose, scv_flags, standardized_share,
)

HEADER = ["#VariationID", "ClinicalSignificance", "DateLastEvaluated", "Description",
          "ExplanationOfInterpretation", "ReviewStatus", "CollectionMethod", "Submitter", "SCV"]


@pytest.mark.parametrize(
    "text,term,expected",
    [
        ("Supporting evidence (PP3) from computational data", "pp3_bp4", True),
        ("bp4 applied", "pp3_bp4", True),
        ("PP30 is not a code", "pp3_bp4", False),
        ("internally defined REVEL score threshold 0.5", "revel", True),
        ("a revelation", "revel", False),
        ("AlphaMissense predicts likely pathogenic", "alphamissense", True),
        ("Alpha-Missense", "alphamissense", True),
        ("SIFT and PolyPhen-2 predict damaging", "legacy_tools", True),
        ("predicted to be deleterious by in silico analysis", "insilico_generic", True),
        ("predicted to result in nonsense-mediated decay", "insilico_generic", False),
        ("highly conserved across vertebrates", "conservation", True),
        ("absent from gnomAD", "population", True),
        ("PM2 and PP3 and BS1", "acmg_code", True),
        ("PVS1_Strong, PP3_Moderate", "pp3_bp4", True),  # glued by an underscore: \\b would have missed it
        ("PS4_VeryStrong; PM2_Supporting", "acmg_code", True),
        ("BP4_Strong", "pp3_bp4", True),
        ("REVEL_score 0.91", "revel", True),
    ],
)
def test_patterns(text, term, expected):
    df = pl.DataFrame({"t": [text]}).select(pl.col("t").str.contains(PATTERNS[term][0]))
    assert df["t"][0] is expected


def test_scv_flags_end_to_end(write_gz):
    crit, noc = "criteria provided, single submitter", "no assertion criteria provided"
    long_text = "x" * 40
    rows = [
        (1, "Pathogenic", "Jun 25, 2024", "REVEL 0.9 supports PP3. " + long_text, "-", crit, "clinical testing", "Lab A", "SCV000001.2"),
        (1, "Benign", "Feb 01, 2020", "-", "predicted benign by in silico analysis " + long_text, noc, "literature only", "Lab B", "SCV000002.1"),
        (2, "Pathogenic", "Jan 01, 2021", "outside universe", "-", crit, "clinical testing", "Lab A", "SCV000003.1"),
    ]
    path = write_gz("submission_summary.txt.gz", HEADER, rows, preamble=["##Overview", "#VariationID: the id"])
    out = scv_flags(path, pl.DataFrame({"VariationID": [1]})).sort("scv")
    assert out["scv"].to_list() == ["SCV000001", "SCV000002"]  # version stripped, variant 2 dropped
    assert out["revel"].to_list() == [True, False]
    assert out["pp3_bp4"].to_list() == [True, False]
    assert out["insilico_generic"].to_list() == [False, True]  # picked up from the explanation column
    assert out["criteria_provided"].to_list() == [True, False]
    assert out["date"][0] == dt.date(2024, 6, 25)
    assert out["label"].to_list() == ["PLP", "BLB"]


def _scv(rows):
    """rows: (submitter, date, revel) -> frame with every flag column present."""
    return pl.DataFrame(
        {
            "submitter": [r[0] for r in rows],
            "date": [r[1] for r in rows],
            "text_len": [100] * len(rows),
            "revel": [r[2] for r in rows],
            **{f: [False] * len(rows) for f in FLAGS if f != "revel"},
        }
    )


def test_aggregate_share_counts_and_wilson():
    d = dt.date(2023, 2, 1)
    scv = _scv([("A", d, True), ("A", d, False), ("B", d, False), ("B", d, False)])
    out = aggregate(scv).filter((pl.col("term") == "revel")).row(0, named=True)
    assert (out["n_text"], out["k"], out["share"]) == (4, 1, 0.25)
    assert out["lo"] < 0.25 < out["hi"] and 0 < out["lo"] and out["hi"] < 1
    assert (out["n_submitters"], out["top_submitter_share"]) == (1, 1.0)
    zero = aggregate(scv).filter(pl.col("term") == "cadd").row(0, named=True)
    assert (zero["k"], zero["share"], zero["n_submitters"]) == (0, 0.0, 0)


def _two_lab_story():
    """Lab A: 10% -> 10%, volume 100 -> 300.  Lab B: 50% -> 90%, volume 100 -> 100."""
    pre, post = dt.date(2021, 6, 1), dt.date(2024, 6, 1)
    rows = []
    for lab, n, rate, d in [("A", 100, 0.1, pre), ("A", 300, 0.1, post), ("B", 100, 0.5, pre), ("B", 100, 0.9, post)]:
        k = round(n * rate)
        rows += [(lab, d, True)] * k + [(lab, d, False)] * (n - k)
    return _scv(rows)


def test_decompose_pooled_change_is_pure_composition_shift():
    r = decompose(_two_lab_story(), "revel", ("2021-01-01", "2022-01-01"), ("2024-01-01", "2025-01-01"))
    assert r["continuing_pre"] == pytest.approx(0.30)
    assert r["continuing_post"] == pytest.approx(0.30)
    assert r["within_lab"] == pytest.approx(0.15)  # B really did start mentioning more
    assert r["between_lab"] == pytest.approx(-0.15)  # but A (low rate) grew
    assert r["within_lab"] + r["between_lab"] == pytest.approx(r["continuing_post"] - r["continuing_pre"])
    assert r["coverage_pre"] == r["coverage_post"] == 1.0
    assert r["pooled_delta"] == pytest.approx(0.0)
    assert r["entry_exit"] == pytest.approx(0.0)


def test_decompose_attributes_a_new_lab_to_entry():
    story = _two_lab_story()
    post = dt.date(2024, 6, 1)
    entrant = _scv([("C", post, True)] * 100)  # lab C only appears after the milestone
    r = decompose(pl.concat([story, entrant]), "revel", ("2021-01-01", "2022-01-01"), ("2024-01-01", "2025-01-01"))
    assert r["pooled_post"] == pytest.approx((30 + 90 + 100) / 500)
    assert r["pooled_delta"] == pytest.approx(0.14)
    assert r["entry_exit"] == pytest.approx(0.14)  # continuing labs still net to zero
    assert r["within_lab"] + r["between_lab"] + r["entry_exit"] == pytest.approx(r["pooled_delta"])
    assert r["coverage_post"] == pytest.approx(400 / 500)


def test_standardized_share_holds_lab_mix_fixed():
    out = standardized_share(
        _two_lab_story(), "revel", ref=("2021-01-01", "2022-01-01"), min_n=10
    ).sort("quarter")
    pre, post = out["share_std"].to_list()
    assert pre == pytest.approx(0.30)  # 50/50 mix of 10% and 50%
    assert post == pytest.approx(0.50)  # same 50/50 mix of 10% and 90%
    assert out["coverage"].to_list() == [1.0, 1.0]


def test_composition_groups_small_labs():
    d = dt.date(2023, 2, 1)
    scv = _scv([("Big Lab, Inc.", d, False)] * 6 + [("Small 1", d, False), ("Small 2", d, False)])
    out = composition(scv, top=1).sort("group")
    assert out["group"].to_list() == ["All other labs", "Big Lab"]
    assert out["share"].to_list() == [0.25, 0.75]
