import polars as pl
import pytest

from rdv.clinvar import (
    BLB, CONFLICT, OTHER, PLP, VUS,
    UNIVERSE_COLUMNS, classify_significance, iter_tsv_chunks, missense_universe,
    parse_date, read_tsv, review_stars, snapshot_month,
)


@pytest.mark.parametrize(
    "text,expected",
    [
        ("Pathogenic", PLP),
        ("Likely pathogenic", PLP),
        ("Pathogenic/Likely pathogenic", PLP),
        ("Pathogenic, low penetrance", PLP),
        ("Pathogenic; risk factor", PLP),
        ("Pathogenic/Likely pathogenic/Pathogenic, low penetrance", PLP),
        ("Benign", BLB),
        ("Likely benign", BLB),
        ("Benign/Likely benign", BLB),
        ("Uncertain significance", VUS),
        ("Conflicting classifications of pathogenicity", CONFLICT),
        ("Conflicting interpretations of pathogenicity", CONFLICT),
        ("Pathogenic/Benign", CONFLICT),
        ("not provided", OTHER),
        ("risk factor", OTHER),
        ("Likely risk allele", OTHER),
        ("drug response", OTHER),
        ("", OTHER),
        (None, OTHER),
    ],
)
def test_classify_significance(text, expected):
    assert classify_significance(text) == expected


@pytest.mark.parametrize(
    "status,stars",
    [
        ("practice guideline", 4),
        ("reviewed by expert panel", 3),
        ("criteria provided, multiple submitters, no conflicts", 2),
        ("criteria provided, single submitter", 1),
        ("criteria provided, conflicting classifications", 1),
        ("no assertion criteria provided", 0),
        ("no classification provided", 0),
        (None, 0),
    ],
)
def test_review_stars(status, stars):
    assert review_stars(status) == stars


def test_snapshot_month():
    assert snapshot_month("/x/variant_summary_2023-04.txt.gz") == "2023-04"
    with pytest.raises(ValueError):
        snapshot_month("variant_summary.txt.gz")


def test_parse_date():
    df = pl.DataFrame({"d": ["Jun 29, 2010", "-", None, "garbage"]})
    out = df.select(parse_date("d"))["d"].to_list()
    assert str(out[0]) == "2010-06-29"
    assert out[1:] == [None, None, None]


def test_reader_skips_prose_selects_columns_and_fills_missing(write_gz):
    path = write_gz(
        "submission_summary_2019-07.txt.gz",
        ["#VariationID", "ClinicalSignificance", "Description"],
        [(1, "Pathogenic", "a"), (2, "Benign", "b"), (3, "Uncertain significance", "c")],
        preamble=["##Overview of interpretation", "#VariationID: the identifier assigned by ClinVar"],
    )
    chunks = list(iter_tsv_chunks(path, ["VariationID", "Description", "NotInThisRelease"], chunk_lines=2))
    assert [c.height for c in chunks] == [2, 1]
    df = pl.concat(chunks)
    assert df.columns == ["VariationID", "Description", "NotInThisRelease"]
    assert df["Description"].to_list() == ["a", "b", "c"]
    assert df["NotInThisRelease"].null_count() == 3


def test_read_tsv_applies_chunk_filter(write_gz):
    path = write_gz("t.txt.gz", ["#A", "B"], [(i, "x") for i in range(10)])
    out = read_tsv(path, ["A"], lambda d: d.filter(pl.col("A").cast(pl.Int64) % 2 == 0), chunk_lines=3)
    assert out["A"].to_list() == ["0", "2", "4", "6", "8"]


def _vs_row(**kw):
    base = {c: "-" for c in UNIVERSE_COLUMNS}
    base.update(Type="single nucleotide variant", Assembly="GRCh38", ReviewStatus="criteria provided, single submitter",
              ClinicalSignificance="Uncertain significance", Chromosome="1", PositionVCF="100",
              ReferenceAlleleVCF="C", AlternateAlleleVCF="T")
    base.update(kw)
    return [base[c] for c in UNIVERSE_COLUMNS]


def test_missense_universe_keeps_only_grch38_single_aa_substitutions(write_gz):
    rows = [
        _vs_row(AlleleID=1, VariationID=11, GeneSymbol="TP53", Name="NM_000546.6(TP53):c.817C>T (p.Arg273Cys)",
                ClinicalSignificance="Pathogenic"),
        _vs_row(AlleleID=1, VariationID=11, GeneSymbol="TP53", Assembly="GRCh37",
                Name="NM_000546.5(TP53):c.817C>T (p.Arg273Cys)"),  # other assembly: dropped
        _vs_row(AlleleID=2, VariationID=12, Name="NM_1.1(G2):c.1A>G (p.Met1Val)"),  # start-loss
        _vs_row(AlleleID=3, VariationID=13, Name="NM_1.1(G3):c.5C>A (p.Ser2Ter)"),  # nonsense
        _vs_row(AlleleID=4, VariationID=14, Name="NM_1.1(G4):c.6C>T (p.Ser2=)"),  # synonymous
        _vs_row(AlleleID=5, VariationID=15, Name="NM_1.1(G5):c.8_9del (p.Ser3fs)", Type="Deletion"),
        _vs_row(AlleleID=6, VariationID=16, Name="NM_1.1(G6):c.9A>G (p.Lys3Arg)",
                ClinicalSignificance="Benign", ReviewStatus="reviewed by expert panel"),
    ]
    path = write_gz("variant_summary.txt.gz", ["#" + UNIVERSE_COLUMNS[0], *UNIVERSE_COLUMNS[1:]], rows)
    u = missense_universe(path)
    assert u["VariationID"].to_list() == [11, 16]
    r = u.row(0, named=True)
    assert (r["gene"], r["ref_aa"], r["aa_pos"], r["alt_aa"], r["label"]) == ("TP53", "Arg", 273, "Cys", PLP)
    assert u.filter(pl.col("VariationID") == 16).row(0, named=True)["stars"] == 3
