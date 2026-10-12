import polars as pl

from rdv.history import derive, entry_route, extract_snapshot, quarter_of, snapshot_path

HEADER = ["#AlleleID", "ClinicalSignificance", "ReviewStatus", "LastEvaluated"]
CRIT = "criteria provided, single submitter"

# allele: label at 2017-01, 2017-04, 2017-07, 2017-10
STORY = {
    1: ["Pathogenic"] * 4,  # confident from the start: left-censored
    2: ["Uncertain significance", "Uncertain significance", "Likely pathogenic", "Likely pathogenic"],
    3: [None, None, "Benign", "Benign"],  # enters ClinVar already benign
    4: ["Pathogenic", "Conflicting interpretations of pathogenicity", "Benign", "Benign"],
    5: ["Uncertain significance"] * 4,  # never confident: excluded
    6: ["Conflicting interpretations of pathogenicity", "Conflicting interpretations of pathogenicity",
        "Likely benign", "Uncertain significance"],  # resolved, then downgraded
    7: [None, "Pathogenic", "Pathogenic", None],  # later withdrawn from ClinVar
}
MONTHS = ["2017-01", "2017-04", "2017-07", "2017-10"]


def _build(write_gz, tmp_path):
    ids = pl.DataFrame({"allele_id": list(STORY) + [999]})  # 999 never appears
    out = tmp_path / "hist"
    out.mkdir()
    for i, month in enumerate(MONTHS):
        rows = [
            (a, labels[i], CRIT, "Jan 05, 2016")
            for a, labels in STORY.items()
            if labels[i] is not None
        ]
        rows.append((12345, "Pathogenic", CRIT, "-"))  # outside the universe: dropped
        rows.append((1, "Pathogenic", CRIT, "Jan 05, 2016"))  # GRCh37 duplicate row
        src = write_gz(f"variant_summary_{month}.txt.gz", HEADER, rows)
        extract_snapshot(src, ids).write_parquet(snapshot_path(out, month))
    return out


def test_extract_snapshot_filters_dedupes_and_types(write_gz, tmp_path):
    out = _build(write_gz, tmp_path)
    s = pl.read_parquet(snapshot_path(out, "2017-01"))
    assert s["allele_id"].to_list() == [1, 2, 4, 5, 6]
    assert s["snapshot"].unique().to_list() == ["2017-01"]
    assert s["label"].to_list() == ["PLP", "VUS", "PLP", "VUS", "CONFLICT"]
    assert s["stars"].unique().to_list() == [1]
    assert str(s["last_evaluated"][0]) == "2016-01-05"


def test_derive_first_confident(write_gz, tmp_path):
    first, trans = derive(_build(write_gz, tmp_path))
    f = entry_route(first).sort("allele_id")
    assert f["allele_id"].to_list() == [1, 2, 3, 4, 6, 7]  # 5 never confident
    got = {r["allele_id"]: r for r in f.iter_rows(named=True)}
    assert (got[1]["snapshot"], got[1]["route"], got[1]["left_censored"]) == ("2017-01", "left_censored", True)
    assert (got[2]["snapshot"], got[2]["label"], got[2]["route"]) == ("2017-07", "PLP", "from_vus")
    assert (got[3]["snapshot"], got[3]["label"], got[3]["route"]) == ("2017-07", "BLB", "new")
    assert got[3]["first_seen"] == "2017-07"
    # 4 was PLP first (left-censored), later flipped to BLB: first confident stays the first
    assert (got[4]["snapshot"], got[4]["label"]) == ("2017-01", "PLP")
    assert (got[6]["snapshot"], got[6]["label"], got[6]["route"]) == ("2017-07", "BLB", "from_conflict")
    assert got[6]["current_label"] == "VUS"
    assert got[7]["current_label"] is None  # withdrawn before the latest snapshot


def test_derive_transitions(write_gz, tmp_path):
    _, trans = derive(_build(write_gz, tmp_path))
    pairs = {(r["allele_id"], r["snapshot"], r["prior_label"], r["label"]) for r in trans.iter_rows(named=True)}
    assert pairs == {
        (2, "2017-07", "VUS", "PLP"),
        (4, "2017-04", "PLP", "CONFLICT"),
        (4, "2017-07", "CONFLICT", "BLB"),
        (6, "2017-07", "CONFLICT", "BLB"),
        (6, "2017-10", "BLB", "VUS"),
    }


def test_quarter_of():
    assert [quarter_of(m) for m in ["2017-01", "2017-04", "2017-07", "2017-10"]] == [
        "2017Q1", "2017Q2", "2017Q3", "2017Q4"]
