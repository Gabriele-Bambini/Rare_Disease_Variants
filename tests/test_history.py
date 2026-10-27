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


def test_attach_submitters_takes_first_polarity_matching_scv():
    from rdv.history import attach_submitters

    first = pl.DataFrame({
        "allele_id": [1, 2, 3], "snapshot": ["2022-04"] * 3,
        "label": pl.Series(["PLP", "BLB", "PLP"], dtype=pl.Enum(["VUS", "PLP", "BLB", "CONFLICT", "OTHER"])),
        "route": ["new", "from_vus", "new"],
    })
    universe = pl.DataFrame({"AlleleID": [1, 2, 3], "VariationID": [11, 12, 13], "gene": ["A", "B", "C"]})
    scv = pl.DataFrame({
        "variation_id": [11, 11, 11, 12, 13],
        "scv": ["SCV000500", "SCV000020", "SCV000900", "SCV000030", "SCV000040"],
        "label": ["PLP", "PLP", "BLB", "VUS", "BLB"],  # 11: two PLP SCVs; 12: no BLB SCV; 13: only the opposite polarity
        "submitter": ["Lab late", "Lab first", "Lab other", "Lab vus", "Lab opposite"],
    })
    out = attach_submitters(first, universe, scv).sort("allele_id")
    assert out["submitter"].to_list() == ["Lab first", None, None]
    assert out["scv"].to_list() == ["SCV000020", None, None]
    assert out["gene"].to_list() == ["A", "B", "C"]
