import numpy as np
import polars as pl
import pytest

from rdv.rq1 import (
    P1, P2, Comparison, add_weights, build_cohort, event_study, half_year_bins, holm, run_did, shift_snapshot,
    shifted, two_sided_p,
)


def test_shift_snapshot_and_shifted_windows():
    assert shift_snapshot("2021-01", -24) == "2019-01"
    assert shift_snapshot("2022-10", -36) == "2019-10"
    assert shift_snapshot("2023-01", -1) == "2022-12"
    s = shifted(P1, -24)
    assert (s.pre, s.post, s.name) == (("2019-01", "2020-10"), ("2021-01", "2022-10"), "P1_pejaver_2022-12_shift-24m")


def test_primary_tests_are_as_frozen():
    assert (P1.treated, P1.control, P1.pre, P1.post) == ("revel", "esm1b", ("2021-01", "2022-10"), ("2023-01", "2024-10"))
    assert (P2.treated, P2.control, P2.pre, P2.post) == ("am", "esm1b", ("2022-01", "2023-07"), ("2023-10", "2025-01"))


def test_half_year_bins():
    b = half_year_bins("2019-01", "2020-10")
    assert b == {"2019H1": ("2019-01", "2019-04"), "2019H2": ("2019-07", "2019-10"),
                 "2020H1": ("2020-01", "2020-04"), "2020H2": ("2020-07", "2020-10")}


def test_build_cohort_orients_scores_and_fills_unknown_submitter():
    first = pl.DataFrame({
        "allele_id": [1, 2, 3], "gene": ["A", "A", "B"], "snapshot": ["2022-04"] * 3,
        "label": pl.Series(["PLP", "BLB", "CONFLICT"], dtype=pl.Enum(["VUS", "PLP", "BLB", "CONFLICT", "OTHER"])),
        "submitter": ["Lab", None, "Lab"], "route": ["new"] * 3, "stars": [1, 1, 1],
    })
    scores = pl.DataFrame({"allele_id": [1, 2, 3], "revel": [0.9, 0.1, 0.5], "am_score": [0.8, 0.2, 0.5],
                           "esm1b_llr": [-9.0, -1.0, -3.0], "grantham": [100.0, 20.0, 50.0], "blosum62": [-3, 2, 0]})
    c = build_cohort(first, scores).sort("allele_id")
    assert c["allele_id"].to_list() == [1, 2]  # the conflict is not a confident entry
    assert c["y"].to_list() == [1, 0] and c["submitter"].to_list() == ["Lab", "unknown"]
    assert c["esm1b"].to_list() == [9.0, 1.0] and c["blosum62"].to_list() == [3, -2]  # larger = more pathogenic
    assert c["revel"].to_list() == [0.9, 0.1] and c["am"].to_list() == [0.8, 0.2]


def test_add_weights_caps_each_submitter_quarter_label_cell():
    rows = [("2022-04", "Big", "BLB")] * 1000 + [("2022-04", "Small", "BLB")] * 10 + [("2022-04", "Big", "PLP")] * 400
    c = pl.DataFrame(rows, schema=["snapshot", "submitter", "label"], orient="row")
    w = add_weights(c, 500)
    by = w.group_by("submitter", "label").agg(pl.col("w").first().alias("w"), pl.len().alias("n"))
    got = {(r["submitter"], r["label"]): r["w"] for r in by.iter_rows(named=True)}
    assert got == {("Big", "BLB"): 0.5, ("Small", "BLB"): 1.0, ("Big", "PLP"): 1.0}
    assert add_weights(c, None)["w"].unique().to_list() == [1.0]  # no cap


def _synthetic_cohort(gain, n_genes=200, per_gene=24, seed=5):
    """T = revel gains separation after 2022; C = esm1b does not. Two big labs and a long tail."""
    rng = np.random.default_rng(seed)
    rows = []
    for g in range(n_genes):
        for snap in ("2021-04", "2021-10", "2022-04", "2023-04", "2023-10", "2024-04"):
            post = snap >= "2023-01"
            for _ in range(per_gene // 6):
                y = int(rng.random() < 0.5)
                lab = "Big" if rng.random() < 0.6 else f"Small{int(rng.integers(0, 30))}"
                rows.append((f"G{g}", snap, y, lab, "new", 1,
                             float(rng.normal() + (1.0 + (gain if post else 0.0)) * y),
                             float(rng.normal() + 1.0 * y),
                             float(rng.normal() + 1.0 * y)))
    return pl.DataFrame(rows, schema=["gene", "snapshot", "y", "submitter", "route", "stars", "revel", "esm1b", "am"], orient="row") \
        .with_columns(pl.col("y").cast(pl.Int8), pl.Series("allele_id", range(len(rows))),
                      pl.when(pl.col("y") == 1).then(pl.lit("PLP")).otherwise(pl.lit("BLB")).alias("label"))


def test_run_did_detects_a_planted_gain_and_reports_the_frozen_fields():
    r = run_did(_synthetic_cohort(0.6), P1, cap=500, replicates=300, seed=1)
    assert r["estimate"] > 0.05 and r["lo"] > 0
    for k in ("se", "mde", "pairs_pre", "pairs_post", "genes_pre", "genes_post", "treated_pre", "control_post"):
        assert r[k] is not None
    null = run_did(_synthetic_cohort(0.0, seed=6), P1, cap=500, replicates=300, seed=1)
    assert null["lo"] < 0 < null["hi"]
    # excluding the big lab changes the sample but not the logic
    r2 = run_did(_synthetic_cohort(0.6), P1, exclude_submitters=("Big",), replicates=100, seed=1)
    assert r2["variants_pre"] < r["variants_pre"]


def test_run_did_uses_common_support_only():
    c = _synthetic_cohort(0.0).with_columns(
        pl.when(pl.col("allele_id") % 2 == 0).then(None).otherwise(pl.col("esm1b")).alias("esm1b"))
    r = run_did(c, P1, replicates=50, seed=1)
    full = run_did(_synthetic_cohort(0.0), P1, replicates=50, seed=1)
    assert r["variants_pre"] < full["variants_pre"]


def test_event_study_returns_auc_and_differences_per_bin():
    c = _synthetic_cohort(0.6)
    es = event_study(c, ["revel", "esm1b"], diffs=(("revel", "esm1b"),), replicates=50, seed=1)
    assert set(es["kind"]) == {"auc", "diff"}
    assert set(es.filter(pl.col("kind") == "diff")["pred"]) == {"revel-esm1b"}
    assert {"2021H1", "2021H2", "2022H1", "2023H1", "2023H2", "2024H1"} <= set(es["bin"])
    d = es.filter((pl.col("kind") == "diff") & (pl.col("bin") == "2023H1")).row(0, named=True)
    assert d["lo"] <= d["value"] <= d["hi"] and d["value"] > 0


def test_two_sided_p_and_holm():
    import numpy as np

    assert two_sided_p(np.linspace(0.1, 1, 1000)) == pytest.approx(1 / 1001)  # all above zero: floor
    assert two_sided_p(np.linspace(-1, 1, 1001)) == pytest.approx(1.0, abs=0.01)  # centred on zero
    assert holm([0.01, 0.04]) == pytest.approx([0.02, 0.04])
    assert holm([0.04, 0.01]) == pytest.approx([0.04, 0.02])  # order preserved
    assert holm([0.3, 0.2]) == pytest.approx([0.4, 0.4])  # monotone after adjustment


def test_run_did_can_restrict_to_one_lab():
    c = _synthetic_cohort(0.0)
    big = run_did(c, P1, include_submitters=("Big",), replicates=50, seed=1)
    allr = run_did(c, P1, replicates=50, seed=1)
    assert 0 < big["variants_pre"] < allr["variants_pre"]


def test_clean_boundary_drops_the_first_post_snapshot():
    c = __import__("rdv.rq1", fromlist=["clean_boundary"]).clean_boundary
    assert c(P1).post == ("2023-04", "2024-10") and c(P1).pre == P1.pre
    assert c(P2).post == ("2024-01", "2025-01") and c(P2).name == "P2_alphamissense_2023-09_clean"


def test_group_contrast_sees_a_planted_group_difference_and_ignores_shared_case_mix():
    from rdv.rq1 import group_contrast

    rng = np.random.default_rng(11)
    rows = []
    for g in range(250):
        for k in range(24):
            cites = k % 2 == 0
            y = int(rng.random() < 0.5)
            shared = 0.3 if cites else 0.0  # easier case mix for the citing group: every predictor gains
            rows.append((f"G{g}", "2022-04", y, "Lab", "new", 1, 80, cites,
                         float(rng.normal() + (1.0 + shared + (0.6 if cites else 0.0)) * y),  # revel: extra gain when citing
                         float(rng.normal() + (1.0 + shared) * y)))                          # esm1b: only the shared gain
    df = pl.DataFrame(rows, schema=["gene", "snapshot", "y", "submitter", "route", "stars", "text_len", "cites", "revel", "esm1b"],
                      orient="row").with_columns(pl.col("y").cast(pl.Int8), pl.lit("PLP").alias("label"),
                                                 pl.Series("allele_id", range(250 * 24)))
    out = group_contrast(df, "cites", replicates=300, seed=1)
    assert out["gap_difference"] > 0.05 and out["lo"] > 0  # planted extra separation of 0.6 for revel
    assert out["auc_A_esm1b"] > out["auc_B_esm1b"]  # shared case-mix shift is visible in the control too
    assert out["variants_A"] == out["variants_B"] == 3000
