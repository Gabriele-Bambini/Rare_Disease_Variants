import numpy as np
import polars as pl
import pytest

from rdv import promotion as pr


def _write_history(tmp_path, states):
    """states: {snapshot: {allele_id: label}} -> history directory of snapshot parquets."""
    d = tmp_path / "hist"
    d.mkdir()
    for snap, st in states.items():
        pl.DataFrame({"allele_id": list(st), "label": list(st.values())}, schema={"allele_id": pl.Int64, "label": pl.Utf8}) \
            .write_parquet(d / f"snapshot={snap}.parquet")
    return d


def test_interval_events_classifies_every_rule(tmp_path):
    states = {
        "2021-01": {1: "VUS", 2: "VUS", 3: "VUS", 4: "VUS", 5: "VUS", 6: "VUS", 7: "VUS", 8: "VUS"},
        "2021-04": {1: "PLP", 2: "BLB", 3: "VUS", 4: "CONFLICT", 5: "VUS", 6: "VUS", 8: "PLP"},  # 7 withdrawn
    }
    hist = _write_history(tmp_path, states)
    first = pl.DataFrame({"allele_id": [1, 2, 5, 6, 8], "first_snapshot": ["2021-04", "2021-04", "2021-10", "2018-01", "2019-04"]})
    out = pr.interval_events(hist, "2021-01", "2021-04", first).sort("allele_id")
    got = dict(zip(out["allele_id"], out["outcome"]))
    assert got == {
        1: "up",        # VUS -> P/LP, first confident label
        2: "down",      # VUS -> B/LB, first confident label
        3: "control",   # stays VUS, never confident
        4: "other",     # moved to conflicting: neither case nor control
        5: "control",   # stays VUS, will be promoted later (first confident after the interval)
        6: "other",     # was confident in 2018, downgraded to VUS: not a clean control
        8: "other",     # P/LP now but first confident in 2019: a re-promotion, not a first event
    }                    # 7 was withdrawn and does not appear


def test_case_weights_cap_and_waves():
    ev = pl.DataFrame({"s": ["2022-04"] * 6 + ["2022-07"], "submitter": ["A"] * 4 + ["B"] * 2 + ["A"], "outcome": ["up"] * 7,
                       "allele_id": range(7), "gene": ["G"] * 7})
    w = pr.case_weights(ev, 2).sort("allele_id")["w"].to_list()
    assert w == [0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0]  # a cell of 4 is capped at 2
    assert pr.case_weights(ev, None)["w"].unique().to_list() == [1.0]
    assert pr.wave_intervals(ev, share=0.5) == ["2022-04"]


def _world(tmp_path, slope_pre, slope_post, seed=0, n_genes=150, per_gene=40):
    """Variants with a covariate; each interval a VUS is promoted with odds ~ exp(slope * revel)."""
    rng = np.random.default_rng(seed)
    n = n_genes * per_gene
    ids = np.arange(1, n + 1)
    gene = np.array([f"G{i // per_gene}" for i in range(n)])
    X = rng.normal(size=(n, 4))
    snaps = ["2021-01", "2021-04", "2021-07", "2021-10", "2022-01", "2022-04", "2022-07", "2022-10", "2023-01"]
    label = np.array(["VUS"] * n, dtype=object)
    first = {}
    states = {snaps[0]: dict(zip(ids, label))}
    for t, s in enumerate(snaps[1:], start=1):
        slope = slope_pre if s <= "2021-10" else slope_post
        p = 0.06 * np.exp(slope * X[:, 0]) / np.exp(slope**2 / 2)
        up = (label == "VUS") & (rng.random(n) < np.clip(p, 0, 0.9))
        label[up] = "PLP"
        for i in ids[up]:
            first[int(i)] = s
        states[s] = dict(zip(ids, label.copy()))
    hist = _write_history(tmp_path, states)
    scores = pl.DataFrame({"allele_id": ids, "VariationID": ids, "revel": X[:, 0], "am_score": X[:, 1], "esm1b_llr": -X[:, 2],
                           "grantham": X[:, 3]})
    af = pl.DataFrame({"VariationID": ids, "af_proxy": [None] * n}, schema={"VariationID": pl.Int64, "af_proxy": pl.Float64})
    universe = pl.DataFrame({"allele_id": ids, "gene": gene})
    first_snap = pl.DataFrame({"allele_id": list(first), "first_snapshot": list(first.values())},
                              schema={"allele_id": pl.Int64, "first_snapshot": pl.Utf8})
    return hist, first_snap, universe, pr.covariate_table(scores, af, universe)


def test_covariate_table_orients_and_standardises(tmp_path):
    _, _, _, cov = _world(tmp_path, 0, 0, n_genes=20)
    assert cov["revel"].mean() == pytest.approx(0, abs=1e-9) and cov["revel"].std() == pytest.approx(1, abs=1e-9)
    assert cov["esm1b"].mean() == pytest.approx(0, abs=1e-9)  # esm1b_llr was negated into a pathogenicity score
    assert cov["log_af"].to_list()[:3] == [0.0, 0.0, 0.0] and cov["af_seen"].sum() == 0


def test_end_to_end_recovers_a_planted_change_in_the_covariate_effect(tmp_path):
    hist, first, universe, cov = _world(tmp_path, slope_pre=0.2, slope_post=0.9, seed=1)
    ev = pr.all_events(hist, first, universe, pl.DataFrame({"allele_id": [], "submitter": []}, schema={"allele_id": pl.Int64, "submitter": pl.Utf8}))
    st = pr.build_strata(hist, first, ev, universe, cov, "up", m=10, seed=3)
    assert st.X.shape[1] == 11 and st.valid[:, 0].all()
    win = pr.Window("synthetic", "revel", ("2021-04", "2021-10"), ("2022-04", "2022-10"))
    r = pr.fit_window(st, win, ev, cap=None)
    assert r["converged"] and r["n_pre"] > 100 and r["n_post"] > 100
    assert abs(r["pre_revel"] - 0.2) < 3 * r["se_revel"] + 0.05
    assert r["delta_revel"] == pytest.approx(0.7, abs=0.25) and r["lo_revel"] > 0  # planted change of 0.7
    assert abs(r["delta_esm1b"]) < 3.5 * r["se_esm1b"]  # a covariate with no planted change
    mde = pr.pool_information(hist, first, ev, universe, cov, "up", win, cap=None)
    assert 0 < mde["mde_revel"] < 0.8


def test_benign_direction_flips_the_score_orientation(tmp_path):
    hist, first, universe, cov = _world(tmp_path, 0.5, 0.5, seed=2, n_genes=40)
    ev = pr.all_events(hist, first, universe, pl.DataFrame({"allele_id": [], "submitter": []}, schema={"allele_id": pl.Int64, "submitter": pl.Utf8}))
    # relabel the synthetic promotions as benign events: the first four covariates must come out negated
    ev_down = ev.with_columns(pl.lit("down").alias("outcome"))
    up = pr.build_strata(hist, first, ev, universe, cov, "up", m=5, seed=3)
    down = pr.build_strata(hist, first, ev_down, universe, cov, "down", m=5, seed=3)
    assert up.X[0, 0, :4] == pytest.approx(-down.X[0, 0, :4])
    assert up.X[0, 0, 4:] == pytest.approx(down.X[0, 0, 4:])  # frequency covariates are not flipped
