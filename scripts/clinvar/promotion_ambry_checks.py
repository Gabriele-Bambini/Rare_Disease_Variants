"""Two post-hoc checks on the Ambry benign-direction lead, on data already seen (see the Deviations log).

    python scripts/clinvar/promotion_ambry_checks.py --data-dir $DATA_DIR [--results results]

1. Specificity: do AlphaMissense or ESM1b move with REVEL inside Ambry?
2. Portfolio: with controls restricted to variants Ambry itself has a submission on, does the effect remain?
"""

import argparse
from pathlib import Path

import numpy as np
import polars as pl

from rdv import promotion as pr
from promotion_prefit import load

AMB = "Ambry Genetics"


def show(label, r):
    print(f"{label}  strata {r['n_pre']}/{r['n_post']}")
    for name in ("revel", "am", "esm1b", "grantham", "log_af", "af_seen"):
        print(f"    {name:<9} {r['pre_' + name]:+.2f} -> {r['post_' + name]:+.2f}   change {r['delta_' + name]:+.2f} (SE {r['se_' + name]:.2f})")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", required=True, type=Path)
    ap.add_argument("--results", type=Path, default=Path("results"))
    a = ap.parse_args()
    hist, fc, gene_of, sub, cov = load(a.data_dir)
    ev = pr.all_events(hist, fc, gene_of, sub)
    per = ev.group_by("s").len()
    waves = tuple(per.filter(pl.col("len") > 5 * per["len"].median())["s"].to_list())

    st = pr.build_strata(hist, fc, ev, gene_of, cov, "down")
    others = tuple(x for x in st.meta["submitter"].unique().to_list() if x != AMB)
    r1 = pr.fit_window(st, pr.Q1, ev, exclude_s=waves, exclude_submitters=others)
    print("1. Ambry only, benign direction, Q1 windows, controls drawn from all VUS in the gene")
    show("   ", r1)

    scv = pl.read_parquet(a.data_dir / "derived/scv_flags.parquet", columns=["variation_id", "submitter"])
    uni = pl.read_parquet(a.data_dir / "derived/missense_universe.parquet", columns=["AlleleID", "VariationID"])
    ids = uni.join(scv.filter(pl.col("submitter") == AMB).select("variation_id").unique(), left_on="VariationID", right_on="variation_id")
    amb_alleles = np.sort(ids["AlleleID"].to_numpy().astype(np.int64))
    print(f"\nAmbry has a submission on {len(amb_alleles):,} missense variants")
    matched = pr.build_strata(hist, fc, ev, gene_of, cov, "down",
                              allowed_controls=lambda row: amb_alleles if row["submitter"] == AMB else np.empty(0, dtype=np.int64))
    r2 = pr.fit_window(matched, pr.Q1, ev, exclude_s=waves, exclude_submitters=others)
    print("2. Ambry only, benign direction, Q1 windows, controls restricted to variants Ambry also holds")
    show("   ", r2)
    pl.DataFrame([{"check": "specificity_all_controls", **r1}, {"check": "portfolio_lab_matched", **r2}], infer_schema_length=None) \
        .write_csv(a.results / "promotion_ambry_checks.csv")


if __name__ == "__main__":
    main()
