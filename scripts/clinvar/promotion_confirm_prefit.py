"""Outcome-blind preparation for the confirmation of the Ambry benign-direction lead.

    python scripts/clinvar/promotion_confirm_prefit.py --data-dir $DATA_DIR [--results results]

Counts cases and computes the minimum detectable lab difference from risk sets alone. Never relates a
score to a promotion. See docs/protocol_promotion_confirm.md.
"""

import argparse
from pathlib import Path

import polars as pl

from rdv import promotion as pr
from promotion_prefit import load

LAB = "Ambry Genetics"
FIRST_FRESH = "2025-07"  # first interval not used when the lead was found (the discovery windows end at 2025-01)


def fresh_intervals(ev: pl.DataFrame) -> list[str]:
    """Non-bulk intervals from FIRST_FRESH on, with the same rule as the discovery: not above 5x the median."""
    per = ev.group_by("s").len()
    waves = set(per.filter(pl.col("len") > 5 * per["len"].median())["s"].to_list())
    return sorted(s for s in per["s"].to_list() if s >= FIRST_FRESH and s not in waves)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", required=True, type=Path)
    ap.add_argument("--results", type=Path, default=Path("results"))
    a = ap.parse_args()
    hist, fc, gene_of, sub, cov = load(a.data_dir)
    ev = pr.all_events(hist, fc, gene_of, sub)
    fresh = fresh_intervals(ev)
    cases = ev.filter((pl.col("outcome") == "down") & pl.col("s").is_in(pl.Series(fresh, dtype=pl.Utf8).implode()))
    n_lab, n_other = cases.filter(pl.col("submitter") == LAB).height, cases.filter(pl.col("submitter") != LAB).height
    mde = pr.lab_contrast_mde(hist, fc, ev, gene_of, cov, "down", set(fresh), LAB, cap=pr.CAP)
    print("fresh intervals:", fresh)
    print(f"benign promotions in them: {LAB} {n_lab}, all other labs {n_other}")
    print("minimum detectable lab difference (model-based lower bound, per SD): " +
          "  ".join(f"{k[4:]} {v:.3f}" for k, v in mde.items() if k.startswith("mde_")))
    pl.DataFrame([{"intervals": ",".join(fresh), "n_lab": n_lab, "n_other": n_other, **mde}]).write_csv(a.results / "promotion_confirm_prefit.csv")


if __name__ == "__main__":
    main()
