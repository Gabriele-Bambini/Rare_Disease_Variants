"""Confirmation of the Ambry benign-direction lead on intervals not used to find it.

    python scripts/clinvar/promotion_confirm.py --data-dir $DATA_DIR [--results results]

Implements docs/protocol_promotion_confirm.md. Writes promotion_confirm.csv.
"""

import argparse
import math
from pathlib import Path

import numpy as np
import polars as pl

from rdv import promotion as pr
from promotion_prefit import load
from promotion_confirm_prefit import LAB, fresh_intervals


def one_sided_p(diff: float, se: float) -> float:
    return 0.5 * math.erfc((diff / se) / math.sqrt(2))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", required=True, type=Path)
    ap.add_argument("--results", type=Path, default=Path("results"))
    a = ap.parse_args()
    hist, fc, gene_of, sub, cov = load(a.data_dir)
    ev = pr.all_events(hist, fc, gene_of, sub)
    fresh = fresh_intervals(ev)
    print("intervals:", fresh, flush=True)

    scv = pl.read_parquet(a.data_dir / "derived/scv_flags.parquet", columns=["variation_id", "submitter"])
    uni = pl.read_parquet(a.data_dir / "derived/missense_universe.parquet", columns=["AlleleID", "VariationID"])
    catalogue = {}

    def own_catalogue(row):
        lab = row["submitter"]
        if lab == "unknown":
            return np.empty(0, dtype=np.int64)
        if lab not in catalogue:
            ids = uni.join(scv.filter(pl.col("submitter") == lab).select("variation_id").unique(), left_on="VariationID", right_on="variation_id")
            catalogue[lab] = np.sort(ids["AlleleID"].to_numpy().astype(np.int64))
        return catalogue[lab]

    rows = []
    for label, kw in (("primary_lab_matched", {"allowed_controls": own_catalogue}), ("secondary_unmatched", {})):
        st = pr.build_strata(hist, fc, ev, gene_of, cov, "down", intervals=fresh, **kw)
        r = pr.fit_lab_contrast(st, ev, set(fresh), LAB)
        r["analysis"], r["p_one_sided_revel"] = label, one_sided_p(r["diff_revel"], r["se_revel"])
        rows.append(r)
        print(f"{label}: Ambry cases {r['n_lab']}, other labs {r['n_other']}")
        for name in ("revel", "am", "esm1b", "grantham"):
            print(f"    {name:<9} others {r['other_' + name]:+.2f}   Ambry {r['lab_' + name]:+.2f}   difference {r['diff_' + name]:+.2f} (SE {r['se_' + name]:.2f})")
        print(f"    one-sided p for a larger REVEL coefficient at Ambry: {r['p_one_sided_revel']:.4f}", flush=True)
    pl.DataFrame(rows, infer_schema_length=None).write_csv(a.results / "promotion_confirm.csv")
    prim, sec = rows
    confirmed = prim["p_one_sided_revel"] < 0.05 and sec["diff_revel"] > 0
    print("\nVERDICT by the pre-specified rule:", "CONFIRMED" if confirmed else "NOT CONFIRMED")


if __name__ == "__main__":
    main()
