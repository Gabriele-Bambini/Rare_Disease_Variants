"""VUS-promotion test: fit every analysis listed in docs/protocol_promotion.md.

    python scripts/clinvar/promotion_fit.py --data-dir $DATA_DIR [--results results]

Needs promotion_prefit.py's inputs. Writes promotion_primary.csv and promotion_secondary.csv.
"""

import argparse
import math
import time
from pathlib import Path

import numpy as np
import polars as pl

from rdv import promotion as pr
from rdv.rq1 import holm
from promotion_prefit import load


def pvalue(delta, se):
    return math.erfc(abs(delta / se) / math.sqrt(2))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", required=True, type=Path)
    ap.add_argument("--results", type=Path, default=Path("results"))
    a = ap.parse_args()
    hist, fc, gene_of, sub, cov = load(a.data_dir)
    ev = pr.all_events(hist, fc, gene_of, sub)
    waves = tuple(pr.wave_intervals(ev))
    top2 = tuple(ev.filter(pl.col("submitter") != "unknown").group_by("submitter").len().sort("len", descending=True).head(2)["submitter"])
    print("waves:", waves, "\ntop-2 labs:", top2, flush=True)

    t = time.time()
    strata = {d: pr.build_strata(hist, fc, ev, gene_of, cov, d) for d in ("up", "down")}
    print(f"strata: up {strata['up'].X.shape[0]:,}, down {strata['down'].X.shape[0]:,} ({time.time() - t:.0f}s)", flush=True)

    def run(label, direction, w, **kw):
        r = {"analysis": label, "direction": direction, "focus": w.focus,
             **pr.fit_window(strata[direction] if "strata" not in kw else kw.pop("strata"), w, ev, **kw)}
        d, se = r[f"delta_{w.focus}"], r[f"se_{w.focus}"]
        r["delta"], r["se"], r["lo"], r["hi"], r["p"], r["mde"] = d, se, d - 1.96 * se, d + 1.96 * se, pvalue(d, se), 2.8 * se
        print(f"  {label:<22} {w.name:<32} {direction:<5} d{w.focus} {d:+.3f} [{r['lo']:+.3f}, {r['hi']:+.3f}] "
              f"p={r['p']:.3f} MDE {r['mde']:.3f}  strata {r['n_pre']}/{r['n_post']}", flush=True)
        return r

    primary = [run("primary", "up", w) for w in pr.PRIMARY]
    adj = holm([r["p"] for r in primary])
    for r, q in zip(primary, adj):
        r["p_holm"] = q
    pl.DataFrame(primary, infer_schema_length=None).write_csv(a.results / "promotion_primary.csv")

    sec = []
    for w in pr.PRIMARY:
        sec.append(run("down_direction", "down", w))
        for direction in ("up", "down"):
            sec.append(run("waves_out", direction, w, exclude_s=waves))
            for cap in (250, 1000, None):
                sec.append(run(f"cap_{cap or 'none'}", direction, w, cap=cap))
            sec.append(run("exclude_top2_labs", direction, w, exclude_submitters=top2))
        for months in (-24, -36):
            sec.append(run(f"placebo_{months:+d}m", "up", pr.shifted(w, months)))

    # lab-matched controls: the case's own lab must also have a submission on the control
    scv = pl.read_parquet(a.data_dir / "derived/scv_flags.parquet", columns=["variation_id", "submitter"])
    uni = pl.read_parquet(a.data_dir / "derived/missense_universe.parquet", columns=["AlleleID", "VariationID"])
    lab_alleles = {}

    def allowed(row):
        lab = row["submitter"]
        if lab == "unknown":
            return np.empty(0, dtype=np.int64)
        if lab not in lab_alleles:
            ids = uni.join(scv.filter(pl.col("submitter") == lab).select("variation_id").unique(), left_on="VariationID", right_on="variation_id")
            lab_alleles[lab] = np.sort(ids["AlleleID"].to_numpy().astype(np.int64))
        return lab_alleles[lab]

    matched = pr.build_strata(hist, fc, ev, gene_of, cov, "up", allowed_controls=allowed)
    for w in pr.PRIMARY:
        sec.append(run("lab_matched_controls", "up", w, strata=matched))
    pl.DataFrame(sec, infer_schema_length=None).write_csv(a.results / "promotion_secondary.csv")
    print("done")


if __name__ == "__main__":
    main()
