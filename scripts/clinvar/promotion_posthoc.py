"""Post-hoc analyses for the VUS-promotion test, added after the primary results were seen.

    python scripts/clinvar/promotion_posthoc.py --data-dir $DATA_DIR [--results results]

Logged in the Deviations section of docs/protocol_promotion.md. Writes promotion_posthoc.csv.
"""

import argparse
import math
from pathlib import Path

import polars as pl

from rdv import promotion as pr
from promotion_prefit import load


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", required=True, type=Path)
    ap.add_argument("--results", type=Path, default=Path("results"))
    a = ap.parse_args()
    hist, fc, gene_of, sub, cov = load(a.data_dir)
    ev = pr.all_events(hist, fc, gene_of, sub)
    per = ev.group_by("s").len()
    median = per["len"].median()
    # symmetric rule, from counts alone: any interval with more than 5x the median number of promotions
    waves = tuple(per.filter(pl.col("len") > 5 * median)["s"].to_list())
    print(f"median promotions per interval {median}; waves (> 5x): {sorted(waves)}", flush=True)
    strata = {d: pr.build_strata(hist, fc, ev, gene_of, cov, d) for d in ("up", "down")}
    rows = []
    for w in pr.PRIMARY:
        for direction in ("up", "down"):
            r = {"analysis": "posthoc_symmetric_waves_out", "direction": direction, "focus": w.focus,
                 **pr.fit_window(strata[direction], w, ev, exclude_s=waves)}
            d, se = r[f"delta_{w.focus}"], r[f"se_{w.focus}"]
            r.update(delta=d, se=se, lo=d - 1.96 * se, hi=d + 1.96 * se, p=math.erfc(abs(d / se) / math.sqrt(2)), mde=2.8 * se)
            print(f"  {w.name:<26} {direction:<5} d{w.focus} {d:+.3f} [{r['lo']:+.3f}, {r['hi']:+.3f}] p={r['p']:.3f} "
                  f"MDE {r['mde']:.3f} strata {r['n_pre']}/{r['n_post']}", flush=True)
            rows.append(r)
    pl.DataFrame(rows, infer_schema_length=None).write_csv(a.results / "promotion_posthoc.csv")

    # Q1 benign direction in detail: every coefficient, then who supplies the cases
    q1 = next(r for r in rows if r["test"] == pr.Q1.name and r["direction"] == "down")
    print("\nQ1 down, symmetric waves out: coefficient per SD, pre -> post (change, SE)")
    for name in pr.COVARIATES:
        print(f"  {name:<9} {q1['pre_' + name]:+.3f} -> {q1['post_' + name]:+.3f}   change {q1['delta_' + name]:+.3f} (SE {q1['se_' + name]:.3f})")
    print(f"  sum of the three ML predictors: change {q1['delta_sum_ml']:+.3f} (SE {q1['se_sum_ml']:.3f})")
    meta = strata["down"].meta.filter(~pl.col("s").is_in(pl.Series(list(waves), dtype=pl.Utf8).implode()))
    w = pr.Q1
    meta = meta.with_columns(
        pl.when((pl.col("s") >= w.pre[0]) & (pl.col("s") <= w.pre[1])).then(pl.lit("pre"))
        .when((pl.col("s") >= w.post[0]) & (pl.col("s") <= w.post[1])).then(pl.lit("post")).otherwise(None).alias("period")
    ).drop_nulls("period")
    comp = meta.group_by("period", "submitter").len().sort("period", "len", descending=[False, True])
    print("\nlargest case suppliers per period (benign direction, waves removed):")
    for per_ in ("pre", "post"):
        print(" ", per_, [(r["submitter"][:22], r["len"]) for r in comp.filter(pl.col("period") == per_).head(5).iter_rows(named=True)])

    # leave one lab out, and lab-matched controls, for the same estimate
    big = comp.group_by("submitter").agg(pl.col("len").sum()).sort("len", descending=True).head(5)["submitter"].to_list()
    loo = []
    for lab in big:
        r = pr.fit_window(strata["down"], pr.Q1, ev, exclude_s=waves, exclude_submitters=(lab,))
        loo.append({"analysis": f"posthoc_leave_out:{lab[:30]}", "direction": "down", **r})
        print(f"  without {lab[:34]:<34} change in REVEL {r['delta_revel']:+.3f} (SE {r['se_revel']:.3f}), strata {r['n_pre']}/{r['n_post']}", flush=True)
    pl.DataFrame(loo, infer_schema_length=None).write_csv(a.results / "promotion_posthoc_leave_out.csv")

    # one lab at a time: is the change inside a lab, or in who supplies the cases?
    all_labs = strata["down"].meta["submitter"].unique().to_list()
    lab_rows = []
    for direction in ("down", "up"):
        top = (strata[direction].meta.filter(pl.col("submitter") != "unknown").group_by("submitter").len()
               .sort("len", descending=True).head(6)["submitter"].to_list())
        for lab in top:
            others = tuple(x for x in strata[direction].meta["submitter"].unique().to_list() if x != lab)
            for scope, excl in (("this lab only", others), ("everyone except this lab", (lab,))):
                try:
                    r = pr.fit_window(strata[direction], pr.Q1, ev, exclude_s=waves, exclude_submitters=excl)
                except Exception:
                    continue
                lab_rows.append({"direction": direction, "lab": lab, "scope": scope, **{k: r[k] for k in (
                    "n_pre", "n_post", "pre_revel", "post_revel", "delta_revel", "se_revel", "pre_am", "post_am", "delta_am", "se_am")}})
    labs_df = pl.DataFrame(lab_rows, infer_schema_length=None)
    labs_df.write_csv(a.results / "promotion_posthoc_by_lab.csv")
    print("\nQ1 window, waves removed, REVEL coefficient pre -> post, by lab:")
    for r in labs_df.filter(pl.col("n_pre") + pl.col("n_post") >= 150).iter_rows(named=True):
        print(f"  {r['direction']:<5} {r['lab'][:24]:<24} {r['scope']:<24} {r['pre_revel']:+.2f} -> {r['post_revel']:+.2f}  change {r['delta_revel']:+.2f} (SE {r['se_revel']:.2f})  n {r['n_pre']}/{r['n_post']}")

    # rolling placebo: the same pre/gap/post geometry slid across every position in the history
    snaps = [s_ for s_ in pr.snapshots(hist)][1:]  # event snapshots
    grid = []
    for direction in ("down", "up"):
        for i in range(len(snaps) - 19):
            win = pr.Window(f"roll_{snaps[i]}", "revel", (snaps[i], snaps[i + 7]), (snaps[i + 12], snaps[i + 19]))
            try:
                r = pr.fit_window(strata[direction], win, ev, exclude_s=waves)
            except Exception as e:  # a window with no strata on one side
                continue
            grid.append({"direction": direction, "pre_start": snaps[i], "post_start": snaps[i + 12],
                         "n_pre": r["n_pre"], "n_post": r["n_post"], "delta_revel": r["delta_revel"], "se_revel": r["se_revel"],
                         "delta_am": r["delta_am"], "se_am": r["se_am"]})
    g = pl.DataFrame(grid)
    g.write_csv(a.results / "promotion_rolling.csv")
    for direction in ("down", "up"):
        d = g.filter((pl.col("direction") == direction) & (pl.col("n_pre") >= 100) & (pl.col("n_post") >= 100))
        z = d["delta_revel"] / d["se_revel"]
        print(f"\nrolling windows, {direction}, symmetric waves out: {d.height} positions with 100+ strata per side; "
              f"change in REVEL: median {d['delta_revel'].median():+.3f}, 5th-95th pct {d['delta_revel'].quantile(0.05):+.3f} to "
              f"{d['delta_revel'].quantile(0.95):+.3f}; |z|>2 in {(z.abs() > 2).sum()} of {d.height}")


if __name__ == "__main__":
    main()
