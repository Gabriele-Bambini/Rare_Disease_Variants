"""RQ1: did endorsed predictors gain on new ClinVar labels? Implements docs/protocol.md
(primary and secondary analyses) plus the post-hoc analyses logged in its Deviations section.

    python scripts/clinvar/rq1_did.py --data-dir $DATA_DIR [--results results] [--replicates 2000]

Needs derived/scores.parquet (join_scores.py) and derived/first_confident_submitter.parquet
(attach_submitters.py). Writes rq1_primary.csv, rq1_secondary.csv, rq1_posthoc.csv, rq1_positive_control.csv,
rq1_event_study.csv.
"""

import argparse
import time
from pathlib import Path

import polars as pl

from rdv import rq1
from rdv.rq1 import Comparison, P1, P2, PRIMARY


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", required=True, type=Path)
    ap.add_argument("--results", type=Path, default=Path("results"))
    ap.add_argument("--replicates", type=int, default=rq1.REPLICATES)
    a = ap.parse_args()
    d = a.data_dir / "derived"
    B = a.replicates
    cohort = rq1.build_cohort(pl.read_parquet(d / "first_confident_submitter.parquet"), pl.read_parquet(d / "scores.parquet"))
    top2 = tuple(
        cohort.filter(pl.col("submitter") != "unknown").group_by("submitter").len()
        .sort("len", descending=True).head(2)["submitter"].to_list()
    )
    print(f"cohort: {cohort.height:,} variants; two largest labs: {top2}", flush=True)
    a.results.mkdir(parents=True, exist_ok=True)

    def run(label, test, **kw):
        t = time.time()
        r = {"analysis": label, **rq1.run_did(cohort, test, replicates=B, **kw)}
        e = r.get("estimate")
        print(f"  {label:<26} {test.name:<34} " + (f"DiD {e:+.4f} [{r['lo']:+.4f}, {r['hi']:+.4f}] "
              f"pairs {int(r['pairs_pre']):,}/{int(r['pairs_post']):,} ({time.time() - t:.0f}s)" if e is not None else "no data"),
              flush=True)
        return r

    primary = []
    for test in PRIMARY:
        primary.append(run("primary", test))
        for cap in rq1.CAP_SENSITIVITY:
            primary.append(run(f"cap_{cap or 'none'}", test, cap=cap))
        primary.append(run("exclude_top2_labs", test, exclude_submitters=top2))
        for months in rq1.PLACEBO_SHIFTS:
            primary.append(run(f"placebo_{months:+d}m", rq1.shifted(test, months)))
    df = pl.DataFrame(primary, infer_schema_length=None)
    prim = df.filter(pl.col("analysis") == "primary")
    adj = rq1.holm(prim["p_boot"].to_list())
    df = df.join(pl.DataFrame({"test": prim["test"], "p_holm": adj}), on="test", how="left").with_columns(
        pl.when(pl.col("analysis") == "primary").then(pl.col("p_holm")).otherwise(None).alias("p_holm")
    )
    df.write_csv(a.results / "rq1_primary.csv")

    secondary = []
    for test in PRIMARY:
        for control in ("grantham", "blosum62"):
            secondary.append(run(f"control_{control}", Comparison(test.name, test.treated, control, test.pre, test.post)))
        secondary.append(run("stars_2plus", test, min_stars=2))
        secondary.append(run("from_vus_only", test, routes=("from_vus",)))
    secondary.append(run("revel_as_comparator", Comparison(P2.name, "revel", "esm1b", P2.pre, P2.post)))
    ext = {"pre": ("2023-10", "2025-01"), "post": ("2025-04", "2026-10")}
    for treated in ("am", "esm1b"):
        secondary.append(run("ext_2025", Comparison(f"ext2025_{treated}", treated, "grantham", **ext)))
    for lab in cohort.filter(pl.col("submitter") != "unknown").group_by("submitter").len().sort("len", descending=True).head(6)["submitter"]:
        for test in PRIMARY:
            secondary.append(run(f"lab:{lab[:28]}", test, include_submitters=(lab,)))
    pl.DataFrame(secondary, infer_schema_length=None).write_csv(a.results / "rq1_secondary.csv")

    # Post-hoc analyses, added after the primary results were seen (protocol Deviations log).
    posthoc = [run("posthoc_clean_boundary", rq1.clean_boundary(t)) for t in PRIMARY]
    pl.DataFrame(posthoc, infer_schema_length=None).write_csv(a.results / "rq1_posthoc.csv")
    scv = pl.read_parquet(d / "scv_flags.parquet")
    flagged = rq1.attach_labelling_scv_flags(
        cohort, pl.read_parquet(d / "first_confident_submitter.parquet"), scv)
    contrasts = []
    for flag in ("pp3_bp4", "insilico_generic"):
        for treated in ("revel", "am"):
            r = rq1.group_contrast(flagged, flag, treated=treated, replicates=B)
            print(f"  positive control {flag:<17} {treated} minus esm1b: {r['gap_difference']:+.4f} "
                  f"[{r['lo']:+.4f}, {r['hi']:+.4f}] ({r['variants_A']:,} citing vs {r['variants_B']:,} other)", flush=True)
            contrasts.append(r)
    pl.DataFrame(contrasts, infer_schema_length=None).write_csv(a.results / "rq1_positive_control.csv")

    es = rq1.event_study(cohort, ["revel", "am", "esm1b", "grantham", "blosum62"], replicates=B)
    es.write_csv(a.results / "rq1_event_study.csv")
    print("wrote rq1_primary, rq1_secondary, rq1_posthoc, rq1_positive_control and rq1_event_study tables")


if __name__ == "__main__":
    main()
