"""Outcome-blind preparation for the VUS-promotion test: event counts and minimum detectable effects.

    python scripts/clinvar/promotion_prefit.py --data-dir $DATA_DIR [--results results]

Uses labels (who was promoted, when, by whom) and the covariates of the unpromoted VUS only. It never
relates a score to a promotion. Run and committed before promotion_fit.py.
"""

import argparse
from pathlib import Path

import polars as pl

from rdv import promotion as pr


def load(d: Path):
    der = d / "derived"
    fc = pl.read_parquet(der / "first_confident.parquet", columns=["allele_id", "snapshot"]).rename({"snapshot": "first_snapshot"})
    uni = pl.read_parquet(der / "missense_universe.parquet", columns=["AlleleID", "VariationID", "gene"]).rename({"AlleleID": "allele_id"})
    sub = pl.read_parquet(der / "first_confident_submitter.parquet", columns=["allele_id", "submitter"])
    scores = pl.read_parquet(der / "scores.parquet")
    af = pl.read_parquet(der / "af_proxy.parquet")
    cov = pr.covariate_table(scores, af, uni)
    gene_of = uni.select("allele_id", "gene")
    return der / "history", fc, gene_of, sub, cov


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", required=True, type=Path)
    ap.add_argument("--results", type=Path, default=Path("results"))
    a = ap.parse_args()
    hist, fc, gene_of, sub, cov = load(a.data_dir)
    ev = pr.all_events(hist, fc, gene_of, sub)
    a.results.mkdir(parents=True, exist_ok=True)

    per = ev.group_by("s", "outcome").len().pivot(on="outcome", index="s", values="len").fill_null(0).sort("s")
    # promotions to a conflicting label are outside the design; count them so the limit is visible
    snaps = pr.snapshots(hist)
    conflicts = []
    for p, s in zip(snaps[:-1], snaps[1:]):
        a_ = pl.read_parquet(hist / f"snapshot={p}.parquet", columns=["allele_id", "label"]).filter(pl.col("label") == "VUS")
        b_ = pl.read_parquet(hist / f"snapshot={s}.parquet", columns=["allele_id", "label"]).rename({"label": "ls"})
        conflicts.append((s, a_.join(b_, on="allele_id").filter(pl.col("ls") == "CONFLICT").height))
    per = per.join(pl.DataFrame(conflicts, schema=["s", "vus_to_conflict"], orient="row"), on="s", how="left")
    per.write_csv(a.results / "promotion_events_by_interval.csv")
    waves = pr.wave_intervals(ev)
    print("wave intervals (top 10% by promotions):", sorted(waves))

    rows = []
    for w in pr.PRIMARY:
        for direction in ("up", "down"):
            inwin = lambda lo, hi: (pl.col("s") >= pl.lit(lo)) & (pl.col("s") <= pl.lit(hi))  # noqa: E731
            d_ = pl.col("outcome") == direction
            n_pre = ev.filter(inwin(*w.pre) & d_).height
            n_post = ev.filter(inwin(*w.post) & d_).height
            mde = pr.pool_information(hist, fc, ev, gene_of, cov, direction, w, cap=pr.CAP)
            rows.append({"test": w.name, "direction": direction, "focus": w.focus, "events_pre": n_pre, "events_post": n_post, **mde})
            print(f"  {w.name:<26} {direction:<5} events {n_pre:>5}/{n_post:<5}  MDE(change, per SD) "
                  f"revel {mde['mde_revel']:.3f}  am {mde['mde_am']:.3f}  esm1b {mde['mde_esm1b']:.3f}  grantham {mde['mde_grantham']:.3f}", flush=True)
    pl.DataFrame(rows).write_csv(a.results / "promotion_prefit.csv")


if __name__ == "__main__":
    main()
