"""Build the missense universe and the per-snapshot label history (RQ1, week 2).

    python scripts/clinvar/build_history.py --data-dir $DATA_DIR [--workers 3]

Reads  $DATA_DIR/clinvar/current/variant_summary.txt.gz
       $DATA_DIR/clinvar/archive/variant_summary_YYYY-MM.txt.gz   (download_all.sh clinvar_archives)
Writes $DATA_DIR/derived/missense_universe.parquet
       $DATA_DIR/derived/history/snapshot=YYYY-MM.parquet           (skipped if present)
       $DATA_DIR/derived/first_confident.parquet, transitions.parquet
"""

import argparse
import multiprocessing
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import polars as pl

from rdv.clinvar import missense_universe, snapshot_month
from rdv.history import derive, entry_route, extract_snapshot, snapshot_path


def _extract(args):
    src, out, universe_path = args
    t = time.time()
    ids = pl.read_parquet(universe_path, columns=["AlleleID"]).rename({"AlleleID": "allele_id"})
    df = extract_snapshot(src, ids)
    df.write_parquet(out, compression="zstd")
    return src.name, df.height, time.time() - t


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--data-dir", required=True, type=Path)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--only", help="comma-separated snapshot months, for testing")
    a = ap.parse_args()

    derived = a.data_dir / "derived"
    hist = derived / "history"
    hist.mkdir(parents=True, exist_ok=True)

    universe_path = derived / "missense_universe.parquet"
    if not universe_path.exists():
        t = time.time()
        u = missense_universe(a.data_dir / "clinvar/current/variant_summary.txt.gz")
        u.write_parquet(universe_path, compression="zstd")
        print(f"universe: {u.height:,} missense variants ({time.time() - t:.0f}s)", flush=True)

    only = set(a.only.split(",")) if a.only else None
    jobs = []
    for src in sorted((a.data_dir / "clinvar/archive").glob("variant_summary_*.txt.gz")):
        month = snapshot_month(src)
        out = snapshot_path(hist, month)
        if (only and month not in only) or out.exists():
            continue
        jobs.append((src, out, universe_path))
    print(f"{len(jobs)} snapshots to extract", flush=True)
    # spawn, not fork: polars keeps a thread pool, and forking after it has been used deadlocks.
    with ProcessPoolExecutor(a.workers, mp_context=multiprocessing.get_context("spawn")) as ex:
        for f in as_completed([ex.submit(_extract, j) for j in jobs]):
            name, n, secs = f.result()
            print(f"  {name}: {n:,} universe variants ({secs:.0f}s)", flush=True)

    if only:
        return
    first, trans = derive(hist)
    entry_route(first).write_parquet(derived / "first_confident.parquet", compression="zstd")
    trans.write_parquet(derived / "transitions.parquet", compression="zstd")
    print(f"first_confident: {first.height:,} variants; transitions: {trans.height:,}")


if __name__ == "__main__":
    main()
