"""Download AlphaFold DB models for the UniProt accessions of a gene list.

Uses the prediction API to resolve file URLs, so it keeps working across
AlphaFold DB model versions. Long proteins come as several fragments (F1, F2…).
"""

import argparse
import csv
import gzip
import json
import os
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor

API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"


def accessions_for_genes(uniprot_tsv, genes):
    opener = gzip.open if uniprot_tsv.endswith(".gz") else open
    accs = {}
    with opener(uniprot_tsv, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene = (row.get("Gene Names (primary)") or "").split(";")[0].strip()
            if gene in genes:
                accs[row["Entry"]] = gene
    return accs


def _get(url, retries=4):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                return resp.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if attempt == retries - 1:
                raise
        except Exception:
            if attempt == retries - 1:
                raise
        time.sleep(2 ** attempt)


def fetch(acc, out_dir, api=API):
    marker = os.path.join(out_dir, f"{acc}.done")
    if os.path.exists(marker):
        return acc, "skipped"
    body = _get(api.format(acc=acc))
    if body is None:
        return acc, "not in AlphaFold DB"
    entries = json.loads(body)
    for entry in entries:
        url = entry.get("pdbUrl") or entry.get("cifUrl")
        data = _get(url)
        if data is None:
            return acc, f"missing file {url}"
        with open(os.path.join(out_dir, url.rsplit("/", 1)[-1]), "wb") as fh:
            fh.write(data)
    with open(marker, "w") as fh:
        json.dump([{k: e.get(k) for k in ("entryId", "latestVersion", "uniprotStart", "uniprotEnd")}
                   for e in entries], fh)
    return acc, f"{len(entries)} model(s)"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genes", required=True, help="TSV with an hgnc_symbol column")
    ap.add_argument("--uniprot", required=True, help="UniProt TSV with Entry and 'Gene Names (primary)'")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--api", default=API, help="prediction API URL with an {acc} placeholder")
    args = ap.parse_args()

    with open(args.genes) as fh:
        genes = {r["hgnc_symbol"] for r in csv.DictReader(fh, delimiter="\t")}
    accs = accessions_for_genes(args.uniprot, genes)
    print(f"{len(accs)} UniProt accessions for {len(genes)} genes")
    os.makedirs(args.out_dir, exist_ok=True)

    failed = 0
    with ThreadPoolExecutor(args.workers) as pool:
        for i, (acc, status) in enumerate(pool.map(lambda a: fetch(a, args.out_dir, args.api), sorted(accs)), 1):
            if status.startswith(("not", "missing")):
                failed += 1
                print(f"  {acc} ({accs[acc]}): {status}")
            if i % 500 == 0:
                print(f"  {i}/{len(accs)} done", flush=True)
    print(f"finished; {failed} accessions without a model")


if __name__ == "__main__":
    main()
