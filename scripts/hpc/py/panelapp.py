"""Download green (high-evidence) genes from every public Genomics England PanelApp panel.

Writes one row per (gene, panel) with mode of inheritance, so the gene set and
the AD/AR annotation used later come from the same snapshot.
"""

import argparse
import csv
import json
import os
import time
import urllib.request

DEFAULT_BASE = "https://panelapp.genomicsengland.co.uk/api/v1/"
FIELDS = ["hgnc_symbol", "hgnc_id", "panel_id", "panel_name", "panel_version",
          "confidence_level", "mode_of_inheritance", "phenotypes"]


def get_json(url, retries=5):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                return json.load(resp)
        except Exception:
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)


def iter_panels(base):
    url = base + "panels/?page=1"
    while url:
        page = get_json(url)
        yield from page["results"]
        url = page.get("next")


def green_genes(panel_detail):
    for gene in panel_detail.get("genes", []):
        if str(gene.get("confidence_level")) != "3":
            continue
        data = gene.get("gene_data") or {}
        yield {
            "hgnc_symbol": data.get("hgnc_symbol") or gene.get("entity_name"),
            "hgnc_id": data.get("hgnc_id", ""),
            "panel_id": panel_detail["id"],
            "panel_name": panel_detail["name"],
            "panel_version": panel_detail.get("version", ""),
            "confidence_level": gene["confidence_level"],
            "mode_of_inheritance": gene.get("mode_of_inheritance", ""),
            "phenotypes": ";".join(gene.get("phenotypes") or []),
        }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True)
    ap.add_argument("--base", default=DEFAULT_BASE)
    args = ap.parse_args()
    base = args.base if args.base.endswith("/") else args.base + "/"

    rows = []
    panels = list(iter_panels(base))
    for i, panel in enumerate(panels, 1):
        detail = get_json(f"{base}panels/{panel['id']}/")
        rows.extend(green_genes(detail))
        if i % 50 == 0:
            print(f"  {i}/{len(panels)} panels", flush=True)

    tmp = args.out + ".tmp"
    with open(tmp, "w", newline="") as fh:
        writer = csv.DictWriter(fh, FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, args.out)
    genes = {r["hgnc_symbol"] for r in rows}
    print(f"{len(rows)} green gene-panel rows, {len(genes)} unique genes, {len(panels)} panels")


if __name__ == "__main__":
    main()
