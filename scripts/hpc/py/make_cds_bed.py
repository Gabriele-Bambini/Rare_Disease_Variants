"""Build a BED of MANE Select coding exons for a gene list.

Used to slice gnomAD remotely: only coding positions are needed for missense
variants, which keeps the slices to a few GB instead of whole gene bodies.
"""

import argparse
import csv
import gzip
import re
from collections import defaultdict

_ATTR = re.compile(r'(\S+) "([^"]*)"')


def read_genes(path):
    with open(path) as fh:
        return {row["hgnc_symbol"] for row in csv.DictReader(fh, delimiter="\t") if row["hgnc_symbol"]}


def mane_select_transcripts(summary_path, genes):
    """Map versioned Ensembl transcript id -> gene symbol for MANE Select only."""
    out = {}
    with gzip.open(summary_path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["MANE_status"] == "MANE Select" and row["symbol"] in genes:
                out[row["Ensembl_nuc"]] = row["symbol"]
    return out


def cds_intervals(gtf_path, transcripts, pad):
    by_chrom = defaultdict(list)
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            attrs = dict(_ATTR.findall(f[8]))
            symbol = transcripts.get(attrs.get("transcript_id"))
            if symbol is None:
                continue
            chrom = f[0] if f[0].startswith("chr") else "chr" + f[0]
            start = max(0, int(f[3]) - 1 - pad)   # GTF 1-based inclusive -> BED 0-based
            end = int(f[4]) + pad
            by_chrom[chrom].append((start, end, symbol))
    return by_chrom


def merge(intervals):
    merged = []
    for start, end, sym in sorted(intervals):
        if merged and start <= merged[-1][1]:
            prev = merged[-1]
            syms = prev[2] if sym in prev[2].split(",") else prev[2] + "," + sym
            merged[-1] = (prev[0], max(prev[1], end), syms)
        else:
            merged.append((start, end, sym))
    return merged


def _chrom_key(c):
    c = c[3:]
    return (0, int(c)) if c.isdigit() else (1, c)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genes", required=True, help="TSV with an hgnc_symbol column (panelapp.py output)")
    ap.add_argument("--mane-summary", required=True)
    ap.add_argument("--mane-gtf", required=True, help="MANE Ensembl genomic GTF")
    ap.add_argument("--pad", type=int, default=2)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    genes = read_genes(args.genes)
    transcripts = mane_select_transcripts(args.mane_summary, genes)
    by_chrom = cds_intervals(args.mane_gtf, transcripts, args.pad)

    n, bp = 0, 0
    with open(args.out, "w") as out:
        for chrom in sorted(by_chrom, key=_chrom_key):
            for start, end, syms in merge(by_chrom[chrom]):
                out.write(f"{chrom}\t{start}\t{end}\t{syms}\n")
                n += 1
                bp += end - start
    found = set(transcripts.values())
    print(f"{len(found)}/{len(genes)} genes with a MANE Select transcript; "
          f"{n} intervals, {bp / 1e6:.2f} Mb")
    missing = sorted(genes - found)
    if missing:
        print(f"no MANE Select for {len(missing)} genes, e.g. {', '.join(missing[:10])}")


if __name__ == "__main__":
    main()
