"""Extract gnomAD v4.1 allele frequencies for BED regions via remote tabix.

The full joint release is over a terabyte; this reads only the requested
regions over HTTPS and writes one gzipped TSV per chromosome (SNVs only).
Finished chromosomes are skipped on rerun.
"""

import argparse
import gzip
import os
from collections import defaultdict
from multiprocessing import Pool

ANCESTRIES = ["afr", "amr", "asj", "eas", "fin", "mid", "nfe", "sas", "remaining"]
INFO_FIELDS = (
    ["AC_joint", "AN_joint", "AF_joint", "nhomalt_joint",
     "grpmax_joint", "AF_grpmax_joint",
     "fafmax_faf95_max_joint", "fafmax_faf95_max_gen_anc_joint", "faf95_joint",
     "exomes_filters", "genomes_filters"]
    + [f"AF_joint_{a}" for a in ANCESTRIES]
    + [f"nhomalt_joint_{a}" for a in ANCESTRIES]
)


def _value(v):
    if v is None:
        return ""
    if isinstance(v, tuple):  # Number=A fields; sites VCFs are biallelic per row
        return ",".join("" if x is None else str(x) for x in v)
    if isinstance(v, bool):
        return "1" if v else "0"
    return str(v)


def slice_chrom(job):
    import pysam

    chrom, regions, template, out_dir = job
    out_path = os.path.join(out_dir, f"{chrom}.tsv.gz")
    if os.path.exists(out_path + ".done"):
        return chrom, "skipped"
    vcf = pysam.VariantFile(template.format(chrom=chrom))
    fields = [f for f in INFO_FIELDS if f in vcf.header.info]
    seen = set()
    n = 0
    tmp = out_path + ".tmp"
    with gzip.open(tmp, "wt") as out:
        out.write("\t".join(["chrom", "pos", "ref", "alt", "filter"] + fields) + "\n")
        for start, end in regions:
            for rec in vcf.fetch(chrom, start, end):
                alt = rec.alts[0] if rec.alts else ""
                if len(rec.ref) != 1 or len(alt) != 1:
                    continue
                key = (rec.pos, rec.ref, alt)
                if key in seen:
                    continue
                seen.add(key)
                info = rec.info
                row = [chrom, str(rec.pos), rec.ref, alt, ";".join(rec.filter.keys()) or "PASS"]
                row += [_value(info.get(f)) for f in fields]
                out.write("\t".join(row) + "\n")
                n += 1
    os.replace(tmp, out_path)
    open(out_path + ".done", "w").close()
    return chrom, f"{n} SNVs"


def read_bed(path):
    regions = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.strip() and not line.startswith(("#", "track")):
                f = line.split("\t")
                regions[f[0]].append((int(f[1]), int(f[2])))
    return regions


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--template", required=True, help="VCF URL with a {chrom} placeholder")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    regions = read_bed(args.bed)
    jobs = [(c, sorted(r), args.template, args.out_dir) for c, r in regions.items()]
    with Pool(min(args.workers, len(jobs) or 1)) as pool:
        for chrom, status in pool.imap_unordered(slice_chrom, jobs):
            print(f"  {chrom}: {status}", flush=True)


if __name__ == "__main__":
    main()
