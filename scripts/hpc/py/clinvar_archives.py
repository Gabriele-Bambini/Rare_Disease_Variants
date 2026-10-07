"""List ClinVar monthly archive URLs (variant_summary / submission_summary).

The archive directory keeps recent files at the top level and may move older
ones into per-year subdirectories, so both levels are crawled.
"""

import argparse
import re
import urllib.parse

from listing import list_links

DEFAULT_BASE = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/"


def archive_urls(base, kinds, year_from, months):
    pattern = re.compile(rf"^({'|'.join(map(re.escape, kinds))})_(\d{{4}})-(\d{{2}})\.txt\.gz$")
    found = {}

    def scan(url, depth):
        for href in list_links(url):
            name = href.rstrip("/").rsplit("/", 1)[-1]
            m = pattern.match(name)
            if m:
                kind, year, month = m.group(1), int(m.group(2)), m.group(3)
                if year >= year_from and month in months:
                    # Top-level and year-folder copies of the same file collapse here.
                    found.setdefault((kind, year, month), urllib.parse.urljoin(url, href))
            elif depth == 0 and re.fullmatch(r"\d{4}/?", href):
                if int(href[:4]) >= year_from:
                    scan(urllib.parse.urljoin(url, href.rstrip("/") + "/"), depth + 1)

    scan(base if base.endswith("/") else base + "/", 0)
    return [found[k] for k in sorted(found)]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base", default=DEFAULT_BASE)
    ap.add_argument("--kinds", default="variant_summary,submission_summary")
    ap.add_argument("--from-year", type=int, default=2017)
    ap.add_argument("--months", default="01,04,07,10")
    args = ap.parse_args()
    urls = archive_urls(args.base, args.kinds.split(","), args.from_year,
                        set(args.months.split(",")))
    if not urls:
        raise SystemExit(f"no archive files found under {args.base}")
    print("\n".join(urls))


if __name__ == "__main__":
    main()
