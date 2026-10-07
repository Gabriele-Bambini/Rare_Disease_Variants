"""Discover files in HTTP directory listings (NCBI FTP over HTTPS and similar)."""

import argparse
import re
import sys
import urllib.parse
import urllib.request

_HREF = re.compile(r'href="([^"?#]+)"', re.IGNORECASE)


def list_links(url, timeout=60):
    """Return the hrefs of an HTTP directory listing, relative to ``url``."""
    if not url.endswith("/"):
        url += "/"
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        html = resp.read().decode("utf-8", errors="replace")
    links = []
    for href in _HREF.findall(html):
        href = urllib.parse.unquote(href)
        if href.startswith(("http://", "https://", "/", "..")):
            # Absolute or parent links are not children of this directory,
            # except absolute links that point inside it.
            absolute = urllib.parse.urljoin(url, href)
            if not absolute.startswith(url) or absolute == url:
                continue
            href = absolute[len(url):]
        links.append(href)
    return links


def _version_key(name):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", name)]


def latest(url, pattern):
    """Full URL of the highest-versioned file in ``url`` matching ``pattern``."""
    regex = re.compile(pattern)
    matches = [h for h in list_links(url) if regex.fullmatch(h)]
    if not matches:
        raise SystemExit(f"no file matching {pattern!r} in {url}")
    return urllib.parse.urljoin(url if url.endswith("/") else url + "/",
                                max(matches, key=_version_key))


def main():
    ap = argparse.ArgumentParser(description="Print the latest file matching a regex in a listing.")
    ap.add_argument("url")
    ap.add_argument("pattern", help="regex matched against the whole file name")
    args = ap.parse_args()
    print(latest(args.url, args.pattern))


if __name__ == "__main__":
    sys.exit(main())
