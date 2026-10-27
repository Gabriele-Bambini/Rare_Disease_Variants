"""Shared HTTP opener. Some servers (PanelApp's CDN, for one) answer 403 to Python's default
User-Agent, so every request here identifies the project instead."""

import urllib.request

USER_AGENT = (
    "rare-disease-variants/0.1 (academic research; "
    "https://github.com/Gabriele-Bambini/Rare_Disease_Variants)"
)


def open_url(url, timeout=60):
    return urllib.request.urlopen(
        urllib.request.Request(url, headers={"User-Agent": USER_AGENT}), timeout=timeout
    )
