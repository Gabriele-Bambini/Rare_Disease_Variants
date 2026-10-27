"""Reading ACMG/AMP evidence codes out of submitter text, and asking which labels hang on computational evidence.

A submission that lists its criteria, e.g. "PM3:Strong, PM2, PS3:Supporting, BP4", exposes the evidence
behind its classification. From the list we derive the class two ways (Tavtigian points and the 2015
combining rules) and ask a counterfactual: would the class survive without PP3 (or BP4)? A label whose
class does not survive is *decisive on a computational prediction*.

Parsing is heuristic (free text), so every function here is tested on hand-checked cases, and the
analysis only trusts submissions whose stated class agrees with the class derived from their own list.
"""

from __future__ import annotations

import re
from typing import NamedTuple

# Lookarounds, not \b: a code is often glued to its strength with an underscore ("PP3_Moderate"), and
# underscore is a word character, so \b would see no boundary after the code and skip it.
CODE_RE = re.compile(
    r"(?<![A-Za-z0-9])(PVS1|PS[1-4]|PM[1-6]|PP[1-5]|BA1|BS[1-4]|BP[1-7])(?![A-Za-z0-9])"
    r"(?:[\s:_\-\(\[]+(very[\s-]?strong|strong|moderate|supporting|stand[\s-]?alone|vs|mod|sup)(?![A-Za-z])"
    r"|_(s|m|p)(?![A-Za-z0-9]))?",
    re.I,
)
# Negation is scoped to the item (text between commas, semicolons, full stops, bars, newlines) holding the
# code, so "PP3 not met, PM2" does not also negate PM2.
_DELIMS = ",;|.\n"
_NEG_BEFORE = re.compile(r"\b(?:not|no|cannot|can't|without|fail\w*|lack\w*|unmet|insufficient|n/a)\b", re.I)
_NEG_AFTER = re.compile(
    r"not\s+(?:met|applied|applicable|used|satisfied|fulfilled|assigned|invoked|reached)|\bwas\s+not\b|\bis\s+not\b|"
    r"\bn/a\b|\bunmet\b|\bfail\w*", re.I)

POINTS = {"stand_alone": 8, "very_strong": 8, "strong": 4, "moderate": 2, "supporting": 1}
COMPUTATIONAL = {"PP3": "P", "BP4": "B"}


class Criterion(NamedTuple):
    code: str       # 'PM2'
    side: str       # 'P' (pathogenic) or 'B' (benign)
    strength: str   # very_strong | strong | moderate | supporting | stand_alone
    explicit: bool  # the text states the strength; otherwise it is the code's default


def _default(code: str) -> str:
    if code == "PVS1":
        return "very_strong"
    if code == "BA1":
        return "stand_alone"
    return {"PS": "strong", "PM": "moderate", "PP": "supporting", "BS": "strong", "BP": "supporting"}[code[:2]]


def parse_criteria(text: str) -> list[Criterion]:
    """Distinct, non-negated codes in `text`; a repeated code keeps its strongest stated strength."""
    best: dict[str, Criterion] = {}
    for m in CODE_RE.finditer(text):
        code = m.group(1).upper()
        left = max(text.rfind(d, 0, m.start()) for d in _DELIMS) + 1
        right = min((i for i in (text.find(d, m.end()) for d in _DELIMS) if i != -1), default=len(text))
        if _NEG_BEFORE.search(text[left:m.start()]) or _NEG_AFTER.search(text[m.end():right]):
            continue
        side = "B" if code[0] == "B" else "P"
        strength, explicit = _default(code), False
        word = (m.group(2) or "").lower()
        if word:
            word = re.sub(r"[\s-]+", "_", word)
            strength = {"very_strong": "very_strong", "verystrong": "very_strong", "vs": "very_strong", "strong": "strong", "moderate": "moderate",
                        "mod": "moderate", "supporting": "supporting", "sup": "supporting", "stand_alone": "stand_alone",
                        "standalone": "stand_alone"}[word]
            explicit = True
        elif m.group(3):
            strength, explicit = {"s": "strong", "m": "moderate", "p": "supporting"}[m.group(3).lower()], True
        if side == "P" and strength == "stand_alone":  # not a pathogenic strength; keep the default
            strength, explicit = _default(code), False
        new = Criterion(code, side, strength, explicit)
        if code not in best or POINTS[strength] > POINTS[best[code].strength]:
            best[code] = new
    return list(best.values())


def _effective(c: Criterion, pm2_supporting: bool) -> str:
    # ClinGen SVI (2020) downgraded PM2 to supporting; an unqualified PM2 may mean either
    return "supporting" if pm2_supporting and c.code == "PM2" and not c.explicit else c.strength


def points(criteria: list[Criterion], pm2_supporting: bool = False) -> int:
    """Tavtigian et al. 2020: pathogenic evidence adds, benign subtracts (8 / 4 / 2 / 1)."""
    return sum((1 if c.side == "P" else -1) * POINTS[_effective(c, pm2_supporting)] for c in criteria)


def tavtigian_class(p: int) -> str:
    return "P" if p >= 10 else "LP" if p >= 6 else "VUS" if p >= 0 else "LB" if p >= -6 else "B"


def richards_class(criteria: list[Criterion], pm2_supporting: bool = False) -> str:
    """Richards et al. 2015, Table 5 (combining rules). Conflicting pathogenic and benign evidence gives VUS."""
    n = {"P": {"very_strong": 0, "strong": 0, "moderate": 0, "supporting": 0},
         "B": {"stand_alone": 0, "very_strong": 0, "strong": 0, "moderate": 0, "supporting": 0}}
    for c in criteria:
        n[c.side][_effective(c, pm2_supporting)] += 1
    vs, s, m, p = (n["P"][k] for k in ("very_strong", "strong", "moderate", "supporting"))
    ba = n["B"]["stand_alone"] + n["B"]["very_strong"]
    bs, bp = n["B"]["strong"], n["B"]["moderate"] + n["B"]["supporting"]
    path = (vs >= 1 and (s >= 1 or m >= 2 or (m == 1 and p >= 1) or p >= 2)) or s >= 2 or (
        s == 1 and (m >= 3 or (m == 2 and p >= 2) or (m == 1 and p >= 4)))
    lpath = (vs >= 1 and m >= 1) or (s == 1 and m in (1, 2)) or (s == 1 and p >= 2) or m >= 3 or (
        m == 2 and p >= 2) or (m == 1 and p >= 4)
    ben = ba >= 1 or bs >= 2
    lben = (bs == 1 and bp >= 1) or bp >= 2
    if (path or lpath) and (ben or lben):
        return "VUS"
    return "P" if path else "LP" if lpath else "B" if ben else "LB" if lben else "VUS"


def family(cls: str) -> str:
    return "PLP" if cls in ("P", "LP") else "BLB" if cls in ("B", "LB") else "VUS"


def derive(criteria: list[Criterion], rule: str, pm2_supporting: bool = False) -> str:
    """Class family ('PLP', 'BLB' or 'VUS') implied by a criteria list under 'tavtigian' or 'richards'."""
    if rule == "tavtigian":
        return family(tavtigian_class(points(criteria, pm2_supporting)))
    if rule == "richards":
        return family(richards_class(criteria, pm2_supporting))
    raise ValueError(rule)


def stated_family(significance: str | None) -> str:
    s = (significance or "").strip().lower()
    if s.startswith(("pathogenic", "likely pathogenic")):
        return "PLP"
    if s.startswith(("benign", "likely benign")):
        return "BLB"
    return "VUS" if s.startswith("uncertain") else "other"


def is_decisive(criteria: list[Criterion], code: str, rule: str, stated: str, pm2_supporting: bool = False) -> bool:
    """The stated class agrees with the listed criteria, but loses its family if `code` is removed."""
    if stated not in ("PLP", "BLB") or derive(criteria, rule, pm2_supporting) != stated:
        return False
    if code not in {c.code for c in criteria}:
        return False
    rest = [c for c in criteria if c.code != code]
    return derive(rest, rule, pm2_supporting) != stated


def decisive_codes(criteria: list[Criterion], rule: str, stated: str, pm2_supporting: bool = False) -> list[str]:
    """Every listed code whose removal would change the class family (computational or not)."""
    return [c.code for c in criteria if is_decisive(criteria, c.code, rule, stated, pm2_supporting)]
