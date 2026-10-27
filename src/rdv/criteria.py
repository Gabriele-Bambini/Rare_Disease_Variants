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

from rdv.clinvar import BLB, PLP, VUS, classify_significance

# Lookarounds, not \b: a code is often glued to its strength with an underscore ("PP3_Moderate"), and
# underscore is a word character, so \b would see no boundary after the code and skip it.
_CODE = r"PVS1|PS[1-4]|PM[1-6]|PP[1-5]|BA1|BS[1-4]|BP[1-7]"
CODE_RE = re.compile(rf"(?<![A-Za-z0-9])({_CODE})(?![0-9])", re.I)
# Strength words as submitters write them: spelled out, abbreviated (a survey of 116k submissions found
# _STR, _SUPP, _VSTR, _SUP, _MOD in the thousands), or mistyped ("moderated", "supporti").
_WORD = (r"very[\s-]?strong|vstr|vs|strong|str|stand[\s-]?alone|moderated?|mod|medium|"
         r"supporting|supporti|supp|sup")
_STRENGTH = {"verystrong": "very_strong", "vstr": "very_strong", "vs": "very_strong", "strong": "strong", "str": "strong",
             "standalone": "stand_alone", "moderate": "moderate", "moderated": "moderate", "mod": "moderate",
             "medium": "moderate", "supporting": "supporting", "supporti": "supporting", "supp": "supporting",
             "sup": "supporting"}
# What may follow a code: "PP3_Moderate", "PP3 (supporting)", "PM2sup", "PVS1_S" / "PM2-P", or a stated change
# of strength, "PM1 downgraded in strength to Supporting". 'vs' needs a bracket or underscore, since " vs " is
# also "versus".
_AFTER = re.compile(
    rf"(?P<sep>[\s:_\-\(\[]*)(?P<word>{_WORD})(?![A-Za-z0-9])"
    rf"|[_\-](?P<letter>[smp])(?![A-Za-z0-9])"
    rf"|\s+(?:(?:was|is|been)\s+)?(?:down|up)graded(?:\s+in\s+strength)?\s+(?:to|as|at)\s+(?P<moved>{_WORD})(?![A-Za-z0-9])",
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
    """Distinct, non-negated codes in `text`. A repeated code keeps a stated strength over a default one,
    and the strongest among stated ones: a bare "PM2" later in the text carries no strength information."""
    best: dict[str, Criterion] = {}
    for m in CODE_RE.finditer(text):
        code = m.group(1).upper()
        end, word, letter = m.end(), None, None
        mod = _AFTER.match(text, end)
        if mod:
            end = mod.end()
            word = mod.group("word") or mod.group("moved")
            letter = mod.group("letter")
            if (mod.group("word") or "").lower() == "vs" and not set(mod.group("sep")) & set(":_-(["):
                word, end = None, m.end()  # " vs " is "versus"
        if end == m.end() and text[end:end + 1].isalpha():
            continue  # a code glued to an unknown word ("PS4M", "PM2ndry") is not read as a code
        left = max(text.rfind(d, 0, m.start()) for d in _DELIMS) + 1
        right = min((i for i in (text.find(d, end) for d in _DELIMS) if i != -1), default=len(text))
        if _NEG_BEFORE.search(text[left:m.start()]) or _NEG_AFTER.search(text[end:right]):
            continue
        side = "B" if code[0] == "B" else "P"
        strength, explicit = _default(code), False
        if word:
            strength, explicit = _STRENGTH[re.sub(r"[\s-]+", "", word.lower())], True
        elif letter:
            strength, explicit = {"s": "strong", "m": "moderate", "p": "supporting"}[letter.lower()], True
        if side == "P" and strength == "stand_alone":  # not a pathogenic strength; keep the default
            strength, explicit = _default(code), False
        new = Criterion(code, side, strength, explicit)
        if code not in best or (explicit, POINTS[strength]) > (best[code].explicit, POINTS[best[code].strength]):
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
    """'PLP', 'BLB' or 'VUS' for a submitter's classification, 'other' for anything else (risk allele, conflicting...).

    Delegates to the normaliser that labels the variant universe, so a submission and the variant it sits on
    are read with one definition.
    """
    c = classify_significance(significance)
    return c if c in (PLP, BLB, VUS) else "other"


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


_TAGS = (("t", "tavtigian", False), ("t_pm2s", "tavtigian", True), ("r", "richards", False), ("r_pm2s", "richards", True))


def describe(text: str, significance: str | None) -> dict | None:
    """Everything the audit needs from one submission, or None if it lists no (non-negated) code.

    For each of four reference derivations (Tavtigian points or Richards 2015 rules, with an unqualified PM2
    read as moderate or as supporting) it records the derived class family, whether PP3 / BP4 is decisive,
    whether a non-computational code is decisive, and which codes are decisive. A single listed code is kept:
    BA1 alone is a complete benign justification, and dropping it would inflate BP4's share of benign labels.
    """
    crit = parse_criteria(text)
    if not crit:
        return None
    stated = stated_family(significance)
    pp3 = next((c for c in crit if c.code == "PP3"), None)
    bp4 = next((c for c in crit if c.code == "BP4"), None)
    pm2 = next((c for c in crit if c.code == "PM2"), None)
    out = {
        "stated": stated, "n_codes": len(crit), "points_t": points(crit),
        "codes": ",".join(f"{c.code}:{c.strength}{'' if c.explicit else '*'}" for c in crit),
        "has_pp3": pp3 is not None, "pp3_strength": pp3.strength if pp3 else None, "pp3_explicit": bool(pp3 and pp3.explicit),
        "has_bp4": bp4 is not None, "bp4_strength": bp4.strength if bp4 else None, "bp4_explicit": bool(bp4 and bp4.explicit),
        "has_pm2": pm2 is not None, "pm2_explicit": bool(pm2 and pm2.explicit),
    }
    for tag, rule, pm2_sup in _TAGS:
        dec = decisive_codes(crit, rule, stated, pm2_sup)
        out[f"derived_{tag}"] = derive(crit, rule, pm2_sup)
        out[f"dec_pp3_{tag}"] = "PP3" in dec
        out[f"dec_bp4_{tag}"] = "BP4" in dec
        out[f"dec_other_{tag}"] = any(c not in COMPUTATIONAL for c in dec)
        out[f"decisive_{tag}"] = ",".join(dec)
    return out
