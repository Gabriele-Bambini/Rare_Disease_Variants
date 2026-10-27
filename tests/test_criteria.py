import pytest

from rdv.criteria import (
    Criterion, decisive_codes, derive, is_decisive, parse_criteria, points, richards_class, stated_family, tavtigian_class,
)


def codes(text):
    return {c.code: (c.side, c.strength, c.explicit) for c in parse_criteria(text)}


def test_parse_a_real_cegat_style_list():
    got = codes("HFE: PM3:Strong, PM2, PS3:Supporting, BP4")
    assert got == {"PM3": ("P", "strong", True), "PM2": ("P", "moderate", False),
                   "PS3": ("P", "supporting", True), "BP4": ("B", "supporting", False)}


@pytest.mark.parametrize("text,expected", [
    ("PVS1_Strong, PM2_Supporting, PP3_Moderate", {"PVS1": "strong", "PM2": "supporting", "PP3": "moderate"}),
    ("PS4_VeryStrong; PMIDs: 31939074", {"PS4": "very_strong"}),
    ("PVS1_S, PS3_M, PM2_P", {"PVS1": "strong", "PS3": "moderate", "PM2": "supporting"}),
    ("PP3 (supporting) and PM1 moderate", {"PP3": "supporting", "PM1": "moderate"}),
    ("BA1, BP4", {"BA1": "stand_alone", "BP4": "supporting"}),
    ("BP4_Strong", {"BP4": "strong"}),
])
def test_parse_strengths(text, expected):
    assert {c.code: c.strength for c in parse_criteria(text)} == expected


@pytest.mark.parametrize("text,kept", [
    ("PM2 supporting, PP3 not met", {"PM2"}),
    ("BP4 is not applicable. PM2", {"PM2"}),
    ("no PP3 applied; PS3 met", {"PS3"}),
    ("PP3 was not applied, PM2, PS4", {"PM2", "PS4"}),
    ("PP3, PM2", {"PP3", "PM2"}),
])
def test_negated_codes_are_dropped(text, kept):
    assert set(codes(text)) == kept


def test_repeated_code_prefers_a_stated_strength_then_the_strongest_and_pathogenic_stand_alone_is_ignored():
    assert codes("PM2, PM2_Supporting")["PM2"] == ("P", "supporting", True)  # a bare mention carries no strength
    assert codes("PM2_Supporting, PM2")["PM2"] == ("P", "supporting", True)
    assert codes("PM2_Supporting, PM2_Moderate")["PM2"][1] == "moderate"
    assert codes("PS3 stand-alone")["PS3"] == ("P", "strong", False)  # not a pathogenic strength: keep the default


@pytest.mark.parametrize("text,expected", [
    ("PM1_STR, PM2_SUP, PP3_MOD, PM3_VSTR", {"PM1": "strong", "PM2": "supporting", "PP3": "moderate", "PM3": "very_strong"}),
    ("PS4_supp, PP1_Str", {"PS4": "supporting", "PP1": "strong"}),
    ("PM2 moderated, PS4 supporting", {"PM2": "moderate", "PS4": "supporting"}),
    ("PP3 (medium pathogenic)", {"PP3": "moderate"}),
    ("assigned evidence codes ['BS3_Supporti", {"BS3": "supporting"}),
    ("PM2sup, PM1mod", {"PM2": "supporting", "PM1": "moderate"}),
    ("PM1 downgraded in strength to Supporting", {"PM1": "supporting"}),
    ("PS2 downgraded to moderate", {"PS2": "moderate"}),
    ("PP1 upgraded in strength to Strong", {"PP1": "strong"}),
    ("PM1-Supporting => PM1 downgraded in strength to Supporting", {"PM1": "supporting"}),
    ("PS3-Very Strong => PS3 upgraded in strength to Very Strong", {"PS3": "very_strong"}),
    ("PP3 vs BP4", {"PP3": "supporting", "BP4": "supporting"}),  # " vs " is "versus", not very strong
    ("PVS1_VS, PS4-P, PM2 (VS)", {"PVS1": "very_strong", "PS4": "supporting", "PM2": "very_strong"}),
])
def test_parse_observed_strength_notation(text, expected):
    assert {c.code: c.strength for c in parse_criteria(text)} == expected


@pytest.mark.parametrize("text,kept", [
    ("ACMG codes: PS2; PS4M; PM2", {"PS2", "PM2"}),  # glued to an unknown suffix: not read as a code
    ("PM2ndry findings", set()),
    ("PP31 and BA10", set()),
])
def test_codes_glued_to_unknown_text_are_not_codes(text, kept):
    assert set(codes(text)) == kept


@pytest.mark.parametrize("p,cls", [(12, "P"), (10, "P"), (9, "LP"), (6, "LP"), (5, "VUS"), (0, "VUS"), (-1, "LB"), (-6, "LB"), (-7, "B")])
def test_tavtigian_boundaries(p, cls):
    assert tavtigian_class(p) == cls


def test_points_and_derived_family_for_the_cegat_example():
    c = parse_criteria("HFE: PM3:Strong, PM2, PS3:Supporting, BP4")
    assert points(c) == 6 and derive(c, "tavtigian") == "PLP" and derive(c, "richards") == "PLP"
    assert [x for x in decisive_codes(c, "tavtigian", "PLP")] == ["PM3", "PM2", "PS3"]  # BP4 removal leaves LP


def _c(*items):
    out = []
    for code, strength in items:
        side = "B" if code[0] == "B" else "P"
        out.append(Criterion(code, side, strength, True))
    return out


@pytest.mark.parametrize("items,cls", [
    ([("PS1", "strong"), ("PS3", "strong")], "P"),
    ([("PVS1", "very_strong"), ("PM2", "moderate")], "LP"),
    ([("PM1", "moderate"), ("PM2", "moderate"), ("PM4", "moderate")], "LP"),
    ([("PVS1", "very_strong"), ("PS3", "strong")], "P"),
    ([("BS1", "strong"), ("BP4", "supporting")], "LB"),
    ([("BP4", "supporting"), ("BP6", "supporting")], "LB"),
    ([("BA1", "stand_alone")], "B"),
    ([("BS1", "strong"), ("BS2", "strong")], "B"),
    ([("PM2", "moderate"), ("PP3", "supporting")], "VUS"),
    ([("PS3", "strong"), ("PM2", "moderate"), ("BP4", "supporting"), ("BP6", "supporting")], "VUS"),  # conflict
])
def test_richards_2015_combining_rules(items, cls):
    assert richards_class(_c(*items)) == cls


def test_pp3_decisive_under_both_rules():
    c = parse_criteria("PVS1_Strong, PM2_Supporting, PP3_Moderate")
    assert points(c) == 7
    assert is_decisive(c, "PP3", "tavtigian", "PLP") and is_decisive(c, "PP3", "richards", "PLP")
    assert not is_decisive(c, "PP3", "tavtigian", "BLB")  # the stated class must agree with the list


def test_pp3_not_decisive_when_other_evidence_is_enough():
    c = parse_criteria("PS3, PM2, PP3")  # 4 + 2 + 1 = 7; without PP3 still 6
    assert derive(c, "tavtigian") == "PLP" and not is_decisive(c, "PP3", "tavtigian", "PLP")


def test_unqualified_pm2_downgrade_changes_decisiveness():
    c = parse_criteria("PS3, PM2, PP3")
    assert is_decisive(c, "PP3", "tavtigian", "PLP", pm2_supporting=True)  # 4 + 1 + 1 = 6; without PP3 it is 5
    assert is_decisive(c, "PP3", "richards", "PLP", pm2_supporting=True)
    explicit = parse_criteria("PS3, PM2_Moderate, PP3")
    assert not is_decisive(explicit, "PP3", "tavtigian", "PLP", pm2_supporting=True)  # an explicit strength is never downgraded


def test_bp4_decisive_only_when_it_tips_benign_to_uncertain():
    c = parse_criteria("BS1, BP4")  # -4 - 1 = -5: LB. Without BP4 -4: still LB under points, but VUS-free? -4 is LB
    assert not is_decisive(c, "BP4", "tavtigian", "BLB")
    c2 = parse_criteria("BP1, BP4")  # -2: LB; without BP4 -1: still LB
    assert not is_decisive(c2, "BP4", "tavtigian", "BLB")
    c3 = parse_criteria("BP4_Moderate")  # a single code: -2 LB; without it 0: VUS
    assert is_decisive(c3, "BP4", "tavtigian", "BLB")


def test_inconsistent_stated_class_is_never_decisive():
    c = parse_criteria("PM2, PP3")  # 3 points: VUS under the list
    assert not is_decisive(c, "PP3", "tavtigian", "PLP")


@pytest.mark.parametrize("sig,fam", [
    ("Likely pathogenic", "PLP"), ("Pathogenic/Likely pathogenic", "PLP"), ("Benign/Likely benign", "BLB"),
    ("Uncertain significance", "VUS"), ("Conflicting classifications of pathogenicity", "other"), (None, "other"),
])
def test_stated_family(sig, fam):
    assert stated_family(sig) == fam


def test_describe_summarises_one_submission():
    from rdv.criteria import describe

    d = describe("HFE: PM3:Strong, PM2, PS3:Supporting, BP4", "Likely pathogenic")
    assert d["stated"] == "PLP" and d["n_codes"] == 4 and d["points_t"] == 6
    assert d["derived_t"] == "PLP" and d["derived_r"] == "PLP"
    assert d["has_bp4"] and not d["has_pp3"] and not d["dec_bp4_t"] and d["dec_other_t"]
    assert d["codes"] == "PM3:strong,PM2:moderate*,PS3:supporting,BP4:supporting*"
    solo = describe("PP3 only", "Pathogenic")  # one code is kept, because BA1 alone is a complete benign justification
    assert solo["n_codes"] == 1 and solo["derived_r"] == "VUS" and not solo["dec_pp3_r"]
    ba1 = describe("BA1", "Benign")
    assert ba1["derived_r"] == "BLB" and ba1["derived_t"] == "BLB" and ba1["dec_other_r"] and not ba1["dec_bp4_r"]
    assert ba1["decisive_r"] == "BA1"
    assert describe("no evidence listed here", "Pathogenic") is None
    assert describe("PP3 not met", "Pathogenic") is None
    e = describe("PVS1_Strong, PM2_Supporting, PP3_Moderate", "Pathogenic")
    assert e["dec_pp3_t"] and e["dec_pp3_r"] and e["pp3_explicit"] and e["pp3_strength"] == "moderate"
