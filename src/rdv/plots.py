"""Figures for the label-history and SCV-text analyses.

Style follows the dataviz method: one axis per panel (small multiples, never dual axes),
thin 1.5 pt lines, hairline recessive grid, categorical hues in fixed slot order, neutral
gray for controls/context, text in ink colours rather than series colours. Light surface only.
"""

from __future__ import annotations

import datetime as dt
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.ticker as mticker  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
import polars as pl  # noqa: E402

from .scv_text import MILESTONES  # noqa: E402

SURFACE, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
GRID, AXIS = "#e7e6e2", "#c9c8c3"
SLOT = {"blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a"}  # categorical slots 1-3
NEUTRALS = ["#52514e", "#8b8a85", "#b3b2ab"]  # dark -> light, for context shares
OTHER_LABS = "#d9d8d2"

TITLES = {
    "pp3_bp4": "Cites PP3 / BP4",
    "revel": "Names REVEL",
    "alphamissense": "Names AlphaMissense",
    "insilico_generic": "Says \"in silico\" (any tool)",
    "legacy_tools": "Names SIFT / PolyPhen / similar",
}
COLOURS = {"pp3_bp4": SLOT["blue"], "revel": SLOT["orange"], "alphamissense": SLOT["aqua"],
           "insilico_generic": INK2, "legacy_tools": INK2}


def _style(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=8, length=3, color=AXIS)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.xaxis.set_major_locator(mdates.YearLocator(2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))


def _milestones(ax):
    for day in MILESTONES:
        ax.axvline(dt.date.fromisoformat(day), color=INK2, linewidth=0.8, alpha=0.5, zorder=1)


MILESTONE_NOTE = "Vertical lines: " + "; ".join(f"{t} ({d[:7]})" for d, t in MILESTONES.items()) + "."


def fig_scv_mentions(
    agg: pl.DataFrame,
    standardized: dict[str, pl.DataFrame],
    comp: pl.DataFrame,
    out: str | Path,
    min_n: int = 1000,
    min_coverage: float = 0.6,
    ref_label: str = "2021-22",
) -> None:
    """2x3 grid: five term panels (pooled + lab-mix-fixed) and the lab composition."""
    fig, axes = plt.subplots(2, 3, figsize=(11.5, 6.6), dpi=200, facecolor=SURFACE)
    terms = list(TITLES)
    for ax, term in zip(axes.flat, terms):
        _style(ax)
        d = agg.filter((pl.col("term") == term) & (pl.col("n_text") >= min_n)).sort("quarter")
        ax.plot(d["quarter"], d["share"] * 100, color=COLOURS[term], linewidth=1.8, zorder=3)
        s = standardized[term].filter(pl.col("coverage") >= min_coverage)
        ax.plot(s["quarter"], s["share_std"] * 100, color=INK2, linewidth=1.2, alpha=0.9, zorder=2)
        _milestones(ax)
        title = TITLES[term] + ("  (note the scale)" if term == "alphamissense" else "")
        ax.set_title(title, loc="left", fontsize=10, color=INK, pad=6)
        ax.set_ylim(bottom=0)
        ax.yaxis.set_major_formatter(lambda v, _: f"{v:g}%")
    # composition panel
    ax = axes.flat[5]
    _style(ax)
    groups = comp.group_by("group").agg(pl.col("n").sum()).sort("n", descending=True)["group"].to_list()
    named = [g for g in groups if g != "All other labs"]
    wide = comp.pivot(on="group", index="quarter", values="share").sort("quarter").fill_null(0)
    wide = wide.join(comp.group_by("quarter").agg(pl.col("n").sum().alias("n")), on="quarter")
    wide = wide.filter(pl.col("n") >= min_n)
    order = named + ["All other labs"]
    shades = NEUTRALS[: len(named)] + [OTHER_LABS]
    ax.stackplot(wide["quarter"], *[wide[g] * 100 for g in order], colors=shades,
                 edgecolor=SURFACE, linewidth=1.2, zorder=3)
    _milestones(ax)
    ax.set_ylim(0, 100)
    ax.yaxis.set_major_formatter(lambda v, _: f"{v:g}%")
    ax.set_title("Who writes the text (share of SCVs)", loc="left", fontsize=10, color=INK, pad=6)
    ax.legend(
        handles=[Patch(facecolor=c, edgecolor=SURFACE, label="Other labs" if g == "All other labs" else g.split()[0])
                 for g, c in zip(order, shades)][::-1],
        loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=len(order), frameon=False, fontsize=7.5,
        labelcolor=INK2, handlelength=1.2, columnspacing=1.2,
    )
    for ax in axes.flat:
        ax.set_xlim(dt.date(2016, 1, 1), dt.date(2026, 10, 1))
    handles = [Line2D([], [], color=INK, linewidth=1.8, label="share of text-bearing SCVs, all labs pooled"),
               Line2D([], [], color=INK2, linewidth=1.2, alpha=0.9,
                      label=f"same labs held at their {ref_label} volume weights")]
    fig.legend(handles=handles, loc="upper left", ncol=2, frameon=False, fontsize=8,
               labelcolor=INK2, bbox_to_anchor=(0.01, 1.0))
    fig.text(0.01, 0.005,
             "Text-bearing SCVs on ClinVar missense variants from clinical-testing labs with assertion criteria; "
             "quarter = DateLastEvaluated.\nLab-mix line drawn only for quarters where those labs cover >=60% of "
             "the reference weight. " + MILESTONE_NOTE,
             fontsize=7, color=INK2, va="bottom")
    fig.tight_layout(rect=(0, 0.05, 1, 0.95), h_pad=2.2, w_pad=2.0)
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)


def fig_label_flow(first_confident: pl.DataFrame, out: str | Path) -> None:
    """New confident labels per snapshot quarter, by how the variant got there."""
    routes = [("new", "enters ClinVar already labelled", SLOT["blue"]),
              ("from_vus", "resolved from VUS", SLOT["orange"]),
              ("from_conflict", "resolved from a conflict", SLOT["aqua"])]
    d = first_confident.filter(~pl.col("left_censored")).with_columns(
        pl.col("snapshot").str.to_date("%Y-%m").alias("quarter")
    )
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 3.9), dpi=200, facecolor=SURFACE)
    for ax, (lab, title) in zip(axes, [("PLP", "Pathogenic / likely pathogenic"),
                                       ("BLB", "Benign / likely benign")]):
        _style(ax)
        sub = d.filter(pl.col("label") == lab)
        quarters = sorted(sub["quarter"].unique().to_list())
        bottom = [0] * len(quarters)
        for route, name, colour in routes:
            counts = {r["quarter"]: r["n"] for r in
                      sub.filter(pl.col("route") == route).group_by("quarter").agg(pl.len().alias("n")).iter_rows(named=True)}
            vals = [counts.get(q, 0) for q in quarters]
            ax.bar(quarters, vals, bottom=bottom, width=70, color=colour, edgecolor=SURFACE,
                   linewidth=1.0, label=name, zorder=3)
            bottom = [b + v for b, v in zip(bottom, vals)]
        _milestones(ax)
        for i in sorted(range(len(quarters)), key=lambda i: -bottom[i])[:3]:
            ax.annotate(f"{bottom[i]:,}", (quarters[i], bottom[i]), xytext=(0, 3), textcoords="offset points",
                        ha="center", fontsize=7, color=INK, zorder=5,
                        bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.6))
        ax.set_title(title, loc="left", fontsize=10, color=INK, pad=6)
        ax.yaxis.set_major_formatter(lambda v, _: f"{v:,.0f}")
        ax.set_xlim(dt.date(2017, 1, 1), dt.date(2026, 12, 1))
    axes[0].legend(frameon=False, fontsize=8, labelcolor=INK2, loc="upper left")
    fig.text(0.01, 0.005,
             "Variants first seen with a confident label in each quarterly ClinVar snapshot (2017-04 onward; "
             "labels already present in 2017-01 are excluded). Missense variants only.\n" + MILESTONE_NOTE,
             fontsize=7, color=INK2, va="bottom")
    fig.tight_layout(rect=(0, 0.07, 1, 1), w_pad=2.5)
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)


PRED_STYLE = {  # fixed identity per predictor across every figure
    "revel": ("REVEL", SLOT["orange"]),
    "am": ("AlphaMissense", SLOT["aqua"]),
    "esm1b": ("ESM1b (control)", SLOT["blue"]),
    "grantham": ("Grantham (control)", INK2),
}


def _bin_date(label: str) -> dt.date:
    """'2023H1' -> 2023-01-01 (entries at the Jan and Apr snapshots), 'H2' -> 2023-07-01."""
    return dt.date(int(label[:4]), 1 if label.endswith("H1") else 7, 1)


def fig_rq1(event_study: pl.DataFrame, primary: pl.DataFrame, posthoc: pl.DataFrame, out: str | Path) -> None:
    """Three panels: within-gene AUC by half-year; predictor minus ESM1b; forest of every primary variant."""
    fig = plt.figure(figsize=(13, 5.2), dpi=200, facecolor=SURFACE)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1.1, 1.0], wspace=0.55)
    ax1, ax2, ax3 = (fig.add_subplot(gs[0, i]) for i in range(3))
    for ax in (ax1, ax2, ax3):
        _style(ax)
    for ax in (ax1, ax2):
        _milestones(ax)

    auc = event_study.filter(pl.col("kind") == "auc")
    for pred, (name, colour) in PRED_STYLE.items():
        d = auc.filter(pl.col("pred") == pred).sort("bin")
        if d.is_empty():
            continue
        ax1.plot([_bin_date(b) for b in d["bin"]], d["value"], color=colour, linewidth=1.6, label=name, zorder=3)
    ax1.set_ylim(0.55, 1.0)
    ax1.set_title("Within-gene AUC of new labels", loc="left", fontsize=10, color=INK, pad=6)
    ax1.legend(frameon=False, fontsize=7.5, labelcolor=INK2, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2)

    diff = event_study.filter(pl.col("kind") == "diff")
    for key, pred in (("revel-esm1b", "revel"), ("am-esm1b", "am")):
        d = diff.filter(pl.col("pred") == key).sort("bin")
        x = [_bin_date(b) for b in d["bin"]]
        colour = PRED_STYLE[pred][1]
        ax2.fill_between(x, d["lo"], d["hi"], color=colour, alpha=0.14, linewidth=0, zorder=2)
        ax2.plot(x, d["value"], color=colour, linewidth=1.6, label=f"{PRED_STYLE[pred][0]} minus ESM1b", zorder=3)
    ax2.axhline(0, color=INK2, linewidth=0.8, zorder=2)
    ax2.set_ylim(-0.15, 0.30)
    ax2.set_title("Advantage over the control, with 95% CI", loc="left", fontsize=10, color=INK, pad=6)
    ax2.legend(frameon=False, fontsize=7.5, labelcolor=INK2, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=1)
    for ax in (ax1, ax2):
        ax.set_xlim(dt.date(2019, 1, 1), dt.date(2026, 10, 1))

    rows = []
    names = {"primary": "primary (cap 500)", "cap_250": "cap 250", "cap_1000": "cap 1000", "cap_none": "no cap",
             "exclude_top2_labs": "without top 2 labs", "posthoc_clean_boundary": "clean boundary*",
             "placebo_-24m": "placebo, 24 mo earlier", "placebo_-36m": "placebo, 36 mo earlier"}
    both = pl.concat([primary.select(posthoc.columns), posthoc], how="vertical_relaxed") if set(posthoc.columns) <= set(primary.columns) else primary
    for test, colour, title in (("P1", SLOT["orange"], "P1  REVEL vs ESM1b, Dec 2022"),
                                ("P2", SLOT["aqua"], "P2  AlphaMissense vs ESM1b, Sep 2023")):
        rows.append((title, None, None, None, None))
        for key, label in names.items():
            r = both.filter(pl.col("test").str.starts_with(test) & (pl.col("analysis") == key))
            if r.is_empty():
                continue
            r = r.row(0, named=True)
            rows.append((label, r["estimate"], r["lo"], r["hi"], INK2 if key.startswith("placebo") else colour))
    y = len(rows)
    for label, est, lo, hi, colour in rows:
        y -= 1
        if est is None:
            ax3.text(0.01, y, label, fontsize=8.5, color=INK, weight="bold", va="center", ha="left",
                     transform=ax3.get_yaxis_transform())
            continue
        ax3.plot([lo, hi], [y, y], color=colour, linewidth=1.6, zorder=3)
        ax3.plot([est], [y], marker="o", markersize=5, color=colour, markeredgecolor=SURFACE, markeredgewidth=1.2, zorder=4)
        ax3.text(-0.02, y, label, fontsize=7.5, color=INK2, va="center", ha="right", transform=ax3.get_yaxis_transform())
    ax3.axvline(0, color=INK2, linewidth=0.8, zorder=1)
    ax3.set_yticks([])
    ax3.set_ylim(-0.8, len(rows) + 0.2)
    ax3.set_xlim(-0.08, 0.08)
    ax3.xaxis.set_major_locator(mticker.MultipleLocator(0.04))  # replaces the date locator set by _style
    ax3.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: "0" if abs(v) < 1e-9 else f"{v:+.2f}"))
    ax3.spines["left"].set_visible(False)
    ax3.grid(axis="y", visible=False)
    ax3.grid(axis="x", color=GRID, linewidth=0.8)
    ax3.set_title("Difference-in-differences, 95% CI", loc="left", fontsize=10, color=INK, pad=6)
    ax3.set_xlabel("change in (predictor minus control) AUC", fontsize=8, color=INK2)
    fig.text(0.01, 0.005,
             "Within-gene weighted AUC of variants as they enter ClinVar as P/LP or B/LB (lab-quarter weights, cap 500), half-year bins of "
             "entry snapshot; gene-cluster bootstrap.\n*Post window starting one snapshot later (post hoc). " + MILESTONE_NOTE,
             fontsize=7, color=INK2, va="bottom")
    fig.subplots_adjust(left=0.05, right=0.985, top=0.9, bottom=0.2)
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)
