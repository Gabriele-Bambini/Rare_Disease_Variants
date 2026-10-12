"""Readers and normalisers for ClinVar tab-delimited releases.

ClinVar files are large (the current submission_summary is ~3.7 GB uncompressed)
and their columns drift between releases, so everything here streams in chunks,
selects columns by name, and fills columns an older release lacks with nulls.
"""

from __future__ import annotations

import gzip
import io
import re
import zipfile
from pathlib import Path
from typing import Callable, Iterator

import polars as pl

# --- reading -----------------------------------------------------------------


def _open(path: str | Path):
    """Binary stream of a plain, gzipped or zipped (first member) file."""
    name = str(path)
    if name.endswith(".gz"):
        return gzip.open(path, "rb")
    if name.endswith(".zip"):
        zf = zipfile.ZipFile(path)
        return zf.open(zf.namelist()[0])
    return open(path, "rb")


def _is_header(line: bytes) -> bool:
    # variant_summary starts with the header; submission_summary has "##" prose
    # and "#Column: explanation" lines first. Only the real header contains tabs.
    return line.startswith(b"#") and not line.startswith(b"##") and b"\t" in line


def iter_tsv_chunks(
    path: str | Path, columns: list[str], chunk_lines: int = 400_000, separator: str = "\t",
    commented_header: bool = True,
) -> Iterator[pl.DataFrame]:
    """Yield DataFrames of `columns` (all Utf8) from a ClinVar-style table, chunk by chunk.

    Columns missing from this release come back as all-null. With `commented_header` (ClinVar,
    AlphaMissense) the header is the first '#' line containing a tab; otherwise (REVEL) it is the
    first line. `separator` is passed to the CSV parser.
    """
    with _open(path) as fh:
        for line in fh:
            if (_is_header(line) if commented_header else True):
                header = line.lstrip(b"#").rstrip(b"\r\n") if commented_header else line.rstrip(b"\r\n")
                break
        else:
            raise ValueError(f"no header line found in {path}")
        present_names = header.decode().split(separator)
        present = [c for c in columns if c in present_names]
        missing = [c for c in columns if c not in present_names]

        def parse(buf: list[bytes]) -> pl.DataFrame:
            df = pl.read_csv(
                io.BytesIO(header + b"\n" + b"".join(buf)),
                separator=separator,
                quote_char=None,
                infer_schema=False,
                columns=present,
                truncate_ragged_lines=True,
                null_values=None,
            )
            if missing:
                df = df.with_columns(pl.lit(None, dtype=pl.Utf8).alias(c) for c in missing)
            return df.select(columns)

        buf: list[bytes] = []
        for line in fh:
            if commented_header and line.startswith(b"#"):
                continue
            buf.append(line if line.endswith(b"\n") else line + b"\n")
            if len(buf) >= chunk_lines:
                yield parse(buf)
                buf = []
        if buf:
            yield parse(buf)


def read_tsv(
    path: str | Path,
    columns: list[str],
    chunk_filter: Callable[[pl.DataFrame], pl.DataFrame] | None = None,
    chunk_lines: int = 400_000,
    separator: str = "\t",
    commented_header: bool = True,
) -> pl.DataFrame:
    """Read selected columns, applying `chunk_filter` to each chunk to bound memory."""
    parts = []
    for df in iter_tsv_chunks(path, columns, chunk_lines, separator, commented_header):
        parts.append(chunk_filter(df) if chunk_filter else df)
    if not parts:
        return pl.DataFrame(schema={c: pl.Utf8 for c in columns})
    return pl.concat(parts, how="vertical_relaxed")


SNAPSHOT_RE = re.compile(r"(\d{4})-(\d{2})")


def snapshot_month(path: str | Path) -> str:
    """'variant_summary_2023-04.txt.gz' -> '2023-04'."""
    m = SNAPSHOT_RE.search(Path(path).name)
    if not m:
        raise ValueError(f"no YYYY-MM in file name {path}")
    return f"{m.group(1)}-{m.group(2)}"


# --- normalising labels --------------------------------------------------------

PLP, BLB, VUS, CONFLICT, OTHER = "PLP", "BLB", "VUS", "CONFLICT", "OTHER"
_PLP_TERMS = {"pathogenic", "likely pathogenic"}
_BLB_TERMS = {"benign", "likely benign"}


def classify_significance(value: str | None) -> str:
    """Collapse ClinVar's free-form ClinicalSignificance into PLP/BLB/VUS/CONFLICT/OTHER.

    Only the germline pathogenicity part counts: secondary terms after ';' or ','
    ("Pathogenic; risk factor", "Likely pathogenic, low penetrance") are dropped.
    Mixed P and B terms in one aggregate (rare, old releases) count as CONFLICT.
    """
    if value is None:
        return OTHER
    v = value.strip().lower()
    if "conflicting" in v:
        return CONFLICT
    head = re.split(r"[;,]", v, maxsplit=1)[0].strip()
    parts = {p.strip() for p in head.split("/") if p.strip()}
    if not parts:
        return OTHER
    if parts <= _PLP_TERMS:
        return PLP
    if parts <= _BLB_TERMS:
        return BLB
    if parts <= _PLP_TERMS | _BLB_TERMS:
        return CONFLICT
    if parts == {"uncertain significance"}:
        return VUS
    return OTHER


def review_stars(status: str | None) -> int:
    """ClinVar review status -> gold stars (0-4)."""
    if status is None:
        return 0
    s = status.strip().lower()
    if s.startswith("practice guideline"):
        return 4
    if s.startswith("reviewed by expert panel"):
        return 3
    if s.startswith("criteria provided, multiple submitters"):
        return 2
    if s.startswith("criteria provided"):  # single submitter, or conflicting
        return 1
    return 0


def map_unique(col: str, fn: Callable[[str | None], object], dtype) -> pl.Expr:
    """Apply a pure Python function via a lookup over unique values (fast: few uniques)."""
    return pl.col(col).fill_null("").map_batches(
        lambda s: s.replace_strict(
            {v: fn(v) for v in s.unique().to_list()}, return_dtype=dtype
        ),
        return_dtype=dtype,
    )


def parse_date(col: str) -> pl.Expr:
    """ClinVar dates look like 'Jun 29, 2010'; '-' and blanks become null."""
    return pl.col(col).str.strptime(pl.Date, "%b %d, %Y", strict=False)


# --- missense universe -----------------------------------------------------------

# "NM_000546.6(TP53):c.817C>T (p.Arg273Cys)"
_NAME_RE = r"^(?<transcript>N[MR]_\d+(?:\.\d+)?)\((?<gene>[^)]+)\):c\.[^ ]+ \(p\.(?<ref_aa>[A-Z][a-z]{2})(?<aa_pos>\d+)(?<alt_aa>[A-Z][a-z]{2})\)$"

UNIVERSE_COLUMNS = [
    "AlleleID", "Type", "Name", "GeneSymbol", "HGNC_ID", "Assembly", "Chromosome",
    "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF", "VariationID",
    "ClinicalSignificance", "ReviewStatus", "LastEvaluated", "NumberSubmitters",
]


def missense_universe(variant_summary: str | Path) -> pl.DataFrame:
    """All GRCh38 single-nucleotide missense variants in a (current) variant_summary.

    Missense = SNV whose HGVS protein change is a single amino-acid substitution.
    Excluded: synonymous, nonsense (Ter), stop-loss, and initiator Met1 changes
    (ACMG treats start-loss like loss of function, not via PP3/BP4).
    """

    def keep(df: pl.DataFrame) -> pl.DataFrame:
        return df.filter(
            (pl.col("Assembly") == "GRCh38") & (pl.col("Type") == "single nucleotide variant")
        )

    df = read_tsv(variant_summary, UNIVERSE_COLUMNS, keep)
    parsed = pl.col("Name").str.extract_groups(_NAME_RE)
    df = (
        df.with_columns(parsed.alias("p"))
        .unnest("p")
        .filter(
            pl.col("ref_aa").is_not_null()
            & (pl.col("ref_aa") != pl.col("alt_aa"))
            & (pl.col("ref_aa") != "Ter")
            & (pl.col("alt_aa") != "Ter")
            & ~((pl.col("ref_aa") == "Met") & (pl.col("aa_pos") == "1"))
        )
    )
    return (
        df.select(
            pl.col("VariationID").cast(pl.Int64),
            pl.col("AlleleID").cast(pl.Int64),
            "gene",  # from the HGVS name: one gene, unlike GeneSymbol ("A;B")
            pl.col("HGNC_ID").alias("hgnc_id"),
            "transcript",
            pl.col("Chromosome").alias("chrom"),
            pl.col("PositionVCF").cast(pl.Int64).alias("pos"),
            pl.col("ReferenceAlleleVCF").alias("ref"),
            pl.col("AlternateAlleleVCF").alias("alt"),
            "ref_aa",
            pl.col("aa_pos").cast(pl.Int64),
            "alt_aa",
            map_unique("ClinicalSignificance", classify_significance, pl.Utf8).alias("label"),
            map_unique("ReviewStatus", review_stars, pl.Int8).alias("stars"),
        )
        .unique("VariationID", keep="first")
        .sort("VariationID")
    )
