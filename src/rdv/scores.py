"""Attach precomputed predictor scores to the ClinVar missense universe.

Identifier hygiene (docs/brainstorm.md section 6): variants are matched on GRCh38 genomic
coordinates (chrom, pos, ref, alt), never on gene names or ClinVar's transcript numbering. The
protein-level coordinates used by protein-keyed resources (ESM1b, structures) come from
AlphaMissense's own genomic-to-protein mapping, restricted to the MANE transcript of the variant's
gene whenever AlphaMissense has it. Every join reports what it could not match.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from .aa import THREE_TO_ONE
from .clinvar import read_tsv

AM_COLUMNS = ["CHROM", "POS", "REF", "ALT", "uniprot_id", "transcript_id", "protein_variant",
              "am_pathogenicity", "am_class"]
_STATUS_RANK = {"MANE Select": 0, "MANE Plus Clinical": 1}


def mane_table(summary: str | Path) -> pl.DataFrame:
    """MANE v1.x summary: one row per MANE transcript, with versionless Ensembl id."""
    df = read_tsv(summary, ["HGNC_ID", "symbol", "RefSeq_nuc", "Ensembl_nuc", "Ensembl_prot", "MANE_status"])
    return df.select(
        pl.col("HGNC_ID").alias("hgnc_id"),
        pl.col("symbol").alias("mane_symbol"),
        pl.col("RefSeq_nuc").alias("refseq_nuc"),
        pl.col("Ensembl_nuc").alias("enst"),
        pl.col("Ensembl_nuc").str.split(".").list.first().alias("enst_base"),
        pl.col("Ensembl_prot").alias("ensp"),
        pl.col("MANE_status").alias("mane_status"),
    )


def alphamissense(universe: pl.DataFrame, path: str | Path, mane: pl.DataFrame) -> pl.DataFrame:
    """One row per universe variant that AlphaMissense scores.

    Columns: allele_id, am_score, am_class, uniprot_id, transcript, protein_variant (one-letter,
    e.g. 'R273C'), am_n_rows (transcripts AlphaMissense lists for this change), am_mane (the chosen
    row is on the MANE transcript of the variant's gene), protein_agrees (AlphaMissense's protein
    change equals the one in ClinVar's HGVS name).
    """
    keys = universe.select(
        pl.col("AlleleID").alias("allele_id"), "hgnc_id", "chrom", "pos", "ref", "alt",
        pl.col("ref_aa"), pl.col("aa_pos"), pl.col("alt_aa"),
    )
    join_keys = keys.select("allele_id", "chrom", "pos", "ref", "alt")

    def keep(df: pl.DataFrame) -> pl.DataFrame:
        df = df.with_columns(
            pl.col("CHROM").str.strip_prefix("chr").alias("chrom"),
            pl.col("POS").cast(pl.Int64).alias("pos"),
        )
        return df.join(join_keys, left_on=["chrom", "pos", "REF", "ALT"],
                       right_on=["chrom", "pos", "ref", "alt"], how="inner").select(
            "allele_id", "uniprot_id", "transcript_id", "protein_variant", "am_pathogenicity", "am_class"
        )

    rows = read_tsv(path, AM_COLUMNS, keep, chunk_lines=2_000_000)
    n_rows = rows.group_by("allele_id").agg(pl.len().alias("am_n_rows"))
    gene_of = mane.select("enst_base", pl.col("hgnc_id").alias("mane_hgnc"), "mane_status")
    ranked = (
        rows.with_columns(pl.col("transcript_id").str.split(".").list.first().alias("enst_base"))
        .join(gene_of, on="enst_base", how="left")
        .join(keys.select("allele_id", "hgnc_id", "ref_aa", "aa_pos", "alt_aa"), on="allele_id")
        .with_columns(
            (pl.col("mane_hgnc") == pl.col("hgnc_id")).fill_null(False).alias("am_mane"),
            pl.col("mane_status").replace_strict(_STATUS_RANK, default=2, return_dtype=pl.Int8).alias("rank"),
            (pl.col("ref_aa").replace_strict(THREE_TO_ONE, default="?") + pl.col("aa_pos").cast(pl.Utf8)
             + pl.col("alt_aa").replace_strict(THREE_TO_ONE, default="?")).alias("clinvar_pv"),
        )
        .with_columns((pl.col("protein_variant") == pl.col("clinvar_pv")).alias("protein_agrees"))
        # prefer the MANE transcript of the right gene, then ClinVar's own protein change
        .sort("allele_id", pl.col("am_mane"), pl.col("protein_agrees"), pl.col("rank"),
              descending=[False, True, True, False])
        .group_by("allele_id", maintain_order=True)
        .first()
    )
    return (
        ranked.join(n_rows, on="allele_id")
        .select(
            "allele_id",
            pl.col("am_pathogenicity").cast(pl.Float64).alias("am_score"),
            pl.col("am_class"),
            "uniprot_id",
            pl.col("transcript_id").alias("transcript"),
            "protein_variant",
            "am_n_rows",
            "am_mane",
            "protein_agrees",
        )
        .sort("allele_id")
    )


REVEL_COLUMNS = ["chr", "grch38_pos", "ref", "alt", "aaref", "aaalt", "REVEL", "Ensembl_transcriptid"]


def revel(universe: pl.DataFrame, path: str | Path, am: pl.DataFrame | None = None) -> pl.DataFrame:
    """REVEL v1.3 score for each universe variant it covers.

    A genomic change can have several REVEL rows when overlapping transcripts translate it
    differently. The row is chosen deterministically: first the one whose amino-acid change equals
    ClinVar's, then the one listing the transcript AlphaMissense used (`am.transcript`), then the
    highest score. Columns: allele_id, revel, revel_n_rows, revel_aa_agrees.
    """
    keys = universe.select(
        pl.col("AlleleID").alias("allele_id"), "chrom", "pos", "ref", "alt",
        pl.col("ref_aa").replace_strict(THREE_TO_ONE, default="?").alias("ref1"),
        pl.col("alt_aa").replace_strict(THREE_TO_ONE, default="?").alias("alt1"),
    )

    def keep(df: pl.DataFrame) -> pl.DataFrame:
        df = df.with_columns(pl.col("grch38_pos").cast(pl.Int64, strict=False).alias("pos"))
        return df.join(keys, left_on=["chr", "pos", "ref", "alt"], right_on=["chrom", "pos", "ref", "alt"],
                       how="inner").select("allele_id", "aaref", "aaalt", "REVEL", "Ensembl_transcriptid",
                                           "ref1", "alt1")

    rows = read_tsv(path, REVEL_COLUMNS, keep, chunk_lines=2_000_000, separator=",", commented_header=False)
    rows = rows.with_columns(pl.col("REVEL").cast(pl.Float64, strict=False)).filter(pl.col("REVEL").is_not_null())
    n_rows = rows.group_by("allele_id").agg(pl.len().alias("revel_n_rows"))
    rows = rows.with_columns(((pl.col("aaref") == pl.col("ref1")) & (pl.col("aaalt") == pl.col("alt1"))).alias("revel_aa_agrees"))
    if am is not None:
        rows = rows.join(
            am.select("allele_id", pl.col("transcript").str.split(".").list.first().alias("am_enst")),
            on="allele_id", how="left",
        ).with_columns(
            pl.col("Ensembl_transcriptid").str.contains(pl.col("am_enst"), literal=True).fill_null(False).alias("on_am")
        )
    else:
        rows = rows.with_columns(pl.lit(False).alias("on_am"))
    best = (
        rows.sort("allele_id", "revel_aa_agrees", "on_am", "REVEL", descending=[False, True, True, True])
        .group_by("allele_id", maintain_order=True)
        .first()
    )
    return (
        best.join(n_rows, on="allele_id")
        .select("allele_id", pl.col("REVEL").alias("revel"), "revel_n_rows", "revel_aa_agrees")
        .sort("allele_id")
    )


def esm1b(am: pl.DataFrame, zip_path: str | Path) -> pl.DataFrame:
    """ESM1b log-likelihood ratio (Brandes et al. 2023) for variants with an AlphaMissense mapping.

    The zip holds one 20 x L matrix per UniProt isoform: rows are the substituted amino acid, columns
    are 'WT position' (e.g. 'R 273'). Variants are looked up by the UniProt accession and one-letter
    protein change that AlphaMissense reports. The column's wild-type residue must equal the
    variant's, otherwise `esm1b_llr` stays null (a numbering mismatch, never silently scored).
    Columns: allele_id, esm1b_llr (negative = more damaging), esm1b_ref_ok.
    """
    import zipfile

    import numpy as np

    todo = am.filter(pl.col("uniprot_id").is_not_null()).select(
        "allele_id", "uniprot_id",
        pl.col("protein_variant").str.slice(0, 1).alias("wt"),
        pl.col("protein_variant").str.slice(1).str.head(-1).alias("pos"),
        pl.col("protein_variant").str.tail(1).alias("alt"),
    )
    out_ids, out_llr, out_ok = [], [], []
    with zipfile.ZipFile(zip_path) as zf:
        members = {n.rsplit("/", 1)[-1]: n for n in zf.namelist() if n.endswith("_LLR.csv")}
        for (uid,), g in todo.group_by("uniprot_id"):
            member = members.get(f"{uid}_LLR.csv")
            if member is None:
                continue
            # Hand parser: these matrices are 20 rows by up to 35,000 columns, where a general CSV
            # reader spends its time on column bookkeeping.
            lines = [ln for ln in zf.read(member).decode().split("\n") if ln]
            col_of = {c: j for j, c in enumerate(lines[0].split(",")[1:])}
            alt_row, rows = {}, []
            for k, ln in enumerate(lines[1:]):
                head, _, rest = ln.partition(",")
                alt_row[head] = k
                rows.append(np.array(rest.split(","), dtype=float))
            values = np.vstack(rows)
            for aid, wt, pos, alt in zip(g["allele_id"], g["wt"], g["pos"], g["alt"]):
                j = col_of.get(f"{wt} {pos}")
                i = alt_row.get(alt)
                out_ids.append(aid)
                out_ok.append(j is not None)
                out_llr.append(float(values[i, j]) if (j is not None and i is not None) else None)
    return pl.DataFrame(
        {"allele_id": out_ids, "esm1b_llr": out_llr, "esm1b_ref_ok": out_ok},
        schema={"allele_id": pl.Int64, "esm1b_llr": pl.Float64, "esm1b_ref_ok": pl.Boolean},
    ).sort("allele_id")
