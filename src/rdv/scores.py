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
