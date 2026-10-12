import polars as pl

from rdv.scores import alphamissense

AM_HEADER = ["#CHROM", "POS", "REF", "ALT", "genome", "uniprot_id", "transcript_id", "protein_variant",
             "am_pathogenicity", "am_class"]


def _mane():
    rows = [("HGNC:1", "GENE1", "NM_1.1", "ENST00000000001.2", "ENST00000000001", "ENSP1", "MANE Select"),
            ("HGNC:2", "GENE2", "NM_2.1", "ENST00000000002.5", "ENST00000000002", "ENSP2", "MANE Select")]
    return pl.DataFrame(rows, schema=["hgnc_id", "mane_symbol", "refseq_nuc", "enst", "enst_base", "ensp", "mane_status"], orient="row")


def _universe():
    cols = ["AlleleID", "hgnc_id", "chrom", "pos", "ref", "alt", "ref_aa", "aa_pos", "alt_aa"]
    rows = [
        (1, "HGNC:1", "1", 100, "G", "A", "Arg", 10, "Cys"),  # scored on MANE and a non-MANE transcript
        (2, "HGNC:1", "1", 200, "C", "T", "Ala", 20, "Val"),  # only a non-MANE transcript
        (3, "HGNC:1", "1", 300, "C", "T", "Ala", 30, "Val"),  # not in AlphaMissense
        (4, "HGNC:1", "1", 400, "T", "C", "Leu", 40, "Pro"),  # only the MANE transcript of another gene
        (5, "HGNC:1", "1", 500, "A", "G", "Lys", 50, "Glu"),  # MANE row exists but protein change differs
    ]
    return pl.DataFrame(rows, schema=cols, orient="row")


def test_alphamissense_prefers_mane_of_the_right_gene(write_gz):
    rows = [
        ("chr1", 100, "G", "A", "hg38", "P1", "ENST00000999999.1", "R12C", 0.50, "ambiguous"),
        ("chr1", 100, "G", "A", "hg38", "P1", "ENST00000000001.2", "R10C", 0.90, "likely_pathogenic"),
        ("chr1", 200, "C", "T", "hg38", "P1", "ENST00000999999.1", "A20V", 0.10, "likely_benign"),
        ("chr1", 400, "T", "C", "hg38", "P2", "ENST00000000002.5", "L40P", 0.95, "likely_pathogenic"),
        ("chr1", 500, "A", "G", "hg38", "P1", "ENST00000000001.2", "K51E", 0.40, "ambiguous"),
        ("chr1", 600, "A", "G", "hg38", "P1", "ENST00000000001.2", "K1E", 0.40, "ambiguous"),  # not in universe
    ]
    path = write_gz("AlphaMissense_hg38.tsv.gz", AM_HEADER, rows, preamble=["# Copyright 2023", "#"])
    out = alphamissense(_universe(), path, _mane()).sort("allele_id")
    got = {r["allele_id"]: r for r in out.iter_rows(named=True)}
    assert sorted(got) == [1, 2, 4, 5]  # variant 3 and the extra AlphaMissense row do not appear

    assert (got[1]["am_score"], got[1]["protein_variant"], got[1]["am_mane"], got[1]["am_n_rows"]) == (0.90, "R10C", True, 2)
    assert got[1]["protein_agrees"] is True and got[1]["uniprot_id"] == "P1"
    assert (got[2]["am_mane"], got[2]["am_n_rows"]) == (False, 1)
    assert got[4]["am_mane"] is False  # MANE transcript of GENE2 does not count for a GENE1 variant
    assert (got[5]["am_mane"], got[5]["protein_agrees"]) == (True, False)  # numbering differs from ClinVar's
