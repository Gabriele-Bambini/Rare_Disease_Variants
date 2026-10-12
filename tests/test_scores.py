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


def test_revel_picks_matching_amino_acid_then_am_transcript_then_highest(tmp_path):
    import zipfile

    from rdv.scores import revel

    header = "chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,Ensembl_transcriptid"
    rows = [
        "1,90,100,G,A,R,C,0.30,ENST00000000011",  # matches ClinVar (Arg>Cys) but a different transcript
        "1,90,100,G,A,R,C,0.40,ENST00000000001;ENST00000000012",  # matches, and lists AlphaMissense's transcript
        "1,90,100,G,A,R,W,0.99,ENST00000000013",  # different amino acid: never preferred
        "1,190,200,C,T,A,V,0.10,ENST00000000014",
        "1,190,200,C,T,A,V,0.60,ENST00000000015",  # same aa, no AM info: highest wins
        "1,290,.,C,T,A,V,0.50,ENST00000000016",  # no GRCh38 position
        "1,390,400,T,C,L,P,.,ENST00000000017",  # no score
        "2,490,500,A,G,K,E,0.70,ENST00000000018",  # variant not in the universe
    ]
    path = tmp_path / "revel.zip"
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("revel_with_transcript_ids", "\n".join([header, *rows]) + "\n")
    am = pl.DataFrame({"allele_id": [1, 2], "transcript": ["ENST00000000001.2", "ENST00000999999.1"]})
    out = revel(_universe().filter(pl.col("AlleleID").is_in([1, 2, 3, 4])).with_columns(
        pl.when(pl.col("AlleleID") == 3).then(pl.lit(400)).otherwise(pl.col("pos")).alias("pos"),
        pl.when(pl.col("AlleleID") == 3).then(pl.lit("T")).otherwise(pl.col("ref")).alias("ref"),
        pl.when(pl.col("AlleleID") == 3).then(pl.lit("C")).otherwise(pl.col("alt")).alias("alt"),
    ), path, am).sort("allele_id")
    got = {r["allele_id"]: r for r in out.iter_rows(named=True)}
    assert sorted(got) == [1, 2]  # variant 3 has no score, 4 is absent from REVEL
    assert (got[1]["revel"], got[1]["revel_n_rows"], got[1]["revel_aa_agrees"]) == (0.40, 3, True)
    assert (got[2]["revel"], got[2]["revel_n_rows"]) == (0.60, 2)


def test_esm1b_lookup_checks_the_wild_type_residue(tmp_path):
    import zipfile

    from rdv.scores import esm1b

    prot = ",M 1,R 2,A 3\nR,-1.0,0.000,-2.5\nC,-3.0,-7.5,-4.0\nA,-0.5,-6.0,0.000\n"
    path = tmp_path / "esm.zip"
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("content/ALL_hum_isoforms_ESM1b_LLR/", "")
        zf.writestr("content/ALL_hum_isoforms_ESM1b_LLR/P1_LLR.csv", prot)
    am = pl.DataFrame({
        "allele_id": [1, 2, 3, 4, 5],
        "uniprot_id": ["P1", "P1", "P1", "P1", "P9"],
        "protein_variant": ["R2C", "A3R", "R3C", "M1A", "R2C"],  # 3rd: residue 3 is A, not R; 5th: protein not in zip
    })
    out = esm1b(am, path).sort("allele_id")
    got = {r["allele_id"]: (r["esm1b_llr"], r["esm1b_ref_ok"]) for r in out.iter_rows(named=True)}
    assert got == {1: (-7.5, True), 2: (-2.5, True), 3: (None, False), 4: (-0.5, True)}
