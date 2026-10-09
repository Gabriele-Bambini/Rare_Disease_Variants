# Downloading the data on an HPC cluster

One command fetches every dataset listed in `docs/brainstorm.md` §6 into `$DATA_DIR`.
It is resumable, idempotent and logs to `$DATA_DIR/logs/`.

## Bocconi HPC notes

What public sources say. The official pages could not be opened from where this was
written, so check what follows against your account email or with IT.

- **Hardware.** SLURM cluster of about 24 nodes, 752 cores, InfiniBand, several GPU types,
  about 250 TB of storage.
- **Access.** Aimed at faculty and PhD students. MSc DSBA/AI students get an account
  through the "Introduction to Linux for HPC" course. **If you don't have an account yet,
  ask the BSML board** (BSML ran an HPC tutorial in April 2025) or Bocconi's Technology
  Office.
- **Off campus.** SSH to the login host on `unibocconi.it` needs the Bocconi VPN.

**First step on the cluster:** `bash scripts/hpc/check_cluster.sh --compute`. It reports:
- the partitions and accounts to put in `download_all.sbatch`;
- where there is space for `DATA_DIR`;
- whether each data host is reachable from the login node *and* from a compute node.

## Setup (once)

```bash
git clone <this repo> && cd Rare_Disease_Variants
module load python/3.11            # or whatever your cluster provides
python3 -m venv venv && source venv/bin/activate
pip install -r scripts/hpc/requirements.txt
export DATA_DIR=/scratch/$USER/rare_disease_variants/data   # default: $SCRATCH or $HOME
```

## Run

```bash
# Option A: SLURM (if compute or transfer nodes have internet)
sbatch scripts/hpc/download_all.sbatch

# Option B: login / data-transfer node, inside tmux
tmux new -s dl
bash scripts/hpc/download_all.sh

# Only some stages
bash scripts/hpc/download_all.sh clinvar_current clinvar_archives
```

When the run ends, failed stages are listed together with the command that reruns only
those. Every downloaded file is recorded with its URL and UTC time in
`$DATA_DIR/MANIFEST.tsv`. This doubles as the data snapshot record for the paper.

## Stages

| Stage | What | Size (approx.) | Depends on |
|---|---|---|---|
| `clinvar_current` | `variant_summary`, `submission_summary` (MD5-checked), GRCh38 VCF | ~1–2 GB | |
| `clinvar_archives` | Quarterly `variant_summary` + `submission_summary` snapshots since 2017, for the label history (RQ1) | ~15–25 GB | |
| `mane` | Latest MANE summary + Ensembl genomic GTF | <100 MB | |
| `panelapp` | Green genes of all public PanelApp panels, with inheritance mode | <10 MB | |
| `uniprot` | Human Swiss-Prot: accession, gene, Ensembl xrefs, sequence | ~30 MB | |
| `alphamissense` | `AlphaMissense_hg38`, `aa_substitutions`, gene-level | ~1.9 GB | |
| `gnomad_constraint` | gnomAD v4.1 constraint metrics | ~0.1 GB | |
| `hpo` | `genes_to_phenotype.txt`, `phenotype.hpoa`, `hp.obo` | ~70 MB | |
| `proteingym` | DMS reference file + DMS substitutions (unzipped) | ~1–2 GB | |
| `conservation` | phyloP 241-mammal and 100-vertebrate bigWigs | ~15–25 GB | |
| `gnomad_slices` | gnomAD v4.1 **joint** AFs at MANE Select coding positions of the green genes, read by remote tabix (SNVs, per-ancestry AF, grpmax, FAF95, nhomalt) | a few GB | `panelapp`, `mane` |
| `alphafold` | AlphaFold DB models for those genes' UniProt accessions | ~2–4 GB | `panelapp`, `uniprot` |
| `esm` | ESM-2 650M weights, so GPU jobs can run offline | ~2.6 GB | |

Plan for roughly **50–60 GB**. The full gnomAD release (well over 1 TB) is deliberately
**not** downloaded. The slices contain only the positions the project uses. To slice more
genes later, edit `$DATA_DIR/panelapp/green_genes.tsv` (or point `make_cds_bed.py` at
another gene list), delete `gnomad/joint_v4.1_slices/`, and rerun `gnomad_slices`.

## What was tested

These parts were run end-to-end against the real servers:
- the `gnomad_slices` stage, on the gnomAD v4.1 joint VCFs on Google Cloud Storage,
  including skip-on-rerun;
- the `gnomad_constraint` and `hpo` stages.

The AlphaMissense and ProteinGym reference URLs were also verified to be reachable.

NCBI (ClinVar, MANE), PanelApp, UniProt, AlphaFold DB, UCSC and Hugging Face could not be
reached from the environment where this was written. For those, the parsing and crawling
logic was tested against mock servers, but the live URLs were not tested. If one of them
fails on the cluster, check the log first: the URLs most likely to drift are the phyloP
bigWigs and the ProteinGym zip, and both can be overridden in `config.env`.

**Not included yet:**
- REVEL, CADD and BayesDel scores. They will come from the myvariant.info batch API once the
  variant list exists, which avoids the dbNSFP download.
- Precomputed ESM1b scores (Brandes et al. 2023). Their hosting location still needs to be
  confirmed.
