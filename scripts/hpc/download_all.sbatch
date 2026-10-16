#!/bin/bash
#SBATCH --job-name=rdv-download
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=rdv-download-%j.log
# Adapt to your cluster, e.g.:
##SBATCH --partition=transfer
##SBATCH --account=your_account
#
# Submit from the repository root:
#   sbatch scripts/hpc/download_all.sbatch                  # all stages
#   sbatch scripts/hpc/download_all.sbatch gnomad_slices    # selected stages
#
# Many clusters give compute nodes no internet access. If so, use a
# data-transfer partition, or run download_all.sh on a login node inside tmux.

set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

# module load python/3.11   # if your cluster uses environment modules
# source venv/bin/activate  # environment with scripts/hpc/requirements.txt

export WORKERS="${SLURM_CPUS_PER_TASK:-8}"
bash scripts/hpc/download_all.sh "$@"
