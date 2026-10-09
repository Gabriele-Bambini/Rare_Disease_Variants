#!/usr/bin/env bash
# Run once on the cluster login node before downloading anything:
#   bash scripts/hpc/check_cluster.sh            # login node only
#   bash scripts/hpc/check_cluster.sh --compute  # also test a compute node via srun
# Prints what download_all.sh needs to know: partitions, storage, Python,
# and which data hosts are reachable from login and compute nodes.

HOSTS=(
  https://ftp.ncbi.nlm.nih.gov/pub/clinvar/
  https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv
  https://storage.googleapis.com/dm_alphamissense/AlphaMissense_gene_hg38.tsv.gz
  https://panelapp.genomicsengland.co.uk/api/v1/panels/
  https://rest.uniprot.org/uniprotkb/P38398.fasta
  https://alphafold.ebi.ac.uk/api/prediction/P38398
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/
  https://github.com/obophenotype/human-phenotype-ontology/releases
  https://marks.hms.harvard.edu/proteingym/
  https://huggingface.co/facebook/esm2_t33_650M_UR50D
)

probe_hosts() {
  for u in "${HOSTS[@]}"; do
    code="$(curl -sS -o /dev/null -m 20 -L -w '%{http_code}' "$u" 2>/dev/null)"
    printf '  %-4s %s\n' "${code:-ERR}" "$u"
  done
}

section() { printf '\n== %s\n' "$*"; }

section "Host"
echo "  $(hostname)  user=$USER"

section "Scheduler"
if command -v sinfo > /dev/null; then
  sinfo -o '  %-15P %-6a %-12l %-6D %-20G %N' 2>/dev/null | head -30
  echo "  accounts/QoS:"; sacctmgr -nP show assoc user="$USER" format=account,partition,qos 2>/dev/null | sed 's/^/    /' | head -10
else
  echo "  sinfo not found (not SLURM, or not on a login node)"
fi

section "Storage (candidate DATA_DIR locations)"
for d in "$HOME" "${SCRATCH:-}" "${WORK:-}" /scratch/"$USER" /work/"$USER" /data/"$USER" /mnt/scratch/"$USER"; do
  [[ -n "$d" && -d "$d" ]] && df -h "$d" 2>/dev/null | awk -v d="$d" 'NR==2{printf "  %-35s free %-8s of %s\n", d, $4, $2}'
done
command -v quota > /dev/null && { echo "  quota:"; quota -s 2>/dev/null | sed 's/^/    /' | head -8; }

section "Python and tools"
for t in python3 module conda mamba curl wget unzip md5sum tmux screen; do
  if command -v "$t" > /dev/null || type "$t" > /dev/null 2>&1; then echo "  ok   $t $(command -v "$t" 2>/dev/null)"; else echo "  --   $t"; fi
done
python3 --version 2>/dev/null | sed 's/^/  /'
type module > /dev/null 2>&1 && { echo "  python modules:"; module -t avail 2>&1 | grep -iE '^(python|anaconda|miniconda|miniforge)' | head -10 | sed 's/^/    /'; }

section "Internet from the login node"
probe_hosts

if [[ "${1:-}" == "--compute" ]] && command -v srun > /dev/null; then
  section "Internet from a compute node (srun, 2 min)"
  srun -n1 -t 2 bash -c "$(declare -f probe_hosts); HOSTS=(${HOSTS[*]}); probe_hosts" 2>&1 | tail -n +1
fi

section "Reading the result"
cat <<'EOF'
  200/301/302 = reachable. 000/ERR = blocked.
  - If compute nodes are blocked but the login node is not: run download_all.sh on the
    login node inside tmux, not with sbatch.
  - Set DATA_DIR to the location with the most free space (≥ 60 GB), then:
      export DATA_DIR=/path/with/space/rare_disease_variants/data
  - Copy the partition/account names above into download_all.sbatch if you use sbatch.
EOF
