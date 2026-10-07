#!/usr/bin/env bash
# Download every dataset in docs/brainstorm.md §6 into $DATA_DIR.
#
#   bash scripts/hpc/download_all.sh                 # all stages, in order
#   bash scripts/hpc/download_all.sh clinvar_current alphamissense   # selected stages
#
# Safe to rerun: finished files carry a .done marker and are skipped, partial
# downloads resume. A failing stage is reported and the others still run.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.env
source "$HERE/config.env"
PY="$HERE/py"

ALL_STAGES=(clinvar_current clinvar_archives mane panelapp uniprot alphamissense
            gnomad_constraint hpo proteingym conservation gnomad_slices alphafold esm)
STAGES=("$@")
[[ ${#STAGES[@]} -eq 0 ]] && STAGES=("${ALL_STAGES[@]}")

mkdir -p "$DATA_DIR/logs"
MANIFEST="$DATA_DIR/MANIFEST.tsv"
[[ -f "$MANIFEST" ]] || printf "downloaded_utc\tpath\turl\n" > "$MANIFEST"

log() { printf '[%s] %s\n' "$(date -u +%H:%M:%S)" "$*"; }

# fetch URL DIR [NAME]: resumable download, skipped once a .done marker exists.
fetch() {
  local url="$1" dir="$2" name="${3:-$(basename "${1%%\?*}")}"
  local dest="$dir/$name"
  mkdir -p "$dir"
  if [[ -f "$dest.done" ]]; then log "  have $name"; return 0; fi
  log "  get  $name"
  curl -fL --retry 5 --retry-delay 15 --retry-connrefused -C - -sS -o "$dest" "$url" || return 1
  touch "$dest.done"
  printf "%s\t%s\t%s\n" "$(date -u +%FT%TZ)" "${dest#"$DATA_DIR"/}" "$url" >> "$MANIFEST"
}

# check_md5 FILE MD5FILE: compare against NCBI's published checksum.
check_md5() {
  local expected actual
  expected="$(awk '{print $1; exit}' "$2")"
  actual="$(md5sum "$1" | awk '{print $1}')"
  if [[ "$expected" != "$actual" ]]; then
    log "  MD5 MISMATCH for $1; deleting so the next run re-downloads"
    rm -f "$1" "$1.done"
    return 1
  fi
}

CLINVAR=https://ftp.ncbi.nlm.nih.gov/pub/clinvar

stage_clinvar_current() {
  local d="$DATA_DIR/clinvar/current"
  for f in variant_summary submission_summary; do
    fetch "$CLINVAR/tab_delimited/$f.txt.gz" "$d" || return 1
    fetch "$CLINVAR/tab_delimited/$f.txt.gz.md5" "$d" || return 1
    check_md5 "$d/$f.txt.gz" "$d/$f.txt.gz.md5" || return 1
  done
  for f in clinvar.vcf.gz clinvar.vcf.gz.tbi clinvar.vcf.gz.md5; do
    fetch "$CLINVAR/vcf_GRCh38/$f" "$d" || return 1
  done
  check_md5 "$d/clinvar.vcf.gz" "$d/clinvar.vcf.gz.md5"
}

stage_clinvar_archives() {
  local urls
  urls="$("$PYTHON" "$PY/clinvar_archives.py" --from-year "$CLINVAR_ARCHIVE_FROM" \
            --months "$CLINVAR_ARCHIVE_MONTHS" --kinds "$CLINVAR_ARCHIVE_KINDS")" || return 1
  log "  $(wc -l <<< "$urls") archive files"
  local rc=0
  while read -r url; do
    fetch "$url" "$DATA_DIR/clinvar/archive" || rc=1
  done <<< "$urls"
  return $rc
}

stage_mane() {
  local base=https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/ url
  for pat in 'MANE\.GRCh38\.v[0-9.]+\.summary\.txt\.gz' 'MANE\.GRCh38\.v[0-9.]+\.ensembl_genomic\.gtf\.gz'; do
    url="$("$PYTHON" "$PY/listing.py" "$base" "$pat")" || return 1
    fetch "$url" "$DATA_DIR/mane" || return 1
  done
}

stage_panelapp() {
  local out="$DATA_DIR/panelapp/green_genes.tsv"
  mkdir -p "$DATA_DIR/panelapp"
  if [[ -f "$out" ]]; then log "  have green_genes.tsv"; return 0; fi
  "$PYTHON" "$PY/panelapp.py" --out "$out"
}

stage_uniprot() {
  # The stream endpoint cannot resume, so drop any partial file first.
  [[ -f "$DATA_DIR/uniprot/human_reviewed.tsv.gz.done" ]] || rm -f "$DATA_DIR/uniprot/human_reviewed.tsv.gz"
  local q='(organism_id:9606)%20AND%20(reviewed:true)'
  fetch "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=tsv&fields=accession,gene_primary,gene_names,length,xref_ensembl,sequence&query=$q" \
        "$DATA_DIR/uniprot" human_reviewed.tsv.gz
}

stage_alphamissense() {
  local base=https://storage.googleapis.com/dm_alphamissense
  for f in AlphaMissense_hg38.tsv.gz AlphaMissense_aa_substitutions.tsv.gz AlphaMissense_gene_hg38.tsv.gz; do
    fetch "$base/$f" "$DATA_DIR/alphamissense" || return 1
  done
}

stage_gnomad_constraint() {
  fetch https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv \
        "$DATA_DIR/gnomad"
}

stage_hpo() {
  local base=https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download
  for f in genes_to_phenotype.txt phenotype.hpoa hp.obo; do
    fetch "$base/$f" "$DATA_DIR/hpo" || return 1
  done
}

stage_proteingym() {
  local d="$DATA_DIR/proteingym"
  fetch https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/main/reference_files/DMS_substitutions.csv "$d" || return 1
  fetch "$PROTEINGYM_DMS_URL" "$d" || return 1
  local zip="$d/$(basename "$PROTEINGYM_DMS_URL")"
  if [[ ! -f "$d/unzipped.done" ]]; then
    unzip -q -o "$zip" -d "$d/dms" && touch "$d/unzipped.done"
  fi
}

stage_conservation() {
  fetch "$PHYLOP241_URL" "$DATA_DIR/conservation" || return 1
  fetch "$PHYLOP100_URL" "$DATA_DIR/conservation"
}

stage_gnomad_slices() {
  local genes="$DATA_DIR/panelapp/green_genes.tsv" bed="$DATA_DIR/gnomad/mane_cds_green_genes.bed"
  local summary gtf
  summary="$(ls "$DATA_DIR"/mane/MANE.GRCh38.*.summary.txt.gz 2>/dev/null | tail -1)"
  gtf="$(ls "$DATA_DIR"/mane/MANE.GRCh38.*.ensembl_genomic.gtf.gz 2>/dev/null | tail -1)"
  if [[ ! -f "$genes" || -z "$summary" || -z "$gtf" ]]; then
    log "  needs the panelapp and mane stages first"; return 1
  fi
  mkdir -p "$DATA_DIR/gnomad"
  "$PYTHON" "$PY/make_cds_bed.py" --genes "$genes" --mane-summary "$summary" \
            --mane-gtf "$gtf" --out "$bed" || return 1
  "$PYTHON" "$PY/gnomad_slices.py" --bed "$bed" --out-dir "$DATA_DIR/gnomad/joint_v4.1_slices" \
            --template "$GNOMAD_VCF_TEMPLATE" --workers "$WORKERS"
}

stage_alphafold() {
  local genes="$DATA_DIR/panelapp/green_genes.tsv" uniprot="$DATA_DIR/uniprot/human_reviewed.tsv.gz"
  if [[ ! -f "$genes" || ! -f "$uniprot" ]]; then
    log "  needs the panelapp and uniprot stages first"; return 1
  fi
  "$PYTHON" "$PY/alphafold.py" --genes "$genes" --uniprot "$uniprot" \
            --out-dir "$DATA_DIR/alphafold" --workers "$WORKERS"
}

stage_esm() {
  local out="$DATA_DIR/models/$(basename "$ESM_MODEL")"
  if [[ -f "$out/.done" ]]; then log "  have $ESM_MODEL"; return 0; fi
  "$PYTHON" "$PY/esm_weights.py" --model "$ESM_MODEL" --out-dir "$out" && touch "$out/.done"
}

log "DATA_DIR=$DATA_DIR"
FAILED=()
for s in "${STAGES[@]}"; do
  if ! declare -F "stage_$s" > /dev/null; then
    log "unknown stage '$s' (choose from: ${ALL_STAGES[*]})"; FAILED+=("$s"); continue
  fi
  log "== $s"
  if "stage_$s" 2>&1 | tee -a "$DATA_DIR/logs/$s.log"; then
    log "== $s ok"
  else
    log "== $s FAILED (see $DATA_DIR/logs/$s.log)"; FAILED+=("$s")
  fi
done

du -sh "$DATA_DIR" 2>/dev/null
if [[ ${#FAILED[@]} -gt 0 ]]; then
  log "failed stages: ${FAILED[*]}  (rerun: bash $0 ${FAILED[*]})"
  exit 1
fi
log "all stages done"
