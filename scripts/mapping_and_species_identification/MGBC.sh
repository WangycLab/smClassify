# Microbial scRNA-seq pipeline using the MGBC reference genome catalog
#
# Stage 1 (per sample):
#   - Prepare per-sample working directory
#   - Link / decompress FASTQ files
#   - Run Kraken2/Bracken classification
#
# Stage 2 (per sample):
#   - Build MGBC-based reference (FASTA + GTF) from classified species list
#   - Generate STAR index
#   - Run STARsolo for UMI-aware gene quantification
#
# All paths below are placeholders and should be customized by the user.
###############################################################################

set -euo pipefail

# -----------------------------
# 0. User-configurable paths
# -----------------------------

# Project root (this script can live here)
PROJECT_DIR="/path/to/project"

# Raw FASTQ root (each sample has its own subfolder under this directory)
# Expecting: ${RAW_DIR}/${sample}/clean_data/${sample}_1.fq(.gz), ${sample}_2.fq(.gz)
RAW_DIR="/path/to/raw_data/03.Add_cecum"

# Output root for MGBC-based workflow (per-sample output under this folder)
OUT_DIR="${PROJECT_DIR}/04.MGBC"

# Python and tools
PYTHON_BIN="python"   # or absolute path, e.g. /opt/conda/bin/python
STAR_BIN="STAR"       # or absolute path, e.g. /usr/local/bin/STAR
PIGZ_BIN="pigz"       # for (re)compression of FASTQ

# Kranken2_smClassification / classification pipeline
Kranken2_smClassification_SCRIPT="/path/to/Kraken2-based classification pipline.py"
NUBEAM_DEDUP_BIN="/path/to/nubeam-dedup"
KRAKEN_BIN="/path/to/kraken2"
BRAKEN_BIN="/path/to/est_abundance_bracken"
KRAKEN_DB_DIR="/path/to/kraken2_MGBC_db"

# MGBC genome resources
MGBC_GENOME_MAP="/path/to/MGBC/genomes/genome_map.tsv"
MGBC_GENOME_DIR="/path/to/MGBC/genomes"

# Script to build sample-specific reference from MGBC species list
BUILD_REF_SCRIPT="/path/to/build_mgbc_reference.py"

# -----------------------------
# 1. SLURM threads & sample list
# -----------------------------

THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Sample names (update to your own IDs if needed)
samples=(
  "Colon-DB-1"   "Cecum-DB-1"
  "Colon-DB-2"   "Cecum-DB-2"
  "Colon-WT-1"   "Cecum-WT-1"
  "Colon-WT-2"   "Cecum-WT-2"
  "Rectum-DB-1"  "Rectum-DB-2"
  "Rectum-WT-1"  "Rectum-WT-2"
)

# -----------------------------
# 2. Stage 1: per-sample preprocessing + classification
# -----------------------------
process_sample() {
    local sample="$1"

    # Per-sample working directory under OUT_DIR
    local sample_dir="${OUT_DIR}/${sample}"
    local clean_dir="${sample_dir}/clean_data"
    mkdir -p "${clean_dir}"

    # Raw data location for this sample
    local raw_clean_dir="${RAW_DIR}/${sample}/clean_data"

    # Decompress any gzipped FASTQ in raw folder (if needed)
    if compgen -G "${raw_clean_dir}/*.gz" > /dev/null; then
        echo "[INFO] Decompressing raw FASTQ for ${sample} ..."
        "${PIGZ_BIN}" -d "${raw_clean_dir}"/*.gz
    fi

    # Link FASTQ files into the working clean_data directory
    ln -sf "${raw_clean_dir}/${sample}_1.fq" "${clean_dir}/"
    ln -sf "${raw_clean_dir}/${sample}_2.fq" "${clean_dir}/"

    cd "${sample_dir}"

    # ----------------- Kranken2_smClassification configuration -----------------
    echo "[INFO] Start Kranken2_smClassification for ${sample}: $(date)"

    cat > config.ini <<EOF
sample=${sample}
outdir=${sample_dir}
process=${THREADS}
nubeam_dedup=${NUBEAM_DEDUP_BIN}
kraken=${KRAKEN_BIN}
braken=${BRAKEN_BIN}
filter_threshold=0
krakenDb=${KRAKEN_DB_DIR}
# mode=m20  # optional: barcode mode, e.g. 20 bp cell barcode
EOF

    # Run Kranken2_smClassification pipeline
    "${PYTHON_BIN}" "${Kranken2_smClassification_SCRIPT}" --cfg config.ini

    echo "[INFO] Kranken2_smClassification finished for ${sample}: $(date)"
    echo
}

# -----------------------------
# 3. Stage 2: mapping & quantification with MGBC
# -----------------------------
process_mapping() {
    local sample="$1"

    cd "${OUT_DIR}"

    local sample_dir="${OUT_DIR}/${sample}"
    local result_dir="${sample_dir}/Result"
    local ref_dir="${sample_dir}/ref"

    mkdir -p "${ref_dir}"

    # 3.1 Generate species list from Bracken/kraken taxonomy report
    local tax_report="${result_dir}/${sample}_sc_taxonomy.report"
    local species_list="${sample_dir}/species.txt"

    if [[ ! -f "${tax_report}" ]]; then
        echo "[ERROR] Taxonomy report not found for ${sample}: ${tax_report}" >&2
        return 1
    fi

    tail -n +2 "${tax_report}" \
      | awk -F '\t' '$7 >= 0 {print $2}' \
      | sort -u > "${species_list}"

    # 3.2 Build MGBC-based reference (FASTA + GTF)
    local genome_fa="${ref_dir}/genome.fna"
    local genome_gtf="${ref_dir}/genome.gtf"
    local acc_list="${sample_dir}/acc.txt"

    touch "${genome_fa}" "${genome_gtf}"

    "${PYTHON_BIN}" "${BUILD_REF_SCRIPT}" \
        "${MGBC_GENOME_MAP}" \
        "${species_list}" \
        "${acc_list}" \
        "${MGBC_GENOME_DIR}" \
        "${genome_fa}" \
        "${genome_gtf}"

    # 3.3 Generate STAR genome index
    local star_index="${ref_dir}/STAR_index"
    mkdir -p "${star_index}"

    "${STAR_BIN}" \
      --runThreadN "${THREADS}" \
      --runMode genomeGenerate \
      --genomeDir "${star_index}" \
      --genomeFastaFiles "${genome_fa}" \
      --sjdbGTFfile "${genome_gtf}" \
      --sjdbOverhang 122 \
      --sjdbGTFfeatureExon gene \
      --sjdbGTFtagExonParentTranscript gene_id \
      --sjdbGTFtagExonParentGene gene_id \
      --genomeSAindexNbases 10 \
      --limitGenomeGenerateRAM 500000000000

    # 3.4 Run STARsolo (CB/UMI layout: 20 bp CB + 8 bp UMI on read1)
    local clean_dir="${sample_dir}/clean_data"

    "${STAR_BIN}" \
      --runThreadN "${THREADS}" \
      --soloType CB_UMI_Simple \
      --soloCBwhitelist None \
      --genomeDir "${star_index}" \
      --soloCBstart 1 \
      --soloCBlen 20 \
      --soloUMIstart 21 \
      --soloUMIlen 8 \
      --readFilesIn "${clean_dir}/${sample}_2.fq" "${clean_dir}/${sample}_1.fq" \
      --outSAMtype BAM Unsorted \
      --outMultimapperOrder Random \
      --runRNGseed 1 \
      --outSAMattributes NH HI AS CR UR GX GN \
      --alignSJoverhangMin 1000 \
      --soloFeatures GeneFull \
      --outFileNamePrefix "${sample_dir}/${sample}_" \
      --soloStrand Reverse \
      --limitBAMsortRAM 40000000000 \
      --soloCellFilter TopCells 15000 \
      --outFilterScoreMinOverLread 0.5

    # 3.5 Downstream matrix filtering + QC plotting
    local solo_dir="${sample_dir}/${sample}_Solo.out/GeneFull"
    local raw_dir="${solo_dir}/raw"
    local raw_sorted_dir="${solo_dir}/rawSorted"

    mkdir -p "${raw_sorted_dir}"

    "${FILTER_MTX_BIN}" \
      -i "${raw_dir}" \
      -o "${raw_sorted_dir}"

    "${PYTHON_BIN}" "${PLOT_SCATTER_SCRIPT}" "${sample}"
}

# -----------------------------
# 4. Cleanup utilities
# -----------------------------
cleanup_sample() {
    local sample="$1"

    local sample_dir="${OUT_DIR}/${sample}"
    local clean_dir="${sample_dir}/clean_data"
    local result_dir="${sample_dir}/Result"

    # Remove large intermediate FASTQ
    rm -f "${clean_dir}/${sample}.dedup.fq" \
          "${clean_dir}/${sample}.prededup.fq" \
          "${clean_dir}/${sample}_final_1.fq" \
          "${clean_dir}/${sample}_final_2.fq"

    # Remove intermediate Kraken classified/unclassified reads and output
    rm -f "${result_dir}/output_classified.fq" \
          "${result_dir}/output_unclassified.fq" \
          "${result_dir}/${sample}_kraken.output"

    # Re-compress original raw FASTQ files (optional)
    local raw_clean_dir="${RAW_DIR}/${sample}/clean_data"
    if compgen -G "${raw_clean_dir}/${sample}_1.fq" > /dev/null; then
        "${PIGZ_BIN}" -p 8 "${raw_clean_dir}/${sample}_1.fq"
    fi
    if compgen -G "${raw_clean_dir}/${sample}_2.fq" > /dev/null; then
        "${PIGZ_BIN}" -p 8 "${raw_clean_dir}/${sample}_2.fq"
    fi

    # Compress STARsolo matrices (optional)
    local solo_dir="${sample_dir}/${sample}_Solo.out/GeneFull"
    if [ -d "${solo_dir}/filtered" ]; then
        gzip -f "${solo_dir}/filtered"/*
    fi
    if [ -d "${solo_dir}/raw" ]; then
        gzip -f "${solo_dir}/raw"/*
    fi
}

# -----------------------------
# 5. Main loop
# -----------------------------
mkdir -p "${OUT_DIR}"

for sample in "${samples[@]}"; do
    echo "============================="
    echo "[INFO] Processing sample: ${sample}"
    echo "============================="

    # Stage 1: FASTQ preprocessing + classification (MGBC Kraken2/Bracken)
    process_sample "${sample}"

    # Stage 2: MGBC-based mapping + quantification
    process_mapping "${sample}"

    # Cleanup (optional but recommended for large datasets)
    cleanup_sample "${sample}"

    echo "[INFO] Finished sample: ${sample}"
    echo
done

echo "[INFO] All tasks completed!"
