# Microbial scRNA-seq pipeline using a RefSeq-based reference catalog
#
# Stage 1 (per sample):
#   - Symlink host-filtered FASTQ from a shared "clean_data" directory
#   - Run Kraken2_smClassification (Nubeam-dedup + Kraken2/Bracken) with a RefSeq Kraken2 DB
#
# Stage 2 (per sample):
#   - Build sample-specific reference (FASTA + GTF) from RefSeq genomes
#   - Generate STAR genome index
#   - Run STARsolo for UMI-aware gene quantification
#   - Filter & sort count matrices and generate QC scatter plot
#
# All paths below are placeholders and should be customized by the user.
###############################################################################

set -euo pipefail

# -----------------------------
# 0. User-configurable paths
# -----------------------------

# Working directory where this script is launched and where per-sample folders
# will be created (e.g. RefSeq project root)
PROJECT_ROOT="/path/to/RefSeq_project"

# Directory containing host-filtered read pairs for each sample:
#   ${CLEAN_SOURCE_ROOT}/${sample}/clean_data/${sample}_1.fq(.gz)
#   ${CLEAN_SOURCE_ROOT}/${sample}/clean_data/${sample}_2.fq(.gz)
CLEAN_SOURCE_ROOT="/path/to/host_filtered_clean_data"

# Kraken2_smClassification / classification (RefSeq DB)
PYTHON_BIN="python"
Kraken2_smClassification_SCRIPT="/path/to/microbiomePipe/Kraken2-based classification pipline.py"
NUBEAM_DEDUP_BIN="/path/to/nubeam-dedup"
KRAKEN_BIN="/path/to/kraken2"
BRAKEN_BIN="/path/to/est_abundance_bracken"
REFSEQ_KRAKEN_DB="/path/to/kraken2_refseq_db"

# Build RefSeq-based reference (selectGenomes.py)
SELECT_GENOMES_SCRIPT="/path/to/selectGenomes.py"
REFSEQ_MAP_TSV="/path/to/refseq_mapping.tsv"
REFSEQ_RAW_GENOME_DIR="/path/to/refseq_raw_genomes"

# STAR & utilities
STAR_BIN="STAR"
PIGZ_BIN="pigz"
FILTER_MTX_BIN="/path/to/filtermtx"       # filtermtx.o
SCATTER_PLOT_SCRIPT="/path/to/plotScatter-deliver.py"

# Optional: extra libraries
# export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/path/to/htslib:/path/to/extra_libs/*"

# -----------------------------
# 1. SLURM threads & sample list
# -----------------------------

THREADS="${SLURM_CPUS_PER_TASK:-24}"

# Example sample list (12 microbial scRNA-seq libraries)
samples=(
  "Colon-DB-1"   "Cecum-DB-1"
  "Colon-DB-2"   "Cecum-DB-2"
  "Colon-WT-1"   "Cecum-WT-1"
  "Colon-WT-2"   "Cecum-WT-2"
  "Rectum-DB-1"  "Rectum-DB-2"
  "Rectum-WT-1"  "Rectum-WT-2"
)

cd "${PROJECT_ROOT}"

# -----------------------------
# 2. Stage 1: sequence filtering / Kraken2_smClassification
# -----------------------------
process_sample() {
    local sample="$1"
    echo "[INFO] Stage 1 – Kraken2_smClassification for sample: ${sample}"

    local sample_dir="${PROJECT_ROOT}/${sample}"
    local clean_dir="${sample_dir}/clean_data"

    # Create working and clean_data dirs if not present
    if [[ ! -d "${sample_dir}" ]]; then
        mkdir -p "${clean_dir}"

        # Decompress source clean_data if needed
        local src_clean="${CLEAN_SOURCE_ROOT}/${sample}/clean_data"
        for file in "${src_clean}"/*; do
            if [[ -f "$file" && "$file" == *.gz ]]; then
                echo "[INFO] Decompressing $file ..."
                "${PIGZ_BIN}" -d -p 8 "$file"
            fi
        done

        # Symlink input FASTQs (host-filtered, CB+UMI layout)
        ln -s "${src_clean}/${sample}_1.fq" "${clean_dir}/${sample}_1.fq"
        ln -s "${src_clean}/${sample}_2.fq" "${clean_dir}/${sample}_2.fq"
    fi

    cd "${sample_dir}"

    # Run Kraken2_smClassification (Nubeam-dedup + Kraken2/Bracken on RefSeq DB)
    echo "[INFO] Start Kraken2_smClassification for ${sample}: $(date)"

    cat > config.ini <<EOF
sample=${sample}
outdir=${PROJECT_ROOT}/${sample}
process=${THREADS}
nubeam_dedup=${NUBEAM_DEDUP_BIN}
kraken=${KRAKEN_BIN}
braken=${BRAKEN_BIN}
filter_threshold=0
krakenDb=${REFSEQ_KRAKEN_DB}
EOF

    "${PYTHON_BIN}" "${Kraken2_smClassification_SCRIPT}" --cfg config.ini

    echo "[INFO] Done Kraken2_smClassification for ${sample}: $(date)"
    echo
}

# -----------------------------
# 3. Stage 2: mapping & quantification
# -----------------------------
process_mapping() {
    local sample="$1"
    echo "[INFO] Stage 2 – mapping & quantification for sample: ${sample}"

    cd "${PROJECT_ROOT}"

    local sample_dir="${PROJECT_ROOT}/${sample}"
    local result_dir="${sample_dir}/Result"
    local clean_dir="${sample_dir}/clean_data"

    # 3.1 Build species list from Bracken report
    local tax_report="${result_dir}/${sample}_sc_taxonomy.report"
    local species_list="${sample_dir}/species.txt"

    if [[ ! -f "${tax_report}" ]]; then
        echo "[ERROR] Taxonomy report not found: ${tax_report}" >&2
        return 1
    fi

    # Use column 3 (NCBI taxid) or adjust as needed
    tail -n +2 "${tax_report}" \
      | cut -f 3 \
      | sort -n \
      | uniq > "${species_list}"

    # 3.2 Extract RefSeq genomes and build FASTA + GTF
    local ref_dir="${sample_dir}/ref"
    mkdir -p "${ref_dir}"

    local genome_fa="${ref_dir}/genome.fna"
    local genome_gtf="${ref_dir}/genome.gtf"

    touch "${genome_fa}" "${genome_gtf}"

    echo "[INFO] Start extracting RefSeq genomes for ${sample}: $(date)"
    "${PYTHON_BIN}" "${SELECT_GENOMES_SCRIPT}" \
        -i "${species_list}" \
        -r "${REFSEQ_MAP_TSV}" \
        -d "${REFSEQ_RAW_GENOME_DIR}" \
        -o "${ref_dir}/genome"
    echo "[INFO] Done extracting RefSeq genomes for ${sample}: $(date)"

    # 3.3 Generate STAR index
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
        --limitGenomeGenerateRAM 495864700000

    # 3.4 Run STARsolo on the per-sample RefSeq reference
    "${STAR_BIN}" \
        --runThreadN "${THREADS}" \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist None \
        --genomeDir "${star_index}" \
        --soloCBstart 1 --soloCBlen 20 \
        --soloUMIstart 21 --soloUMIlen 8 \
        --readFilesIn "${clean_dir}/${sample}_2.fq" "${clean_dir}/${sample}_1.fq" \
        --outSAMtype BAM Unsorted \
        --outMultimapperOrder Random \
        --runRNGseed 1 \
        --outSAMattributes NH HI AS CR UR GX GN \
        --alignSJoverhangMin 1000 \
        --soloFeatures GeneFull \
        --outFileNamePrefix "${sample_dir}/${sample}_" \
        --soloStrand Reverse \
        --limitBAMsortRAM 41143265264 \
        --soloCellFilter TopCells 15000 \
        --outFilterScoreMinOverLread 0.5

    # 3.5 Downstream matrix filtering and QC plot
    local solo_dir="${sample_dir}/${sample}_Solo.out/GeneFull"
    local raw_dir="${solo_dir}/raw"
    local raw_sorted_dir="${solo_dir}/rawSorted"

    mkdir -p "${raw_sorted_dir}"

    "${FILTER_MTX_BIN}" \
        -i "${raw_dir}" \
        -o "${raw_sorted_dir}"

    "${PYTHON_BIN}" "${SCATTER_PLOT_SCRIPT}" "${sample}"
}

# -----------------------------
# 4. Cleanup (optional but recommended)
# -----------------------------
cleanup_sample() {
    local sample="$1"
    echo "[INFO] Cleanup intermediate files for sample: ${sample}"

    local sample_dir="${PROJECT_ROOT}/${sample}"
    local clean_dir="${sample_dir}/clean_data"
    local result_dir="${sample_dir}/Result"

    # Large intermediate FASTQs
    rm -f "${clean_dir}/${sample}.dedup.fq" \
          "${clean_dir}/${sample}.prededup.fq" \
          "${clean_dir}/${sample}_final_1.fq" \
          "${clean_dir}/${sample}_final_2.fq"

    # Optional: pre-host-filter FASTQs from source directory can also be gzipped
    if compgen -G "${clean_dir}/*pre_*.fq" > /dev/null; then
        "${PIGZ_BIN}" -p 8 "${clean_dir}"/"*pre_"*.fq
    fi

    # Kraken classified/unclassified reads
    rm -f "${result_dir}/output_classified.fq" \
          "${result_dir}/output_unclassified.fq" \
          "${result_dir}/${sample}_kraken.output"

    # Compress original clean_data FASTQ again if desired
    if [[ -f "${CLEAN_SOURCE_ROOT}/${sample}/clean_data/${sample}_1.fq" ]]; then
        "${PIGZ_BIN}" -p 8 "${CLEAN_SOURCE_ROOT}/${sample}/clean_data/${sample}_1.fq"
    fi
    if [[ -f "${CLEAN_SOURCE_ROOT}/${sample}/clean_data/${sample}_2.fq" ]]; then
        "${PIGZ_BIN}" -p 8 "${CLEAN_SOURCE_ROOT}/${sample}/clean_data/${sample}_2.fq"
    fi

    # Compress STARsolo matrices
    local solo_dir="${sample_dir}/${sample}_Solo.out/GeneFull"
    if [[ -d "${solo_dir}/filtered" ]]; then
        gzip -f "${solo_dir}/filtered"/*
    fi
    if [[ -d "${solo_dir}/raw" ]]; then
        gzip -f "${solo_dir}/raw"/*
    fi
}

# -----------------------------
# 5. Main loop
# -----------------------------

mkdir -p "${PROJECT_ROOT}"

for sample in "${samples[@]}"; do
    echo "========================================="
    echo "[INFO] Processing sample: ${sample}"
    echo "========================================="

    process_sample "${sample}"
    process_mapping "${sample}"
    cleanup_sample "${sample}"

    echo "[INFO] Finished sample: ${sample}"
    echo
done

echo "[INFO] All tasks completed!"
