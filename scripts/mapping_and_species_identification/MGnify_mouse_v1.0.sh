# Microbial scRNA-seq pipeline using the m-MGnify Mouse Gut v1.0 catalog
#
# Stage 1 (per sample):
#   - Merge / preprocess raw FASTQ with anchorbcm
#   - Remove host-aligned reads (mouse genome) with STAR
#   - Run Kraken2_smClassificaiton (Nubeam-dedup + Kraken2/Bracken) using m-MGnify mouse-gut DB
#
# Stage 2 (per sample):
#   - Build m-MGnify-based reference (FASTA + GTF) from Bracken species list
#   - Generate STAR genome index
#   - Run STARsolo for UMI-aware gene quantification
#   - Filter & sort count matrices
#
# All paths below are placeholders and should be customized by the user.
###############################################################################

set -euo pipefail

# -----------------------------
# 0. User-configurable paths
# -----------------------------

# Project root
PROJECT_DIR="/path/to/project"

# Raw FASTQ root, each sample in its own subfolder:
#   ${RAW_DIR}/${sample}/.../Sample_*-${sample}-*/...combined_R1.fastq.gz
RAW_DIR_1="/path/to/raw_run_1"   
RAW_DIR_2="/path/to/raw_run_2"  

# Intermediate output: per-sample clean_data under this directory
ADD_CECUM_DIR="${PROJECT_DIR}/03.Add_cecum"

# Host (mouse) STAR index for host-read removal
MOUSE_STAR_INDEX="/path/to/mouse_index"

# Python & tools
PYTHON_BIN="python"
STAR_BIN="STAR"
PIGZ_BIN="pigz"
ANCHORBCM_BIN="/path/to/anchorbcm"
ANCHORBCM_PREINDEX="/path/to/prev6v2.msr"

# Kraken2_smClassificaiton / classification (using m-MGnify Mouse Gut v1.0 Kraken DB)
Kraken2_smClassificaiton_SCRIPT="/path/to/microbiomePipe/Kraken2_smClassificaiton.py"
NUBEAM_DEDUP_BIN="/path/to/nubeam-dedup"
KRAKEN_BIN="/path/to/kraken2"
BRAKEN_BIN="/path/to/est_abundance_bracken"
MGNIFY_KRAKEN_DB_DIR="/path/to/mouse-gut/kraken"

# m-MGnify genome resources for building sample-specific reference
BUILD_REF_SCRIPT="/path/to/build_mgnify_reference.py"
MGNIFY_GENOME_MAP="/path/to/mgnify_mouse_v1/genome_map.tsv"
MGNIFY_GENOME_DIR="/path/to/mgnify_mouse_v1/genomes"


# -----------------------------
# 1. SLURM threads & sample list
# -----------------------------

THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Use the same 12 samples as in the MGBC script
samples=(
  "Colon-DB-1"   "Cecum-DB-1"
  "Colon-DB-2"   "Cecum-DB-2"
  "Colon-WT-1"   "Cecum-WT-1"
  "Colon-WT-2"   "Cecum-WT-2"
  "Rectum-DB-1"  "Rectum-DB-2"
  "Rectum-WT-1"  "Rectum-WT-2"
)

# -----------------------------
# 2. Stage 1: sequence filtering & host removal
# -----------------------------
process_sample() {
    local sample="$1"
    echo "[INFO] Stage 1 – processing sample: ${sample}"

    # Per-sample clean_data folder
    local sample_root="${ADD_CECUM_DIR}/${sample}"
    local clean_dir="${sample_root}/clean_data"
    mkdir -p "${clean_dir}"

    # Build possible input paths (two runs / lanes)
    # Pattern: .../Sample_*-${sample}-*/...combined_R1/2.fastq.gz
    local R1_FILES=(
      "${RAW_DIR_1}"/Sample_*-"${sample}"-*/*combined_R1.fastq.gz
      "${RAW_DIR_2}"/Sample_*-"${sample}"-*/*combined_R1.fastq.gz
    )
    local R2_FILES=(
      "${RAW_DIR_1}"/Sample_*-"${sample}"-*/*combined_R2.fastq.gz
      "${RAW_DIR_2}"/Sample_*-"${sample}"-*/*combined_R2.fastq.gz
    )

    # Check if any files exist
    shopt -s nullglob
    local r1_count=${#R1_FILES[@]}
    local r2_count=${#R2_FILES[@]}
    shopt -u nullglob

    if [[ "${r1_count}" -eq 0 || "${r2_count}" -eq 0 ]]; then
        echo "[WARN] Input FASTQ not found for ${sample}; skipping."
        return 1
    fi

    # 2.1 anchorbcm: merge + barcode / UMI processing
    "${ANCHORBCM_BIN}" \
      -1 <(for f in "${R1_FILES[@]}"; do zcat "$f" 2>/dev/null; done) \
      -2 <(for f in "${R2_FILES[@]}"; do zcat "$f" 2>/dev/null; done) \
      -b "${ANCHORBCM_PREINDEX}" \
      -o "${clean_dir}/${sample}_"

    cd "${sample_root}"

    # 2.2 Remove host-aligned reads with STAR (mouse reference)
    "${STAR_BIN}" \
      --runThreadN "${THREADS}" \
      --soloType CB_UMI_Simple \
      --soloCBwhitelist None \
      --genomeDir "${MOUSE_STAR_INDEX}" \
      --soloCBstart 1 --soloCBlen 20 \
      --soloUMIstart 21 --soloUMIlen 8 \
      --readFilesIn "clean_data/${sample}_2.fq" "clean_data/${sample}_1.fq" \
      --outSAMtype BAM Unsorted \
      --outMultimapperOrder Random \
      --runRNGseed 1 \
      --outSAMattributes NH HI AS CR UR GX GN \
      --alignSJoverhangMin 1000 \
      --soloFeatures GeneFull \
      --outFileNamePrefix "${sample}_filt_" \
      --soloStrand Reverse \
      --limitBAMsortRAM 40000000000 \
      --outReadsUnmapped Fastx \
      --outFilterScoreMinOverLread 0.5

    # Rename: keep unmapped microbial reads as new input
    mv "clean_data/${sample}_2.fq" "clean_data/${sample}_pre_2.fq"
    mv "clean_data/${sample}_1.fq" "clean_data/${sample}_pre_1.fq"
    mv "${sample}_filt_Unmapped.out.mate2" "clean_data/${sample}_1.fq"
    mv "${sample}_filt_Unmapped.out.mate1" "clean_data/${sample}_2.fq"

    # 2.3 Kraken2_smClassificaiton: dedup + Kraken2/Bracken classification (m-MGnify mouse-gut DB)
    echo "[INFO] Start Kraken2_smClassificaiton for ${sample}: $(date)"

    cat > config.ini <<EOF
sample=${sample}
outdir=${sample_root}
process=${THREADS}
nubeam_dedup=${NUBEAM_DEDUP_BIN}
kraken=${KRAKEN_BIN}
braken=${BRAKEN_BIN}
filter_threshold=0
krakenDb=${MGNIFY_KRAKEN_DB_DIR}
EOF

    "${PYTHON_BIN}" "${Kraken2_smClassificaiton_SCRIPT}" --cfg config.ini

    echo "[INFO] Done Kraken2_smClassificaiton for ${sample}: $(date)"
    echo
}

# -----------------------------
# 3. Stage 2: mapping & quantification (m-MGnify)
# -----------------------------
process_mapping() {
    local sample="$1"
    echo "[INFO] Stage 2 – mapping & quantification: ${sample}"

    local sample_root="${ADD_CECUM_DIR}/${sample}"
    cd "${ADD_CECUM_DIR}"

    # 3.1 Species list from Bracken report
    local tax_report="${sample_root}/Result/${sample}_sc_taxonomy.report"
    local species_list="${sample_root}/species.txt"

    if [[ ! -f "${tax_report}" ]]; then
        echo "[ERROR] Taxonomy report not found: ${tax_report}" >&2
        return 1
    fi

    tail -n +2 "${tax_report}" \
      | awk -F '\t' '$7 >= 0 {print $2}' \
      | sort -u > "${species_list}"

    # 3.2 Build m-MGnify-based reference (FASTA + GTF)
    local ref_dir="${sample_root}/ref"
    mkdir -p "${ref_dir}"
    local genome_fa="${ref_dir}/genome.fna"
    local genome_gtf="${ref_dir}/genome.gtf"
    local acc_list="${sample_root}/acc.txt"

    touch "${genome_fa}" "${genome_gtf}"

    "${PYTHON_BIN}" "${BUILD_REF_SCRIPT}" \
      "${MGNIFY_GENOME_MAP}" \
      "${species_list}" \
      "${acc_list}" \
      "${MGNIFY_GENOME_DIR}" \
      "${genome_fa}" \
      "${genome_gtf}"

    # 3.3 STAR genome index for m-MGnify subset
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
      --sjdbGTFtagExonParentTranscript ID \
      --sjdbGTFtagExonParentGene ID \
      --genomeSAindexNbases 10 \
      --limitGenomeGenerateRAM 500000000000 \
      --limitSjdbInsertNsj 10000000

    # 3.4 Ensure clean_data FASTQ are uncompressed
    local clean_dir="${sample_root}/clean_data"
    if compgen -G "${clean_dir}/*.gz" > /dev/null; then
        echo "[INFO] Decompressing clean_data FASTQ for ${sample} ..."
        "${PIGZ_BIN}" -d "${clean_dir}"/*.gz
    fi

    # 3.5 STARsolo (microbial reference)
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
      --outFileNamePrefix "${sample_root}/${sample}_combine_name_" \
      --soloStrand Reverse \
      --limitBAMsortRAM 40000000000 \
      --soloCellFilter TopCells 15000 \
      --outFilterScoreMinOverLread 0.5

    # 3.6 Downstream matrix filtering
    local solo_dir="${sample_root}/${sample}_Solo.out/GeneFull"
    local raw_dir="${solo_dir}/raw"
    local raw_sorted_dir="${solo_dir}/rawSorted"

    mkdir -p "${raw_sorted_dir}"

    "${FILTER_MTX_BIN}" \
      -i "${raw_dir}" \
      -o "${raw_sorted_dir}"
}

# -----------------------------
# 4. Cleanup (optional)
# -----------------------------
cleanup_sample() {
    local sample="$1"

    local sample_root="${ADD_CECUM_DIR}/${sample}"
    local clean_dir="${sample_root}/clean_data"
    local result_dir="${sample_root}/Result"

    # Large intermediate FASTQ
    rm -f "${clean_dir}/${sample}.dedup.fq" \
          "${clean_dir}/${sample}.prededup.fq" \
          "${clean_dir}/${sample}_final_1.fq" \
          "${clean_dir}/${sample}_final_2.fq"

    # Pre-host-filter FASTQ (optional to keep)
    if [[ -f "${clean_dir}/${sample}_pre_1.fq" ]]; then
        "${PIGZ_BIN}" -p 8 "${clean_dir}/${sample}_pre_1.fq"
    fi
    if [[ -f "${clean_dir}/${sample}_pre_2.fq" ]]; then
        "${PIGZ_BIN}" -p 8 "${clean_dir}/${sample}_pre_2.fq"
    fi

    # Kraken classified/unclassified reads
    rm -f "${result_dir}/output_classified.fq" \
          "${result_dir}/output_unclassified.fq" \
          "${result_dir}/${sample}_kraken.output"

    # Compress final clean_data FASTQ (optional)
    if [[ -f "${clean_dir}/${sample}_1.fq" ]]; then
        "${PIGZ_BIN}" -p 8 "${clean_dir}/${sample}_1.fq"
    fi
    if [[ -f "${clean_dir}/${sample}_2.fq" ]]; then
        "${PIGZ_BIN}" -p 8 "${clean_dir}/${sample}_2.fq"
    fi

    # Compress STARsolo matrices (if exist)
    local solo_dir="${sample_root}/${sample}_Solo.out/GeneFull"
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

mkdir -p "${ADD_CECUM_DIR}"

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
