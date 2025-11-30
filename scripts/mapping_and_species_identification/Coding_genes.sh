# STARsolo preprocessing pipeline
# This script:
#   1) Processes GTF annotations
#   2) Masks ncRNA regions in the genome
#   3) Builds a STAR genome index
#   4) Runs STARsolo for multiple samples
#   5) Compresses STARsolo output matrices
# ---------------------------

# If not running under SLURM, fall back to 12 threads
THREADS=${SLURM_CPUS_PER_TASK:-12}

# Project root directory (edit this for your environment)
PROJECT_DIR=/path/to/project

# Reference files (relative to PROJECT_DIR or absolute paths)
GENOME_FASTA="${PROJECT_DIR}/reference/M1_genome.fna"
NCRNA_GTF="${PROJECT_DIR}/reference/M1_genome_all_ncRNA.gtf"
CDS_GTF="${PROJECT_DIR}/reference/M1_genome_all.cds.gtf"
GENOME_DIR="${PROJECT_DIR}/STAR_index/M1_all"

# Script and tool locations
PROCESS_GTF_SCRIPT="${PROJECT_DIR}/scripts/process_gtf.py"

STAR_BIN=STAR           # or /path/to/STAR
BEDTOOLS_BIN=bedtools   # or /path/to/bedtools

# Working directory where sample folders are located
WORK_DIR="${PROJECT_DIR}/samples"

# Sample list (folder names under WORK_DIR)
samples=(
  "Colon-DB-1" "Cecum-DB-1" "Colon-DB-2" "Cecum-DB-2"
  "Colon-WT-1" "Cecum-WT-1" "Colon-WT-2" "Cecum-WT-2"
  "Rectum-DB-1" "Rectum-DB-2" "Rectum-WT-1" "Rectum-WT-2"
)

cd "${PROJECT_DIR}"

# ---------------------------
# 1. Process GTF
# ---------------------------
# This script should generate:
#   - M1_genome_all.cds.gtf
#   - M1_genome_all_ncRNA.gtf
python "${PROCESS_GTF_SCRIPT}"

# ---------------------------
# 2. Mask ncRNA regions in the genome
# ---------------------------
${BEDTOOLS_BIN} maskfasta \
  -fi "${GENOME_FASTA}" \
  -bed "${NCRNA_GTF}" \
  -fo "${PROJECT_DIR}/reference/M1_ncRNA_masked.fasta"

# ---------------------------
# 3. Build STAR genome index
# ---------------------------
mkdir -p "${GENOME_DIR}"

${STAR_BIN} \
  --runThreadN "${THREADS}" \
  --runMode genomeGenerate \
  --genomeDir "${GENOME_DIR}" \
  --genomeFastaFiles "${PROJECT_DIR}/reference/M1_ncRNA_masked.fasta" \
  --sjdbGTFfile "${CDS_GTF}" \
  --sjdbOverhang 122 \
  --sjdbGTFfeatureExon gene \
  --sjdbGTFtagExonParentTranscript ID \
  --sjdbGTFtagExonParentGene ID \
  --genomeSAindexNbases 10 \
  --limitGenomeGenerateRAM=286034755850

# ---------------------------
# 4. Run STARsolo for each sample
# ---------------------------

run_star_solo() {
    local sample=$1

    # Each sample is assumed to be under ${WORK_DIR}/${sample}
    local sample_dir="${WORK_DIR}/${sample}"

    ${STAR_BIN} \
      --runThreadN "${THREADS}" \
      --soloType CB_UMI_Simple \
      --soloCBwhitelist None \
      --genomeDir "${GENOME_DIR}" \
      --soloCBstart 1 \
      --soloCBlen 20 \
      --soloUMIstart 21 \
      --soloUMIlen 8 \
      --readFilesIn "${sample_dir}/clean_data/${sample}_2.fq.gz" \
                    "${sample_dir}/clean_data/${sample}_1.fq.gz" \
      --outSAMtype BAM SortedByCoordinate \
      --outMultimapperOrder Random \
      --runRNGseed 1 \
      --outSAMattributes NH HI AS CR UR GX GN \
      --alignSJoverhangMin 1000 \
      --soloFeatures GeneFull \
      --outFileNamePrefix "${sample_dir}/${sample}_nocdhit_" \
      --soloStrand Reverse \
      --soloCellFilter TopCells 15000 \
      --outFilterScoreMinOverLread 0.5 \
      --readFilesCommand zcat
}

compress_files() {
    local sample=$1
    local sample_dir="${WORK_DIR}/${sample}"
    local solo_out="${sample_dir}/${sample}_nocdhit_Solo.out/GeneFull"

    gzip "${solo_out}/filtered/"* 
    gzip "${solo_out}/raw/"*
}

# ---------------------------
# 5. Loop over samples
# ---------------------------
for sample in "${samples[@]}"; do
    run_star_solo "${sample}"
    compress_files "${sample}"
done

wait

echo "All tasks completed!"
