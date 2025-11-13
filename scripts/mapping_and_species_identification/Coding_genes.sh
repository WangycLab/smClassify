#!/bin/bash
#SBATCH -J coding
#SBATCH -p fat
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12

THREADS=$SLURM_CPUS_PER_TASK
samples=("241224-XGR-D71JC" "241224-XGR-D71MC" "241224-XGR-D72JC" "241224-XGR-D72MC" "241224-XGR-W71JC" "241224-XGR-W71MC" \
 "241224-XGR-W72JC" "241224-XGR-W72MC" "41217N01" "41217N02" "41217N03" "41217N04")

# Process gtf file to split files
python process_gtf.py

# mask ncRNA region of genome
bedtools maskfasta -fi M1_genome.fna -bed M1_genome_all_ncRNA.gtf -fo M1_ncRNA_masked.fasta


/public/home/wangycgroup/public/software/bin/STAR --runThreadN $THREADS --runMode genomeGenerate --genomeDir ./M1_all \
 --genomeFastaFiles ./M1_ncRNA_masked.fasta --sjdbGTFfile M1_genome_all.cds.gtf --sjdbOverhang 122 --sjdbGTFfeatureExon gene \
 --sjdbGTFtagExonParentTranscript ID --sjdbGTFtagExonParentGene ID --genomeSAindexNbases 10 --limitGenomeGenerateRAM=286034755850


# Run STAR Solo
run_star_solo() {
    local sample=$1

    /public/home/wangycgroup/public/software/bin/STAR --runThreadN $THREADS --soloType CB_UMI_Simple \
        --soloCBwhitelist None --genomeDir ./M1_all \
        --soloCBstart 1 --soloCBlen 20 --soloUMIstart 21 --soloUMIlen 8 \
        --readFilesIn ${sample}/clean_data/${sample}_2.fq.gz ${sample}/clean_data/${sample}_1.fq.gz \
        --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --runRNGseed 1 \
        --outSAMattributes NH HI AS CR UR GX GN --alignSJoverhangMin 1000 \
        --soloFeatures GeneFull --outFileNamePrefix ${sample}/${sample}_nocdhit_ \
        --soloStrand Reverse --soloCellFilter TopCells 15000 \
        --outFilterScoreMinOverLread 0.5 --readFilesCommand zcat
}

# Compress files
compress_files() {
    local sample=$1
    gzip /public/home/wangycgroup/xinlong/mouse_colon/03.Add_cecum/${sample}/${sample}_nocdhit_Solo.out/GeneFull/filtered/*
    gzip /public/home/wangycgroup/xinlong/mouse_colon/03.Add_cecum/${sample}/${sample}_nocdhit_Solo.out/GeneFull/raw/* 

}

for sample in "${samples[@]}"; do
    run_star_solo "$sample"
    compress_files "$sample"
done

wait

echo "All tasks completed!"
