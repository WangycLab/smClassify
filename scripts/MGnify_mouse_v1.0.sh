#!/bin/bash
#SBATCH -J MGnify_mouse_v1.0
#SBATCH -p cu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/public/home/wangycgroup/public/panhuaran/htslib:/public/home/wangycgroup/public/software/lib2/*

# Select samples to process
samples=("241224-XGR-D71JC" "241224-XGR-D71MC" "241224-XGR-D72JC" "241224-XGR-D72MC" "241224-XGR-W71JC" "241224-XGR-W71MC" "241224-XGR-W72JC" "241224-XGR-W72MC" "41217N01" "41217N02" "41217N03" "41217N04")

THREADS=$SLURM_CPUS_PER_TASK

# Step 1: Sequence filtering and processing
process_sample() {
    local sample=$1
    echo "Processing sample: $sample"
    
    # Create clean data directory
    mkdir -p /public/home/wangycgroup/xinlong/mouse_colon/03.Add_cecum/${sample}/clean_data
    
    # Build all possible input file paths
    R1_FILES=(/public/home/wangycgroup/LZ_public/LZ_Raw2025/241224/Sample_*-${sample}-*/*combined_R1.fastq.gz /public/home/wangycgroup/LZ_public/LZ_Raw2025/241224add/Sample_*-${sample}-*/*combined_R1.fastq.gz)
    R2_FILES=(/public/home/wangycgroup/LZ_public/LZ_Raw2025/241224/Sample_*-${sample}-*/*combined_R2.fastq.gz /public/home/wangycgroup/LZ_public/LZ_Raw2025/241224add/Sample_*-${sample}-*/*combined_R2.fastq.gz)
    
    # Check if files are found
    if [ ${#R1_FILES[@]} -eq 0 ] || [ ${#R2_FILES[@]} -eq 0 ]; then
        echo "Warning: Input files for sample ${sample} not found"
        return 1
    fi
    
    # anchorbcm processing, using process substitution to merge multiple input files
    /public/home/wangycgroup/public/software/anchorbcm.o \
        -1 <(for f in "${R1_FILES[@]}"; do zcat "$f" 2>/dev/null; done) \
        -2 <(for f in "${R2_FILES[@]}"; do zcat "$f" 2>/dev/null; done) \
        -b /public/home/wangycgroup/public/Sendru/prev6v2.msr \
        -o /public/home/wangycgroup/xinlong/mouse_colon/03.Add_cecum/${sample}/clean_data/${sample}_
    
    cd /public/home/wangycgroup/xinlong/mouse_colon/03.Add_cecum/${sample}
    
    # STAR remove reads aligned to host genome
    /public/home/wangycgroup/public/software/bin/STAR --runThreadN $THREADS --soloType CB_UMI_Simple \
        --soloCBwhitelist None --genomeDir /public/home/wangycgroup/xinlong/mouse_colon/03.Add_cecum/mouse_index \
        --soloCBstart 1 --soloCBlen 20 --soloUMIstart 21 --soloUMIlen 8 \
        --readFilesIn clean_data/${sample}_2.fq clean_data/${sample}_1.fq \
        --outSAMtype BAM Unsorted --outMultimapperOrder Random --runRNGseed 1 \
        --outSAMattributes NH HI AS CR UR GX GN --alignSJoverhangMin 1000 \
        --soloFeatures GeneFull --outFileNamePrefix ${sample}_filt_ \
        --soloStrand Reverse --limitBAMsortRAM 41143265264 \
        --outReadsUnmapped Fastx --outFilterScoreMinOverLread 0.5
    
    # Rename files
    mv clean_data/${sample}_2.fq clean_data/${sample}_pre_2.fq
    mv clean_data/${sample}_1.fq clean_data/${sample}_pre_1.fq
    mv ${sample}_filt_Unmapped.out.mate2 clean_data/${sample}_1.fq
    mv ${sample}_filt_Unmapped.out.mate1 clean_data/${sample}_2.fq
    
    # Run scMeta for species classification
    echo "# Start scMeta for ${sample}: $(date)"
    touch config.ini
    echo "sample=${sample}" > config.ini
	echo "outdir=/public/home/wangycgroup/xinlong/mouse_colon/03.Add_cecum/${sample}" >> config.ini
	echo "process=$THREADS" >> config.ini
	echo "nubeam_dedup=/public/home/wangycgroup/public/software/nubeamdedup-master/Linux/nubeam-dedup" >> config.ini
	echo "kraken=/public/home/wangycgroup/public/software/kraken2/kraken2" >> config.ini
	echo "braken=/public/home/wangycgroup/public/01_Pipeline/microbiomePipe/est_abundance_240613.o" >> config.ini
	echo "filter_threshold=0" >> config.ini
	echo "krakenDb=/public/home/wangycgroup/public/00_Genome_ref/genomesForPipeline/mouse-gut/kraken" >> config.ini
    /public/home/wangycgroup/public/software/bin/python /public/home/wangycgroup/public/01_Pipeline/microbiomePipe/scMeta_xzy.py --cfg config.ini
    echo "# Done: $(date)"
    echo ""
}

# Step 2: Mapping and quantification
process_mapping() {
    local sample=$1
    cd /public/home/wangycgroup/xinlong/mouse_colon/03.Add_cecum/
    # Generate species list
    tail -n +2 ${sample}/Result/${sample}_sc_taxonomy.report | awk -F '\t' '$7 >= 0 {print $2}' | sort -n | uniq > ${sample}/species.txt
    
    # Generate genome index
    mkdir -p ${sample}/ref
    touch ${sample}/ref/genome.fna
    touch ${sample}/ref/genome.gtf
    
    /public/home/wangycgroup/public/software/bin/python /public/home/wangycgroup/xinlong/mouse_colon/02.M1/build.py \
        /public/home/wangycgroup/xinlong/mouse_colon/02.M1/genome_map.tsv \
        ${sample}/species.txt ${sample}/acc.txt \
        /public/home/wangycgroup/xinlong/mouse_colon/genomes \
        ${sample}/ref/genome.fna ${sample}/ref/genome.gtf

    # Generate STAR index
    /public/home/wangycgroup/public/software/bin/STAR --runThreadN $THREADS --runMode genomeGenerate \
        --genomeDir ${sample}/ref/STAR_index \
        --genomeFastaFiles ${sample}/ref/genome.fna \
        --sjdbGTFfile ${sample}/ref/genome.gtf \
        --sjdbOverhang 122 --sjdbGTFfeatureExon gene \
        --sjdbGTFtagExonParentTranscript ID \
        --sjdbGTFtagExonParentGene ID \
        --genomeSAindexNbases 10 --limitGenomeGenerateRAM 495864700000 --limitSjdbInsertNsj 10000000
    
    # Run STAR Solo
    for file in ${sample}/clean_data/*; do
    # Check if the file is a .gz file
		if [[ "$file" == *.gz ]]; then
    # Decompress the file
		echo "Decompressing $file..."
		/public/home/wangycgroup/public/software/pigz-2.7/unpigz "$file"
		fi
	done
	
    /public/home/wangycgroup/public/software/bin/STAR --runThreadN $THREADS --soloType CB_UMI_Simple \
        --soloCBwhitelist None --genomeDir ${sample}/ref/STAR_index \
        --soloCBstart 1 --soloCBlen 20 --soloUMIstart 21 --soloUMIlen 8 \
        --readFilesIn ${sample}/clean_data/${sample}_2.fq ${sample}/clean_data/${sample}_1.fq \
        --outSAMtype BAM Unsorted --outMultimapperOrder Random --runRNGseed 1 \
        --outSAMattributes NH HI AS CR UR GX GN --alignSJoverhangMin 1000\
        --soloFeatures GeneFull --outFileNamePrefix ${sample}/${sample}_combine_name_ \
        --soloStrand Reverse --limitBAMsortRAM 41143265264 --soloCellFilter TopCells 15000 \
        --outFilterScoreMinOverLread 0.5
    
    # Downstream processing
    mkdir -p $sample/${sample}_Solo.out/GeneFull/rawSorted
    /public/home/wangycgroup/public/software/filtermtx.o \
        -i ${sample}/${sample}_Solo.out/GeneFull/raw \
        -o ${sample}/${sample}_Solo.out/GeneFull/rawSorted
    
}

# Main process
for sample in "${samples[@]}"; do
    # First stage: Sequence filtering
    process_sample "$sample"
    
    # Second stage: Mapping and quantification
    process_mapping "$sample"
    
    # Clean temporary files
    rm ${sample}/clean_data/${sample}.dedup.fq
    rm ${sample}/clean_data/${sample}.prededup.fq
    rm ${sample}/clean_data/${sample}_final_1.fq
    rm ${sample}/clean_data/${sample}_final_2.fq
    /public/home/wangycgroup/public/software/pigz-2.7/pigz -p 8 ${sample}/clean_data/${sample}_pre_1.fq
    /public/home/wangycgroup/public/software/pigz-2.7/pigz -p 8 ${sample}/clean_data/${sample}_pre_2.fq
    rm ${sample}/Result/output_classified.fq
    rm ${sample}/Result/output_unclassified.fq
    rm ${sample}/Result/${sample}_kraken.output
    /public/home/wangycgroup/public/software/pigz-2.7/pigz -p 8 ${sample}/clean_data/${sample}_1.fq
    /public/home/wangycgroup/public/software/pigz-2.7/pigz -p 8 ${sample}/clean_data/${sample}_2.fq
    /public/home/wangycgroup/public/software/pigz-2.7/pigz /public/home/wangycgroup/xinlong/mouse_colon/03.Add_cecum/${sample}/${sample}_combine_name_Solo.out/GeneFull/filtered/*
    /public/home/wangycgroup/public/software/pigz-2.7/pigz /public/home/wangycgroup/xinlong/mouse_colon/03.Add_cecum/${sample}/${sample}_combine_name_Solo.out/GeneFull/raw/*  
done

echo "All tasks completed!"

