#!/bin/bash
#SBATCH -J p1
#SBATCH -p fat
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/public/home/wangycgroup/public/panhuaran/htslib:/public/home/wangycgroup/public/software/lib2/*

# Define sample array
samples=("241224-XGR-D71JC" "241224-XGR-D71MC" "241224-XGR-D72JC" "241224-XGR-D72MC" "241224-XGR-W71JC" "241224-XGR-W71MC" "241224-XGR-W72JC" "241224-XGR-W72MC" "41217N01" "41217N02" "41217N03" "41217N04")
THREADS=$SLURM_CPUS_PER_TASK
p=`pwd`
# Stage1 : Sequence filtering processing
process_sample() {
    local sample=$1
    if [ ! -d "$sample" ]; then
        mkdir $p/$sample
        mkdir $p/$sample/clean_data

        for file in /public/home/wangycgroup/xinlong/mouse_colon/03.Add_cecum/${sample}/clean_data/*; do
            # Check if the file is a .gz file
            if [[ "$file" == *.gz ]]; then
                # Decompress the file
                echo "Decompressing $file..."
                /public/home/wangycgroup/public/software/pigz-2.7/unpigz -p 8 "$file"
            fi
        done
        ln -s /public/home/wangycgroup/xinlong/mouse_colon/03.Add_cecum/$sample/clean_data/${sample}_1.fq ./${sample}/clean_data/
        ln -s /public/home/wangycgroup/xinlong/mouse_colon/03.Add_cecum/$sample/clean_data/${sample}_2.fq ./${sample}/clean_data/
    fi
    cd ${sample}
    # Run scMeta
    echo "# Start scMeta for ${sample}: $(date)"
    touch config.ini
    echo "sample=${sample}" > config.ini
    echo "outdir=/public/home/wangycgroup/xinlong/mouse_colon/05.Refseq/${sample}" >> config.ini
    echo "process=$THREADS" >> config.ini
    echo "nubeam_dedup=/public/home/wangycgroup/public/software/nubeamdedup-master/Linux/nubeam-dedup" >> config.ini
    echo "kraken=/public/home/wangycgroup/public/software/kraken2/kraken2" >> config.ini
    echo "braken=/public/home/wangycgroup/public/01_Pipeline/microbiomePipe/est_abundance_240613.o" >> config.ini
    echo "filter_threshold=0" >> config.ini
    echo "krakenDb=/public/home/wangycgroup/public/Database/Microbiome/kraken2" >> config.ini
    /public/home/wangycgroup/public/software/bin/python /public/home/wangycgroup/public/01_Pipeline/microbiomePipe/scMeta_xzy.py --cfg config.ini
    echo "# Done: $(date)"
    echo ""
}

# Stage2 : Mapping and quantification processing
process_mapping() {
    local sample=$1
    cd /public/home/wangycgroup/xinlong/mouse_colon/05.Refseq/
    # Generate species list
    tail -n +2 $sample/Result/${sample}_sc_taxonomy.report | cut -f 3 | sort -n | uniq > $sample/species.txt
    # Generate genome index
    mkdir -p ${sample}/ref
    touch ${sample}/ref/genome.fna
    touch ${sample}/ref/genome.gtf
    
    echo '# Start extracting genome info for ${sample}: '`date`
	/public/home/wangycgroup/public/software/bin/python /public/home/wangycgroup/public/01_Pipeline/microbiomePipe/selectGenomes.py -i ${sample}/species.txt -r /public/home/wangycgroup/public/00_Genome_ref/genomesForPipeline/Bacteria/refseq_mapping.tsv -d /public/home/wangycgroup/public/00_Genome_ref/genomesForPipeline/Bacteria/raw -o ${sample}/ref/genome
	echo '# Done:                                  '`date`
    
    # Generate STAR index
    /public/home/wangycgroup/public/software/bin/STAR --runThreadN $THREADS --runMode genomeGenerate \
        --genomeDir ${sample}/ref/STAR_index \
        --genomeFastaFiles ${sample}/ref/genome.fna \
        --sjdbGTFfile ${sample}/ref/genome.gtf \
        --sjdbOverhang 122 --sjdbGTFfeatureExon gene \
        --sjdbGTFtagExonParentTranscript gene_id \
        --sjdbGTFtagExonParentGene gene_id \
        --genomeSAindexNbases 10 --limitGenomeGenerateRAM 495864700000
    
    # Run STAR Solo
    /public/home/wangycgroup/public/software/bin/STAR --runThreadN $THREADS --soloType CB_UMI_Simple \
        --soloCBwhitelist None --genomeDir ${sample}/ref/STAR_index \
        --soloCBstart 1 --soloCBlen 20 --soloUMIstart 21 --soloUMIlen 8 \
        --readFilesIn ${sample}/clean_data/${sample}_2.fq ${sample}/clean_data/${sample}_1.fq \
        --outSAMtype BAM Unsorted --outMultimapperOrder Random --runRNGseed 1 \
        --outSAMattributes NH HI AS CR UR GX GN --alignSJoverhangMin 1000 \
        --soloFeatures GeneFull --outFileNamePrefix ${sample}/${sample}_ \
        --soloStrand Reverse --limitBAMsortRAM 41143265264 --soloCellFilter TopCells 15000 \
        --outFilterScoreMinOverLread 0.5
    
    # Downstream processing
    mkdir -p $sample/${sample}_Solo.out/GeneFull/rawSorted
    /public/home/wangycgroup/public/software/filtermtx.o \
        -i ${sample}/${sample}_Solo.out/GeneFull/raw \
        -o ${sample}/${sample}_Solo.out/GeneFull/rawSorted
    
    /public/home/wangycgroup/public/software/bin/python \
        /public/home/wangycgroup/public/01_Pipeline/starSoloPipe-preindex/plotScatter-deliver.py ${sample}
}

    # Main process
for sample in "${samples[@]}"; do
    # Stage1 : Sequence filtering
    process_sample "$sample"
    
    # Stage2 : Mapping and quantification
    process_mapping "$sample"
    
    # Clean temporary files
    rm ${sample}/clean_data/${sample}.dedup.fq
    rm ${sample}/clean_data/${sample}.prededup.fq
    rm ${sample}/clean_data/${sample}_final_1.fq
    rm ${sample}/clean_data/${sample}_final_2.fq
    rm ${sample}/Result/output_classified.fq
    rm ${sample}/Result/output_unclassified.fq
    rm ${sample}/Result/${sample}_kraken.output
    /public/home/wangycgroup/public/software/pigz-2.7/pigz -p 8 /public/home/wangycgroup/xinlong/mouse_colon/03.Add_cecum/${sample}/clean_data/${sample}_1.fq
    /public/home/wangycgroup/public/software/pigz-2.7/pigz -p 8 /public/home/wangycgroup/xinlong/mouse_colon/03.Add_cecum/${sample}/clean_data/${sample}_2.fq
    /public/home/wangycgroup/public/software/pigz-2.7/pigz -p 8 /public/home/wangycgroup/xinlong/mouse_colon/03.Add_cecum/${sample}/clean_data/*pre_*.fq
    /public/home/wangycgroup/public/software/pigz-2.7/pigz -p 8 /public/home/wangycgroup/xinlong/mouse_colon/05.Refseq/${sample}/${sample}_Solo.out/GeneFull/filtered/*
    /public/home/wangycgroup/public/software/pigz-2.7/pigz -p 8 /public/home/wangycgroup/xinlong/mouse_colon/05.Refseq/${sample}/${sample}_Solo.out/GeneFull/raw/*  
done

echo "All tasks completed!"

