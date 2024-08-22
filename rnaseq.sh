#!/bin/bash

SECONDS=0

# Set variables
DATA="/home/cnelsonlab/Downloads/RawData/OMRF2035_RNAseq"
TRIMMED_DATA="/home/cnelsonlab/Downloads/RawData/trimmed"
ALIGNED_DATA="/home/cnelsonlab/Downloads/RawData/Mouse_grcm_aligned_RNAseq"
HISAT2_INDEX="/home/cnelsonlab/rnaseq/HISAT2/grcm38/genome"
LOG_FILE="alignment_log.txt"

# Generate an array of samples
SAMPLES=()
for i in {1..17}; do
    SAMPLES+=("Sample_$i")
done

# Create directories
mkdir -p "$TRIMMED_DATA" "$ALIGNED_DATA"

# Run FastQC
for sample in "${SAMPLES[@]}"; do
    fastqc "$DATA/${sample}_R1_001.fastq.gz" "$DATA/${sample}_R2_001.fastq.gz" -o "$TRIMMED_DATA/"
    echo "FastQC done for $sample!"
done
echo "FastQC processing completed for all samples!"

# Run Trimmomatic
for sample in "${SAMPLES[@]}"; do
    R1_FILE=$(ls "$DATA/${sample}_R1_001.fastq.gz")
    R2_FILE=$(ls "$DATA/${sample}_R2_001.fastq.gz")
    trimmomatic PE -phred33 \
        "$R1_FILE" "$R2_FILE" \
        "$TRIMMED_DATA/${sample}_paired1.fq.gz" "$TRIMMED_DATA/${sample}_unpaired1.fq.gz" \
        "$TRIMMED_DATA/${sample}_paired2.fq.gz" "$TRIMMED_DATA/${sample}_unpaired2.fq.gz" \
        TRAILING:10
    echo "Trimmomatic is done for $sample!"
done
echo "Trimmomatic processing completed for all samples!"

# Run HISAT2, samtools, and featureCounts for each sample
for sample in "${SAMPLES[@]}"; do
    echo "Processing sample: $sample"

    # Run HISAT2 and samtools
    hisat2 -q --rna-strandness R -x "$HISAT2_INDEX" \
           -1 "${TRIMMED_DATA}/${sample}_paired1.fq.gz" \
           -2 "${TRIMMED_DATA}/${sample}_paired2.fq.gz" | \
           samtools sort -o "${ALIGNED_DATA}/${sample}.bam" || { echo "Error in HISAT2 for $sample"; exit 1; }

    # Index the BAM file
    samtools index "${ALIGNED_DATA}/${sample}.bam"
    
    # Run featureCounts
    featureCounts -p -s 2 \
      -a /home/cnelsonlab/rnaseq/HISAT2/Mus_musculus.GRCm39.110.gtf \
      -o "${ALIGNED_DATA}/${sample}_readcount.txt" \
      "${ALIGNED_DATA}/${sample}.bam" || { echo "Error in featureCounts for $sample"; exit 1; }

    echo "Alignment, indexing, and feature counting for $sample completed" | tee -a "$LOG_FILE"
done

echo "All samples processed in $SECONDS seconds"
