**General Information for RNA-Seq Pipeline**
This document provides a general overview of an RNA-Seq pipeline, guiding you from raw FASTQ files to quantification, differential expression analysis, and gene ontology analysis.

**Pipeline Overview:**
FASTQ files → Quality Control → Trimming → Mapping to Reference Genome/Transcriptome → Quantification → Differential Analysis → Gene Ontology.

**FASTQ Files**
FASTQ files are the output from sequencing platforms. Typically, these files are compressed as .fastq.gz. Most available tools can work directly with .fastq.gz files, so there's no need to decompress them before analysis.

**Quality Control**
FastQC is the most widely used tool for RNA-Seq quality control. It generates an HTML report that displays key metrics such as per-base sequence quality (good quality if >Q20 or Q30), per-base N content, per-sequence GC content (check for poly-G tails), short read content, and adapter contamination. Based on the FastQC report, you can adjust trimming parameters to improve read quality.

**Trimming**
Trimmomatic and Cutadapt are the most commonly used tools for trimming. In this pipeline, we use Trimmomatic to remove low-quality bases and adapter sequences. After trimming, it is recommended to rerun FastQC to verify that the trimming has improved read quality.

**Mapping**
After trimming, align the reads to a reference genome or transcriptome. For genomic alignment, popular tools include STAR, HISAT2, and TopHat2 (though TopHat2 is generally slower and less commonly used now). For transcriptomic alignment, tools like Salmon (more popular) and Kallisto are preferred.

In the provided bash script, HISAT2 is used for alignment, as it offers pre-built indexed reference genomes available for download at https://daehwankimlab.github.io/hisat2/download/. While indexing the reference genome is straightforward, it requires a good amount of time and computational resources.

**Quantification**
The final output from genome alignment tools, such as STAR or HISAT2, are BAM files. These files contain the aligned reads, and to get quantification, a tool like featureCounts is required to count the number of reads mapped to each gene or feature. Salmon and Kallisto skip the need for BAM files by directly quantifying transcript abundances from the RNA-Seq reads. They provide raw read counts and TPM values as their final output, so no additional quantification tool is needed.
https://github.com/harumiswari/RNAseq_analysis/blob/master/rnaseq.sh

----------end of bash script code --------------


----------start of R code-------------------

**Differential expression**
https://github.com/harumiswari/RNAseq_analysis/blob/master/DESvolcano_final.R
A key step in RNA-Seq data analysis, used to identify genes or transcripts that show statistically significant differences in expression levels between different conditions or groups. The two most popular R packages are DESeq2 and edgeR. The raw read count output from the bash script will be the input for this differential analysis. These tools will normalize the counts, model the data, and identify genes that are significantly differentially expressed between the conditions. In the provided code, it will also demonstrate how to generate PCA and volcano plots, as well as how to add gene annotations to the plots.

**Gene ontology analysis**
https://github.com/harumiswari/RNAseq_analysis/blob/master/geneontology.R
GO analysis is commonly performed after differential expression analysis to identify which GO terms are overrepresented among the differentially expressed genes. This helps to understand the biological significance of the gene expression changes in the context of broader biological functions, processes, and locations. Several R packages for GO analysis are topGO and clusterprofiler. In the provided code, we used clusterProfiler.  
