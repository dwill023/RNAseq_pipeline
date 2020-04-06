# RNAseq_pipeline
RNAseq pipeline for analysis of transcriptional changes.

High quality RNA should be used in library prep before sequencing. Below is what a bioanalyzer graph of eukaryotic ribosomal RNA should look like.

![Bioanalyzer](https://github.com/dwill023/RNAseq_pipeline/blob/master/RNA-seq%20Data%20Analysis_files/Image.png)

Alternatively, running the extracted RNA on an agarose gel (below) will be sufficient to assess RNA quality.
![Agarose gel](https://github.com/dwill023/RNAseq_pipeline/blob/master/RNA-seq%20Data%20Analysis_files/Image%20%5B1%5D.png)

Below is a schematic of a typical RNAseq analyses. 
![schematic](https://github.com/dwill023/RNAseq_pipeline/blob/master/RNA-seq%20Data%20Analysis_files/Untitled%20Diagram.jpg)

The majority of the analysis (with the exception of differential expression) is performed in UCR's high performance commputing cluster (HPCC) which contains the software indicated.

Once samples have been sequenced at a core facility, the reads are automatically trimmed of adapter sequences. However, low quality sequcences may exist (usually at the 3'end) which should be removed.

## Trimming low quality bases
Using the package [Trimmomatic](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) in the HPCC as shown below. 
```Shell
module load trimmomatic/0.36

trimmomatic SE -threads 4 <input fastq> <output fastq> LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:36
```
This will perform the following in order:
- Remove leading low quality or N bases (below quality 3)
- Remove trailing low quality or N bases (below quality 3)
- Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 25
- Drop reads which are less than 36 bases long after these steps

After trimming the low quality reads, run fastQC to check the quality stats.
```Shell
module load fastqc/0.11.7

fastqc <trimmed.fastq> 
```
This will output the quality stats into an html file. Check that the low quality bases, read length and 'N' sequences have been removed.

## Align reads to the genome
Use a short read aligner such as [HISAT2](https://ccb.jhu.edu/software/hisat2/manual.shtml), which has similar algorithms to other aligners like Bowtie2 and Tophat2, but outperfoms them at much faster speeds.

Download the Genome sequence, primary assembly (fasta file) for:
Human: https://www.gencodegenes.org/human/
Mouse: https://www.gencodegenes.org/mouse/

First an index of the genome you plan to align against must be created.
```Shell
module load hisat2/2.1.0
hisat2-build <genome_file.fa> <path to where you want your index stored>
```
Run the alignment.
hisat2 [options] -x <hisat2 index> -U <trimmed.fastq> -S <name of outfile> -k 1 --new-summary <name of summary file>
- -k: searches for at most 1 distinct, primary alinment for each read. The default is 5.
```Shell
hisat2 -p 10 -x <hisat2 index> -U <trimmed.fastq> -S <name of outfile> -k 1 --new-summary <name of summary file>
```
For multiple files use a shell script (hisat2_align.sh) to run the multiple fastq files in succession.

```Shell
# edit the hisat2_align.sh file with the appropriate index path and fastq file paths.
bash hisat2_align.sh
```

Run samtools to convert to bam format and sort the bam file
```Shell
module load parallel
module load samtools

ls *.sam | parallel --eta --verbose "samtools view -@ 4 -Sb {} | samtools sort > {.}.sorted.bam"
```

## Generating counts
This estimates the abundance of the reads that align over an interval. This is called read quantification and is required for gene expression analysis. The output is a count table that gives the number of reads within a gene feature. This is done with the tool featureCounts within the [subread package](http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf)

This will require a gtf file which can be downloaded:
Human: https://www.gencodegenes.org/human/
Mouse: https://www.gencodegenes.org/mouse/

```Shell
module load subread/1.6.2
GTF= <path to .gtf file>
# -T specifies the number of cores, -g to specify list generated with gene names instead of the default IDs.

featureCounts -T 4 -a $GTF -g gene_name -o counts.txt <path to *.bam>

# Simplify the counts to have only gene names and counts
cat counts.txt | cut -f 1,6-20 > simple_counts.txt
```
#### Getting RPKM from raw counts
In R
```R
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(DESeq2)
  library(edgeR)
})

# Getting the RPKM from the raw read counts
counts <- read_delim("simplecounts.txt", delim = "\t") # load count matrices

counts <- as.data.frame(counts)
# Remove duplicated rows in Geneid column
counts <- counts %>% distinct(Geneid, .keep_all = TRUE)

# remove the first column from counts
countdata <- counts[,-(1)]

# Store Geneid as rownames
rownames(countdata) <- counts[,1]
genes <- countdata[,-(1)]
head(countdata)
colnames(genes)

# obtain RPKM
RPKM <- rpkm(genes, gene.length = countdata$Length)

# Filter the RPKM for lowly expressed microRNAs.
# which RPKMs are > 0.1. A threshold that keeps genes with a raw count ~1
thresh <- RPKM > 0.1
head(thresh)

# Summary of how many TRUEs are in each row
table(rowSums(thresh))
# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2

# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- RPKM[keep,]
summary(keep)
dim(counts.keep)

# export the logcounts into a file
write.csv(counts.keep, "RPKM.csv")
```
## Differential gene exression (DEG) analysis
Using R package Deseq2 and following script in file DESeq2.R.

Before beginning, create a metadata table in .csv format like the one below.

| id            | treatment |
|     :---:     |   :---:   |
| control1      | control   |
| control2      | control   |
| control3      | control   |
| Experiment1   | experiment|
| Experiment2   | experiment|
| Experiment3   | experiment|

- PCA plot (Principal component analysis): A form of dimension reduction that converts the correlations (or lack of) among all samples an plots them onto a 2D plot. Samples that are highly correlated cluster together. Differences along the first axis (PC1) are more important than differences along the PC2 axis. A measure of the variance between samples. Heatmaps generated using the rlog transformation for samples < 30. For samples > 30 use the VST normalization. 
- MA-plot: plot visualizes the gene expression differences between two samples. 
- Volcano plot: Generates visualization similar to MA plot but with more details (pvalues, gene names.)
- Gene Set Enrichment Analysis with `ClusterProfiler`

For the results if you see p-value & p.adjust values of NA it indicates outliers detected by Cook’s distance. NA only for p.adjust means the gene is filtered by automatic independent filtering for having a low mean normalized count. 
