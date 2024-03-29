# RNAseq_pipeline
RNAseq pipeline for analysis of transcriptional changes.

High quality RNA should be used in library prep before sequencing. Below is what a bioanalyzer graph of eukaryotic ribosomal RNA should look like.

![Bioanalyzer](https://github.com/dwill023/RNAseq_pipeline/blob/master/RNA-seq%20Data%20Analysis_files/Image.png)

Alternatively, running the extracted RNA on an agarose gel (below) will be sufficient to assess RNA quality.
![Agarose gel](https://github.com/dwill023/RNAseq_pipeline/blob/master/RNA-seq%20Data%20Analysis_files/Image%20%5B1%5D.png)

Below is a schematic of a typical RNAseq analyses. 
![schematic](https://github.com/dwill023/RNAseq_pipeline/blob/master/RNA-seq%20Data%20Analysis_files/Untitled%20Diagram.jpg)

# Workflow Option

An automatic way to process the raw reads from quality analysis to read counting can be performed using a [shell script](https://github.com/dwill023/RNAseq_pipeline/blob/master/rna-seq-hisat2-workflow.sh) that can take a list of your raw files and process them in parallel. The script takes a ids.txt file containing the root names of the raw sequencing file. For example if your files are sample_1.fastq.gz, sample_2.fastq.gz and so on, the ids.txt file should have each file name on a new line like so:

```
sample_1
sample_2
```
The script then takes the ids.txt file and runs the following steps:

1. Generate read quality files using fastqc.
2. Read alignment with HISAT2 using a separate makefile [hisat2.mk](https://github.com/dwill023/RNAseq_pipeline/blob/master/hisat2.mk) that contains parameters to align the file. 
  - NCPU : The number of CPUs to use, default is 4
  - MODE : SE or PE for single-end or paired-end reads
  - REF : Path to the genomic reference (.fa) file in a directory called refs/
  - IDX : Path to the indexed genome file stored in a directory the refs/idx, if available. If not the script will automatically create an index.
3.  Generates counts of the reads using Subread's featurecounts function.
4. Runs MultiQC to output the metrics of the read quality, alignment and read counts.

Before you run the shell script make sure you have:
1. The raw reads in a directory called reads/
2. A reference genome and a genome annotation file (.gtf) in a directory called refs/
3. The ids.txt containing the names of the raw sequencing files.
4. Amend the script with the name of your genome and annotation files on:
  - line 9 : REF=refs/your_genome_file.fa
  - line 12 : GTF=refs/your_annotation_file.gtf
5. If you have paired end reads amend the MODE with PE in the script at line 21.   

Run the script in the same directory where the hisat2.mk file is.
```
bash rna-seq-hisat2-workflow.sh
```

After the counts file is generated follow the section [Differential gene exression (DEG) analysis](https://github.com/dwill023/RNAseq_pipeline/tree/master#differential-gene-exression-deg-analysis).


# Step-by-Step Option

The below outlines the individual steps taken in the workflow. These steps go into detail about what is occurring in the automatic workflow and can also be followed to get the same results.

## Trimming low quality bases (optional, not part of workflow)
If trimming is needed the package [Trimmomatic](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) can be run as shown below:
```Shell
module load trimmomatic/0.36

trimmomatic SE -threads 4 <input fastq> <output fastq> LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
This will perform the following in order:
- Remove leading low quality or N bases (below quality 3)
- Remove trailing low quality or N bases (below quality 3)
- Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 25
- Drop reads which are less than 36 bases long after these steps

After trimming the low quality reads, re-run MultiQC to check the quality stats.

This will output the quality stats into an html file. Check that the low quality bases, read length and 'N' sequences have been removed.

## Align reads to the genome
Use a short read aligner such as [HISAT2](https://ccb.jhu.edu/software/hisat2/manual.shtml), which has similar algorithms to other aligners like Bowtie2 and Tophat2, but outperfoms them at much faster speeds.

Download the Genome sequence, primary assembly (fasta file) for:
Human: https://www.gencodegenes.org/human/
Mouse: https://www.gencodegenes.org/mouse/

First an index of the genome you plan to align against must be created.
```Shell
module load hisat2/2.1.0
hisat2-build <genome_file.fa> <basename of the index file>
```
Run the alignment.
hisat2 [options] -k 1 -x <hisat2 index basename> -U <trimmed.fastq> -S <name of outfile.sam> --summary-file <name of summary file.txt>
- -k: searches for at most 1 distinct, primary alinment for each read. The default is 5.
```Shell
hisat2 -p 10 -k 1 -x hg38 -U control.fastq.gz -S control.sam --summary-file control.txt
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
Normalizing counts by RPKM/CPM/TPM is good for comparing genes within a sample. It is not good for comparing **between** samples. This is because these normalizations assume RNA abundance and distributions are similar across compared samples. See this [article](https://rnajournal.cshlp.org/content/early/2020/04/13/rna.074922.120) discussing misuse of RPKM/TPM.

The below code will provide RPKM normalization if needed.
In R, load count file which must have a column for the gene lengths.
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
#### To get normalized counts
To obtain normalized counts for comparisons between samples, DESeq2 has a function that normalizes counts divided by sample-specific size factors. These size factors are determined by median ratio of gene counts relative to geometric mean per gene. See the [DESeq2.R](https://github.com/dwill023/RNAseq_pipeline/blob/master/DESeq2.R) file.


# Differential gene exression (DEG) analysis
Using R package Deseq2 and following script in [DESeq2.R](https://github.com/dwill023/RNAseq_pipeline/blob/master/DESeq2.R).

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
