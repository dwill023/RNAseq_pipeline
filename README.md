# RNAseq_pipeline
RNAseq pipeline for analysis of transcriptional changes.

High quality RNA should be used in library prep before sequencing. Below is what a bioanalyzer graph of eukaryotic ribosomal RNA should look like.
![Bioanalyzer](https://github.com/dwill023/RNAseq_pipeline/blob/master/RNA-seq%20Data%20Analysis_files/Image.png)

Alternatively, running the extracted RNA on an agarose gel (below) will be sufficient to assess RNA quality.
![Agarose gel](https://github.com/dwill023/RNAseq_pipeline/blob/master/RNA-seq%20Data%20Analysis_files/Image%20%5B1%5D.png)

The below schematic is an overview of the RNAseq analyses performed in our lab. 
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




