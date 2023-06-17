# 
# Workflow for RNA-seq alignment and counting
#

# bash strict mode
set -uex

# genome reference
REF=refs/GRCm39_genome.fa

# genome annotation file
GTF=refs/GRCm39_annotation.gtf

# Run fastqc
cat ids.txt | parallel fastqc reads/{}.fastq.gz

# generate the hisat2 index
make -f hisat2.mk index NCPU=8 REF=${REF}

# Run the alignments.
cat ids.txt | parallel --eta --verbose make -f hisat2.mk align MODE=SE SRR={} REF=${REF} BAM=bam/{}.bam R1=reads/{}.fastq.gz

# Count the features
cat ids.txt | parallel -k echo bam/{}.bam | parallel -u --xargs featureCounts -T 8 -a ${GTF} -g gene_name -o counts.txt {}

# simplify the counts
cat counts.txt | cut -f 1,6-20 > simple_counts.txt             

# run multiqc
multiqc .



