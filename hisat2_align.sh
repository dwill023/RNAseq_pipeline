#!/bin/bash

# Immediately stop on errors
 
 set -uex pipefail


# The reference genome path
REF=Homo_sapiens.GRCh38.fa 


module load hisat2/2.1.0

hisat2 -p 8 -k 1 -x $REF -U /rawdata/3_S13_R1_001.fastq.gz -S /aligned/D4_WT.sam --summary-file /aligned/D4_WT.txt

hisat2 -p 8 -k 1 -x $REF -U /rawdata/6_S14_R1_001.fastq.gz -S /aligned/D4_mir690s.sam --summary-file /aligned/D4_mir690s.txt

hisat2 -p 8 -k 1 -x $REF -U /rawdata/13_S1_R1_001.fastq.gz -S /aligned/D3_1.sam --summary-file /aligned/D3_1.txt
hisat2 -p 8 -k 1 -x $REF -U /rawdata/14_S2_R1_001.fastq.gz -S /aligned/D3_2.sam --summary-file /aligned/D3_2.txt
hisat2 -p 8 -k 1 -x $REF -U /rawdata/15_S3_R1_001.fastq.gz -S /aligned/D3_3.sam --summary-file /aligned/D3_3.txt
