#
# Generate alignments with hisat2
#

# A root to derive output default names from.
SRR ?= SRR1553425

# Number of CPUS
NCPU ?= 4

# Additional flags to pass to HISAT.
#HISAT2_FLAGS ?= --sensitive --rg-id ${ID} --rg SM:${SM} --rg LB:${LB} --rg PL:${PL}

# Sam filter flags to filter the BAM file before sorting.
#SAM_FLAGS ?=

# Default alignment mode is single end.
MODE ?= SE

# FASTQ read pair.
R1 ?= reads/${SRR}_1.fastq
R2 ?= reads/${SRR}_2.fastq

# The reference genome.
REF ?= refs/AF086833.fa

# The indexed reference genome.
IDX ?= $(dir ${REF})/idx/$(notdir ${REF})

# A file in the index directory.
IDX_FILE ?= ${IDX}.1.ht2

# The alignment file.
BAM ?= bam/${SRR}.hisat2.bam

# The unsorted BAM file.
BAM_TMP ?= $(basename ${BAM}).unsorted.bam


# Makefile customizations.
.RECIPEPREFIX = >
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Print usage information.
usage::
> @echo "#"
> @echo "# hisat2.mk: aligns reads using HISAT2"
> @echo "#"
> @echo "# make index align IDX=${IDX}"
> @echo "#"

# Build the index for the reference genome.
${IDX_FILE}:
> mkdir -p $(dir $@)
> hisat2-build --threads ${NCPU} ${REF} ${IDX}

# Create the index.
index: ${IDX_FILE}
> @ls -lh ${IDX_FILE}

# Remove the index.
#index!:
#> rm -rf ${IDX_FILE}

# We do not list the index as a dependency to avoid accidentally triggering the index build.

# Paired end alignment.
ifeq ($(MODE), PE)
${BAM}: ${R1} ${R2}
> @mkdir -p $(dir $@)
> hisat2 --threads ${NCPU} -k 1 -x ${IDX} -1 ${R1} -2 ${R2} | samtools view -Sb > ${BAM_TMP}
> samtools sort -@ ${NCPU} ${BAM_TMP} > ${BAM}
> rm -f ${BAM_TMP}
endif

# Single end alignment.
ifeq ($(MODE), SE)
${BAM}: ${R1}
> @mkdir -p $(dir $@)
> hisat2 --threads ${NCPU} -k 1 -x ${IDX} -U ${R1} | samtools view -Sb > ${BAM_TMP}
> samtools sort -@ ${NCPU} ${BAM_TMP} > ${BAM}
> rm -f ${BAM_TMP}
endif

# Create the BAM index file.
${BAM}.bai: ${BAM}
> samtools index ${BAM}

# Display the BAM file path.
align: ${BAM}.bai
> @ls -lh ${BAM}

# Remove the BAM file.
align!:
> rm -rf ${BAM} ${BAM}.bai

# Targets that are not files.
.PHONY: align install usage index



