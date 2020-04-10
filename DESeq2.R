#RNA-seq workflow: 
# https://www.bioconductor.org/help/workflows/rnaseqGene/#top

# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear console
cat("\014") 
# Clean workspace
rm(list=ls())
# Set working directory
setwd()


## Using DESeq2 

# For analysis of RNA-seq data. Preparing count matrices: Using the featureCounts function (Liao, Smyth, and Shi 2013) in the Rsubread package. The DESeq2 model internally corrects for library size, so transformed or normalized values such as counts scaled by library size should not be used as input.
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(Glimma)
  library(gplots)
  library(RColorBrewer)
})


mycounts <- read.delim("simple_counts.txt", sep = "\t") # load count matrices 
metadata <-  read_csv("metadata.csv") # Information about each sample.


class(mycounts)
class(metadata)

mycounts <- as.data.frame(mycounts)
metadata <- as.data.frame(metadata)

all(names(mycounts)[-1]==metadata$id) # check that the labeling of samples is the same in the count matrix and metadata file

# Remove duplicated rows in Geneid column
mycounts <- mycounts %>% distinct(Geneid, .keep_all = TRUE)

#########  DEG Analysis ###################


# Set design = the variables you want to compare in the metadata file.
dds <- DESeqDataSetFromMatrix(countData=mycounts, 
                              colData=metadata, 
                              design=~treatment, 
                              tidy=TRUE)
dds

dds <- DESeq(dds)
# pre-filtering the dataset. Removes rows that have no or nearly no information about the amount of gene expression by removing rows that have no counts or only a single count across all samples. 
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

# To make all comparisons against a control or specific sample.
resultsNames(dds)
dds$treatment <- relevel(dds$treatment, ref = "control")
dds <- DESeq(dds)

#rlog and variance: transformations for count data that stabilize the variance across the mean: the regularized-logarithm transformation or rlog
# rlog tends to work well on small datasets (n < 30), sometimes outperforming the VST when there is a large range of sequencing depth across samples (an order of magnitude difference). 

rld <- rlog(dds, blind = FALSE) # blind = FALSE, which means that differences between cell lines and treatment (the variables in the design) will not contribute to the expected variance-mean trend of the experiment.
head(assay(rld), 3)

# to take out only specific treatments you want to compare on heatmap
# rld.sub <- rld[ , rld$treatment %in% c("UNT","VPA","CPA","MTX")] 

## Visualization

# principal components analysis (PCA)
plotPCA(rld, intgroup = c("treatment")) + geom_text(aes(label=name),vjust=2)

# Clustering and heat maps
library("genefilter")
library(pheatmap)

topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 100) #a subset of the most highly variable genes.

mat  <- assay(rld)[ topVarGenes, ]
mat  <- mat - rowMeans(mat) # The heatmap becomes more interesting if we do not look at absolute expression strength but rather at the amount by which each gene deviates in a specific sample from the gene's average across all samples. 
anno <- as.data.frame(colData(rld)[, c("id","treatment")])
pheatmap(mat, annotation_col = anno, main = "Top 100 most variable genes")

# Box-plots
# The DESeq function calculates, for every gene and for every sample, a
# diagnostic test for outliers called Cook's distance. Cook's distance is a
# measure of how much a single sample is influencing the fitted coefficients for
# a gene, and a large value of Cook's distance is intended to indicate an
# outlier count. The Cook's distances are stored as a matrix available in
# assays(dds)[["cooks"]].

par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2, ylab="Log10 Cook's distance", main= "Boxplots to access outliers")


#comparing results between two conditions must do this for results!~~~~~~~~~~~~~~~~~~~~~~~~~
# results for a comparison of any two levels of a variable can be extracted using the contrast argument to results.
#Raise the log2 fold change threshold by supplying a value on the log2 scale
#For example, by specifying lfcThreshold = 1, we test for genes that show significant effects of treatment on gene 
#counts more than doubling or less than halving, because 2^1=2
#  for example, a log2 fold change of 1.5 means that the gene's expression is increased by a multiplicative factor of 2^1.5 ??? 2.82.

resultsNames(dds) #shows us whats in the dataframe to know what to compare

#Getting results
res_2D <- results(dds, contrast=c("treatment","2D_T75","control")) #  lfcThreshold=1. The name of the variable, the name of the factor level for the numerator of the log2 ratio, and the name of the factor level for the denominator.
res_60 <- results(dds, contrast=c("treatment","60RPM","control"))
res_0 <- results(dds, contrast=c("treatment","Stat_0RPM","control"))

mcols(res_2D, use.names = TRUE)
summary(res_60)
sum(res_2D$padj < 0.05, na.rm=TRUE)

# subset the results table to significant genes you put thrreshholds on above. We subset the results table to these genes and 
# then sort it by the log2 fold change estimate to get the significant genes with the strongest down-regulation:

resSig1 <- subset(res_2D, padj < 0.05)
ressig2 <- subset(res_60, padj < 0.05)
ressig3 <- subset(res_0, padj < 0.05)

summary(ressig2)
sigres3 <- as.data.frame(ressig3)
write.csv(sigres3, file = "0RPM sig_DEG.csv")

# Alternatively, you can filter the result for both padj < 0.05 and lfc > 1
log2FCcutoff <- 1
BHcutoff <- 0.05
sigGeneList <- subset(res, abs(log2FoldChange)>=log2FCcutoff & padj<=BHcutoff)


# MA-Plot: shows the log2 fold changes attributable to a given variable over the mean of normalized counts
# Points will be colored red if the adjusted p value is less than 0.1. 
# Points which fall out of the window are plotted as open triangles pointing either up or down.
# It is more useful visualize the MA-plot for the shrunken log2 fold changes, which remove the noise 
# associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.

plotMA(res_60,ylim=c(-3,3), cex=.8, alpha = 0.05, main="60 RPM vs Control")
abline(h=c(-1,1), col="dodgerblue", lwd=2)

# Make a awesome volcano plots
library(EnhancedVolcano)
# The default cut-off for log2FC is >|2|; the default cut-off for P value is 10e-6.
EnhancedVolcano(res_0,
                lab = rownames(res_0),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-5, 8),
                subtitle = "0 RPM",selectLab = c('JUN'))

# use selectLab = c('A2m') in above volcano plot to get rid of labels

## Gene ontology Analysis with ClusterProfiler 
# The over-representation analysis (enrichGO) looks at whether a subset of genes that you have already separated
#out associate significantly with certain pathways. The gene set enrichment
#analysis (gseGO) takes differential data from every measured gene and looks for
#pathways displaying significantly co-ordinated shifts in those values.

library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

# use sigres data frame for GO
# sigres1 = 2D_T75
# sigres2 = 60 RPM
# sigres3 = 0 RPM


# we want the log2 fold change 
original_gene_list1 <- sigres1$log2FoldChange
original_gene_list2 <- sigres2$log2FoldChange
original_gene_list3 <- sigres3$log2FoldChange
# name the vector
names(original_gene_list1) <- row.names(sigres1)
names(original_gene_list2) <- row.names(sigres2)
names(original_gene_list3) <- row.names(sigres3)
# omit any NA values 
gene_list1<-na.omit(original_gene_list1)
gene_list2<-na.omit(original_gene_list2)
gene_list3<-na.omit(original_gene_list3)
# sort the list in decreasing order (required for clusterProfiler)
gene_list1 = sort(gene_list1, decreasing = TRUE)
gene_list2 = sort(gene_list2, decreasing = TRUE)
gene_list3 = sort(gene_list3, decreasing = TRUE)
# filter on min log2fold change (log2FoldChange > 1)
genes1 <- names(gene_list1)[abs(gene_list1) > 0]
genes2 <- names(gene_list2)[abs(gene_list2) > 0]
genes3 <- names(gene_list3)[abs(gene_list3) > 0]

mydf <- data.frame(Entrez=names(gene_list3), FC=gene_list3)
mydf <- mydf[abs(mydf$FC) > 0,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"


formula_res <- compareCluster(Entrez~group, data=mydf, fun="enrichGO", ont= "BP", OrgDb= org.Hs.eg.db, keyType = "SYMBOL")

head(as.data.frame(formula_res))
dotplot(formula_res, showCategory = 25, font.size = 16, title = "Zero RPM Enriched Ontologies")

# Get GO list
write.csv(as.data.frame(formula_res), file ="0RPM EGO list.csv")

# Over-representation test
ego <- enrichGO(gene           = genes3,
                OrgDb         = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

dotplot(ego, showCategory = 25, font.size = 18, title = "Gene Enrichment Analysis")
emapplot(ego, showCategory = 25, font.size = 30)
barplot(ego, showCategory = 25)

# Gene set enrichment analysis
gse <- gseGO(geneList=gene_list3, 
             ont ="BP", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")

head(gse)
# visualization
dotplot(gse, showCategory = 25, title = "Gene set enrichment") 

# Ridgeplot: plots the frequency of fold change values per gene within each set. Helps to view up/down-regulated pathways. 
ridgeplot(gse) + labs(x = "enrichment distribution")

# Get GSE list
write.csv(as.data.frame(gse), file ="0 RPM GSE list.csv")

sessionInfo()
