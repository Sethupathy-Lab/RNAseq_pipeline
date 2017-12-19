###################################################################################################
#
#     Prelude: You only need to run this code the first time
#
###################################################################################################

source("https://bioconductor.org/biocLite.R")
biocLite("tximport")
biocLite("readr")
biocLite("tximportData")
biocLite("biomaRt")
biocLite("GenomicFeatures")
biocLite("DESeq2")
install.packages("ggplot2")
install.packages("gplots")
install.packages("ggrepel")


###################################################################################################
#
#     Set up working environment and import data
#
###################################################################################################

## Load libraries (install if necessary)
library(tximport)
library(readr)
library(GenomicFeatures)
library(biomaRt)
library(DESeq2)
library(ggplot2)
library(gplots)
library(ggrepel)

## Change to location of results directory and samples.csv file
setwd("/Users/Matthias/Dropbox/praveen/sequencing_projects/RNAseq_tutorial/")
base_dir = getwd()

# Create directory for output results
dir.create(file.path(getwd(), 'DESeq_output'), showWarnings = FALSE)

## Import sample and condition file
samples = read.csv(file.path(base_dir, "samples.csv"), header = TRUE, stringsAsFactors=FALSE)
samples$condition <- factor(samples$condition)
samples$condition <- relevel()
samples        # Prints the sample / condition list

## Make a TxDb object from transcript annotations
## - available as a GFF3 or GTF file
## - can be downloaded from gencode for mouse / human, ensembl for other species
## - this dataset was mapped against mm10

gtf="gencode.vM10.annotation.gtf"   # Download this file from here; ftp://cbsuftp.tc.cornell.edu/pr46ftp/tutorial_files/DESeq/
                                    # then move to the working directory (where the script is)
txdb=makeTxDbFromGFF(gtf,
                     format="gtf",
                     organism="Mus musculus",
                     taxonomyId="10090")

k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1]
head(tx2gene)

## Import Salmon quant files and create counts table
files <- file.path(base_dir, "results", samples$sample, "quant.sf")
all(file.exists(files))        # Verify names of files match names in samples.csv, should return True
names(files)=samples$sample
txi <- tximport(files, 
                type = "salmon", 
                tx2gene = tx2gene)
head(txi$counts)               # This is the counts table for all of our samples

## Now to import the data into a DESeq Data Set (dds)
## Verify that sample names and colnames are the same
identical(samples$sample,colnames(txi$counts))

## Create a DEseqDataSet from txi count table
dds <- DESeqDataSetFromTximport(txi, samples, ~condition)

###################################################################################################
#
#    EXPLORATORY DATA ANALYSIS
#
###################################################################################################

library(tidyverse)
library("RColorBrewer")
library(pheatmap)

## Set color palette for figures
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

## Perform a rlog transformation on count data (essentially a puts on a log2 scale)
## This helps our data assume a normal distribution and is good to do before these analyses
rld <- rlog(dds, blind=TRUE)

## Setup annotation file to show the conditions on the figures
treat_ann <- samples
rownames(treat_ann) <- treat_ann$sample
treat_ann$sample <- NULL
treat_ann

## SAMPLE TO SAMPLE DISTANCE & CORRELATION HEATMAPS

## Sample correlation heatmap
corr_samps <- cor(as.matrix(assay(rld)))      # Computes pairwise correlations between samples based on gene expression
png(filename="DESeq_output/DESeq_sampleCorr_HM.png", units = 'in', width = 8, height = 8, res = 250)
pheatmap(corr_samps,
         annotation = treat_ann,
         col=colors,
         main="Sample Correlations")
dev.off()

# Sample distance heatmap
sampleDists <- dist(t(assay(rld)))            # Computes Euclidean distance between samples based on gene expression
sampleDistMatrix <- as.matrix(sampleDists)

png(filename="DESeq_output/DESeq_sampleDist_HM.png", units = 'in', width = 8, height = 8, res = 250)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         annotation = treat_ann,
         col=colors,
         main="Sample to Sample Distances")
dev.off()

## Principal Component Analysis
## Separates samples based on variation between sample's gene expression
## Greater variation will affect separation to a greater degree

data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

png('DESeq_output/DESeq_PCA.png', units='in', width=8, height=6, res=250)
ggplot(data, aes(PC1, PC2, color=condition)) +
  geom_point(size=3.5) +
  geom_text_repel(aes(label=name)) +
  scale_colour_manual(values = c("orange", "steelblue", 'red')) +
  theme_bw() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("")
dev.off()

###################################################################################################
#
#   Get table to convert names
#
###################################################################################################

## Convert the ensembl gene ID to gene name
## This will require an active internet connection
## Need to change according to species of interest

mart<- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))

ensembl_2_geneName <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  mart = mart)
head(ensembl_2_geneName)
names(ensembl_2_geneName) <- c("GENEID", "geneName")

###################################################################################################
#
#     DIFFERENTIAL EXPRESSSION ANALYSIS 
#
###################################################################################################
## DESeq = fx to calculate DE
## Combines multiple steps from DESeq
dds <- DESeq(dds)

## Write normalized counts to file
normalized.counts <- as.data.frame(counts(dds, normalized=TRUE ))
rownames(normalized.counts)<-sub("\\.[0-9]*", "", rownames(normalized.counts))
head(normalized.counts)

# Add gene name column
idx <- match( rownames(normalized.counts), ensembl_2_geneName$GENEID )
normalized.counts$geneName <- ensembl_2_geneName$geneName[ idx ]
head(normalized.counts)
write.table(normalized.counts, file = 'DESeq_output/DESeq_normalized_counts.csv', qmethod = NULL, sep = ',')


###################################################################################################
#
#  Getting fold changes from direct comparisons with control
#
###################################################################################################

## Volcano plot function
## Will be used to visualize the differential expression of genes in the next section
volcanoPlot <- function(df, line_val, fc_cut = 0, pv_cut = .05, padj_cut = .1) {
    log2_lim <- 10
    pval_lim <- 20
  
    df_plt <- as.data.frame(df) %>%
    mutate(threshold = ifelse(pvalue <= pv_cut & log2FoldChange < fc_cut, -1,
                              ifelse(pvalue <= pv_cut & log2FoldChange > fc_cut, 1, 0))) %>%
    mutate(threshold = as.factor(threshold)) %>%
    mutate(log_pval = -log10(pvalue)) %>%
    mutate(shape = ifelse(log_pval > pval_lim | abs(log2FoldChange) > log2_lim, 17, 16)) %>%
    mutate(log_pval = ifelse(log_pval > pval_lim, pval_lim * .99, log_pval),
           log2FoldChange = ifelse(abs(log2FoldChange) > log2_lim, log2_lim * .99 * sign(log2FoldChange), log2FoldChange))
  

  
  ##Construct the plot object
  g = ggplot(df_plt, aes(log2FoldChange, y=log_pval, color=threshold)) +
    geom_point(alpha=0.4, size=1.75, shape = df_plt$shape) +
    scale_colour_manual(values = c("blue", "gray", "red")) +
    geom_hline(yintercept = -log10(pv_cut), linetype = 'dashed') +
    annotate('text', label = 'p-value', x = -log2_lim * .98, y = -log10(pv_cut) + .15, vjust = 0, hjust = 0) +
    ggtitle('') +
    scale_x_continuous("log2 fold change",
                       limit = c(-log2_lim, log2_lim),
                       expand = c(0,0)) +
    scale_y_continuous("-log10 p-value",
                       limit = c(0, pval_lim),
                       expand = c(0, 0)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", fill = "white"), 
          text = element_text(size = 24),
          legend.position = "none")
  return(g)
}

###################################################################################################
#
#  Getting fold changes from direct comparisons with control
#
###################################################################################################
library(dplyr)

## Make comparison between conditions
samples
res <- results( dds, contrast = c("condition", "CondA", "Control") )
head(res)

## Add gene name to res file
rownames(res)<-sub("\\.[0-9]*", "", rownames(res))
idx <- match( rownames(res), ensembl_2_geneName$GENEID )
res$geneName <- ensembl_2_geneName$geneName[ idx ]

## Filter to remove genes with a baseMean of 5 or less
## baseMean is the average expression for that gene across all samples
res.5<-res[res$baseMean>5, ]

## Adjust p-value according to Benjamini & Hochberg method (need to do this since we filtered out genes by base mean 5 above)
res.5$padj <- p.adjust(res.5$pvalue, method="BH")

## Remove lines where pvalue is NA
res.5 <- res.5[!is.na(res.5$pvalue),]

## Check number of differentially expressed genes
as.data.frame(res.5) %>% filter(pvalue < .05) %>% dim()                       # Number of sig genes by p-value
as.data.frame(res.5) %>% filter(pvalue < .05 & log2FoldChange < 0) %>% dim()  # Number of sig down genes by p-value
as.data.frame(res.5) %>% filter(pvalue < .05 & log2FoldChange > 0) %>% dim()  # Number of sig up genes by p-value

## Write res.cont DESeq data to output file
write.csv(res.5, file="DESeq_output/DESeq_CONDAvsCONTROL.csv", quote=F, row.names = F)

## Construct gene list for miRhub
as.data.frame(res.5) %>%
  filter(pvalue < .05) %>%                  # p-value must be less than .05
  filter(log2FoldChange > 0) %>%            # log2FoldChange must be positive (aka upregulated genes)
  mutate(Name = toupper(geneName)) %>%      # Capitalize gene names
  pull(Name) %>%                            # Take only gene names column
  as.vector() %>%                           # Set as vector
  write.table(., file = 'DESeq_output/DEG_CONDAvsCONTROL_up_PVAL.05.txt', quote = F, row.name = F, col.names = F)     # Write to file

## Volcano plot
png('DESeq_output/DEG_CONDAvsCONTROL_VolcanoPlot.png', units = 'in', width = 6, height = 6, res = 250)
g <- volcanoPlot(res.5)
g
dev.off()

###########################

# Make comparison between conditions
samples
res <- results( dds, contrast = c("condition", "CondB", "Control") )
head(res)


## Add gene name to res file
rownames(res)<-sub("\\.[0-9]*", "", rownames(res))
idx <- match( rownames(res), ensembl_2_geneName$GENEID )
res$geneName <- ensembl_2_geneName$geneName[ idx ]

## Filter to remove genes with a baseMean of 5 or less
## baseMean is the average expression for that gene across all samples
res.5<-res[res$baseMean>5, ]

## Adjust p-value according to Benjamini & Hochberg method (need to do this since we filtered out by base mean 5 above)
res.5$padj <- p.adjust(res.5$pvalue, method="BH")

## Remove lines where pvalue is NA
res.5 <- res.5[!is.na(res.5$pvalue),]

## Check number of differentially expressed genes
as.data.frame(res.5) %>% filter(pvalue < .05) %>% dim()                       # Number of sig genes by p-value
as.data.frame(res.5) %>% filter(pvalue < .05 & log2FoldChange < 0) %>% dim()  # Number of sig down genes by p-value
as.data.frame(res.5) %>% filter(pvalue < .05 & log2FoldChange > 0) %>% dim()  # Number of sig up genes by p-value

## Write res.cont DESeq data to output file
write.csv(res.5, file="DESeq_output/DESeq_CONDBvsCONTROL.csv", quote=F, row.names = F)

## Construct gene list for miRhub
as.data.frame(res.5) %>%
  filter(pvalue < .05) %>%                  # p-value must be less than .05
  filter(log2FoldChange > 0) %>%            # log2FoldChange must be positive (aka upregulated genes)
  mutate(Name = toupper(geneName)) %>%      # Capitalize gene names
  pull(Name) %>%                            # Take only gene names column
  as.vector() %>%                           # Set as vector
  write.table(., file = 'DESeq_output/DEG_CONDBvsCONTROL_up_PVAL.05.txt', quote = F, row.name = F, col.names = F)     # Write to file

## Volcano plot
png('DESeq_output/DEG_CONDBvsCONTROL_VolcanoPlot.png', units = 'in', width = 6, height = 6, res = 250)
g <- volcanoPlot(res.5)
g
dev.off()
