
# load packages --------------------------------------------------------
library(tidyverse)
library(biomaRt)
library(DESeq2)
library(tximport)
library(here)
library(pheatmap)
library(genefilter)
library(limma)
library(apeglm)
library(IHW)

# load data ---------------------------------------------------------------
load(here("clean_data_RECISTcut.RData"))


# analyze data ------------------------------------------------------------
samples <- data.frame(sample = s2c$file_id, condition = s2c$RECIST_cut, cohort = s2c$cohort) 
dds <- DESeqDataSetFromTximport(txi, samples, design = ~condition+cohort)
# make boxplot of counts
boxplot(log10(counts(dds)+1))

# keep transcripts with at least 5 reads in 3 samples  --------------------
keep <- rowSums(counts(dds)) >= 10 # as depicted in deseq2 tutorial
dds <- dds[keep, ]


# define responder and non-responder --------------------------------------
dds$condition <- relevel(dds$condition, ref = "R")


# initial data prep ---------------------------------------------------------------
boxplot(log10(counts(dds)+1))
dds <- estimateSizeFactors(dds)
boxplot(log10(counts(dds, normalized=TRUE)+1))
vsd <- vst(dds)

plotPCA(vsd, intgroup = "condition")
plotPCA(vsd, intgroup = "cohort")

# remove batch effect
assay(vsd) <- limma::removeBatchEffect(x=assay(vsd), batch = vsd$cohort) # remove batch effect of cohort
plotPCA(vsd, intgroup = "cohort")
plotPCA(vsd, intgroup = "condition")



all_data_df <- assay(vsd)
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 100) # finder mest variable gener, bruges ikke 
mat <- assay(vsd)[ topVarGenes, ]
mat <- mat - rowMeans(mat) # center the gene expressions around 0 - why?
anno <- as.data.frame(colData(vsd)[,c("condition", "cohort")])

# make boxplot of variable genes 
pheatmap(mat, annotation_col=anno)

# differential expression analysis ----------------------------------------
dds<-DESeq(dds)
res <- results(dds, cooksCutoff = TRUE, contrast = c("condition", "R", "NR"), name = "R_vs_NR") 
df <- as.data.frame(res)


# plot results from DEA ---------------------------------------------------
library(data.table)
df <- setDT(df, keep.rownames= TRUE)[] # df = as.dataframe(res)
dfp <- filter(df, padj < 0.01) # filtrering
dfp <- filter(dfp, log2FoldChange > 2 | log2FoldChange < -2 ) # filtrering
dfp <- dfp %>% arrange(desc(abs(log2FoldChange)))


diff_genes <- head(dfp$rn, n=100) 
length(diff_genes)
df_expression_DEA_genes <-all_data_df[rownames(all_data_df) %in% diff_genes,]
pheatmap(df_expression_DEA_genes, annotation_col=anno)
scales_df <- t(scale(t(df_expression_DEA_genes)))



# plot heatmaps -----------------------------------------------------------
tiff(file='results/heatmap_both_cohorts_non-targeted_RECIST_top100.tiff', height = 1000, width = 500)
pheatmap(scales_df, annotation_col=anno, clustering_method = "complete") # default 
dev.off()

# pheatmap(scales_df, annotation_col=anno, clustering_method = "average")
pheatmap(scales_df, annotation_col=anno, clustering_method = "single") # seem to have a good effect
# pheatmap(scales_df, annotation_col=anno, clustering_method = "mcquitty")
# pheatmap(scales_df, annotation_col=anno, clustering_method = "median")
# pheatmap(scales_df, annotation_col=anno, clustering_method = "centroid")
# pheatmap(scales_df, annotation_col=anno, clustering_method = "ward.D")
# pheatmap(scales_df, annotation_col=anno, clustering_method = "ward.D2")


###########################################################################
# do targeted analysis considering immunological genes --------------------
###########################################################################

# read immunogene panel 
IMM_Genes <- read.table(here('genes_750_immunological'), header = TRUE, sep = '\t', stringsAsFactors = FALSE)

dds_2 <- dds[rownames(dds)  %in% IMM_Genes$genes,]

dds_2<-DESeq(dds_2)
res_2 <- results(dds_2, cooksCutoff = TRUE, contrast = c("condition", "R", "NR")) 
df_2 <- as.data.frame(res_2)


df_2<- setDT(df_2, keep.rownames= TRUE)[]
dfp_2 <- filter (df_2, padj < 0.05)
dfp_2 <- filter(dfp_2, log2FoldChange > 2 | log2FoldChange < -2 )

diff_genes_2 <- dfp_2$rn 
length(diff_genes_2)
df_expression_DEA_genes_2 <-all_data_df[rownames(all_data_df) %in% diff_genes_2,]
scales_df <- t(scale(t(df_expression_DEA_genes_2)))

tiff(file='results/heatmap_both_cohorts_targeted_RECIST.tiff', height = 1000, width = 500)
pheatmap(scales_df, annotation_col=anno, clustering_method = "complete")
dev.off()



