
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
load(here("clean_data.RData"))
source("00_functions.R")

# select observations with neoepitope >  median (all datapoints)
top_s2c <- s2c %>% filter(neo_rank2>median(s2c$neo_rank2))

# select observations with neoepitope < median (all datapoints)
#top_s2c <- s2c %>% filter(neo_rank2<median(s2c$neo_rank2))

# make txi match the new s2c
txi_top <- lapply(txi[1:3], function(x) keep_columns(df1 = x, s2c_filtered = top_s2c))
txi_top$countsFromAbundance <- "no"


# standard DEA ------------------------------------------------------------

# change naming
s2c <- top_s2c
txi <- txi_top

# analyze data ------------------------------------------------------------
samples <- data.frame(sample = s2c$file_id, condition = s2c$PFS_cut, 
                      cohort = s2c$cohort, cancer_type = s2c$cancer_type, 
                      T_cell_diversity = s2c$T_cell_diversity, treatment=s2c$treatment,
                      neoepitope = s2c$neo_rank2, treatment_target = s2c$treatment_target)
dds <- DESeqDataSetFromTximport(txi, samples, design = ~cohort + condition)
# make boxplot of counts
boxplot(log10(counts(dds)+1))

# keep transcripts with at least 5 reads in 3 samples  --------------------
keep <- rowSums(counts(dds)) >= 10 # as depicted in deseq2 tutorial
dds <- dds[keep, ]


# define responder and non-responder --------------------------------------
dds$condition <- relevel(dds$condition, ref = "high_PFS")


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
# topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 100) # finder mest variable gener, bruges ikke 
# mat <- assay(vsd)[ topVarGenes, ]
# mat <- mat - rowMeans(mat) # center the gene expressions around 0 - why?
anno <- as.data.frame(colData(vsd)[,c("condition", "cohort", "treatment_target", 
                                      "cancer_type", "T_cell_diversity",
                                      "neoepitope")])

# make boxplot of variable genes 
# pheatmap(mat, annotation_col=anno)

# differential expression analysis ----------------------------------------
dds<-DESeq(dds)
resultsNames(dds)
res <- lfcShrink(dds, coef="condition_low_PFS_vs_high_PFS", type = "apeglm")

#res <- results(dds, cooksCutoff = TRUE, contrast = c("condition", "high_PFS", "low_PFS"), name = "highPFS_vs_lowPFS") 
df <- as.data.frame(res)

library(data.table)
df <- setDT(df, keep.rownames= TRUE)[] # df = as.dataframe(res)
dfp <- filter(df, padj < 0.05) # filtrering
dfp <- filter(dfp, log2FoldChange > 2 | log2FoldChange < -2 ) # filtrering
#dfp <- filter(dfp, log2FoldChange > 2 ) # get genes with positive logf fold change
dfp <- dfp %>% arrange(desc(abs(log2FoldChange)))


diff_genes <- head(dfp$rn, n=50) 
length(diff_genes)
#write_csv(x = as.data.frame(diff_genes), "diff_genes_p0.02_lfcShrinkage.csv")
#write_csv(x = as.data.frame(diff_genes), "diff_genes_p0.01_only_low_neoepi_lfcShrinkage.csv")
df_expression_DEA_genes <-all_data_df[rownames(all_data_df) %in% diff_genes,]
scales_df <- t(scale(t(df_expression_DEA_genes)))
pheatmap(scales_df, annotation_col=anno, clustering_method = "complete") # default 
pheatmap(scales_df, annotation_col=anno, clustering_method = "single") # seem to have a good effect

# export heatmaps -----------------------------------------------------------
# tiff(file= 'results/heatmap_both_cohort_targeted_only_high_neoepitope_padj0.02_lfcShrinkage.tiff', height = 1000, width = 800)
# pheatmap(scales_df, annotation_col=anno, clustering_method = "complete") # default 
# dev.off()


# export heatmaps low epitope sort ----------------------------------------
#tiff(file= 'results/heatmap_both_cohort_non-targeted_only_low_neoepitope_padj0.02_lfcShrinkage.tiff', height = 1000, width = 800)
#pheatmap(scales_df, annotation_col=anno, clustering_method = "single") # default
#dev.off()

###########################################################################
# do targeted analysis considering immunological genes --------------------
###########################################################################

# read immunogene panel 
IMM_Genes <- read.table(here('genes_750_immunological'), header = TRUE, sep = '\t', stringsAsFactors = FALSE)

dds_2 <- dds[rownames(dds)  %in% IMM_Genes$genes,]

dds_2<-DESeq(dds_2)
res_2 <- lfcShrink(dds_2, coef="condition_low_PFS_vs_high_PFS", type = "apeglm")

#res_2 <- results(dds_2, cooksCutoff = TRUE, contrast = c("condition", "high_PFS", "low_PFS")) 
df_2 <- as.data.frame(res_2)


df_2<- setDT(df_2, keep.rownames= TRUE)[]
dfp_2 <- filter (df_2, padj < 0.05)
dfp_2 <- filter(dfp_2, log2FoldChange > 1 | log2FoldChange < -1 )

diff_genes_2 <- dfp_2$rn 
length(diff_genes_2)
df_expression_DEA_genes_2 <-all_data_df[rownames(all_data_df) %in% diff_genes_2,]
scales_df_2 <- t(scale(t(df_expression_DEA_genes_2)))


# plot low ----------------------------------------------------------------
pheatmap(scales_df_2, annotation_col=anno)

# plot high epitopes ------------------------------------------------------
#tiff(file= 'results/heatmap_both_cohort_targeted_only_high_neoepitope_padj0.01.tiff', height = 1000, width = 800)
#pheatmap(scales_df_2, annotation_col=anno)
#dev.off()



