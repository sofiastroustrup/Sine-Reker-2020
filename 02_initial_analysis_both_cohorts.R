
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


# analyze data ------------------------------------------------------------
samples <- data.frame(sample = s2c$file_id, condition = s2c$PFS_cut, 
                      cohort = s2c$cohort, cancer_type = s2c$cancer_type, 
                      T_cell_diversity = s2c$T_cell_diversity, treatment=s2c$treatment,
                      neoepitope = s2c$neo_rank2, treatment_target = s2c$treatment_target)
dds <- DESeqDataSetFromTximport(txi, samples, design = ~condition+cohort)
# make boxplot of counts
boxplot(log10(counts(dds)+1))

# keep transcripts with at least 5 reads in 3 samples  --------------------
#keep <- rowSums(counts(dds) >= 5) >=3  annie version
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
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 100) # finder mest variable gener, bruges ikke 
mat <- assay(vsd)[ topVarGenes, ]
mat <- mat - rowMeans(mat) # center the gene expressions around 0 - why?
anno <- as.data.frame(colData(vsd)[,c("condition", "cohort", "treatment_target", "neoepitope", "T_cell_diversity")])

# make boxplot of variable genes 
#pheatmap(mat, annotation_col=anno)

# differential expression analysis ----------------------------------------
dds<-DESeq(dds)
#res <- results(dds, cooksCutoff = TRUE, contrast = c("condition", "high_PFS", "low_PFS"), name = "highPFS_vs_lowPFS") 
res <- lfcShrink(dds, coef="condition_low_PFS_vs_high_PFS" , type="apeglm")
df <- as.data.frame(res)


# log fold change and shrinkage -------------------------------------------
# not sure whether this is relevant for our approach
# resultsNames(dds) # not sure why newyork_vs_rigshospitalet appears
# resLFC <- lfcShrink(dds, coef="condition_low_PFS_vs_high_PFS" , type="apeglm")
# resLFC
# 
# 
# # independent hypothesis testing ------------------------------------------
# # not sure this is needed ??
# resIHW <- results(dds, filterFun=ihw)
# summary(resIHW)
# sum(resIHW$padj < 0.1, na.rm=TRUE)
# metadata(resIHW)$ihwResult
# 
# 
# # MA plot -----------------------------------------------------------------
# DESeq2::plotMA(res)
# DESeq2::plotMA(resLFC)
# 
# # data transformation and visualization -----------------------------------
# #vsd1 <- vst(dds, blind=FALSE)
# rld <- rlog(dds, blind=FALSE)
# head(assay(vsd), 3)
# 
# 
# # effect of transformation on the variance --------------------------------
# # this gives log2(n + 1)
# ntd <- normTransform(dds)
# library(vsn)
# meanSdPlot(assay(ntd))
# meanSdPlot(assay(vsd))
# meanSdPlot(assay(rld))
# 
# # Data quality assesment by sample clustering and visualization -----------
# select <- order(rowMeans(counts(dds,normalized=TRUE)),
#                 decreasing=TRUE)[1:20]
# df <- as.data.frame(colData(dds)[,c("cohort", "condition")])
# 
# pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#          cluster_cols=FALSE, annotation_col=df)
# 
# # assay = vsd
# pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#          cluster_cols=FALSE, annotation_col=df)
# 
# #assay rld
# pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#          cluster_cols=FALSE, annotation_col=df)
# 
# 
# 
# # heatmap of sample to sample distance ------------------------------------
# # get sample to sample distances
# sampleDists <- dist(t(assay(vsd)))
# 
# library("RColorBrewer")
# sampleDistMatrix <- as.matrix(sampleDists)
# rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
# colnames(sampleDistMatrix) <- NULL
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# pheatmap(sampleDistMatrix,
#          clustering_distance_rows=sampleDists,
#          clustering_distance_cols=sampleDists,
#          col=colors)
# 
# # pca plot of the samples
# plotPCA(vsd, intgroup=c("condition", "cohort"))


library(data.table)
df <- setDT(df, keep.rownames= TRUE)[] # df = as.dataframe(res)
dfp <- filter(df, padj < 0.01) # filtrering
dfp <- filter(dfp, log2FoldChange > 2 | log2FoldChange < -2 ) # filtrering
dfp <- dfp %>% arrange(desc(abs(log2FoldChange)))


diff_genes <- head(dfp$rn, n=50) 
length(diff_genes)
df_expression_DEA_genes <-all_data_df[rownames(all_data_df) %in% diff_genes,]
scales_df <- t(scale(t(df_expression_DEA_genes)))
pheatmap(scales_df, annotation_col=anno)



# plot heatmaps -----------------------------------------------------------
tiff(file= 'results/heatmap_both_cohort_non-targeted_complete_top50_treatment.tiff', height = 1000, width = 500)
pheatmap(scales_df, annotation_col=anno, clustering_method = "complete") # default 
dev.off()

#tiff(file= 'results/heatmap_both_cohort_non-targeted_single_top50.tiff', height = 1000, width = 500)
#pheatmap(scales_df, annotation_col=anno, clustering_method = "single") # seem to have a good effect
#dev.off()

# pheatmap(scales_df, annotation_col=anno, clustering_method = "average")
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
res_2 <- results(dds_2, cooksCutoff = TRUE, contrast = c("condition", "high_PFS", "low_PFS")) 
df_2 <- as.data.frame(res_2)


df_2<- setDT(df_2, keep.rownames= TRUE)[]
dfp_2 <- filter (df_2, padj < 0.05)
dfp_2 <- filter(dfp_2, log2FoldChange > 1 | log2FoldChange < -1 )

diff_genes_2 <- dfp_2$rn 
length(diff_genes_2)
df_expression_DEA_genes_2 <-all_data_df[rownames(all_data_df) %in% diff_genes_2,]
# pheatmap(df_expression_DEA_genes, annotation_col=anno)
scales_df <- t(scale(t(df_expression_DEA_genes_2)))

tiff(file= 'results/heatmap_both_cohort_targeted_single_pval0.05.tiff', height = 1000, width = 500)
pheatmap(scales_df, annotation_col=anno, clustering_method = "single")
dev.off()

tiff(file= 'results/heatmap_both_cohort_targeted_complete_pval0.05.tiff', height = 1000, width = 500)
pheatmap(scales_df, annotation_col=anno, clustering_method = "complete")
dev.off()
