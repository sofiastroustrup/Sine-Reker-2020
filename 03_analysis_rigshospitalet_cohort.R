

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
library(data.table)

# load data ---------------------------------------------------------------
load(here("clean_data.RData"))
source(here("00_functions.R"))


# select pvalue for DEA ---------------------------------------------------
pval <- 0.1

# This doesnt work because I didnt select RH!!
##########################################################################
# comment out if full analysis 
# try to divide data by neoepitope ----------------------------------------
# select observations with neoepitope >  median
#top_s2c <- s2c %>% filter(neo_rank2>median(s2c$neo_rank2))

#top_s2c <- s2c %>% filter(neo_rank2<median(s2c$neo_rank2))

# make txi match the new s2c
#txi_top <- lapply(txi[1:3], function(x) keep_columns(df1 = x, s2c_filtered = top_s2c))
#txi_top$countsFromAbundance <- "no"

# standard DEA ------------------------------------------------------------

# change naming
#s2c <- top_s2c
#txi <- txi_top
##########################################################################

# comment in if standard full analysis 
# select only rigshospitalet cohort ---------------------------------------
clinical_data <- clinical_data %>% filter(cohort == "rigshospitalet")
s2c <- s2c %>% filter(cohort== "rigshospitalet")
txi[1:3] <- lapply(txi[1:3], select_cohort) 


# analyze data ------------------------------------------------------------
samples <- data.frame(sample = s2c$file_id, condition = s2c$PFS_cut, 
                      neoepitope = s2c$neo_rank2, cancer_type = s2c$cancer_type,
                      treatment_target = s2c$treatment_target, 
                      t_cell_diversity = s2c$T_cell_diversity) 
dds <- DESeqDataSetFromTximport(txi, samples, design = ~condition)
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
vsd <- vst(dds) # variance stabilizing transformation?

#pca plot 
plotPCA(vsd, intgroup = "condition")



# consider variable genes (not differential expression analysis) ----------
all_data_df <- assay(vsd)
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 100) # finder mest variable gener, bruges ikke 
mat <- assay(vsd)[ topVarGenes, ]
mat <- mat - rowMeans(mat) # center the gene expressions around 0 - why?

anno <- as.data.frame(colData(vsd)[, c("condition", "neoepitope", "cancer_type",
                                       "treatment_target", "t_cell_diversity")])
#colnames(anno) <- c("condition")
#rownames(anno)<- rownames(colData(vsd))

# make boxplot of variable genes 
# pheatmap(mat, annotation_col=anno)

###########################################################################
# differential expression analysis ----------------------------------------
###########################################################################

dds<-DESeq(dds)
resultsNames(dds)
res <- lfcShrink(dds, coef="condition_low_PFS_vs_high_PFS", type = "apeglm")
#res <- results(dds, cooksCutoff = TRUE, contrast = c("condition", "high_PFS", "low_PFS"), name = "highPFS_vs_lowPFS") 
df <- as.data.frame(res)

df <- setDT(df, keep.rownames= TRUE)[] # df = as.dataframe(res)
dfp <- filter(df, padj < pval) # filtrering
dfp <- filter(dfp, log2FoldChange > 2 | log2FoldChange < -2 ) # filtrering
dfp <- dfp %>% arrange(desc(abs(log2FoldChange)))


diff_genes <- head(dfp$rn, n=50) 
#diff_genes <- dfp$rn 
length(diff_genes)
df_expression_DEA_genes <-all_data_df[rownames(all_data_df) %in% diff_genes,]

# export df
#write.csv(x = df_expression_DEA_genes, file = here("results", paste0("RH_DEA_genes_non-targeted_pval_", pval, ".csv")))

# make heatmap of scaled expression matrix 
scales_df <- t(scale(t(df_expression_DEA_genes)))
pheatmap(scales_df, annotation_col=anno, clustering_method = "complete")

tiff(file= paste0("results/heatmap_RH_non-targeted_only_low_neoepitope_pval_", pval, ".tiff"), height = 1000, width = 1500)
pheatmap(scales_df, annotation_col=anno, clustering_method = "complete") # default is complete 
dev.off()

# plot heatmaps -----------------------------------------------------------
# explore cluserings
#pheatmap(scales_df, annotation_col=anno, clustering_method = "average")
#pheatmap(scales_df, annotation_col=anno, clustering_method = "single") # seem to have a good effect
#pheatmap(scales_df, annotation_col=anno, clustering_method = "mcquitty")
#pheatmap(scales_df, annotation_col=anno, clustering_method = "median")
#pheatmap(scales_df, annotation_col=anno, clustering_method = "centroid")
#pheatmap(scales_df, annotation_col=anno, clustering_method = "ward.D")
#pheatmap(scales_df, annotation_col=anno, clustering_method = "ward.D2")


#png(file= 'results/heatmap_vgs_diff.png', height = 2000, width = 1000)

#dev.off()
#save(df, dfp, scales_df,all_data_df, file = here("results", "mm909_DEA.Rdata") )



# do targeted analysis considering immunological genes --------------------
# read immunogene panel 
IMM_Genes <- read.table(here('genes_750_immunological'), header = TRUE, sep = '\t', stringsAsFactors = FALSE)

dds_2 <- dds[rownames(dds)  %in% IMM_Genes$genes,]

dds_2<-DESeq(dds_2)
res_2 <- lfcShrink(dds_2, coef="condition_low_PFS_vs_high_PFS", type = "apeglm")
#res_2 <- results(dds_2, cooksCutoff = TRUE, contrast = c("condition", "high_PFS", "low_PFS")) 
df_2 <- as.data.frame(res_2)


df_2<- setDT(df_2, keep.rownames= TRUE)[]
dfp_2 <- filter (df_2, padj < pval)
dfp_2 <- filter(dfp_2, log2FoldChange > 2 | log2FoldChange < -2 )

diff_genes_2 <- dfp_2$rn 
length(diff_genes_2)
df_expression_DEA_genes_2 <-all_data_df[rownames(all_data_df) %in% diff_genes_2,]

# export df
#write.csv(x = df_expression_DEA_genes_2, file = here("results", paste0("RH_DEA_genes_targeted_", "pval_", pval, ".csv")))


# pheatmap(df_expression_DEA_genes, annotation_col=anno)
scales_df <- t(scale(t(df_expression_DEA_genes_2)))
pheatmap(scales_df, annotation_col=anno, clustering_method = "complete")

# save plot
#tiff(file= paste0("results/", "heatmap_RH_targeted_DEA_", "pval_", pval, ".tiff"), height = 1000, width = 500)
#pheatmap(scales_df, annotation_col=anno, clustering_method = "complete")
#dev.off()

#save(df, dfp, df_2, scales_df,all_data_df, file = "/home/projects/SRHgroup/projects/MSKCC/data/RNA/Kallisto_sleuth/MSKCC_DEA_NvsNR.Rdata"  )






