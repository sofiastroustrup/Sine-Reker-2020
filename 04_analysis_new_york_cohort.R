
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
source(here("00_functions.R"))

# select pvalue for DEA ---------------------------------------------------
pval <- 0.5

# select only rigshospitalet cohort ---------------------------------------
clinical_data <- clinical_data %>% filter(cohort == "new_york")
s2c <- s2c %>% filter(cohort== "new_york")
txi[1:3] <- lapply(txi[1:3], function(x) select_cohort(df = x, cohort = "new_york")) 


# analyze data ------------------------------------------------------------
samples <- data.frame(sample = s2c$file_id, condition = s2c$PFS_cut) 
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

anno <- as.data.frame(colData(vsd)[, c("condition")])
colnames(anno) <- c("condition")
rownames(anno)<- rownames(colData(vsd))

# make boxplot of variable genes 
pheatmap(mat, annotation_col=anno)

###########################################################################
# differential expression analysis ----------------------------------------
###########################################################################

dds<-DESeq(dds)
res <- results(dds, cooksCutoff = TRUE, contrast = c("condition", "high_PFS", "low_PFS"), name = "highPFS_vs_lowPFS") 
df <- as.data.frame(res)

df <- setDT(df, keep.rownames= TRUE)[] # df = as.dataframe(res)
dfp <- filter(df, padj < 0.05) # filtrering
dfp <- filter(dfp, log2FoldChange > 1 | log2FoldChange < -1 ) # filtrering

diff_genes <- dfp$rn 
length(diff_genes)
df_expression_DEA_genes <-all_data_df[rownames(all_data_df) %in% diff_genes,]

# export df
write.csv(x = df_expression_DEA_genes, file = here("results", paste0("NY_DEA_genes_non-targeted_pval_", pval, ".csv")))

#pheatmap(df_expression_DEA_genes, annotation_col=anno)
scales_df <- t(scale(t(df_expression_DEA_genes)))

# save plot
tiff(file=paste0("results/heatmap_NY_non-targeted_DEA_pval_", pval, ".tiff"), height = 1000, width = 500)
pheatmap(scales_df, annotation_col=anno, clustering_method = "complete") # default, best for this cohort 
dev.off()


# plot heatmaps -----------------------------------------------------------
# explore cluserings
#pheatmap(scales_df, annotation_col=anno, clustering_method = "average")
pheatmap(scales_df, annotation_col=anno, clustering_method = "single") # seem to have a good effect
#pheatmap(scales_df, annotation_col=anno, clustering_method = "mcquitty")
#pheatmap(scales_df, annotation_col=anno, clustering_method = "median")
#pheatmap(scales_df, annotation_col=anno, clustering_method = "centroid")
#pheatmap(scales_df, annotation_col=anno, clustering_method = "ward.D")
#pheatmap(scales_df, annotation_col=anno, clustering_method = "ward.D2")


# do targeted analysis considering immunological genes --------------------
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

# export df
write.csv(x = df_expression_DEA_genes_2, file = here("results", paste0("NY_DEA_genes_targeted_pval_", pval, ".csv")))


# pheatmap(df_expression_DEA_genes, annotation_col=anno)
scales_df <- t(scale(t(df_expression_DEA_genes_2)))

tiff(file=paste0("results/heatmap_NY_targeted_DEA_pval_", pval, ".tiff"), height = 1000, width = 500)
pheatmap(scales_df, annotation_col=anno, clustering_method = "single")
dev.off()

#save(df, dfp, df_2, scales_df,all_data_df, file = "/home/projects/SRHgroup/projects/MSKCC/data/RNA/Kallisto_sleuth/MSKCC_DEA_NvsNR.Rdata"  )


