
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

# separate by cohort and standardize ------------------------------------------------------
#clinical_data_RH <- clinical_data %>% filter(cohort == "rigshospitalet")
#s2c_RH <- s2c %>% filter(cohort== "rigshospitalet")

# standardize rigshospitalet individually
txi_RH <- lapply(txi[1:3], select_cohort) 
txi_RH <- lapply(txi_RH[1:3], scale)
#txi_RH$countsFromAbundance<- txi$countsFromAbundance

# standardize NY cohort individually
txi_NY <-  lapply(txi[1:3], function(x) select_cohort(df = x, cohort = "new_york"))
txi_NY <- lapply(txi_NY[1:3], scale)
#txi_NY$countsFromAbundance<- txi$countsFromAbundance

# combine the two standardized lists
txi_s <- map2(txi_NY, txi_RH, cbind)
txi_s$countsFromAbundance<- txi$countsFromAbundance


# read in txi -------------------------------------------------------------
samples <- data.frame(sample = s2c$file_id, condition = s2c$PFS_cut, cohort = s2c$cohort) 
dds <- DESeqDataSetFromTximport(txi_s, samples, design = ~cohort + condition) 
# deseq returns the fold change that results only from the effect of condition
#  Deseq2 treats the first factor as a co-variate and tries to eliminate the 
# fold change that result because of this co-variate.

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
