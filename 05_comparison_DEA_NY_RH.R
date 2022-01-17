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


# read scripts ------------------------------------------------------------
source(here("00_functions.R"))


# comparison of results from DEA of different cohorts ---------------------

# read data 
pnow <- 0.5
NY_non_t<- read_csv(file = here(paste0("results/NY_DEA_genes_non-targeted_pval_", pnow, ".csv")))
NY_t <- read_csv(file = here(paste0("results/NY_DEA_genes_targeted_pval_", pnow, ".csv")),col_names = TRUE)
RH_non_t <- read_csv(file = here(paste0("results/RH_DEA_genes_non-targeted_pval_", pnow, ".csv")))
RH_t <- read_csv(file = here(paste0("results/RH_DEA_genes_targeted_pval_", pnow, ".csv")))


# compare differentially expressed genes in non targeted analysis ---------
intersect(NY_non_t$X1, RH_non_t$X1)
intersect(NY_t$X1, RH_t$X1)

intersect(NY_non_t$X1, NY_t$X1)
intersect(RH_non_t$X1, RH_t$X1)
