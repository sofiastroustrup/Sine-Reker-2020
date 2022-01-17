
# load packages ----------------------------------------------------------
library(tidyverse)
library(mlr)
library(party)
library(here)
library(magrittr)

# load data ---------------------------------------------------------------
load(here("clean_data.RData"))
source(here("00_functions.R"))
target_genes<- read_delim(here("395_genes.csv"), delim = ";", skip = 1)
rh_tpm<- readRDS("rh_tpm.RDS")

# prepare data  -----------------------------------------------------------

# check that all rh_samples are in clinical 
rh_tpm <- rh_tpm %>% filter(Sample %in% clinical_data$patient) %>% dplyr::select(-c("hugo_symbol"))
# scale data 
rh_tmp<- rh_tpm %>% dplyr::select(-c("Sample"))
rh_tmp <- t(rh_tmp)
rh_scale <- rh_tmp %>% scale() 
rh_scale <- t(rh_scale) %>% as_tibble
rh_scale$sample <- rh_tpm$Sample

# prepare bladder data
bladder_tpm <- readRDS(file=here("tpm_data", "tpm_bladder_cancer_spread2.RDS"))
#scale
bladder_tpm_sample <- bladder_tpm$sample
bladder_tpm %<>% dplyr::select(-c(sample))
bladder_scale <- t(bladder_tpm) %>% scale() 
bladder_scale <- t(bladder_scale) %>% as_tibble() 
bladder_scale$sample <- bladder_tpm_sample

# combine RH and NY
tpm <- rbind(rh_scale, bladder_scale)


# prepare metadata  -------------------------------------------------------
# select only rigshospitalet cohort
clinical_sele <- clinical_data%>% filter(patient %in% tpm$sample) 
clinical_sele <- clinical_sele %>% filter(!is.na(PFS_cut)) 
#clinical_sele <- clinical_sele %>% mutate(PFS, PFS_cut = ifelse(PFS>median(PFS), "high_PFS", "low_PFS"))

metadata <- clinical_sele %>% dplyr::select("PFS", 
                                            "patient",
                                            "T_cell_diversity",
                                            "neo_rank2", "PFS_factor")

# change type of specific columns to character 
metadata[,1:3]%<>% mutate_if(is.character,as.factor)

# add extra columns; neoepitope, cancer_type, treatment, T_cell_diversity
data <- left_join(tpm, metadata, by=c("sample"="patient"))
colnames(data)<- colnames(data) %>% str_replace_all("-", "_") 
data <- as_data_frame(data)

# make categorical data to factor and remove 
data <- data %>% filter(!is.na(PFS))
sample <- data$sample
data$sample <- NULL
data <- as_data_frame(data)


