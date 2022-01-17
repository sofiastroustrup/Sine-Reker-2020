###########################################################################
# make RF for RH cohort using few genes suggested by annie  ---------------
###########################################################################


# data preparation --------------------------------------------------------
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
rh_sample <- rh_tpm$Sample

# prepare data  -----------------------------------------------------------
# select only rigshospitalet cohort
clinical_data <- clinical_data %>% filter(cohort == "rigshospitalet")
clinical_sele <- clinical_data%>% filter(patient %in% rh_tpm$Sample) 
clinical_sele <- clinical_sele %>% filter(!is.na(PFS_cut)) 
clinical_sele <- clinical_sele %>% mutate(PFS, PFS_cut = ifelse(PFS>median(PFS), "high_PFS", "low_PFS"))

metadata <- clinical_sele %>% dplyr::select("PFS_cut", 
                                            "treatment",
                                            "cancer_type",
                                            "patient",
                                            "T_cell_diversity",
                                            "neo_rank2")

# change type of specific columns to character 
metadata[,1:3]%<>% mutate_if(is.character,as.factor)

# add extra columns; neoepitope, cancer_type, treatment, T_cell_diversity
rh_data <- left_join(rh_tpm, metadata, by=c("Sample"="patient"))
colnames(rh_data)<- colnames(rh_data) %>% str_replace_all("-", "_") 
rh_data <- as_data_frame(rh_data)

# make categorical data to factor and remove 
rh_data<- as.data.frame(rh_data)
rh_sample <- rh_data$Sample
rh_data$Sample <- NULL
rh_data <- rh_data %>% filter(!is.na(PFS_cut))


# make rf model -----------------------------------------------------------

# select the 395 genes
rh_select <- rh_data[,colnames(rh_data) %in% c("GZMA", "PRF1", "CD8A", 
                                               "PFS_cut", 
                                               "T_cell_diversity", 
                                               "neo_rank2")]


# generate the task 
task = makeClassifTask(data = rh_select, target = "PFS_cut")

# Generate the learner
classif.lrn = makeLearner("classif.cforest", 
                          predict.type = "prob", 
                          fix.factors.prediction = TRUE, 
                          id = "RH_ranfomforest_all")

# set hyperparameters for random forest 
#mtry <- round(sqrt(n_genes))
mtry<- 2
classif.lrn <- setHyperPars(classif.lrn, ntree = 100, mtry=mtry)

# cross-validatation 
iter <- 5
rdesc = makeResampleDesc("CV", iters = iter)
r_cv <- resample(classif.lrn, task, rdesc, measures = list(acc, fn, fp,
                                                                    mmce, 
                                                                    ber, auc, 
                                                                    timetrain))
cur_auc <- performance(r_cv$pred, auc)
df = generateThreshVsPerfData(r_cv$pred, measures = list(fpr, tpr, mmce))
plotROCCurves(df) + ggtitle(paste0(iter, "-fold CV ROC, AUC:", round(cur_auc,2)," RF annie genes, mtry: ", mtry))
ggsave(filename = here("RF_results", paste0(iter,"-fold_ROC_RF_annie_mtry", mtry, ".tiff")))

# full model and variable importance --------------------------------------
mod = train(classif.lrn, task)
varimp(mod$learner.model) %>% sort() 
