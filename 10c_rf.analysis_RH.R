
###########################################################################
# script for making random forest model of RH cohort ----------------------
###########################################################################

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
                                            #"treatment",
                                            #"cancer_type",
                                            "patient",
                                            "T_cell_diversity",
                                            "neo_rank2")

# change type of specific columns to character 
metadata[,1:3]%<>% mutate_if(is.character,as.factor)

rh_select <- rh_tpm

# add extra columns; neoepitope, cancer_type, treatment, T_cell_diversity
rh_data <- left_join(rh_select, metadata, by=c("Sample"="patient"))
colnames(rh_data)<- colnames(rh_data) %>% str_replace_all("-", "_") 
rh_data <- as_data_frame(rh_data)

# make categorical data to factor and remove 
rh_data<- as.data.frame(rh_data)
rh_sample <- rh_data$sample
rh_data$Sample <- NULL
rh_data <- rh_data %>% filter(!is.na(PFS_cut))

###########################################################################
# random forest with all genes  -------------------------------------------
###########################################################################
# https://mlr.mlr-org.com/articles/tutorial/feature_selection.html

# make random forest model ------------------------------------------------

# generate the task 
task = makeClassifTask(data = rh_data, target = "PFS_cut")

# filter methods for selection of features 
n_genes <- 30
filtered.task = filterFeatures(task, method = "FSelectorRcpp_information.gain", 
                               abs = n_genes,
                               mandatory.feat = c("T_cell_diversity", "neo_rank2"))

# Generate the learner
classif.lrn = makeLearner("classif.cforest", 
                          predict.type = "prob", 
                          fix.factors.prediction = TRUE, 
                          id = "RH_ranfomforest_all")

# set hyperparameters for random forest 
mtry <- round(sqrt(n_genes))
classif.lrn <- setHyperPars(classif.lrn, ntree = 10000, mtry=mtry)

# cross-validatation 
iter <- 5
rdesc = makeResampleDesc("CV", iters = iter)
r_cv <- resample(classif.lrn, filtered.task, rdesc, measures = list(acc, fn, fp,
                                                            mmce, 
                                                            ber, auc, 
                                                            timetrain))
cur_auc <- performance(r_cv$pred, auc)
df = generateThreshVsPerfData(r_cv$pred, measures = list(fpr, tpr, mmce))
plotROCCurves(df) + ggtitle(paste0(iter, "-fold CV ROC, AUC:", round(cur_auc,2)," RF all genes (filtered ", n_genes, ") mtry: ", mtry))
ggsave(filename = here("RF_results", paste0(iter,"-fold_ROC_RF_all_genes_filtered", n_genes, "_mtry", mtry, ".tiff")))

# full model and variable importance --------------------------------------
mod = train(classif.lrn, filtered.task)
varimp(mod$learner.model) %>% sort() 





# LOO - takes a long time, so for now it is outcommented 
# LOO = makeResampleDesc("LOO")
# r = resample(classif.lrn, task, LOO, measures = list(acc, fn, fp, mmce, ber, auc, timetrain))
# pred = getRRPredictions(r)
# df = generateThreshVsPerfData(r$pred, measures = list(fpr, tpr, mmce))
# plotROCCurves(df) + ggtitle("LOO CV ROC curve random forest all genes")




###########################################################################
# if I want to do twolayer CV use this code  ------------------------------
###########################################################################
# classif.lrn = makeLearner("classif.cforest", 
#                           predict.type = "prob", #"response", 
#                           fix.factors.prediction = TRUE, 
#                           id = "RH_ranfomforest_395genes")
# 
# # set hyperparameters for random forest 
# classif.lrn <- setHyperPars(classif.lrn, ntree = 10000, mtry=10)
# 
# # make wrapper for filter 
# lrn = makeFilterWrapper(learner = classif.lrn, fw.abs=1000, "FSelectorRcpp_information.gain")
# 
# rdesc = makeResampleDesc("CV", iters = 10)
# r = resample(learner = lrn, task = task, resampling = rdesc, show.info = FALSE, 
#              models = TRUE, measures = list(acc, fn, fp,
#                                             mmce, 
#                                             ber, auc, 
#                                             timetrain))
# 
# sfeats = sapply(r$models, getFilteredFeatures)
# table(sfeats) %>% as_tibble() %>% arrange(desc(n))
# selected_features <- table(sfeats) %>% as_tibble() %>% arrange(desc(n)) %>% filter(n>=8)

# fit random forest on selected features ----------------------------------

# rh_select <- rh_tpm[,colnames(rh_tpm) %in% c(selected_features$sfeats)]
# rh_select$sample <- rh_tpm$Sample
# 
# # add extra columns; neoepitope, cancer_type, treatment, T_cell_diversity
# rh_data <- left_join(rh_select, metadata, by=c("sample"="patient"))
# colnames(rh_data)<- colnames(rh_data) %>% str_replace_all("-", "_") 
# rh_data <- as_data_frame(rh_data)
# 
# # make categorical data to factor and remove 
# rh_data<- as.data.frame(rh_data)
# rh_sample <- rh_data$sample
# rh_data$sample <- NULL
# rh_data <- rh_data %>% filter(!is.na(PFS_cut))

# generate the ta





