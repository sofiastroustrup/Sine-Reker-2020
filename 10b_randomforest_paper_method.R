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

###########################################################################
# random forest with 395 genes --------------------------------------------
###########################################################################

# select the 395 genes
rh_select <- rh_data[,colnames(rh_data) %in% c(target_genes$Gene, 
                                              "PFS_cut", 
                                              "T_cell_diversity", 
                                              "neo_rank2")]

# make random forest ------------------------------------------------------
# generate the task 
task = makeClassifTask(data = rh_select, target = "PFS_cut")

# Generate the learner
classif.lrn = makeLearner("classif.cforest", 
                          predict.type = "prob",
                          fix.factors.prediction = TRUE, 
                          id = "RH_ranfomforest_395genes")

# set hyperparameters for random forest 
classif.lrn <- setHyperPars(classif.lrn, ntree = 10000, mtry=20)


# cross-validatation ------------------------------------------------------

# 5-fold CV
iter <- 5
rdesc = makeResampleDesc("CV", iters = iter)
r_cv <- resample(classif.lrn, task, rdesc, measures = list(acc, fn, fp, mmce, ber, auc, timetrain))
cur_auc <- performance(r_cv$pred, auc)

df = generateThreshVsPerfData(r_cv$pred, measures = list(fpr, tpr, mmce), gridsize = 100)
plotROCCurves(df) + ggtitle(paste0(iter, "-fold CV ROC, AUC:", round(cur_auc, 2)," random forest 395 genes + TcD + neoepitope"))
ggsave(filename = here("RF_results", paste0(iter,"-fold_ROC_RF_395genes.tiff")))

# full model and variable importance --------------------------------------
mod = train(classif.lrn, task)
varimp(mod$learner.model) %>% sort()


###########################################################################
# random forest with 395 genes and feature selection ----------------------
###########################################################################
# select the 395 genes
rh_select <- rh_data[,colnames(rh_data) %in% c(target_genes$Gene, 
                                               "PFS_cut", 
                                               "T_cell_diversity", 
                                               "neo_rank2")]
# generate the task 
task.395.f = makeClassifTask(data = rh_select, target = "PFS_cut")

# filter task
n_genes <- 30
filtered.task = filterFeatures(task.395.f, method = "FSelectorRcpp_information.gain", abs = n_genes)

# Generate the learner
classif.lrn = makeLearner("classif.cforest", 
                          predict.type = "prob",
                          fix.factors.prediction = TRUE, 
                          id = "RH_ranfomforest_395genes")

# set hyperparameters for random forest 
classif.lrn <- setHyperPars(classif.lrn, ntree = 10000, mtry=20)

# 5-fold CV
iter <- 5
rdesc = makeResampleDesc("CV", iters = iter)
r_cv <- resample(classif.lrn, filtered.task, rdesc, measures = list(acc, fn, fp, mmce, ber, auc, timetrain))

cur_auc <- performance(r_cv$pred, auc)
df = generateThreshVsPerfData(r_cv$pred, measures = list(fpr, tpr, mmce), gridsize = 100)
plotROCCurves(df) + ggtitle(paste0(iter, "-fold ROC, auc=", round(cur_auc, digits = 2)," RF 395 genes filtered ", n_genes, " genes"))
ggsave(filename = here("RF_results", paste0(iter,"-fold_ROC_RF_395genes_filtered_", n_genes, ".tiff")))

###########################################################################
# random forest with M1 signature and peripheral T-cell signature ---------
###########################################################################

M1.sig <- c("CBLB", "CCR7", "CD27", "CD48", "FOXO1", "FYB", "HLA-B", "HLA-G", 
            "IFIH1", "IKZF4", "LAMP3", "NFKBIA", "SAMHD1")
pT.cell.sig <- c("HLA-DOA", "GPR18", "STAT1")

rh_select <- rh_data[,colnames(rh_data) %in% c(M1.sig, pT.cell.sig, "PFS_cut")]



# make random forest ------------------------------------------------------
# generate the task 
task = makeClassifTask(data = rh_select, target = "PFS_cut")

# do oversampling to account for class imbalance
#task = oversample(task, rate = 2)

# Generate the learner
classif.lrn = makeLearner("classif.cforest", 
                          predict.type = "prob", #"response", 
                          fix.factors.prediction = TRUE, 
                          id = "RH_ranfomforest_395genes")

# set hyperparameters for random forest 
classif.lrn <- setHyperPars(classif.lrn, ntree = 10000, mtry=10)


# cross-validatation ------------------------------------------------------

# 5-fold CV
iter <- 5
rdesc = makeResampleDesc("CV", iters = iter)
r_cv <- resample(classif.lrn, task, rdesc, measures = list(acc, fn, fp, mmce, ber, auc, timetrain))
cur_auc <- performance(r_cv$pred, auc)
df = generateThreshVsPerfData(r_cv$pred, measures = list(fpr, tpr, mmce))
plotROCCurves(df) + ggtitle(paste0(iter, "-fold CV ROC, AUC=", round(cur_auc,2), " RF 16 genes"))
ggsave(filename = here("RF_results", paste0(iter,"-fold_ROC_RF_16_genes.tiff")))

# LOO 
#LOO = makeResampleDesc("LOO")
#r = resample(classif.lrn, task, LOO, measures = list(acc, fn, fp, mmce, ber, auc, timetrain))
#pred = getRRPredictions(r)
# df = generateThreshVsPerfData(r_cv$pred, measures = list(fpr, tpr, mmce))
# plotROCCurves(df) + ggtitle("5-fold CV ROC curve random forest 395 genes")

# full model and variable importance --------------------------------------
mod = train(classif.lrn, task)
varimp(mod$learner.model) %>% sort()

###########################################################################
# random forest with M1 signature -----------------------------------------
###########################################################################

M1.sig <- c("CBLB", "CCR7", "CD27", "CD48", "FOXO1", "FYB", "HLA-B", "HLA-G", 
            "IFIH1", "IKZF4", "LAMP3", "NFKBIA", "SAMHD1")

rh_select <- rh_data[,colnames(rh_data) %in% c(M1.sig, "PFS_cut")]



# make random forest ------------------------------------------------------
# generate the task 
task = makeClassifTask(data = rh_select, target = "PFS_cut")

# do oversampling to account for class imbalance
#task = oversample(task, rate = 2)

# Generate the learner
classif.lrn = makeLearner("classif.cforest", 
                          predict.type = "prob", #"response", 
                          fix.factors.prediction = TRUE, 
                          id = "RH_ranfomforest_395genes")

# set hyperparameters for random forest 
classif.lrn <- setHyperPars(classif.lrn, ntree = 10000, mtry=10)


# cross-validatation ------------------------------------------------------

# 5-fold CV
iter <- 5
rdesc = makeResampleDesc("CV", iters = iter)
r_cv <- resample(classif.lrn, task, rdesc, measures = list(acc, fn, fp, mmce, ber, auc, timetrain))
cur_auc <- performance(r_cv$pred, auc)
df = generateThreshVsPerfData(r_cv$pred, measures = list(fpr, tpr, mmce))
plotROCCurves(df) + ggtitle(paste0(iter, "-fold CV ROC, AUC=", round(cur_auc,2), " RF M1 genes"))
ggsave(filename = here("RF_results", paste0(iter,"-fold_ROC_RF_M1_genes.tiff")))

# full model and variable importance --------------------------------------
mod = train(classif.lrn, task)
varimp(mod$learner.model) %>% sort()


###########################################################################
# random forest with T-cell signature -------------------------------------
###########################################################################
pT.cell.sig <- c("HLA-DOA", "GPR18", "STAT1")

rh_select <- rh_data[,colnames(rh_data) %in% c(pT.cell.sig, "PFS_cut")]



# make random forest ------------------------------------------------------
# generate the task 
task = makeClassifTask(data = rh_select, target = "PFS_cut")

# do oversampling to account for class imbalance
#task = oversample(task, rate = 2)

# Generate the learner
classif.lrn = makeLearner("classif.cforest", 
                          predict.type = "prob", #"response", 
                          fix.factors.prediction = TRUE, 
                          id = "RH_ranfomforest_395genes")

# set hyperparameters for random forest 
classif.lrn <- setHyperPars(classif.lrn, ntree = 10000, mtry=3)


# cross-validatation ------------------------------------------------------

# 5-fold CV
iter <- 5
rdesc = makeResampleDesc("CV", iters = iter)
r_cv <- resample(classif.lrn, task, rdesc, measures = list(acc, fn, fp, mmce, ber, auc, timetrain))
cur_auc <- performance(r_cv$pred, auc)
df = generateThreshVsPerfData(r_cv$pred, measures = list(fpr, tpr, mmce))
plotROCCurves(df) + ggtitle(paste0(iter, "-fold CV ROC, AUC=", round(cur_auc,2), " RF peripheral T-cell genes"))
ggsave(filename = here("RF_results", paste0(iter,"-fold_ROC_RF_peripheral_Tcell_genes.tiff")))

# full model and variable importance --------------------------------------
mod = train(classif.lrn, task)
varimp(mod$learner.model) %>% sort()
