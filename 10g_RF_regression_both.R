# REGRESSION TREE
###########################################################################
# Random forest for Rigshospitalet and new york cohort --------------------
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
clinical_sele <- clinical_data%>% filter(patient %in% rh_tpm$Sample) 
clinical_sele <- clinical_sele %>% filter(!is.na(PFS_cut)) 
#clinical_sele <- clinical_sele %>% mutate(PFS, PFS_cut = ifelse(PFS>median(PFS), "high_PFS", "low_PFS"))

metadata <- clinical_sele %>% dplyr::select("PFS", 
                                            "patient",
                                            "T_cell_diversity",
                                            "neo_rank2")

# change type of specific columns to character 
metadata[,1:3]%<>% mutate_if(is.character,as.factor)
#metadata%<>% mutate_if(is.double,as.numeric)
#metadata %<>% mutate(PFS, PFS=as.factor(PFS))

# add extra columns; neoepitope, cancer_type, treatment, T_cell_diversity
data <- left_join(tpm, metadata, by=c("sample"="patient"))
colnames(data)<- colnames(data) %>% str_replace_all("-", "_") 
data <- as_data_frame(data)

# make categorical data to factor and remove 
sample <- data$sample
data$sample <- NULL
data <- data %>% filter(!is.na(PFS))
data <- as_data_frame(data)


###########################################################################
# make random forest with regression output -------------------------------
###########################################################################

# generate the task 
task = makeRegrTask(data = data, target = "PFS")

# filter methods for selection of features 
#n_genes <- 5
fv <- generateFilterValuesData(task,
                               method = "FSelectorRcpp_information.gain", equal=TRUE)
fv$data %>% arrange(desc(value))%>% head(n=100)

filtered.task = filterFeatures(task, method = "FSelectorRcpp_information.gain", 
                               #abs = n_genes,
                               threshold = 0.16,
                               mandatory.feat = c("neo_rank2"),
                               equal=TRUE) #"T_cell_diversity", - from feature importance T-cell diversity doesnt seem to be relevant

# Generate the learner
model <- "regr.nnet" #"regr.randomForestSRC" #"regr.randomForest" #regr.randomForestSRC" # "regr.extraTrees" #"regr.evtree" #"regr.ctree" regr.randomForest regr.ranger 
regr.forest = makeLearner(model, 
                          fix.factors.prediction = TRUE, 
                          id = "RH_ranfomforest_all")

# set hyperparameters for random forest 
mtry <- round(sqrt(filtered.task$task.desc$n.feat[1]))
#regr.forest<- setHyperPars(regr.forest,  ntree = 10000000, mtry=mtry, splitrule="quantile.regr" ) #ntree = 10000,
regr.forest<- setHyperPars(regr.forest, size=2)

# cross-validatation 
#iter <- 10
rdesc = makeResampleDesc("LOO")
r_cv <- resample(regr.forest, filtered.task, rdesc, measures = list(rmse, timetrain))
r_cv$pred %>% print(n=100)

mod = train(regr.forest, filtered.task)
FI<- getFeatureImportance(mod)
FI$res %>% arrange(desc(importance)) %>% print(n=100)

#varimp(mod$learner.model) %>% sort() # varimp works with ctree
# library(cowplot)
# p1 <- plotROCCurves(df) + ggtitle(paste0(iter, "-fold CV ROC, AUC:", round(cur_auc,2)," RF all genes both cohorts (filtered ", n_genes, ") mtry: ", mtry))
# p2<- ggplot(as.data.frame(FI$res), aes(x=variable, y = importance)) + 
#   geom_bar(stat="identity")
# plot_grid(p1, p2, labels = "AUTO")
# ggsave(filename = here("RF_results", paste0(iter,"-fold_ROC_RF_VI_both_cohorts_model", model, "_filtered", n_genes, "_mtry", mtry, ".tiff")))



