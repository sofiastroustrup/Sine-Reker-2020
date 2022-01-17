
###########################################################################
# Random forest for Rigshospitalet and new york cohort --------------------
###########################################################################

# load packages ----------------------------------------------------------
library(tidyverse)
library(mlr)
library(party)
library(here)
library(magrittr)
library(survival)
# load data ---------------------------------------------------------------
load(here("clean_data.RData"))
source(here("00_functions.R"))
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
#NY_spread<-readRDS(file=here("tpm_data", "tpm_bladder_cancer_spread.RDS"))
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
clinical_sele <- clinical_data %>% filter(patient %in% c(tpm$sample)) 
clinical_sele <- clinical_sele %>% filter(!is.na(PFS_cut)) 
clinical_sele <- clinical_sele %>% mutate(PFS, PFS_cut = ifelse(PFS>median(PFS), "high_PFS", "low_PFS"))

metadata <- clinical_sele %>% dplyr::select("PFS_cut", 
                                            "patient",
                                            "T_cell_diversity",
                                            "neo_rank2")

# change type of specific columns to character 
metadata[,1:3]%<>% mutate_if(is.character,as.factor)

# add extra columns; neoepitope, cancer_type, treatment, T_cell_diversity
data <- left_join(tpm, metadata, by=c("sample"="patient"))
colnames(data)<- colnames(data) %>% str_replace_all("-", "_") 
data <- as_data_frame(data)

# make categorical data to factor and remove 
sample <- data$sample
data$sample <- NULL
data <- data %>% filter(!is.na(PFS_cut))
data <- as_data_frame(data)




##########################################################################
# make RF model both with feature selection ------------------------------
##########################################################################

fv_method <- c("FSelectorRcpp_information.gain", # random forest methods very slow
                 "randomForestSRC_var.select", 
                 "randomForestSRC_importance",
                 "randomForest_importance", 
               "FSelectorRcpp_symmetrical.uncertainty")

trees <- c("classif.h2o.randomForest", "classif.randomForest", "classif.randomForestSRC", "classif.ranger ", "classif.cforest")

# generate the task 
task = makeClassifTask(data = data, target = "PFS_cut")

# filter methods for selection of features 
n_genes <- 5
fv <- generateFilterValuesData(task,
                         method = "FSelectorRcpp_information.gain")
fv$data %>% arrange(desc(value))%>% head(n=100)

filtered.task = filterFeatures(task, method = "FSelectorRcpp_information.gain", 
                               abs = n_genes,
                               #threshold = 0.16,
                               mandatory.feat = c("neo_rank2")) #"T_cell_diversity", - from feature importance T-cell diversity doesnt seem to be relevant

# Generate the learner
model <- "classif.randomForestSRC"
classif.lrn = makeLearner(model, 
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
plotROCCurves(df) + ggtitle(paste0(iter, "-fold CV ROC, AUC:", round(cur_auc,2)," RF all genes both cohorts (filtered ", n_genes, ") mtry: ", mtry))
ggsave(filename = here("RF_results", paste0(iter,"-fold_ROC_RF_all_genes_both_cohorts_model", model, "_filtered", n_genes, "_mtry", mtry, ".tiff")))


mod = train(classif.lrn, filtered.task)
FI<- getFeatureImportance(mod)
FI$res
#varimp(mod$learner.model) %>% sort() # varimp works with ctree
library(cowplot)
p1 <- plotROCCurves(df) + ggtitle(paste0(iter, "-fold CV ROC, AUC:", round(cur_auc,2)," RF all genes both cohorts (filtered ", n_genes, ") mtry: ", mtry))
p2<- ggplot(as.data.frame(FI$res), aes(x=variable, y = importance)) + 
  geom_bar(stat="identity")
plot_grid(p1, p2, labels = "AUTO")
ggsave(filename = here("RF_results", paste0(iter,"-fold_ROC_RF_VI_both_cohorts_model", model, "_filtered", n_genes, "_mtry", mtry, ".tiff")))



###########################################################################
# build RF without feature selection --------------------------------------
###########################################################################
# generate the task 
task = makeClassifTask(data = data, target = "PFS_cut")

# Generate the learner
model <- "classif.randomForestSRC"
classif.lrn = makeLearner(model, 
                          predict.type = "prob", 
                          fix.factors.prediction = TRUE, 
                          id = "RH_ranfomforest_all")

# set hyperparameters for random forest 
mtry <- 500
classif.lrn <- setHyperPars(classif.lrn, ntree = 10000, mtry=mtry)

# cross-validatation 
iter <- 5
rdesc = makeResampleDesc("CV", iters = iter)
r_cv <- resample(classif.lrn, task, rdesc, measures = list(acc, fn, fp,
                                                                    mmce, 
                                                                    ber, auc, 
                                                                    timetrain))
cur_auc <- performance(r_cv$pred, auc)
df = generateThreshVsPerfData(r_cv$pred, measures = list(fpr, tpr, mmce))
plotROCCurves(df) + ggtitle(paste0(iter, "-fold CV ROC, AUC:", round(cur_auc,2)," RF all genes both cohorts mtry: ", mtry))
ggsave(filename = here("RF_results", paste0(iter,"-fold_ROC_RF_all_genes_both_cohorts2_model", model, "_not_filtered_mtry", mtry, ".tiff")))


mod = train(classif.lrn, task)
FI<- getFeatureImportance(mod)
FI$res

# boxplot analysis --------------------------------------------------------
library(ggsignif)

p1<- data %>% ggplot(., aes(x=PFS_cut, y=FI$res$variable[1])) +
  geom_violin() + 
  geom_jitter(width = 0.2) + 
  geom_signif(comparisons = list(c('high_PFS', 'low_PFS')),
              na.rm = T,
              data = NULL,
              test = t.test)+ ggtitle("violin plot tpm expression of IFITM2")

p2<- data %>% ggplot(., aes(x=PFS_cut, y= FI$res$variable[2])) +
  geom_violin() + 
  geom_jitter(width = 0.2) + 
  geom_signif(comparisons = list(c('high_PFS', 'low_PFS')),
              na.rm = T,
              data = NULL,
              test = t.test)+ ggtitle("violin plot tpm expression of KAT5")

p3<- data %>% ggplot(., aes(x=PFS_cut, y= FI$res$variable[3])) +
  geom_violin() + 
  geom_jitter(width = 0.2) + 
  geom_signif(comparisons = list(c('high_PFS', 'low_PFS')),
              na.rm = T,
              data = NULL,
              test = t.test)+ ggtitle("violin plot tpm expression of NUMA1")

p4<- data %>% ggplot(., aes(x=PFS_cut, y= FI$res$variable[4])) +
  geom_violin() + 
  geom_jitter(width = 0.2) + 
  geom_signif(comparisons = list(c('high_PFS', 'low_PFS')),
              na.rm = T,
              data = NULL,
              test = t.test)+ ggtitle("violin plot tpm expression of neo_rank2")

p5<- data %>% ggplot(., aes(x=PFS_cut, y= FI$res$variable[5])) +
  geom_violin() + 
  geom_jitter(width = 0.2) + 
  geom_signif(comparisons = list(c('high_PFS', 'low_PFS')),
              na.rm = T,
              data = NULL,
              test = t.test)+ ggtitle("violin plot tpm expression of YIPF2")

library(cowplot)


plot_grid(p1, p2, p3, p4, p5, labels = "AUTO")
ggsave(filename = here("RF_results", paste0("gene_selected_violinplot_", model, "_filtered", n_genes, "_mtry", mtry, ".tiff")))



summary_cohort<- clinical_sele %>% dplyr::select(c(cancer_type, treatment, 
                                  treatment_target, cohort)) %>% 
  group_by(cohort) %>% summarise(no_patients = n(),
                                no.dif.cancer = n_distinct(cancer_type),
                                cancer_types = paste(unique(cancer_type), collapse = ", "), 
                                n_treatment = n_distinct(treatment),
                                treatment = paste(unique(treatment), collapse = ", "),
                                treatment_target = paste(unique(treatment_target), collapse=", "))


write_delim(x = summary_cohort, path=here("summary_cohort.txt"), delim = "\t")

###########################################################################
# make rf model with DE genes  --------------------------------------------
###########################################################################

  




