
###########################################################################
# survival random forest  -------------------------------------------------
###########################################################################
library(tidyverse)
library(here)

# prep data 
source(here("10h_both_cohorts_survival_dataprep.R"))


# make random forest model ------------------------------------------------

# make survival task 
surv.task <- makeSurvTask(data = data, target = c("PFS", "PFS_factor"))

# filter methods for selection of features 
# filter values compatible with survival task: variance, univariate.model.score, ranger_permutation,
# randomForestSRC_var.select, randomForestSRC_importance, praznik_MRMR, permutation.importance, party_cforest.importance, mrmr
fv <- generateFilterValuesData(surv.task,
                               method = "variance")

fv$data %>% arrange(desc(value))%>% head(n=100)

n_genes <- 32
filtered.task = filterFeatures(surv.task, method = "variance", 
                               abs = n_genes,
                               #threshold = 0.16,
                               mandatory.feat = c("neo_rank2", "T_cell_diversity")
                               ) #"T_cell_diversity", - from feature importance T-cell diversity doesnt seem to be relevant

mtry <- round(sqrt(filtered.task$task.desc$n.feat[1]))

# make learner ------------------------------------------------------------
# survival learners: surv.cforest surv.randomForestSRC surv.ranger surv.rpart

# Generate the learner
model <- "surv.randomForestSRC" #"regr.randomForestSRC" #"regr.randomForest" #regr.randomForestSRC" # "regr.extraTrees" #"regr.evtree" #"regr.ctree" regr.randomForest regr.ranger 
surv.forest = makeLearner(model, 
                          predict.type = "response",
                          id = "surv.rf.both.cohort", ntree=1000,
                          mtry=mtry)

# cross-validatation 
iter <- 10
rdesc = makeResampleDesc("CV")
r_cv <- resample(surv.forest, filtered.task, rdesc, measures = list(cindex.uno, iauc.uno, timetrain)) #cindex, cindex.uno, iauc.uno, ibrier, 
r_cv$pred %>% print(n=100)

mod = train(surv.forest, filtered.task)
FI<- getFeatureImportance(mod)
FI$res %>% arrange(desc(importance)) %>% print(n=100)

predictSurvProb(mod, newdata = data, times = data$PFS)

###########################################################################
# try not using mlr -------------------------------------------------------
###########################################################################
library(randomForestSRC)
# select features
sele.genes <- c("B2M", "EEF1A1", "HP", "SNORA73A", "MT_CO3", "T_cell_diversity",
"MIR663AHG", "RNU2_2P", "RNVU1_7",  "RN7SK", "SNORD17", "IGKC", 
"RMRP", "RNA5_8S5", "ITLN1", "RNU4_2", "SNORD3A", "RNU1_28P", "APOA2", 
"RNVU1_18", "RNU1_27P") 

sele.data <- data[,c(sele.genes, "PFS", "PFS_factor")]
# set.seed(123)
# train_ind <- sample(seq_len(nrow(sele.data)), size = nrow(sele.data)/2)
# train <-  data[train_ind,]
# test <- data[-train_ind,] 

require(caret)
folds <- 5
flds <- createFolds(sele.data$B2M, k = folds, list = TRUE, returnTrain = FALSE)


models <- list()
survival <- list()
predictions <- list()
for (i in 1:folds){
  test <- sele.data[flds[[i]],]
  train <- sele.data[-flds[[i]],]
  
  rfsrc <- rfsrc(Surv(time=PFS, event=PFS_factor) ~ .,
                 data = train,
                 type="right", 
                 ntree=10000, ntime = c( 50, 100, 200, 300)) 
  models[[i]] <- rfsrc

  pred <- predict(rfsrc, 
                   newdata = test)

  predictions[[i]] <- pred
  
  cur_sur <- as_tibble(pred$survival) 
  colnames(cur_sur) <- pred$time.interest
  cur_sur$index_data <- flds[[i]]
  cur_sur$cv <- i
  cur_sur$sample <- sample[flds[[i]]]
  cur_meta <- metadata %>% filter(patient %in% cur_sur$sample) %>% dplyr::select(patient, PFS, PFS_factor)
  cur_sur <- left_join(cur_sur, cur_meta, by = c("sample"="patient"))
  survival[[i]] <-cur_sur 
}

library(ggplot2)
input <- survival[[1]]%>% 
  mutate(PFS, PFS_pred = ifelse(PFS>93, "no event before 98", "event before"))

input%>% ggplot(.) + 
  geom_col(aes(sample, `93`,fill=PFS_pred)) + 
  ylim(0,1) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(y="probability of no. event") + 
  ggtitle("Prediction of probability of no event at day 98") 
  


# try holdout method ------------------------------------------------------
library(riskRegression)
library(survival)
percent_train <- 0.7
flds <- createDataPartition(data$B2M,  p = percent_train, list = TRUE)

train <- data[flds[[1]],]
test <- data[-flds[[1]],]

train_forest <- rfsrc(Surv(time=PFS, event=PFS_factor) ~ .,
                      data = train,
                      type="right", 
                      ntree=10000, 
                      ntime = c(1, 10, 50, 100, 150, 200, 300),
                      mtry=round(sqrt(ncol(train))),
                      importance = "none")


ff0 <- Score(list(train_forest), 
             formula=Surv(time=PFS, event=PFS_factor) ~ 1, 
             data=test, 
             times=c(1, 10, 50, 100, 150, 200, 300),
             plots=c("cali", "roc"))

plotROC(ff0, times = 200)

full_forest <- rfsrc(Surv(time=PFS, event=PFS_factor) ~ .,
                      data = data,
                      type="right", 
                      ntree=10000, 
                      ntime = c(1, 10, 50, 100, 150, 200, 300),
                      mtry=round(sqrt(ncol(sele.data))))









