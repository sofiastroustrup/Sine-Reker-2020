
# clean tpm data ----------------------------------------------------------
library(here)
library(tidyr)
library(tidyverse)
library(dplyr)

# prepare RH data  --------------------------------------------------------

# load data 
# fase 
load(here("tpm_data", "RH_data_spread.Rdata")) 
fase <- RNA_dat %>% as_tibble()
fase <- fase %>% mutate(Sample, Sample =ifelse(str_detect(Sample, "^[0-9]"), 
                               paste0("Fase", Sample), Sample)) 
fase$Sample <- fase$Sample%>% gsub(pattern = "T1$", replacement = "")
fase <- fase %>% filter(str_detect(fase$Sample, pattern = "T2$", negate = TRUE))

#fase$Sample <- fase$Sample%>% gsub(pattern = "T2$", replacement = "")

# patien 
#rh_patien <- read.delim(here("tpm_data", "all_rna_RH.txt"), header=TRUE)
load(here("tpm_data", "rna_spread.Rdata"))
patien_data <- RNA_spread %>% as_tibble()

# remove B observations 
# remove "A" form patient id 
RH <- patien_data %>% filter(str_detect(Sample, pattern = "B$", negate = TRUE)) 
RH$Sample <- RH$Sample%>% gsub(pattern = "A$", replacement = "")



# merge RH data and export ------------------------------------------------
# merge
rh_tpm<- rbind(fase, RH)

library(magrittr)
rh_tpm[,2:ncol(rh_tpm)]%<>% mutate_if(is.character,as.numeric)

# export
saveRDS(rh_tpm, file="rh_tpm.RDS")



# bladder cancer ----------------------------------------------------------

# read bladder data
NY<- read_delim(file = here("tpm_data", "all_rna_seq_bladder.txt"), delim="\t")
NY <- NY[-c(1),]
NY_spread <- spread(data = as_data_frame(NY) ,key=hugo_symbol, value=mean_exp)
# check all bladder samples are in clinical 
NY_spread<- NY_spread %>% filter(Sample %in% clinical_data$patient)
NY_sample <- NY_spread$Sample

NY_spread<- NY_spread %>% select(-c(Sample, hugo_symbol))
# make into numeric instead of character 

library(magrittr)
NY_spread %<>% mutate_if(is.character,as.numeric)
NY_spread$sample <- NY_sample
saveRDS(object = NY_spread, file = here("tpm_data", "tpm_bladder_cancer_spread2.RDS"))
