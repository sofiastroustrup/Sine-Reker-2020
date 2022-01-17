# load packages --------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(tximport)
library(here)
library(pheatmap)
library(genefilter)
library(limma)
library(apeglm)
# load data ---------------------------------------------------------------
load(here("old_objects", "analysis_merge_data.Rdata"))
clinical_data <- read_tsv(here("clinical_data.txt"))
source(file = here("00_functions.R"))
load("~/Desktop/22107 Research Immersion/biomarkers/DATA_FRAME.Rdata")
load("~/Desktop/22107 Research Immersion/biomarkers/patien_data_RH_pro.Rdata")

###########################################################################
# clean up data -----------------------------------------------------------
###########################################################################

# add missing data --------------------------------------------------------

# add missing clinical data 
missing_clinical_data <- read_tsv("missing_clinical_data.txt")
missing <- missing_clinical_data %>% dplyr::select(c("patient", "RECIST", 
                                              "PFS", "OS", "OS_factor", 
                                              "PFS_factor",
                                              "sample", "PFS_cut"))

# add missing data, but only the rows with values 
missing <- missing %>% filter(!is.na(RECIST)| !is.na(PFS) | !is.na(OS) | 
                                !is.na(OS_factor) | 
                     !is.na(PFS_factor) | !is.na(sample) | !is.na(PFS_cut))

clinical_data <- rbind(clinical_data, missing)

# remove A from clinical data
clinical_data$patient <- clinical_data$patient %>% gsub(pattern = "A$", replacement = "")

# add disease, treatment, tcell diversity, neoepitopes
# DATA_FRAME cohort
DATA_FRAME$Sample <- ifelse(DATA_FRAME$cancer_type %in% c("Melanoma","Bladder_cancer"), DATA_FRAME$Sample,paste0("Fase",DATA_FRAME$Sample))
NY <- DATA_FRAME %>% filter(str_detect(Sample, pattern = "^M", negate = TRUE))
NY <- NY %>% dplyr::select(c("neo_rank2", "T_cell_diversity", "cancer_type", "Treatment", "Sample")) %>%
  rename("Sample"="sample") %>% rename("Treatment"="treatment")
NY <- NY %>% mutate(treatment_target = simplify_treatment(treatment = treatment))#%>% 

# rigshospitalet cohort
RH <- patien_data %>% filter(str_detect(sample, pattern = "B$", negate = TRUE)) #%>%
RH$sample <- RH$sample%>% gsub(pattern = "A$", replacement = "")
RH <- RH %>% 
  dplyr::select(c("neo_epi_rank2", "T_cell_diversity", "diagnose", "treatment", "sample")) %>% 
  rename("neo_epi_rank2"="neo_rank2")  %>% 
  rename("diagnose"="cancer_type")
RH$sample <- RH$sample %>% tolower()%>% gsub(pattern = "^p", replacement = "P")
RH<-RH %>% mutate(treatment_target = simplify_treatment(treatment = treatment))


# combine the two cohorts with extra info
comb_cohort <- rbind(NY, RH)


# remove samples that ends with B
s2c <- dplyr::slice(s2c, -(s2c$patient %>% grep(pattern = "B$")))
# remove samples that ends with T2
s2c <- dplyr::slice(s2c, -(s2c$patient %>% grep(pattern = "T2$")))
# remove "A" from sample name
s2c$patient <- s2c$patient %>% gsub(pattern = "A$", replacement = "")
# remove T1 from sample name 
s2c$patient <- s2c$patient %>% gsub(pattern = "T1$", replacement = "")

# remove "A" from file_id
s2c$file_id <- s2c$file_id %>% gsub(pattern = "A_", replacement = "_")


# add missing values (extra knowledge obtaiend from annie)
s2c[s2c$patient=="Patien22",]$RECIST <- "PD"
s2c[s2c$patient=="Patien24",]$RECIST <- "PR"
s2c[s2c$patient=="Patien27",]$RECIST <- "PR"
s2c[s2c$patient=="Patien25",]$RECIST  <- "PD"
s2c[s2c$patient=="Patien31",]$RECIST <- "PD"
s2c[s2c$patient=="Patien33",]$RECIST <- "PD"


# add info to clinical data and s2c
clinical_data <- left_join(clinical_data, comb_cohort, by = c("patient"="sample"), keep=FALSE )
s2c <- left_join(s2c, comb_cohort, by = c("patient" = "sample"), keep = FALSE)


# prep txi list sample names (sequencing metadata)
# clean all the dataframe in txi
txi_clean <- lapply(txi[1:3], clean_up_txi)
txi[1:3] <- txi_clean
# get sample names from sequencing data 
abundance_patients <- txi$abundance %>% colnames() 
# remove "_kalisto.txt from sample name
abundance_patients <- sub("_.*", "", abundance_patients) 

# s2c not in clinical 
s2c_not_in_clincal<- s2c[!s2c$patient %in% clinical_data$patient,]
#saveRDS(object = s2c_not_in_clincal, file = "s2c_not_in_clincal2.RDS")
write_csv(x = s2c_not_in_clincal, path = here("s2c_not_in_clincal2.csv"))

# select rows in s2c that are contained in clinical_data file
s2c_filtered <- s2c[s2c$patient %in% clinical_data$patient,]

# check that output is contained in the sequencing files  
s2c_filtered$patient %in% abundance_patients 


# convert months to days 
s2c_filtered <- s2c_filtered %>% mutate(PFS, PFS = ifelse(str_detect(string = patient, pattern = "^[0-9]"), PFS*30.4368499, PFS))
#s2c_filtered <- s2c_filtered %>% mutate(OS, OS = ifelse(is_digit(OS), OS*30.4368499, OS))
s2c_filtered <- s2c_filtered %>% mutate(OS, OS = ifelse(str_detect(string=patient, pattern = "^[0-9]"), OS*30.4368499, OS))


# OBS! this can be changed in we want to use another condition for analysis------------
# remove patients with na in both PFS and RECIST
#s2c_filtered <- s2c_filtered %>% filter(!is.na(PFS) | !is.na(RECIST)) 

# add column specifying cohort  -------------------------------------------
clinical_data <- clinical_data %>% mutate(patient, cohort = ifelse(str_detect(string = patient, pattern = "^[0-9]"), "new_york", "rigshospitalet"))
s2c_filtered <- s2c_filtered %>% mutate(patient, cohort = ifelse(str_detect(string = patient, pattern = "^[0-9]"), "new_york", "rigshospitalet"))


# add recist split --------------------------------------------------------
clinical_data <- clinical_data %>% mutate(RECIST, RECIST_cut = ifelse(RECIST=="PD", "NR", "R"))
s2c_filtered <- s2c_filtered %>% mutate(RECIST, RECIST_cut = ifelse(RECIST=="PD", "NR", "R"))

# make new PFS cut off
s2c_filtered <- s2c_filtered %>% mutate(PFS, PFS_cut = ifelse(PFS < median(PFS, na.rm = TRUE), "low_PFS", "high_PFS"))

# remove the rows with na in PFS column
s2c_filtered_PFS <- s2c_filtered %>% filter(!is.na(PFS))
s2c_filtered_RECIST <- s2c_filtered %>% filter(!is.na(RECIST))

# keep only those columns in txi that are present in s2c_filtered ------------------
txi_filtered_PFS <- lapply(txi[1:3],function(x) keep_columns(df1=x, s2c_filtered = s2c_filtered_PFS))
txi_filtered_RECIST <- lapply(txi[1:3],function(x) keep_columns(df1=x, s2c_filtered = s2c_filtered_RECIST))

txi_filtered_PFS$countsFromAbundance <- txi$countsFromAbundance
txi_filtered_RECIST$countsFromAbundance <- txi$countsFromAbundance



# export all  -------------------------------------------------------------
s2c <- s2c_filtered
save(list = c("s2c", "clinical_data", "txi"), file = "clean_data_nofilter.RData")



# export cleaned data for PFS cut -----------------------------------------------------
s2c <- s2c_filtered_PFS
txi <- txi_filtered_PFS
save(list = c("s2c", "clinical_data", "txi"), file = "clean_data.RData")


# export cleaned data for RECIST cut --------------------------------------
s2c <- s2c_filtered_RECIST
txi <- txi_filtered_RECIST
save(list = c("s2c", "clinical_data", "txi"), file = "clean_data_RECISTcut.RData")

