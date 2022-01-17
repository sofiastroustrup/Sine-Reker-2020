
# load packages -----------------------------------------------------------
library(tidyverse)
library(here)

# define functions to use in other r scripts ------------------------------

clean_up_txi <- function(df, ends = "_kallisto.txt"){
  # remove columns 
  index_remove<- grep(x= df %>% colnames(), pattern = paste0("B", ends, "$"))
  df <- df[,-(index_remove)]
  index_remove <- grep(x= df %>% colnames(), pattern = paste0("T2", ends, "$"))
  df <- df[, -(index_remove)]
  
  # remove T1 and A suffix
  colnames(df) <- df %>% colnames() %>% gsub(pattern = paste0("A", ends, "$"), replacement = ends)
  colnames(df) <- df %>% colnames() %>% gsub(pattern = paste0("T1", ends, "$"), replacement = ends)
  
  return(df)
}

is_digit <- function(x){
  return(str_detect(pattern = "[0-9]", string = x))
}


# filter cohort 

select_cohort <- function(df, cohort="rigshospitalet"){
  out <- NULL
  cnames <- colnames(df)
  if(cohort == "rigshospitalet"){
    keep <- cnames %>% str_detect(pattern = "^[0-9]", negate = TRUE)
    out <- df[,keep]
  }
  if(cohort== "new_york"){
    keep <- cnames %>% str_detect(pattern = "^[0-9]", negate = FALSE)
    out <- df[,keep]
  }
  if(is.null(out)) {print("not a valid cohort")}
  
  return(out)
}


# function for keeping the right columns in txi ---------------------------
keep_columns <- function(df1, s2c_filtered){
  to_keep <- s2c_filtered$file_id
  cnames <- colnames(df1)
  out <- df1[,cnames %in% to_keep]
  return(out)
}


simplify_treatment <- Vectorize(function(treatment){
  if (str_detect(treatment, pattern="atezolizumab")){return("PD-L1")}
  if (str_detect(treatment, pattern="nivolumab" )){return("PD1")}
  if(str_detect(treatment, pattern="prembolizumab")){return("PD1")}
  if(str_detect(treatment, pattern ="ipi+nivo")){return("CTLA-4 + PD1")}
  if(treatment=="ipi+nivo"){return("CTLA-4 + PD1")}
})


