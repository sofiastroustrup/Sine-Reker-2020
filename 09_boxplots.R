library(ggsignif)
library(ggbeeswarm)
library(ggrepel)
library(tidyverse)
library(ggplot2)

diff_genes<- read_csv(file = "diff_genes_p0.02_lfcShrinkage.csv")
load("clean_data.RData")

# select genes to plot  ---------------------------------------------------
cur_data <- assay(vsd)
cur_data <- cur_data[rownames(cur_data) %in% diff_genes$diff_genes,]
#cur_data <- txi$abundance[rownames(txi$abundance) %in% diff_genes$diff_genes,]
patient_id <- colnames(cur_data)
cur_data <- t(cur_data) 
cur_data <- as_tibble(cur_data)
cur_data$file_id<- patient_id

metadata<- s2c %>% dplyr::select(c(file_id, PFS_cut, neo_rank2))
data <- left_join(cur_data, metadata)
#data2plot_high.neo<- data %>% filter(neo_rank2>median(data$neo_rank2))



# Make boxplot for selected genes  ----------------------------------------

data %>% ggplot(., aes(x=PFS_cut, y= KHDC1L)) +
  geom_boxplot() + 
  geom_jitter(width = 0.2) + 
  geom_signif(comparisons = list(c('high_PFS', 'low_PFS')),
              na.rm = T,
              data = NULL,
              test = t.test)

# differentially expressed genes 
# "MAGEA12"   "CSAG2"     "KHDC1L"    "MAGEA3"    "CSAG3"     "MT1A"      "TCL1B"    
# [8] "KRT16P6"   "FXYD2"     "FOXI3"     "XAGE1A"    "PAEP"      "RASGEF1C"  "CSAG1"    
# [15] "SRRM4"     "BCHE"      "AFAP1-AS1" "ZIM2"      "SAA1"      "ZNF560"    "HAS1"     
# [22] "ADGRD1"    "IL6" 




peptide_scrrened_unik %>%
  ggplot(., aes( x = neo_response  , y = Mut_MHCrank_EL)) +
  geom_quasirandom(aes(color = neo_response)) +
  geom_boxplot(alpha = 0.5)  +
  scale_y_log10() +
  theme_bw() +
  #labs(x = "Clonal muatation", y = "Gene Expression level", color = "Cancer driver gene") +
  # geom_label_repel(data=peptide_scrrened_unik %>% filter(neo_response=="yes", Expression_Level>100,Cancer_Driver_Gene=="Yes"),
  #                  aes(label=paste(Gene_Symbol,identity)),
  #                  box.padding = unit(0.35, "lines"),
  #                  point.padding = unit(0.5, "lines"),
  #                  size = 3) +
  geom_signif(comparisons = list(c('yes', 'no')),
              na.rm = T,
              map_signif_level = T,
              data = peptide_scrrened_unik ,
              test = t.test)


require(gridExtra)
par(mfrow = c(3, 2))
plots <- list()
for (i in 1:2){
  print(i)
  print(diff_genes$diff_genes[[i]])
  p <- data_2plot %>% ggplot(., aes(x=PFS_cut, y=diff_genes$diff_genes[[i]])) +
    geom_boxplot(alpha = 0.5) + 
    geom_jitter(width = 0.2) +
    geom_signif(comparisons = list(c('high_PFS', 'low_PFS')),
                na.rm = T,
                map_signif_level = T,
                data = NULL,
                test = t.test)
  plots[[i]] <- p
}

multiplot(plotlist= plots)
grid.arrange(plots[[1]], plots[[2]], ncol=2)
do.call(grid.arrange,plots)

diff_genes$diff_genes[[i]]