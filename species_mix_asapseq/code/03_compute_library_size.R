library(dplyr)
library(data.table)
library(BuenColors)
library(ggrastr)
library(SummarizedExperiment)
library(Matrix)

source("../../global_functions/estimateLibraryComplexity.R")

importQCquick <- function(qc_sc, name, n = 1000){
  
  dt <- fread(qc_sc) %>%
    filter(cell_id != "None") %>%
    mutate(s1 = substr(cell_id, 1,1)) %>% mutate(species = ifelse(s1 == "m", "Mouse", ifelse(s1 == "G", "Human", "d"))) %>% 
    mutate(Experiment = name) %>%
    mutate(DNaseProp = DNase_sensitive_region_fragments/passed_filters) %>%
    mutate(MitoProp = mitochondrial/total) %>% 
    mutate(duplicateProp = duplicate / total) %>%
    mutate(TSSProp = TSS_fragments/passed_filters) %>% data.frame()
  dt$LibraryComplexity <- sapply(1:dim(dt)[1], function(i){
    estimateLibrarySize(dt[i,"total"],dt[i,"passed_filters"])})
  
  top_1k <- dt %>% group_by(species) %>% top_n(n, LibraryComplexity)
  
  return(top_1k)
}

df <- rbind(importQCquick("../data/singlecell1.csv.gz", "apre-SPRI"),
            importQCquick("../data/singlecell2.csv.gz", "post-SPRI")
)


p1 <- ggplot(df, aes(x = species, y = log10(LibraryComplexity), color = Experiment)) +
  geom_boxplot(outlier.shape = NA) + 
  pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = c("dodgerblue4", "green4")) +
  labs(x = "", y = "Chromatin Library Complexity", color = "")

p2 <- ggplot(df, aes(x = species, y = TSSProp * 100, color = Experiment)) +
  geom_boxplot(outlier.shape = NA)  +
  pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = c("dodgerblue4", "green4")) +
  labs(x = "", y = "% Reads near TSS") + theme(legend.position = "none")

ggsave(cowplot::plot_grid(p1, p2, nrow =1, rel_widths = c(1.5, 1)), 
       file = "../plots/qc_SPRI_comparison.pdf", width = 5, height = 1.8)


