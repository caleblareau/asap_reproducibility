library(data.table)
library(dplyr)
source("../../global_functions/estimateLibraryComplexity.R")

pbmc <- fread("../data/broad_pbmcs_aggr_singlecell.csv.gz") %>%
  filter(cell_id != "None")
pbmc$channel <- substr(pbmc$barcode, 18, 18)
pbmc <- pbmc %>% filter(channel %in% c("1", "2")) %>% mutate(total = duplicate + passed_filters) %>% data.frame()
pbmc$LibraryComplexity <- sapply(1:dim(pbmc)[1], function(i){
  estimateLibrarySize(pbmc[i,"total"],pbmc[i,"passed_filters"])})

top_1k <- pbmc  %>% group_by(channel) %>% top_n(5000, LibraryComplexity)
top_1k %>% group_by(channel) %>% summarize(median(LibraryComplexity))
p1 <- ggplot(top_1k, aes(x = channel, y = log10(LibraryComplexity), color = channel)) +
  geom_boxplot(outlier.shape = NA) + 
  pretty_plot(fontsize = 8) + L_border() +
  scale_y_continuous(limits = c(3.5, 5.0)) +
  scale_color_manual(values = c("dodgerblue2", "purple2")) +
  labs(x = "", y = "ATAC Complexity", color = "") +
  theme(legend.position = "bottom")
cowplot::ggsave2(p1, file = "../output/atac_complexity_broad.pdf", 
                 width = 1.5, height = 2)

