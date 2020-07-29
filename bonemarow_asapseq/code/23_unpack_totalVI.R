library(data.table)
library(BuenColors)
library(dplyr)

dat <- fread("../output/totalVI_output_bonemarrow_asap.csv.gz")

ggplot(dat %>% arrange((total_tags)), aes(x = UMAP1, y = UMAP2, color = scVI1)) + 
  geom_point() 
