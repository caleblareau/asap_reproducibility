library(data.table)
library(BuenColors)

dat <- fread("../output/13june2020_citeseq_totalVI_output.csv.gz")
dat$batch <- substr(dat$barcode, 18,18)
dat$total_tags <- rowSums(data.matrix(dat[,seq(28, 472, 2)-1]))
  
ggplot(dat %>% arrange((total_tags)), aes(x = UMAP1, y = UMAP2, color = total_tags > 25000)) + 
  geom_point() 

ggplot(dat %>% arrange((total_tags)), aes(x = UMAP1, y = UMAP2, color = batch)) + 
  geom_point() 
