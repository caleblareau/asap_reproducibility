library(data.table)
library(dplyr)
library(BuenColors)
source('../../global_functions/estimateLibraryComplexity.R')

# Function to process the per-barcode summary statistics
est_comp_tag_sm <- function(dt){
  vec <- c(rep("m", 46), rep("h", 46))
  names(vec) <- as.character(0:91)
  dff <- dt %>%
    mutate(Ab = vec[as.character(set)], umi_ab = paste0(Ab,"_",umi)) %>%
    group_by(barcode,Ab) %>%
    summarise(count = n_distinct(umi), n_total = n()) %>%
    mutate(n_duplicated = n_total - count) 
  dff$complexity <- sapply(1:(dim(dff)[1]), function(i){
    estimateLibrarySize(dff$n_total[i],dff$count[i])
  })
  dff
}

assign_organism <- function(df){
  df$pct_human_atac <- (df$peak_region_fragments_GRCh38)/(df$peak_region_fragments_GRCh38 + df$peak_region_fragments_mm10) * 100

  pha <- df$pct_human_atac
  dha <- df$peak_region_fragments_GRCh38 > 100
  dma <- df$peak_region_fragments_mm10 > 100

  df$assignment_atac <- case_when(
    dha & (pha > 95) ~ "human",
    dma & (pha < 5) ~ "mouse",
    dha & dma ~ "mixed",
    TRUE ~ "none"
  )
  
  df
}

#Import barcodes
cells_df1 <- fread("../data/singlecell1.csv.gz") %>% filter(cell_id != "None")
cells_df2 <- fread("../data/singlecell2.csv.gz") %>% filter(cell_id != "None")

cells_df1$barcode <- gsub("-1","",cells_df1[["barcode"]])
cells_df2$barcode <- gsub("-1","",cells_df2[["barcode"]])

# Assign organism
assign_df1 <- assign_organism(cells_df1)
assign_df2 <- assign_organism(cells_df2)

adt1 <- fread("../../../asap_large_data_files/species_mix/preSPRI_busPT.tsv",col.names = c("barcode", "umi", "set", "idx")) %>%
  dplyr::filter(barcode %in% cells_df1$barcode)
adt2 <- fread("../../../asap_large_data_files/species_mix/postSPRI_busPT.tsv",col.names = c("barcode", "umi", "set", "idx")) %>%
  dplyr::filter(barcode %in% cells_df2$barcode)

# Get per Ab complexity and pull assignments from ATAC
adt1comp <- est_comp_tag_sm(adt1)
adt2comp <- est_comp_tag_sm(adt2)
adt1comp_h <- adt1comp %>% filter(barcode %in% c(assign_df1 %>% filter(assignment_atac == "human") %>% pull(barcode)) & Ab == "h")
adt1comp_m <- adt1comp %>% filter(barcode %in% c(assign_df1 %>% filter(assignment_atac == "mouse") %>% pull(barcode)) & Ab == "m")
adt2comp_h <- adt2comp %>% filter(barcode %in% c(assign_df2 %>% filter(assignment_atac == "human") %>% pull(barcode)) & Ab == "h")
adt2comp_m <- adt2comp %>% filter(barcode %in% c(assign_df2 %>% filter(assignment_atac == "mouse") %>% pull(barcode)) & Ab == "m")

adt1comp_h$Experiment <- "apre-SPRI"; adt1comp_h$species <- "Human"
adt1comp_m$Experiment <- "apre-SPRI"; adt1comp_m$species <- "Mouse"
adt2comp_h$Experiment <- "post-SPRI"; adt2comp_h$species <- "Human"
adt2comp_m$Experiment <- "post-SPRI"; adt2comp_m$species <- "Mouse"

df <- rbind(adt1comp_h,adt1comp_m,adt2comp_h,adt2comp_m)
p1 <- ggplot(df, aes(x = species, y = log10(complexity), color = Experiment)) +
  geom_boxplot(outlier.shape = NA) + 
  pretty_plot(fontsize = 8) + L_border() +
  scale_y_continuous(limits = c(1,4.5)) +
  scale_color_manual(values = c("dodgerblue4", "green4")) +
  labs(x = "", y = "Tag Complexity", color = "") +
  theme(legend.position = "none")

df%>% group_by(species, Experiment) %>% summarize(median(complexity))

ggsave(p1, 
       file = "../plots/qc_tagComplexity_SPRI_comparison.pdf", width = 1.5, height = 1.5)

#---------------------------------------------------------
# Function that summarizes the number of cell barcodes
# associated with each UMI
make_count_df <- function(df, total = 1400000){
  set.seed(1)
  idx_keep1 <- sample(1:dim(df)[1], total, replace = FALSE)
  count_df <- df[idx_keep1,] %>% group_by(V1, V2) %>%
    summarize(count=n()) %>% ungroup() %>%
    group_by(V2) %>% summarize(n_cb = n())
  return(count_df)
}


adt1 <- adt1 %>% group_by(V1, V2) %>% summarize(count = n())
adt2 <- adt2 %>% group_by(V1, V2) %>% summarize(count = n())

count_df_adt1 <- make_count_df(adt1)
count_df_adt2 <- make_count_df(adt2)

ecdf_df <- data.frame(
  n_barcodes = c(count_df_adt1$n_cb, count_df_adt2$n_cb),
  library = c(rep("apre-SPRI", dim(count_df_adt1)[1]),
              rep("post-SPRI", dim(count_df_adt2)[1]))
)

ggplot(ecdf_df %>% dplyr::filter(n_barcodes < 10000), aes(x = n_barcodes, color = library)) +
  stat_ecdf() + theme(legend.position = "bottom") + scale_x_log10() +
  scale_color_manual(values = c("dodgerblue3", "firebrick")) +
  labs(x = "# Barcodes / UMI", y = "Empirical cumulative distribution") +
  pretty_plot() + L_border()

bar_df <- ecdf_df %>% group_by(library, n_barcodes) %>% summarize(count = n()) %>% ungroup() %>%
  group_by(library) %>% mutate(prop = count / sum(count))

bar_df %>% filter(n_barcodes < 6) %>%
  mutate(n_barcodes = as.character(n_barcodes)) -> bdf2

leftover <- (1 - bdf2 %>% group_by(library) %>% summarize(tt = sum(prop)) %>% pull(tt) %>% rev)

bdf2$library <- as.character(bdf2$library)
bar_plot_df <- rbind(data.frame(bdf2[,c(1,2,4)], stringsAsFactors = FALSE),
                     data.frame(
                       library = c("apre-SPRI", "post-SPRI"),
                       n_barcodes = c(">5", ">5"),
                       prop = leftover, stringsAsFactors = FALSE
                     ))
bar_plot_df$n_barcodes <- factor(bar_plot_df$n_barcodes, levels = c("1", "2", "3", "4", "5", ">5"))

p1 <- ggplot(bar_plot_df, aes(x = n_barcodes, y = prop*100, fill = library)) +
  geom_bar(stat="identity", position = position_dodge2(), color = "black") + theme(legend.position = "bottom") + 
  labs(x = "# Barcodes / UBI", y = "% of UBIs", fill = "") +
  pretty_plot(fontsize = 8) + L_border() +
  theme(legend.position = "none") + scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("dodgerblue4", "green4")) +
  theme(legend.position = c(0.8, 0.8))
cowplot::ggsave2(p1, file = "../plots/UBI_dist.pdf", width = 1.8, height = 1.8)

