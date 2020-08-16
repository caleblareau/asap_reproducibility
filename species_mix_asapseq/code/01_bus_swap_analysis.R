library(data.table)
library(dplyr)
library(BuenColors)

cells_df1 <- fread("singlecell1.csv") %>% filter(cell_id != "None")
cells_df2 <- fread("singlecell2.csv") %>% filter(cell_id != "None")

cells_df1$barcode <- gsub("-1","",cells_df1[["barcode"]])
cells_df2$barcode <- gsub("-1","",cells_df2[["barcode"]])

adt1 <- fread("o1.tsv") %>% dplyr::filter(V1 %in% cells_df1$barcode)
adt2 <- fread("o2.tsv") %>% dplyr::filter(V1 %in% cells_df2$barcode)

adt1 <- adt1 %>% group_by(V1, V2) %>% summarize(count = n())
adt2 <- adt2 %>% group_by(V1, V2) %>% summarize(count = n())

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

count_df_adt1 <- make_count_df(adt1)
count_df_adt2 <- make_count_df(adt2)
#count_df_perm1 <- make_count_df(perm1)
#count_df_perm2 <- make_count_df(perm2)

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
cowplot::ggsave2(p1, file = "plots/UBI_dist.pdf", width = 1.8, height = 1.8)

