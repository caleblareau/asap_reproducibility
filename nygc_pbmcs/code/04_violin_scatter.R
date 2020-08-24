library(data.table)
library(BuenColors)
library(Seurat)
library(dplyr)
library(stringr)
source("00_import_kite.R")

#load("../output/AvsB_LLL-Omni.rda")
sc_LLL <- fread("../data/AvsB/singlecell_LLL.csv.gz") %>% filter(cell_id != "None" & peak_region_cutsites > 2000)
sc_OMNI <- fread("../data/AvsB/singlecell_OMNI.csv.gz")%>% filter(cell_id != "None"& peak_region_cutsites > 2000)
dim(sc_OMNI)
dim(sc_LLL)

bc_LLL <- sc_LLL$barcode %>% as.character()
bc_OMNI <- sc_OMNI$barcode %>% as.character()

CLR_Seurat <- function(mat){
  pbmc <- CreateSeuratObject(counts = mat, assay = "ADT")
  pbmc <- NormalizeData(pbmc, assay = "ADT", normalization.method = "CLR")
  pbmc <- ScaleData(pbmc, assay="ADT")
  return(data.matrix(pbmc@assays$ADT@data))
}

# Split into experiment
LLL_kite_A <- import_kite_counts_NYGC("../data/AvsB/TSA_LLL/featurecounts/", "", FALSE, string_col = "TSA_")
OMNI_kite_A <- import_kite_counts_NYGC("../data/AvsB/TSA_OMNI/featurecounts/", "", FALSE, string_col = "TSA_")

LLL_kite_B <- import_kite_counts_NYGC("../data/AvsB/TSB_LLL/featurecounts/", "", FALSE, string_col = "TSB_")
OMNI_kite_B <- import_kite_counts_NYGC("../data/AvsB/TSB_OMNI/featurecounts/", "", FALSE, string_col = "TSB_")

LLL_mat <- t(CLR_Seurat(rbind(t(LLL_kite_A[bc_LLL,]),t(LLL_kite_B[bc_LLL,]))))
OMNI_mat <- t(CLR_Seurat(rbind(t(OMNI_kite_A[bc_OMNI,]),t(OMNI_kite_B[bc_OMNI,]))))

pal <- colorRampPalette(c("white", "dark green"))
greens_color_ramp <- colorRampPalette(c("white", greens = pal(10)))

plot_marker <- function(mat, marker){
  go_df <- data.frame(
    X = mat[,paste0("TSA-", marker)],
    Y = mat[,paste0("TSB-", marker)]
  )
  plot <- smoothScatter(go_df$X, go_df$Y, nbin = 328, colramp = greens_color_ramp, nrpoints = 0,
                        ret.selection = FALSE, xlab="", ylab="", xlim=c(0,4), ylim=c(0,4), labels = FALSE)
  return(plot)
}

markers <- c("CD3", "CD4", "CD8", "CD19", "CD16", "CD56", "CD11c", "CD14", "CD45")

# PLOT OMNI
par(mfrow=c(3,3))
par(mar=c(1,1,1,1))
for(i in markers){
  (plot_marker(OMNI_mat, i))
}

# PLOT LLL
par(mfrow=c(3,3))
par(mar=c(1,1,1,1))
for(i in markers){
  (plot_marker(LLL_mat, i))
}


# Make violin plots
# Function to clip above the 99th% tile to minimize the impact of extreme expressors
rm_above99 <- function(vec){
  q99 <- quantile(vec, 0.99)
  return(vec[vec < q99])
}

LLL_list <- lapply(colnames(LLL_mat), function(val) rm_above99(LLL_mat[,val])); names(LLL_list) <- paste0(colnames(LLL_mat), "-", "LLL")
OMNI_list <- lapply(colnames(OMNI_mat), function(val) rm_above99(OMNI_mat[,val])); names(OMNI_list) <- paste0(colnames(OMNI_mat), "-", "OMNI")

vals_df <- rbind(
  reshape2::melt(LLL_list),
  reshape2::melt(OMNI_list)
)
plot_df <- data.frame(str_split_fixed(vals_df[[2]], "-", 3), 
                      value = vals_df$value)

pA <- ggplot(plot_df %>% filter(X1 == "TSA"), aes(x = X2,  fill = X3, y = value)) +
  geom_violin(scale = "width", adjust = 1, width = 0.8, color = "black") + pretty_plot(fontsize = 7) + L_border()  +
  labs(x = "Antibody", y = "CLR TSA UBI counts") +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#6f94b9", "#e26868"))

pB <- ggplot(plot_df %>% filter(X1 == "TSB"), aes(x = X2,  fill = X3, y = value)) +
  geom_violin(scale = "width", adjust = 1, width = 0.8, color = "black") + pretty_plot(fontsize = 7) + L_border()  +
  labs(x = "Antibody", y = "CLR TSB UBI counts") +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#6f94b9", "#e26868"))

cowplot::ggsave2(pA, file = "../plots/TSA_LLL-OMNI.pdf", width = 2.8, height = 1.7)
cowplot::ggsave2(pB, file = "../plots/TSB_LLL-OMNI.pdf", width = 3, height = 2)


m_plot_df <- rbind(
  data.frame(
    pct_mito = rm_above99((sc_OMNI$mitochondrial)/(sc_OMNI$total))*100,
    condition = "OMNI"
  ),
  data.frame(
    pct_mito = rm_above99((sc_LLL$mitochondrial)/(sc_LLL$total))*100,
    condition = "LLL"
  )
)

pM <- ggplot(m_plot_df, aes(x = condition,  fill = condition, y = pct_mito)) +
  geom_violin(scale = "width", adjust = 1, width = 0.8, color = "black") + pretty_plot(fontsize = 7) + L_border()  +
  labs(x = "Lysis", y = "% mitochondrial fragments") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("#6f94b9", "#e26868"))
cowplot::ggsave2(pM, file = "../plots/mito_pct_LLL-OMNI.pdf", width = 2, height = 1.7)
  
  