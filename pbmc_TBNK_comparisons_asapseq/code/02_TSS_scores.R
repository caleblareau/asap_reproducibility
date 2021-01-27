library(BuenColors)

dt <- read.table("../data/three_comparison/TSS_scores_3channels.tsv", fill = 0, header = TRUE)
dt$pos <- 1:dim(dt)[1]

tb <- pretty_plot(fontsize = 8) + L_border() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 


p1 <- ggplot(dt, aes(x = pos, y = ASAP_A)) +
  geom_line() + tb + labs(y = "") + 
  scale_y_continuous(limits = c(0,12), expand = c(0,0), breaks = c(0, 5, 10))

p2 <- ggplot(dt, aes(x = pos, y = ASAP_B)) +
  geom_line() + tb + labs(y = "") + 
  scale_y_continuous(limits = c(0,12), expand = c(0,0), breaks = c(0, 5, 10))

p3 <- ggplot(dt, aes(x = pos, y = ASAP_C)) +
  geom_line() + tb + labs(y = "") + 
  scale_y_continuous(limits = c(0,12), expand = c(0,0), breaks = c(0, 5, 10))


ggsave(cowplot::plot_grid(p1, p2, p3, ncol = 1), 
       file = "../plots/TSS3.pdf", width = 1.3, height = 2.5)
