### load packages ###
pacman::p_load(tidyverse, cowplot)

### load theme ###
source("scripts/plotting_theme.R")

### load data ####
cpg_windows_min2_sum <- read.csv(file = "results/wgbs/summary_cpg_meth_min2.csv")
load("results/wgbs/cpg_meth_min2.RData")

### plot ###
# assign colours

clr_4 <- c(clrs[1], clrs[2], clrs[8], clrs[13])

# mean TSS methylation (for y value label)

mean_tss_meth <- mean(cpg_windows_min2$methperc_mean[which(cpg_windows_min2$region == "TSS")], na.rm=T)      

ggplot(cpg_windows_min2_sum, aes(x = window_total, y = mean_meth)) + 
  geom_ribbon(aes(ymin = mean_meth - se_meth, ymax =mean_meth + se_meth), fill = clr_high, alpha = 0.6)+
  geom_line(col = "black", linewidth = 0.8) +
  geom_point(size=3) + labs(y = "Mean CpG methylation %",x  = "Window")+
  annotate("text", label="10 kb upstream", x = 20, y = 66, col = clr_high, size = 6)+
  annotate("text", label="Gene body", x = 62, y = 66, col = clr_high, size = 6)+
  annotate("text", label="10 kb downstream", x = 104, y = 66, col = clr_high, size = 6)+
  annotate("text", label="TSS", x = 51, y = mean_tss_meth, col = clr_gerp, size = 6)+ 
  geom_segment(aes(xend = 43, y = mean_tss_meth, x = 47, yend = mean_tss_meth), arrow = arrow(length = unit(0.2, "cm")))+
  geom_vline(xintercept = c(41.5, 42.5, 83.5), linetype = "dotted", col = "darkred", linewidth=1) -> fig_window

#ggsave(fig_window, file = "plots/fig_mean_cpg_meth.png", width=12, height=8)

# fig mean TSS methylation

tss_min2 <- subset(cpg_windows_min2, region == "TSS")
tss_min2$methperc_n <- as.numeric(tss_min2$methperc_n)

ggplot(tss_min2, aes(x = methperc_mean)) + geom_histogram(fill=clr_high, col = "black") +
  labs(x = "Mean CpG methylation % of TSS", y = "Count") -> fig_hist
#ggsave(fig_hist, file = "plots/fig_hist_TSS_meth.png", width=12, height=8)

### or just the histogram and across the genes

cowplot::plot_grid(fig_window, fig_hist, ncol = 1, align = "h", axis = "l", rel_heights=c(1,0.8),  labels="auto", label_fontface = "plain", label_size = 22) -> fig_s1

ggsave(fig_s1, file = "plots/final/supp/sfig_1_wgbs.png", width=12, height=12)
ggsave(fig_s1, file = "plots/final/supp/sfig_1_wgbs.pdf", width=12, height=12, device=cairo_pdf)
