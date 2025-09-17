### packages ###
pacman::p_load(tidyverse, cowplot, data.table)

#### Theme ####
source("scripts/plotting_theme.R")

summary <- read.csv(file = "results/qc/summary_coverage_meth_unfiltered.csv")

### Plot distributions of mean methylation etc. ####

ggplot(summary, aes(mean_cov)) + geom_histogram(col = "black", fill = clr_grey) + 
  labs(x = "Mean coverage", y="Count") +
  scale_y_continuous(labels = scales::unit_format(unit = "K", scale = 1e-3)) -> hist_mean_cov

ggplot(summary, aes(n)) + geom_histogram(col = "black", fill = clr_grey) + 
  labs(x = "Mean sample size", y="Count") +
  scale_y_continuous(labels = scales::unit_format(unit = "K", scale = 1e-3)) -> hist_n

### PCA ####
load(file = "results/pca/quality_check_pca_meanmeth.RData")

# plot PCs
ggplot(merge_pca, aes(x = pc1, y = pc2)) + geom_point(size=3, aes(col = site)) + 
  labs(x = "PC 1", y = "PC 2", col = "Lek") +
  scale_color_manual(values=clrs_hunting) -> pca_site

merge_pca$prepost <- factor(merge_pca$prepost, levels = c("pre", "post"))
ggplot(merge_pca, aes(x = pc1, y = pc2)) + geom_point(size=3, aes(col = prepost)) + 
  labs(x = "PC 1", y = "PC 2", col = "Time period") +
  scale_color_manual(values=clr_prepost, labels = c("Pre-lekking", "Post-lekking")) -> pca_prepost

ggplot(merge_pca, aes(x = pc1, y = pc2)) + geom_point(size=3, aes(col = as.factor(year))) + 
  labs(x = "PC 1", y = "PC 2", col = "Year") +
  scale_color_manual(values=clrs_related[1:3]) -> pca_year

ggplot(merge_pca, aes(x = pc1, y = pc2)) + geom_point(size=3, aes(col = lib)) + 
  labs(x = "PC 1", y = "PC 2", col = "Library") +
  scale_color_manual(values=c(clrs_hunting, clrs_related[1:5])) -> pca_lib

cowplot::plot_grid(hist_mean_cov, hist_n, pca_lib, pca_site, pca_prepost, pca_year,
                   labels="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> hist_plots

ggsave(hist_plots, file="plots/final/supp/summary_stats_cov_meth.png", width=16, height=18)

cowplot::plot_grid(pca_lib, pca_site, pca_prepost, pca_year,
                   labels="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> pcas

ggsave(pcas, file="plots/final/supp/pca_plots.png", width=14, height=14)
