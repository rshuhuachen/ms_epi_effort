pacman::p_load(dplyr, ggplot2, data.table)
source("scripts/plotting_theme.R")

# attend
load(file="results/modeloutput/effort/out_attend_with_pre.RData")
m_attend_out <- data

# dist
load(file="results/modeloutput/effort/out_dist_with_pre.RData")
m_dist_out <- data

# ms
load(file="results/modeloutput/effort/out_MS_with_pre.RData")
m_MS_out <- data

# blue
load(file="results/modeloutput/nextyear/out_blue_ny.RData")
m_blue_out <- data

# lyre
load(file="results/modeloutput/nextyear/out_lyre_ny.RData")
m_lyre_out <- data

# survival
load(file="results/modeloutput/fitness/out_surv_deltameth_filtered.RData")

### volcano plot pre-meth effects ####
models <- list(m_attend_out, m_dist_out, m_MS_out)

plots <- list()
for (i in 1:length(models)){
  df <- as.data.frame(models[[i]])
  df <- df %>% mutate(sig = case_when(pre_pval < 0.05 ~ "sig", TRUE ~ "nonsig"))
  
  ggplot(df, aes(x = pre_estimate, y = -log10(pre_pval))) + 
    geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
    scale_color_manual(values=c(clr_grey, clr_sig)) +
    scale_fill_manual(values=alpha(c(clr_grey, clr_sig), 0.6)) +
    labs(x = expression(paste(beta, " estimate pre-lekking methylation %")), y = "-log10(q-value)", title = as.character(df$parameter[1])) +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") -> volcano
  plots[[i]] <- volcano        
  
  ggsave(volcano, file= paste0("plots/model_out/pre_lekking/volcano_", as.character(df$parameter[1]), ".png"), width=10, height=10)}

### combine in plot
plot_grid(plots[[1]] + labs(title = "Attendance"), 
          plots[[2]] + labs(title = "Centrality"), 
          plots[[3]] + labs(title = "Mating success"),
          ncol=1, align="hv", axis="lb", labels="auto", label_fontface = "plain", label_size = 22 ) -> all_plots
all_plots
ggsave(all_plots, file = "plots/final/supp_volcano_pre_lekking.png", width = 12, height = 18)
