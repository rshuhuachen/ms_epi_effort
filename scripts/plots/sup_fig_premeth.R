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


load(file = "results/modeloutput/all_sites_deltameth.RData")

delta_meth <- subset(delta_meth, chr_pos %in% changing_cpg$chr_pos)

### Plot pre_meth effect on delta_meth ####

lmer_pre_delta <- lmerTest::lmer(delta_meth ~ methperc_pre + (1|chr_pos)+ (1|id), data = delta_meth)

sum_pre_delta <- as.data.frame(summary(lmer_pre_delta)$coef) # negative relationship

ggplot(delta_meth, aes(methperc_pre, delta_meth)) + 
  geom_pointdensity() + 
  scale_color_viridis_c() + 
  geom_abline(intercept=sum_pre_delta$Estimate[1], slope = sum_pre_delta$Estimate[2], 
              color="red", linewidth=1)+
  labs(x = "Methylation % pre-lekking", y = expression(Delta*" methylation %")) -> cor_pre_delta

ggplot(delta_meth, aes(methperc_pre, abs(delta_meth))) + 
  geom_pointdensity() + 
  scale_color_viridis_c() + #geom_smooth() + 
  labs(x = "Methylation % pre-lekking", y = expression("Absolute "*Delta*" methylation %")) -> cor_pre_delta_abs

cowplot::plot_grid(cor_pre_delta, cor_pre_delta_abs, labels="auto", ncol=1, align="hv", axis="lb", label_fontface = "plain", label_size = 22) -> cor_premeths
cor_premeths
ggsave(cor_premeths, file = "plots/final/supp/pre_vs_delta.png", width=6, height=12)

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
    xlim(-2,2)+
    theme(legend.position="none") -> volcano
  plots[[i]] <- volcano        
  
  ggsave(volcano, file= paste0("plots/model_out/pre_lekking/volcano_", as.character(df$parameter[1]), ".png"), width=10, height=10)}

### combine in plot
plot_grid(plots[[1]] + labs(title = "Attendance"), 
          plots[[2]] + labs(title = "Centrality"), 
          plots[[3]] + labs(title = "Mating success"),
          ncol=1, align="hv", axis="lb", labels="auto", label_fontface = "plain", label_size = 22 ) -> all_plots
all_plots
ggsave(all_plots, file = "plots/final/supp_volcano_pre_lekking.png", width = 8, height = 14)
