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

load(file = "results/modeloutput/all_sites_deltameth.RData")

delta_meth <- subset(delta_meth, chr_pos %in% changing_cpg$chr_pos)

### Fig S2: Plot pre_meth effect on delta_meth ####

lmer_pre_delta <- lmerTest::lmer(delta_meth ~ methperc_pre + (1|chr_pos)+ (1|id), data = delta_meth)

sum_pre_delta <- as.data.frame(summary(lmer_pre_delta)$coef) # negative relationship

ggplot(delta_meth, aes(methperc_pre, delta_meth)) + 
  geom_pointdensity() + 
  scale_color_viridis_c(option="A") + 
  geom_abline(intercept=sum_pre_delta$Estimate[1], slope = sum_pre_delta$Estimate[2], 
              color="red", linewidth=1)+
  labs(x = "Methylation % pre-lekking", y = expression(Delta*" methylation %"), col = "Density") -> cor_pre_delta

ggplot(delta_meth, aes(methperc_pre, abs(delta_meth))) + 
  geom_pointdensity() + 
  scale_color_viridis_c(option="A") + #geom_smooth() + 
  labs(x = "Methylation % pre-lekking", y = expression("Absolute "*Delta*" methylation %"), col = "Density") -> cor_pre_delta_abs

cowplot::plot_grid(cor_pre_delta, cor_pre_delta_abs, labels="auto", ncol=1, align="hv", axis="lb", label_fontface = "plain", label_size = 22) -> cor_premeths
cor_premeths
ggsave(cor_premeths, file = "plots/final/supp/sfig_2_pre_vs_delta.png", width=8, height=12)
ggsave(cor_premeths, file = "plots/final/supp/sfig_2_pre_vs_delta.pdf", width=8, height=12, device=cairo_pdf)

### Fig S4: volcano plots pre-meth ####
models <- list(m_attend_out, m_dist_out, m_MS_out)

# attend
m_attend_out$pre_qval <- p.adjust(m_attend_out$pre_pval, method = "fdr", n = nrow(m_attend_out))
m_attend_out <- m_attend_out %>% mutate(sig = case_when(pre_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))
m_attend_out_pmin <-max(m_attend_out$pre_pval[m_attend_out$pre_qval < 0.05])

ggplot(m_attend_out, aes(x = pre_estimate, y = -log10(pre_pval))) + 
    geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
    scale_color_manual(values=c(clr_grey, clr_sig)) +
    scale_fill_manual(values=alpha(c(clr_grey, clr_sig), 0.6)) +
    labs(x = expression(paste(beta, " estimate pre-lekking methylation %")), y = expression(-log[10]*"("*italic(p*")")), title = "Attendance") +
    geom_hline(yintercept = -log10(m_attend_out_pmin), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
    xlim(-2,2)+
    theme(legend.position="none") -> volcano_attend
    
# dist
m_dist_out$pre_qval <- p.adjust(m_dist_out$pre_pval, method = "fdr", n = nrow(m_dist_out))
m_dist_out <- m_dist_out %>% mutate(sig = case_when(pre_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))
m_dist_out_pmin <-max(m_dist_out$pre_pval[m_dist_out$pre_qval < 0.05])

ggplot(m_dist_out, aes(x = pre_estimate, y = -log10(pre_pval))) + 
  geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
  scale_color_manual(values=c(clr_grey, clr_sig)) +
  scale_fill_manual(values=alpha(c(clr_grey, clr_sig), 0.6)) +
  labs(x = expression(paste(beta, " estimate pre-lekking methylation %")), y = expression(-log[10]*"("*italic(p*")")), title = "Centrality") +
  geom_hline(yintercept = -log10(m_dist_out_pmin), col = "darkred", linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
  xlim(-2,2)+
  theme(legend.position="none") -> volcano_dist

# MS
m_MS_out$pre_qval <- p.adjust(m_MS_out$pre_pval, method = "fdr", n = nrow(m_MS_out))
m_MS_out <- m_MS_out %>% mutate(sig = case_when(pre_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))
m_MS_out_pmin <-max(m_MS_out$pre_pval[m_MS_out$pre_qval < 0.05])

ggplot(m_MS_out, aes(x = pre_estimate, y = -log10(pre_pval))) + 
  geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
  scale_color_manual(values=c(clr_grey, clr_sig)) +
  scale_fill_manual(values=alpha(c(clr_grey, clr_sig), 0.6)) +
  labs(x = expression(paste(beta, " estimate pre-lekking methylation %")), y = expression(-log[10]*"("*italic(p*")")), title = "Mating success") +
  geom_hline(yintercept = -log10(m_MS_out_pmin), col = "darkred", linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
  xlim(-2,2)+
  theme(legend.position="none") -> volcano_MS

### combine in plot
plot_grid(volcano_attend, 
          volcano_dist, 
          volcano_MS,
          ncol=1, align="hv", axis="lb", labels="auto", label_fontface = "plain", label_size = 22 ) -> all_plots
all_plots
ggsave(all_plots, file = "plots/final/supp/sfig_4_volcano_pre_lekking.png", width = 8, height = 14)
ggsave(all_plots, file = "plots/final/supp/sfig_4_volcano_pre_lekking.pdf", width = 8, height = 14, device=cairo_pdf)
