#### Packages #####
pacman::p_load(tidyverse, cowplot, data.table, ggrepel, emmeans, lme4)

#### Theme ####
source("scripts/plotting_theme.R")

#### Reproductive effort ####

## attendance
load(file="results/modeloutput/effort/out_attend_with_pre.RData")
attend <- data

## fight 
load(file="results/modeloutput/effort/out_fight_with_pre.RData")
fight <- data

## dist 
load(file="results/modeloutput/effort/out_dist_with_pre.RData")
dist <- data

rm(data)

##### Attend ####
attend <- attend %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))

ggplot(attend, aes(x = parameter_estimate, y = -log10(parameter_qval))) + 
  geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
  scale_color_manual(values=c(clrs[5], clr_sig)) +
  scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
  labs(x = expression(paste(beta, " estimate lek attendance")), y = expression(-log[10]*"("*italic(q*")"))) +
  geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
  theme(legend.position="none") +
  xlim(c(-0.3,0.3))-> volcano_attend

##### Fight ####
fight <- fight %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))

ggplot(fight, aes(x = parameter_estimate, y = -log10(parameter_qval))) + 
  geom_point(size=6, alpha=0.5, aes(col = sig, fill = sig)) +
  labs(x = expression(paste(beta, " estimate fighting rate")), y = expression(-log[10]*"("*italic(q*")"))) +
  scale_color_manual(values=c(clrs[5], clr_sig)) +
  scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
  geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
  theme(legend.position="none") +
  xlim(c(-0.3,0.3))-> volcano_fight

volcano_fight

##### Centrality #### 
dist <- dist %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))

ggplot(dist, aes(x = parameter_estimate, y = -log10(parameter_qval))) + 
  geom_point(size=6, alpha=0.5, aes(col = sig, fill = sig)) +
  labs(x = expression(paste(beta, " estimate lek centrality")), y = expression(-log[10]*"("*italic(q*")"))) +
  scale_color_manual(values=c(clrs[5], clr_sig)) +
  scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
  geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
  theme(legend.position="none") +
  xlim(c(-0.3,0.3))-> volcano_dist
volcano_dist

#### Physiological changes ####

## mass
load(file="results/modeloutput/physio/out_mass_dif_with_pre.RData")
mass_dif <- data

## microf 
load(file="results/modeloutput/physio/out_microf_dif_with_pre.RData")
microf_dif <- data

## trypa 
load(file="results/modeloutput/physio/out_trypa_dif_with_pre.RData")
trypa_dif <- data

## hct 
load(file="results/modeloutput/physio/out_hct_dif_with_pre.RData")
hct_dif <- data

## ig 
load(file="results/modeloutput/physio/out_ig_dif_with_pre.RData")
igg_dif <- data

rm(data)

##### Mass ####

mass_dif <- mass_dif %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))

ggplot(mass_dif, aes(x = parameter_estimate, y = -log10(parameter_qval))) + 
  geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
  scale_color_manual(values=c(clrs[5], clr_sig)) +
  scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
  labs(x = expression(beta*" estimate "*Delta*" body mass"), y = expression(-log[10]*"("*italic(q*")"))) +
  geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
  theme(legend.position="none") +
  xlim(c(-0.3,0.3)) -> volcano_mass_dif

volcano_mass_dif

##### Microfilaria ####
microf_dif <- microf_dif %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))

ggplot(microf_dif, aes(x = parameter_estimate, y = -log(parameter_qval, base=exp(20)))) + 
  geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
  scale_color_manual(values=c(clrs[5], clr_sig)) +
  scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
  labs(x = expression(beta*" estimate "*Delta*italic(" Microfilaria spp.")), y = "-log20(q-value)") +
  geom_hline(yintercept = -log(0.05,base=exp(20)), col = "darkred", linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
  theme(legend.position="none") +
  xlim(c(-0.3,0.3))-> volcano_microf_dif

volcano_microf_dif

##### Trypanosoma #### 
trypa_dif <- trypa_dif %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))

ggplot(trypa_dif, aes(x = parameter_estimate, y = -log10(parameter_qval))) + 
  geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
  scale_color_manual(values=c(clrs[5], clr_sig)) +
  scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
  labs(x = expression(beta*" estimate "*Delta*italic("Trypanosoma spp.")), y = expression(-log[10]*"("*italic(q*")"))) +
  geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
  theme(legend.position="none") +
  xlim(c(-0.3,0.3))-> volcano_trypa_dif

##### HCT ####

hct_dif <- hct_dif %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))

ggplot(hct_dif, aes(x = parameter_estimate, y = -log10(parameter_qval))) + 
  geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
  scale_color_manual(values=c(clrs[5], clr_sig)) +
  scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
  labs(x = expression(beta*" estimate "*Delta*" HCT"), y = expression(-log[10]*"("*italic(q*")"))) +
  geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
  theme(legend.position="none") +
  xlim(c(-0.3,0.3))-> volcano_hct_dif

volcano_hct_dif

##### IgG ####

## igg_dif 
igg_dif <- igg_dif %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))

ggplot(igg_dif, aes(x = parameter_estimate, y = -log10(parameter_qval))) + 
  geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
  scale_color_manual(values=c(clrs[5], clr_sig)) +
  scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
  labs(x = expression(beta*" estimate "*Delta*" IgG"), y = expression(-log[10]*"("*italic(q*")"))) +
  geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
  theme(legend.position="none") +
  xlim(c(-0.3,0.3)) -> volcano_igg_dif

volcano_igg_dif


#### Raw data attendance ####
load(file = "data/processed/delta_meth_ls_attend.RData")
sig_attend <- subset(attend, sig=="sig")

raw_attend <- subset(delta_meth_sub_attend, chr_pos == sig_attend$chr_pos)
model <- lmerTest::lmer(delta_meth ~ attend_scl + methperc_pre + (1|site) , data=raw_attend, REML=FALSE)
gr <- ref_grid(model, cov.keep= c('attend_scl'))
predict <- as.data.frame(emmeans(gr, spec="attend_scl", level=0.95))

ggplot(raw_attend, aes(x = attend_scl, y = delta_meth*100)) + 
  geom_ribbon(data= predict, aes(ymin = lower.CL*100, ymax = upper.CL*100, y= NULL), fill= clrs[5], alpha = 0.6) +
  geom_line(data= predict, aes(y = emmean*100), col = "black", linewidth=1.5) + 
  geom_hline(yintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
  geom_point(size = 7, fill = clr_high, alpha = 0.6, col = clr_high) + 
  labs(x = "z-transformed attendance", y = expression(Delta*" methylation %")) -> attend_raw

#### assemble figure

plot_grid(volcano_attend, attend_raw,
          volcano_fight,
          volcano_dist,
          volcano_mass_dif, 
          volcano_microf_dif, 
          volcano_trypa_dif,
          volcano_igg_dif, 
          volcano_hct_dif,
          ncol=3, labels= "auto", label_fontface = "plain", label_size = 22) -> volcanoes

ggsave(volcanoes, file="plots/final/main/fig_all_traits.png", width=18, height=18)


