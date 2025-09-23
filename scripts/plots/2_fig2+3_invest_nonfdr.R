pacman::p_load(dplyr, ggplot2, data.table)
#### Theme ####
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

rm(data)


#### Volcanoes #####

###### Attend #######
m_attend_out <- m_attend_out %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig_q", 
                                                        parameter_pval < 0.05 ~ "sig_p", 
                                                        TRUE ~ "nonsig"))

ggplot(m_attend_out, aes(x = parameter_estimate, y = -log10(as.numeric(parameter_pval)))) + 
  geom_point(size=5, alpha=0.5, aes(col = as.factor(sig), fill = as.factor(sig))) +
  labs(x = expression("Estimate lek attendance"), y = expression(-log[10]*"("*italic(p*")")), title = "Lek attendance") +
  scale_color_manual(values=c(clrs[5], clrs_related[4], clrs_related[5])) +
  scale_fill_manual(values=alpha(c(clrs[5], clrs_related[4], clrs_related[5]), 0.5)) +
  geom_hline(yintercept = -log10(0.05), col = clrs_related[4], linetype = "dotted", linewidth = 0.5) +
  geom_vline(xintercept = 0, col = "black", linetype = "dotted", linewidth = 0.5) +
  ylim(0, -log10(0.000005))+
  theme(legend.position="none") +
  xlim(-0.2,0.2) -> volcano_attend
volcano_attend

###### Centrality #######
m_dist_out <- m_dist_out %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig_q", 
                                                        parameter_pval < 0.05 ~ "sig_p", 
                                                        TRUE ~ "nonsig"))

ggplot(m_dist_out, aes(x = parameter_estimate, y = -log10(as.numeric(parameter_pval)))) + 
  geom_point(size=5, alpha=0.5, aes(col = as.factor(sig), fill = as.factor(sig))) +
  labs(x = expression("Estimate centrality"), y = expression(-log[10]*"("*italic(p*")")), title = "Centrality") +
  scale_color_manual(values=c(clrs[5], clrs_related[4], clrs_related[5])) +
  scale_fill_manual(values=alpha(c(clrs[5], clrs_related[4], clrs_related[5]), 0.5)) +
  geom_hline(yintercept = -log10(0.05), col = clrs_related[4], linetype = "dotted", linewidth = 0.5) +
  geom_vline(xintercept = 0, col = "black", linetype = "dotted", linewidth = 0.5) +
  ylim(0, -log10(0.000005))+theme(legend.position="none") +
  xlim(-0.2,0.2) -> volcano_dist

volcano_dist
###### MS #######
m_MS_out <- m_MS_out %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig_q", 
                                                    parameter_pval < 0.05 ~ "sig_p", 
                                                    TRUE ~ "nonsig"))

ggplot(m_MS_out, aes(x = parameter_estimate, y = -log10(as.numeric(parameter_pval)))) + 
  geom_point(size=5, alpha=0.5, aes(col = as.factor(sig), fill = as.factor(sig))) +
  labs(x = expression("Estimate mating success"), y = expression(-log[10]*"("*italic(p*")")), title = "Mating success") +
  scale_color_manual(values=c(clrs[5], clrs_related[4], clrs_related[5])) +
  scale_fill_manual(values=alpha(c(clrs[5], clrs_related[4], clrs_related[5]), 0.5)) +
  geom_hline(yintercept = -log10(0.05), col = clrs_related[4], linetype = "dotted", linewidth = 0.5) +
  geom_vline(xintercept = 0, col = "black", linetype = "dotted", linewidth = 0.5) +
  ylim(0, -log10(0.000005))+theme(legend.position="none") +
  xlim(-0.2,0.2) -> volcano_ms

###### Survival #######
m_surv_out <- delta_out_surv %>% mutate(sig = case_when(surv_delta_meth_qval < 0.05 ~ "sig_q", 
                                                        surv_delta_meth_pval < 0.05 ~ "sig_p", 
                                                        TRUE ~ "nonsig"))

ggplot(m_surv_out, aes(x = surv_delta_meth_estimate, y = -log10(as.numeric(surv_delta_meth_pval)))) + 
  geom_point(size=5, alpha=0.5, aes(col = as.factor(sig), fill = as.factor(sig))) +
  labs(x = expression("Estimate "*Delta*" methylation"), y = expression(-log[10]*"("*italic(p*")")), title = "Survival") +
  scale_color_manual(values=c(clrs[5], clrs_related[4], clrs_related[5])) +
  scale_fill_manual(values=alpha(c(clrs[5], clrs_related[4], clrs_related[5]), 0.5)) +
  geom_hline(yintercept = -log10(0.05), col = clrs_related[4], linetype = "dotted", linewidth = 0.5) +
  geom_vline(xintercept = 0, col = "black", linetype = "dotted", linewidth = 0.5) +
  ylim(0, -log10(0.000005))+
  xlim(-10,10)+
  theme(legend.position="none") -> volcano_surv

volcano_surv

###### Blue #######
m_blue_out <- m_blue_out %>% mutate(sig = case_when(deltameth_qval < 0.05 ~ "sig_q", 
                                                    deltameth_pval < 0.05 ~ "sig_p", 
                                                        TRUE ~ "nonsig"))

ggplot(m_blue_out, aes(x = deltameth_estimate, y = -log10(as.numeric(deltameth_pval)))) + 
  geom_point(size=5, alpha=0.5, aes(col = as.factor(sig), fill = as.factor(sig))) +
  labs(x = expression("Estimate "*Delta*" methylation"), y = expression(-log[10]*"("*italic(p*")")), title = "Blue chroma next year") +
  scale_color_manual(values=c(clrs[5], clrs_related[4], clrs_related[5])) +
  scale_fill_manual(values=alpha(c(clrs[5], clrs_related[4], clrs_related[5]), 0.5)) +
  geom_hline(yintercept = -log10(0.05), col = clrs_related[4], linetype = "dotted", linewidth = 0.5) +
  geom_vline(xintercept = 0, col = "black", linetype = "dotted", linewidth = 0.5) +
  ylim(0, -log10(0.000005))+
  xlim(-0.015, 0.015)+
  theme(legend.position="none")  -> volcano_blue

volcano_blue

###### Lyre #######
m_lyre_out <- m_lyre_out %>% mutate(sig = case_when(deltameth_qval < 0.05 ~ "sig_q", 
                                                    deltameth_pval < 0.05 ~ "sig_p", 
                                                    TRUE ~ "nonsig"))

ggplot(m_lyre_out, aes(x = deltameth_estimate, y = -log10(as.numeric(deltameth_pval)))) + 
  geom_point(size=5, alpha=0.5, aes(col = as.factor(sig), fill = as.factor(sig))) +
  labs(x = expression("Estimate "*Delta*" methylation"), y = expression(-log[10]*"("*italic(p*")")), title = "Lyre size next year") +
  scale_color_manual(values=c(clrs[5], clrs_related[4], clrs_related[5])) +
  scale_fill_manual(values=alpha(c(clrs[5], clrs_related[4], clrs_related[5]), 0.5)) +
  geom_hline(yintercept = -log10(0.05), col = clrs_related[4], linetype = "dotted", linewidth = 0.5) +
  geom_vline(xintercept = 0, col = "black", linetype = "dotted", linewidth = 0.5) +
  ylim(0, -log10(0.000005))+
  xlim(-10,10)+
  theme(legend.position="none")  -> volcano_lyre

volcano_lyre

### Raw data sig ones ####

###### Attend #######
# select sig cpg
sig_attend <- subset(m_attend_out, parameter_qval < 0.05)
# load data used for the model
load(file = "data/processed/delta_meth_ls_attend.RData")
data_sig_attend <- subset(delta_meth_sub_attend, chr_pos %in% sig_attend$chr_pos)

model_sig_attend <- lmerTest::lmer(delta_meth ~ attend + age + methperc_pre + (1|site), data = data_sig_attend)

# plot

source("scripts/function_effect_plot_custom.R")
effect_plot(model_sig_attend, pred = attend, plot.points=T, interval=T, data = data_sig_attend,
            line.colors = clrs[5], point.size = 5, colors = clrs_related[5]) + 
  labs(x = "Attendance", y = expression(Delta*" methylation")) +
  scale_color_manual(values = clrs_related[5]) -> raw_attend

###### MS #######
# select sig cpg
sig_MS <- subset(m_MS_out, parameter_qval < 0.05)
# load data used for the model
load(file = "data/processed/delta_meth_MS_sub.RData")

data_sig_MS <- subset(delta_meth_sub_MS, chr_pos %in% sig_MS$chr_pos)

model_sig_MS <- lmerTest::lmer(delta_meth ~ MS + age + methperc_pre + (1|site), data = data_sig_MS)

# plot
effect_plot(model_sig_MS, pred = MS, plot.points=T, interval=T, data = data_sig_MS,
            line.colors = clrs[5], point.size = 5, colors = clrs_related[5]) + 
  labs(x = "Mating success", y = expression(Delta*" methylation")) +
  scale_color_manual(values = clrs_related[5]) -> raw_MS

###### Lyre #######

# select sig cpg
sig_lyre <- subset(m_lyre_out, deltameth_qval < 0.05)

# load data used for the model
load(file = "data/processed/delta_meth_sub_lyre_ny.RData")

data_sig_lyre_1 <- subset(delta_meth_sub_lyre_ny, chr_pos %in% sig_lyre$chr_pos[1])

model_sig_lyre_1 <- lmerTest::lmer(lyre_nextyear ~ scale(delta_meth) + age_cat + (1|site), data = data_sig_lyre_1)

# plot
effect_plot(model_sig_lyre_1, pred = delta_meth, plot.points=T, interval=T, data = data_sig_lyre_1,
            line.colors = clrs[5], point.size = 5, colors = clrs_related[5]) + 
  labs(y = "Lyre size next year", x = expression(Delta*" methylation")) +
  scale_color_manual(values = clrs_related[5]) -> raw_lyre_1

raw_lyre_1

data_sig_lyre_2 <- subset(delta_meth_sub_lyre_ny, chr_pos %in% sig_lyre$chr_pos[2])

model_sig_lyre_2 <- lmerTest::lmer(lyre_nextyear ~ scale(delta_meth) + age_cat + (1|site), data = data_sig_lyre_2)

# plot
effect_plot(model_sig_lyre_2, pred = delta_meth, plot.points=T, interval=T, data = data_sig_lyre_2,
            line.colors = clrs[5], point.size = 5, colors = clrs_related[5]) + 
  labs(y = "Lyre size next year", x = expression(Delta*" methylation")) +
  scale_color_manual(values = clrs_related[5]) -> raw_lyre_2

raw_lyre_2

#### Assemble two figures #####

# volcanoes

plot_grid(volcano_attend, volcano_dist, volcano_ms,ncol=3) -> fig2_volcanoes_invest
plot_grid(volcano_surv, volcano_blue, volcano_lyre, ncol=3) -> fig2_volcanoes_cost

plot_grid(fig2_volcanoes_invest, fig2_volcanoes_cost, labels=c("auto"), ncol=1,
          label_fontface = "plain", label_size = 22)-> fig2
fig2
ggsave(fig2, file = "plots/final/main/fig_2_volcanoes_invest_cost.png", width=14, height=10)
ggsave(fig2, file="plots/final/main/fig_2_volcanoes_invest_cost.pdf", width=15, height=20, device = cairo_pdf)


plot_grid(raw_attend, raw_MS, raw_lyre_1, raw_lyre_2,
          ncol=2, labels="auto", 
          label_fontface = "plain", label_size = 22) -> fig3

ggsave(fig3, file = "plots/final/main/fig_3_raw_invest_cost.png", width=16, height=12)
ggsave(fig3, file="plots/final/main/fig_3_raw_invest_cost.pdf", width=16, height=12, device = cairo_pdf)

