#### Packages #####
#extrafont::loadfonts(device="all")

pacman::p_load(tidyverse, cowplot, data.table)

#### Theme ####
source("scripts/plotting_theme.R")

### Function qq plots ####
source("scripts/plots/function_qqplot.R")

### Load raw model outputs ####

#### Plot 1: changing and fitness
### Changing
load(file="results/modeloutput/prepost_modeloutput_glmer_min0.75_raw.RData")
out_changing <- subset(out_glmer_raw, convergence == "boundary (singular) fit: see help('isSingular')" | is.na(convergence))

### AMS
load(file="results/modeloutput/fitness/out_ams_nopre_raw.RData")
out_ams <- subset(delta_out_ams, ams_message == "relative convergence (4)")

### Survival
load(file="results/modeloutput/fitness/out_surv_nopre_raw.RData")
out_surv <- subset(delta_out_surv, surv_message == "relative convergence (4)")

### Plot 2: effort
load(file="results/modeloutput/effort/out_attend_with_pre.RData")
out_attend <- data
load(file="results/modeloutput/effort/out_fight_with_pre.RData")
out_fight <- data
load(file="results/modeloutput/effort/out_dist_with_pre.RData")
out_dist <- data
rm(data)

### Plot 3: physio
load(file="results/modeloutput/physio/out_mass_dif_with_pre.RData")
out_mass <- data
load(file="results/modeloutput/physio/out_microf_dif_with_pre.RData")
out_microf <- data
load(file="results/modeloutput/physio/out_trypa_dif_with_pre.RData")
out_trypa <- data
load(file="results/modeloutput/physio/out_ig_dif_with_pre.RData")
out_igg <- data
load(file="results/modeloutput/physio/out_hct_dif_with_pre.RData")
out_hct <- data
rm(data)

### Plot 1: changing and fitness ####

## changing
# histogram p-values raw
ggplot(out_changing, aes(x = prepost_pval)) + geom_histogram(col = "black", fill = clr_grey) + 
  scale_y_continuous(labels = scales::unit_format(unit = "K", scale = 1e-3)) + 
  labs(x = expression(italic(p)~values), y = "Count",
       title = "Changing CpG sites")-> hist_pvals_changing

# qq plot 
qq_plot(out_changing$prepost_pval, title = "Changing CpG sites", CB=F, thinning = F) -> qqplot_changing

## ams
# histogram p-values raw
ggplot(out_ams, aes(x = ams_delta_meth_pval)) + geom_histogram(col = "black", fill = clr_grey) + 
  labs(x = expression(italic(p)~values), y = "Count",
       title = "AMS")-> hist_pvals_ams

# qq plot 
qq_plot(out_ams$ams_delta_meth_pval, title = "AMS", CB=F, thinning = F) -> qqplot_ams

## surv
# histogram p-values raw
ggplot(out_surv, aes(x = surv_delta_meth_pval)) + geom_histogram(col = "black", fill = clr_grey) + 
  labs(x = expression(italic(p)~values), y = "Count",
       title = "Survival")-> hist_pvals_surv

# qq plot 
qq_plot(out_surv$surv_delta_meth_pval, title = "Survival", CB=F, thinning = F) -> qqplot_surv

#### Arrange fig 1 ####

plot_grid(hist_pvals_changing, qqplot_changing, 
          hist_pvals_ams, qqplot_ams,
          hist_pvals_surv, qqplot_surv,
          labels=c("auto"), align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22)-> qqs_1

ggsave(qqs_1, file = "plots/final/supp/hist_qqplots_changing_fitness.png", width = 14, height = 18)

#### Plot 2: effort ####
## attend
# histogram p-values raw
ggplot(out_attend, aes(x = parameter_pval)) + geom_histogram(col = "black", fill = clr_grey) + 
  labs(x = expression(italic(p)~values), y = "Count",
       title = "Attendance")-> hist_pvals_attend

# qq plot 
qq_plot(out_attend$parameter_pval, title = "Attendance", CB=F, thinning = F) -> qqplot_attend

## fight
# histogram p-values raw
ggplot(out_fight, aes(x = parameter_pval)) + geom_histogram(col = "black", fill = clr_grey) + 
  labs(x = expression(italic(p)~values), y = "Count",
       title = "Fighting")-> hist_pvals_fight

# qq plot 
qq_plot(out_fight$parameter_pval, title = "Fighting", CB=F, thinning = F) -> qqplot_fight

## dist
# histogram p-values raw
ggplot(out_dist, aes(x = parameter_pval)) + geom_histogram(col = "black", fill = clr_grey) + 
  labs(x = expression(italic(p)~values), y = "Count",
       title = "Centrality")-> hist_pvals_dist

# qq plot 
qq_plot(out_dist$parameter_pval, title = "Centrality", CB=F, thinning = F) -> qqplot_dist

plot_grid(hist_pvals_attend, qqplot_attend, 
          hist_pvals_fight, qqplot_fight,
          hist_pvals_dist, qqplot_dist,
          labels="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22)-> qqs_2

ggsave(qqs_2, file = "plots/final/supp/hist_qqplots_effort.png", width = 14, height = 18)

#### Plot 3: physio ####
## mass
# histogram p-values raw
ggplot(out_mass, aes(x = parameter_pval)) + geom_histogram(col = "black", fill = clr_grey) + 
  labs(x = expression(italic(p)~values), y = "Count",
       title = expression(Delta~mass))-> hist_pvals_mass

# qq plot 
qq_plot(out_mass$parameter_pval, title = expression(Delta~mass), CB=F, thinning = F) -> qqplot_mass

## microf
ggplot(out_microf, aes(x = parameter_pval)) + geom_histogram(col = "black", fill = clr_grey) + 
  labs(x = expression(italic(p)~values), y = "Count",
       title = expression(Delta~Microfilaria~spp.))-> hist_pvals_microf

# qq plot 
qq_plot(out_microf$parameter_pval, title = expression(Delta~Microfilaria~spp.), CB=F, thinning = F) -> qqplot_microf

## trypa
ggplot(out_trypa, aes(x = parameter_pval)) + geom_histogram(col = "black", fill = clr_grey) + 
  labs(x = expression(italic(p)~values), y = "Count",
       title = expression(Delta~Trypanosoma~spp.))-> hist_pvals_trypa

# qq plot 
qq_plot(out_trypa$parameter_pval, title = expression(Delta~Trypanosoma~spp.), CB=F, thinning = F) -> qqplot_trypa

# igg
ggplot(out_igg, aes(x = parameter_pval)) + geom_histogram(col = "black", fill = clr_grey) + 
  labs(x = expression(italic(p)~values), y = "Count",
       title = expression(Delta~IgG))-> hist_pvals_igg

# qq plot 
qq_plot(out_igg$parameter_pval, title = expression(Delta~Delta~IgG), CB=F, thinning = F) -> qqplot_igg

# hct
ggplot(out_hct, aes(x = parameter_pval)) + geom_histogram(col = "black", fill = clr_grey) + 
  labs(x = expression(italic(p)~values), y = "Count",
       title = expression(Delta~HCT))-> hist_pvals_hct

# qq plot 
qq_plot(out_hct$parameter_pval, title = expression(Delta~HCT), CB=F, thinning = F) -> qqplot_hct

plot_grid(hist_pvals_mass, qqplot_mass,
          hist_pvals_microf, qqplot_microf,
          hist_pvals_trypa, qqplot_trypa,
          hist_pvals_igg, qqplot_igg,
          hist_pvals_hct, qqplot_hct,
          labels="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22)-> qqs_3

ggsave(qqs_3, file = "plots/final/supp/hist_qqplots_physio.png", width = 14, height = 26)
