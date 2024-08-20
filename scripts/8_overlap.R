
### load packages
pacman::p_load(tidyverse, data.table, cowplot)

source("scripts/plotting_theme.R")

### load data

load(file = "results/modeloutput/subset_sites_sig_deltameth.RData")

### load phenotypic data

load("data/phenotypes/fulldata_complete_epi_withdates.RData")
load("data/phenotypes/pheno_dif_prepost.RData")

#combine with site and fitness data
pheno_pre <- subset(all_pheno_epi, prepost=="pre")

delta_meth <- left_join(delta_meth, unique(pheno_pre[,c("id", "year", "MS", "surv", "site")]), by = c("id", "year"))
         
### is there any overlap with behavioural or physiological traits?

## load ams model results
load(file="results/modeloutput/AMS_deltameth_modeloutput_filtered.RData")
cpg_sig_ams <- subset(delta_out_ams, ams_delta_meth_qval < 0.05) #481

## load reproductive effort model out
load(file="results/modeloutput/effort_deltameth_modeloutput_filtered.RData")
cpg_effort <- subset(delta_out_all, parameter_qval < 0.05 & abs(parameter_estimate) > 0.1)
  
## load physio changes model out  
load(file="results/modeloutput/physio_deltameth_modeloutput_filtered.RData")
cpg_physio <- subset(delta_out_all, parameter_qval < 0.05 & abs(parameter_estimate) > 0.1)
  
nrow(subset(cpg_effort, chr_pos %in% cpg_sig_ams$chr_pos))  # 8
nrow(subset(cpg_physio, chr_pos %in% cpg_sig_ams$chr_pos))  # 18
nrow(subset(cpg_effort, chr_pos %in% cpg_physio$chr_pos))  # 0


### Merge the data of those #####

### effort and AMS
overlap_ams_effort <- inner_join(cpg_effort, cpg_sig_ams, by="chr_pos")

prepost_dif$attend_scl <- scale(prepost_dif$attend)
prepost_dif$fight_scl <- scale(prepost_dif$fight)
prepost_dif$dist_scl <- scale(prepost_dif$dist)

delta_meth <- left_join(delta_meth, unique(prepost_dif[,c("id", "year", "attend_scl", "fight_scl", "dist_scl")], by = c("id", "year")))

### physio and AMS
overlap_ams_physio <- inner_join(cpg_physio, cpg_sig_ams, by="chr_pos")

prepost_dif$mass_dif_scl <- scale(prepost_dif$mass_dif)
prepost_dif$microf_dif_scl <- scale(prepost_dif$microf_dif)
prepost_dif$trypa_dif_scl <- scale(prepost_dif$trypa_dif)
prepost_dif$ig_dif_scl <- scale(prepost_dif$ig_dif)
prepost_dif$hct_dif_scl <- scale(prepost_dif$hct_dif)

delta_meth <- left_join(delta_meth, unique(prepost_dif[,c("id", "year", "mass_dif", "microf_dif", "trypa_dif", "ig_dif", "hct_dif",
                        "mass_dif_scl", "microf_dif_scl", "trypa_dif_scl", "ig_dif_scl", "hct_dif_scl")], by = c("id", "year")))

### Plot raw data ####

#### Effort and AMS ####

### attend
overlap_attend_ams <- subset(overlap_ams_effort, parameter == "attend") #1

list_plot_attend_a <- list()
list_plot_attend_b <- list()

for (i in 1:nrow(overlap_attend_ams)){
    ggplot(subset(delta_meth, chr_pos == overlap_attend_ams$chr_pos[i]), aes(x = attend_scl, y = scale(delta_meth))) + 
        geom_point(fill=clrs_hunting[1], size=3) + labs(x = "z-transformed attendance", y = expression("z-trans "*Delta*" methylation"),
                     title = paste0("Estimate = ", round(overlap_attend_ams$parameter_estimate[i], 2), 
                                    ", q-value = ", round(overlap_attend_ams$parameter_qval[i], 4)), subtitle = overlap_attend_ams$chr_pos[i]) +
                                        geom_abline(intercept=overlap_attend_ams$intercept[i], slope = overlap_attend_ams$parameter_estimate[i], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_a
    list_plot_attend_a[[i]] <- plot_a 

    ggplot(subset(delta_meth, chr_pos == overlap_attend_ams$chr_pos[i]), aes(x = scale(delta_meth), y = MS)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(y = "AMS", x = expression("z-trans "*Delta*" methylation"),
                      title = paste0("Estimate = ", round(overlap_attend_ams$ams_delta_meth_estimate[i], 2),
                                        ", q-value = ", round(overlap_attend_ams$ams_delta_meth_qval[i], 2)), 
                      subtitle = overlap_attend_ams$chr_pos[i]) +
                                        geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_b
    list_plot_attend_b[[i]] <- plot_b}

cowplot::plot_grid(list_plot_attend_a[[1]], list_plot_attend_b[[1]], 
    labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_attend_ams

ggsave(plots_attend_ams, file = "plots/model_out/rawdata_plot_overlap_attend_AMS.png", width=14, height=12)                                

### fight
overlap_fight_ams <- subset(overlap_ams_effort, parameter == "fight") #6

list_plot_fight_a <- list()
list_plot_fight_b <- list()

for (i in 1:nrow(overlap_fight_ams)){
  ggplot(subset(delta_meth, chr_pos == overlap_fight_ams$chr_pos[i]), aes(x = fight_scl, y = scale(delta_meth))) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = "z-transformed fight", y = expression("z-trans "*Delta*" methylation"),
                                                    title = paste0("Estimate = ", round(overlap_fight_ams$parameter_estimate[i], 2), 
                                                                   ", q-value = ", round(overlap_fight_ams$parameter_qval[i], 4)), subtitle = overlap_fight_ams$chr_pos[i]) +
    geom_abline(intercept=overlap_fight_ams$intercept[i], slope = overlap_fight_ams$parameter_estimate[i], 
                color=clrs_hunting[2], linewidth=1)+
    geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_a
  list_plot_fight_a[[i]] <- plot_a 
  
  ggplot(subset(delta_meth, chr_pos == overlap_fight_ams$chr_pos[i]), aes(x = scale(delta_meth), y = MS)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(y = "AMS", x = expression("z-trans "*Delta*" methylation"),
                                                    title = paste0("Estimate = ", round(overlap_fight_ams$ams_delta_meth_estimate[i], 2),
                                                                   ", q-value = ", round(overlap_fight_ams$ams_delta_meth_qval[i], 2)), 
                                                    subtitle = overlap_fight_ams$chr_pos[i]) +
    geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
    geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_b
  list_plot_fight_b[[i]] <- plot_b}

cowplot::plot_grid(list_plot_fight_a[[1]], list_plot_fight_b[[1]], list_plot_fight_a[[2]], list_plot_fight_b[[2]], 
                   list_plot_fight_a[[3]], list_plot_fight_b[[3]], list_plot_fight_a[[4]], list_plot_fight_b[[4]], 
                   list_plot_fight_a[[5]], list_plot_fight_b[[5]], list_plot_fight_a[[6]], list_plot_fight_b[[6]], 
                   labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_fight_ams

ggsave(plots_fight_ams, file = "plots/model_out/rawdata_plot_overlap_fight_AMS.png", width=14, height=30)                                

### dist
overlap_dist_ams <- subset(overlap_ams_effort, parameter == "dist") #1

list_plot_dist_a <- list()
list_plot_dist_b <- list()

for (i in 1:nrow(overlap_dist_ams)){
  ggplot(subset(delta_meth, chr_pos == overlap_dist_ams$chr_pos[i]), aes(x = dist_scl, y = scale(delta_meth))) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = "z-transformed centrality", y = expression("z-trans "*Delta*" methylation"),
                                                    title = paste0("Estimate = ", round(overlap_dist_ams$parameter_estimate[i], 2), 
                                                                   ", q-value = ", round(overlap_dist_ams$parameter_qval[i], 4)), subtitle = overlap_dist_ams$chr_pos[i]) +
    geom_abline(intercept=overlap_dist_ams$intercept[i], slope = overlap_dist_ams$parameter_estimate[i], 
                color=clrs_hunting[2], linewidth=1)+
    geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_a
  list_plot_dist_a[[i]] <- plot_a 
  
  ggplot(subset(delta_meth, chr_pos == overlap_dist_ams$chr_pos[i]), aes(x = scale(delta_meth), y = MS)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(y = "AMS", x = expression("z-trans "*Delta*" methylation"),
                                                    title = paste0("Estimate = ", round(overlap_dist_ams$ams_delta_meth_estimate[i], 2),
                                                                   ", q-value = ", round(overlap_dist_ams$ams_delta_meth_qval[i], 2)), 
                                                    subtitle = overlap_dist_ams$chr_pos[i]) +
    geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
    geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_b
  list_plot_dist_b[[i]] <- plot_b}

cowplot::plot_grid(list_plot_dist_a[[1]], list_plot_dist_b[[1]], 
                   labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_dist_ams

ggsave(plots_dist_ams, file = "plots/model_out/rawdata_plot_overlap_dist_AMS.png", width=14, height=12)                                

#### Physio and AMS ####
### Body mass  

overlap_mass_ams <- subset(overlap_ams_physio, parameter == "Delta body mass")

list_plot_mass_a <- list()
list_plot_mass_b <- list()

for (i in 1:nrow(overlap_mass_ams)){
    ggplot(subset(delta_meth, chr_pos == overlap_mass_ams$chr_pos[i]), aes(x = mass_dif_scl, y = scale(delta_meth))) + 
        geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression("z-trans "*Delta*" body mass"), y = expression("z-trans "*Delta*" methylation"),
                                                        subtitle = overlap_mass_ams$chr_pos[i],
                     title = paste0("Estimate = ", round(overlap_mass_ams$parameter_estimate[i], 2), ", q-value = ", round(overlap_mass_ams$parameter_qval[i], 4))) +
                                        geom_abline(intercept=overlap_mass_ams$intercept[i], slope = overlap_mass_ams$parameter_estimate[i], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_a
    list_plot_mass_a[[i]] <- plot_a 

    ggplot(subset(delta_meth, chr_pos == overlap_mass_ams$chr_pos[i]), aes(x = scale(delta_meth), y = MS)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(y = "AMS", x = expression("z-trans "*Delta*" methylation"),subtitle = overlap_mass_ams$chr_pos[i],
                      title = paste0("Estimate = ", round(overlap_mass_ams$ams_delta_meth_estimate[i], 2),
                                        ", q-value = ", round(overlap_mass_ams$ams_delta_meth_qval[i], 2))) +
                                        geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_b
    list_plot_mass_b[[i]] <- plot_b}

cowplot::plot_grid(list_plot_mass_a[[1]], list_plot_mass_b[[1]], list_plot_mass_a[[2]], list_plot_mass_b[[2]],
    list_plot_mass_a[[3]], list_plot_mass_b[[3]],
    labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_mass_ams

ggsave(plots_mass_ams, file = "plots/model_out/rawdata_plot_overlap_mass_AMS.png", width=14, height=12)                                

### Delta HCT 

overlap_hct <- subset(overlap_ams_physio, parameter == "Delta HCT") #9

list_plot_hct_a <- list()
list_plot_hct_b <- list()

for (i in 1:nrow(overlap_hct)){
    ggplot(subset(delta_meth, chr_pos == overlap_hct$chr_pos[i]), aes(x = hct_dif_scl, y = scale(delta_meth))) + 
        geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression("z-trans "*Delta*" HCT"), y = expression("z-trans "*Delta*" methylation"),
                                                        subtitle = overlap_hct$chr_pos[i],
                     title = paste0("Estimate = ", round(overlap_hct$parameter_estimate[i], 2), ", q-value = ", round(overlap_hct$parameter_qval[i], 4))) +
                                        geom_abline(intercept=overlap_hct$intercept[i], slope = overlap_hct$parameter_estimate[i], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_a
    list_plot_hct_a[[i]] <- plot_a 

    ggplot(subset(delta_meth, chr_pos == overlap_hct$chr_pos[i]), aes(x = scale(delta_meth), y = MS)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(y = "AMS", x = expression("z-trans "*Delta*" methylation"),subtitle = overlap_hct$chr_pos[i],
                      title = paste0("Estimate = ", round(overlap_hct$ams_delta_meth_estimate[i], 2),
                                        ", q-value = ", round(overlap_hct$ams_delta_meth_qval[i], 2))) +
                                        geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_b
    list_plot_hct_b[[i]] <- plot_b}

cowplot::plot_grid(list_plot_hct_a[[1]], list_plot_hct_b[[1]], list_plot_hct_a[[2]], list_plot_hct_b[[2]],
    list_plot_hct_a[[3]], list_plot_hct_b[[3]],list_plot_hct_a[[4]], list_plot_hct_b[[4]],
    list_plot_hct_a[[5]], list_plot_hct_b[[5]],list_plot_hct_a[[6]], list_plot_hct_b[[6]],
    list_plot_hct_a[[7]], list_plot_hct_b[[7]],list_plot_hct_a[[8]], list_plot_hct_b[[8]],
    list_plot_hct_a[[9]], list_plot_hct_b[[9]],
    labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_hct_ams

ggsave(plots_hct_ams, file = "plots/model_out/rawdata_plot_overlap_hct_AMS.png", width=14, height=24)                                

### delta IgG

overlap_igg_ams <- subset(overlap_ams_physio, parameter == "Delta IgG") #1

list_plot_igg_a <- list()
list_plot_igg_b <- list()

for (i in 1:nrow(overlap_igg_ams)){
  ggplot(subset(delta_meth, chr_pos == overlap_igg_ams$chr_pos[i]), aes(x = ig_dif_scl, y = scale(delta_meth))) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression("z-trans "*Delta*" igg"), y = expression("z-trans "*Delta*" methylation"),
                                                    subtitle = overlap_igg_ams$chr_pos[i],
                                                    title = paste0("Estimate = ", round(overlap_igg_ams$parameter_estimate[i], 2), ", q-value = ", round(overlap_igg_ams$parameter_qval[i], 4))) +
    geom_abline(intercept=overlap_igg_ams$intercept[i], slope = overlap_igg_ams$parameter_estimate[i], 
                color=clrs_hunting[2], linewidth=1)+
    geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_a
  list_plot_igg_a[[i]] <- plot_a 
  
  ggplot(subset(delta_meth, chr_pos == overlap_igg_ams$chr_pos[i]), aes(x = scale(delta_meth), y = MS)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(y = "AMS", x = expression("z-trans "*Delta*" methylation"),
                                                    subtitle = overlap_igg_ams$chr_pos[i],
                                                    title = paste0("Estimate = ", round(overlap_igg_ams$ams_delta_meth_estimate[i], 2),
                                                                   ", q-value = ", round(overlap_igg_ams$ams_delta_meth_qval[i], 2))) +
    geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
    geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_b
  list_plot_igg_b[[i]] <- plot_b}

cowplot::plot_grid(list_plot_igg_a[[1]], list_plot_igg_b[[1]], 
                   labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_igg_ams

ggsave(plots_igg_ams, file = "plots/model_out/rawdata_plot_overlap_igg_AMS.png", width=14, height=12)                                

### delta microf

overlap_microf_ams <- subset(overlap_ams_physio, parameter == "Delta Microfilaria spp.") #2

list_plot_microf_a <- list()
list_plot_microf_b <- list()

for (i in 1:nrow(overlap_microf_ams)){
  ggplot(subset(delta_meth, chr_pos == overlap_microf_ams$chr_pos[i]), aes(x = microf_dif_scl, y = scale(delta_meth))) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression("z-trans "*Delta*" microf"), y = expression("z-trans "*Delta*" methylation"),
                                                    subtitle = overlap_microf_ams$chr_pos[i],
                                                    title = paste0("Estimate = ", round(overlap_microf_ams$parameter_estimate[i], 2), ", q-value = ", round(overlap_microf_ams$parameter_qval[i], 4))) +
    geom_abline(intercept=overlap_microf_ams$intercept[i], slope = overlap_microf_ams$parameter_estimate[i], 
                color=clrs_hunting[2], linewidth=1)+
    geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_a
  list_plot_microf_a[[i]] <- plot_a 
  
  ggplot(subset(delta_meth, chr_pos == overlap_microf_ams$chr_pos[i]), aes(x = scale(delta_meth), y = MS)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(y = "AMS", x = expression("z-trans "*Delta*" methylation"),
                                                    subtitle = overlap_microf_ams$chr_pos[i],
                                                    title = paste0("Estimate = ", round(overlap_microf_ams$ams_delta_meth_estimate[i], 2),
                                                                   ", q-value = ", round(overlap_microf_ams$ams_delta_meth_qval[i], 2))) +
    geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
    geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_b
  list_plot_microf_b[[i]] <- plot_b}

cowplot::plot_grid(list_plot_microf_a[[1]], list_plot_microf_b[[1]], list_plot_microf_a[[2]], list_plot_microf_b[[2]], 
                   labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_microf_ams

ggsave(plots_microf_ams, file = "plots/model_out/rawdata_plot_overlap_microf_AMS.png", width=14, height=16)                                

### delta trypa

overlap_trypa_ams <- subset(overlap_ams_physio, parameter == "Delta Trypanosoma spp.") #3

list_plot_trypa_a <- list()
list_plot_trypa_b <- list()

for (i in 1:nrow(overlap_trypa_ams)){
  ggplot(subset(delta_meth, chr_pos == overlap_trypa_ams$chr_pos[i]), aes(x = trypa_dif_scl, y = scale(delta_meth))) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression("z-trans "*Delta*" trypa"), y = expression("z-trans "*Delta*" methylation"),
                                                    subtitle = overlap_trypa_ams$chr_pos[i],
                                                    title = paste0("Estimate = ", round(overlap_trypa_ams$parameter_estimate[i], 2), ", q-value = ", round(overlap_trypa_ams$parameter_qval[i], 4))) +
    geom_abline(intercept=overlap_trypa_ams$intercept[i], slope = overlap_trypa_ams$parameter_estimate[i], 
                color=clrs_hunting[2], linewidth=1)+
    geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_a
  list_plot_trypa_a[[i]] <- plot_a 
  
  ggplot(subset(delta_meth, chr_pos == overlap_trypa_ams$chr_pos[i]), aes(x = scale(delta_meth), y = MS)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(y = "AMS", x = expression("z-trans "*Delta*" methylation"),
                                                    subtitle = overlap_trypa_ams$chr_pos[i],
                                                    title = paste0("Estimate = ", round(overlap_trypa_ams$ams_delta_meth_estimate[i], 2),
                                                                   ", q-value = ", round(overlap_trypa_ams$ams_delta_meth_qval[i], 2))) +
    geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
    geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_b
  list_plot_trypa_b[[i]] <- plot_b}

cowplot::plot_grid(list_plot_trypa_a[[1]], list_plot_trypa_b[[1]], list_plot_trypa_a[[2]], list_plot_trypa_b[[2]], 
                   list_plot_trypa_a[[3]], list_plot_trypa_b[[3]], 
                   labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_trypa_ams

ggsave(plots_trypa_ams, file = "plots/model_out/rawdata_plot_overlap_trypa_AMS.png", width=14, height=20)                                

#### Highlighted CpG sites ####

# ScEsiA3_16759__HRSCAF_19053_458 is affected by delta mass, delta igg and delta hct and also explains variation in ams
# ScEsiA3_16870__HRSCAF_19426_63522502 is affected by delta trypa and delta hct and also explains variation in ams

cpg1 <- subset(overlap_ams_physio, chr_pos == "ScEsiA3_16759__HRSCAF_19053_458")

ggplot(subset(delta_meth, chr_pos == "ScEsiA3_16759__HRSCAF_19053_458"), aes(x = mass_dif_scl, y = scale(delta_meth))) + 
  geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression("z-trans "*Delta*" mass"), y = expression("z-trans "*Delta*" methylation"),
                                                  subtitle = "ScEsiA3_16759__HRSCAF_19053_458",
                                                  title = paste0("Estimate = ", round(cpg1$parameter_estimate[which(cpg1$parameter=="Delta body mass")], 2), 
                                                                 ", q-value = ", round(cpg1$parameter_qval[which(cpg1$parameter=="Delta body mass")], 4))) +
  geom_abline(intercept=cpg1$intercept[which(cpg1$parameter=="Delta body mass")], slope = cpg1$parameter_estimate[which(cpg1$parameter=="Delta body mass")], 
              color=clrs_hunting[2], linewidth=1)+
  geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_cpg1_mass

ggplot(subset(delta_meth, chr_pos == "ScEsiA3_16759__HRSCAF_19053_458"), aes(x = ig_dif_scl, y = scale(delta_meth))) + 
  geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression("z-trans "*Delta*" IgG"), y = expression("z-trans "*Delta*" methylation"),
                                               #   subtitle = "ScEsiA3_16759__HRSCAF_19053_458",
                                                  title = paste0("Estimate = ", round(cpg1$parameter_estimate[which(cpg1$parameter=="Delta IgG")], 2), 
                                                                 ", q-value = ", round(cpg1$parameter_qval[which(cpg1$parameter=="Delta IgG")], 4))) +
  geom_abline(intercept=cpg1$intercept[which(cpg1$parameter=="Delta IgG")], slope = cpg1$parameter_estimate[which(cpg1$parameter=="Delta IgG")], 
              color=clrs_hunting[2], linewidth=1)+
  geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_cpg1_igg

ggplot(subset(delta_meth, chr_pos == "ScEsiA3_16759__HRSCAF_19053_458"), aes(x = hct_dif_scl, y = scale(delta_meth))) + 
  geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression("z-trans "*Delta*" HCT"), y = expression("z-trans "*Delta*" methylation"),
                                               #   subtitle = "ScEsiA3_16759__HRSCAF_19053_458",
                                                  title = paste0("Estimate = ", round(cpg1$parameter_estimate[which(cpg1$parameter=="Delta HCT")], 2), 
                                                                 ", q-value = ", round(cpg1$parameter_qval[which(cpg1$parameter=="Delta HCT")], 4))) +
  geom_abline(intercept=cpg1$intercept[which(cpg1$parameter=="Delta HCT")], slope = cpg1$parameter_estimate[which(cpg1$parameter=="Delta HCT")], 
              color=clrs_hunting[2], linewidth=1)+
  geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_cpg1_hct

ggplot(subset(delta_meth, chr_pos == "ScEsiA3_16759__HRSCAF_19053_458"), aes(x = scale(delta_meth), y = MS)) + 
  geom_point(fill=clrs_hunting[1], size=3) + labs(y = "AMS", x = expression("z-trans "*Delta*" methylation"),
                                               #   subtitle = "ScEsiA3_16759__HRSCAF_19053_458",
                                                  title = paste0("Estimate = ", round(cpg1$ams_delta_meth_estimate[1], 2),
                                                                 ", q-value = ", round(cpg1$ams_delta_meth_qval[1], 2))) +
  geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
  geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_cpg1_ams


cowplot::plot_grid(plot_cpg1_mass, plot_cpg1_hct, plot_cpg1_igg, plot_cpg1_ams, 
                   labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_cpg1

ggsave(plots_cpg1, file = "plots/model_out/rawdata_plot_overlap_highlight_CpG_1.png", width=14, height=20)                                

## cpg 2

cpg2 <- subset(overlap_ams_physio, chr_pos == "ScEsiA3_16870__HRSCAF_19426_63522502")

ggplot(subset(delta_meth, chr_pos == "ScEsiA3_16870__HRSCAF_19426_63522502"), aes(x = trypa_dif_scl, y = scale(delta_meth))) + 
  geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression("z-trans "*Delta*" Trypanosoma spp."), y = expression("z-trans "*Delta*" methylation"),
                                                #  subtitle = "ScEsiA3_16870__HRSCAF_19426_63522502",
                                                  title = paste0("Estimate = ", round(cpg2$parameter_estimate[which(cpg2$parameter=="Delta Trypanosoma spp.")], 2), 
                                                                 ", q-value = ", round(cpg2$parameter_qval[which(cpg2$parameter=="Delta Trypanosoma spp.")], 4))) +
  geom_abline(intercept=cpg2$intercept[which(cpg2$parameter=="Delta Trypanosoma spp.")], slope = cpg2$parameter_estimate[which(cpg2$parameter=="Delta Trypanosoma spp.")], 
              color=clrs_hunting[2], linewidth=1)+
  geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_cpg2_trypa

ggplot(subset(delta_meth, chr_pos == "ScEsiA3_16870__HRSCAF_19426_63522502"), aes(x = hct_dif_scl, y = scale(delta_meth))) + 
  geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression("z-trans "**Delta*" HCT"), y = expression("z-trans "*Delta*" methylation"),
                                                #  subtitle = "ScEsiA3_16870__HRSCAF_19426_63522502",
                                                  title = paste0("Estimate = ", round(cpg2$parameter_estimate[which(cpg2$parameter=="Delta HCT")], 2), 
                                                                 ", q-value = ", round(cpg2$parameter_qval[which(cpg2$parameter=="Delta HCT")], 4))) +
  geom_abline(intercept=cpg2$intercept[which(cpg2$parameter=="Delta HCT")], slope = cpg2$parameter_estimate[which(cpg2$parameter=="Delta HCT")], 
              color=clrs_hunting[2], linewidth=1)+
  geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_cpg2_hct

ggplot(subset(delta_meth, chr_pos == "ScEsiA3_16870__HRSCAF_19426_63522502"), aes(x = scale(delta_meth), y = MS)) + 
  geom_point(fill=clrs_hunting[1], size=3) + labs(y = "AMS", x = expression("z-trans "*Delta*" methylation"),
                                                 # subtitle = "ScEsiA3_16870__HRSCAF_19426_63522502",
                                                  title = paste0("Estimate = ", round(cpg2$ams_delta_meth_estimate[1], 2),
                                                                 ", q-value = ", round(cpg2$ams_delta_meth_qval[1], 2))) +
  geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
  geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_cpg2_ams

cowplot::plot_grid(plot_cpg2_trypa, plot_cpg2_hct, plot_cpg2_ams, 
                   labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_cpg2

ggsave(plots_cpg2, file = "plots/model_out/rawdata_plot_overlap_highlight_CpG_2.png", width=14, height=20)                                

