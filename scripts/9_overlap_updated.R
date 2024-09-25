
### load packages
pacman::p_load(tidyverse, data.table, cowplot)

source("scripts/plotting_theme.R")

### load data

load(file="results/modeloutput/changing/changing_sites_glmer.RData")

### methylation difference

load(file = "results/modeloutput/all_sites_deltameth.RData")

delta_meth <- subset(delta_meth, chr_pos %in% changing_cpg$chr_pos)

### load phenotypic data

load("data/phenotypes/fulldata_complete_epi_withdates.RData")
load("data/phenotypes/pheno_dif_prepost.RData")

# combine with site and fitness data
pheno_pre <- subset(all_pheno_epi, prepost=="pre")

delta_meth <- left_join(delta_meth, unique(pheno_pre[,c("id", "year", "MS", "surv", "site")]), by = c("id", "year"))

### combine delta methylation data with site and behaviour info
effort <- all_pheno_epi %>% dplyr::select(c("id", "year", "attend", "fight", "dist", "MS")) %>% filter(!is.na(attend)) %>% unique()

effort <- subset(effort, id %in% delta_meth$id)

effort$attend_scl <- scale(effort$attend)
effort$fight_scl <- scale(effort$fight)
effort$dist_scl <- scale(effort$dist)

### combine reproductive effort data with methylation data

delta_meth <- left_join(delta_meth, effort[,c("id", "year", "attend", "fight", "dist", "attend_scl", "fight_scl", "dist_scl")], by = c("id", "year"))

### combine physio data with methylation data 

prepost_dif$mass_dif_scl <- scale(prepost_dif$mass_dif)
prepost_dif$microf_dif_scl <- scale(prepost_dif$microf_dif)
prepost_dif$trypa_dif_scl <- scale(prepost_dif$trypa_dif)
prepost_dif$ig_dif_scl <- scale(prepost_dif$ig_dif)
prepost_dif$hct_dif_scl <- scale(prepost_dif$hct_dif)

delta_meth <- left_join(delta_meth, unique(prepost_dif[,c("id", "year", "mass_dif", "microf_dif", "trypa_dif", "ig_dif", "hct_dif",
                        "mass_dif_scl", "microf_dif_scl", "trypa_dif_scl", "ig_dif_scl", "hct_dif_scl")], by = c("id", "year")))

##### Overlaps #####

## load ams model results
load(file="results/modeloutput/fitness/out_ams_deltameth_filtered.RData")
cpg_sig_ams <- subset(delta_out_ams, ams_delta_meth_qval < 0.05) #362
cpg_sig_ams <- cpg_sig_ams %>% select(c(chr_pos, intercept_ams, ams_delta_meth_estimate, ams_delta_meth_qval))
names(cpg_sig_ams) <- c("chr_pos", "intercept", "parameter_estimate", "parameter_qval")
cpg_sig_ams <- cpg_sig_ams %>% mutate(parameter = "ams", .before=parameter_qval)
cpg_sig_ams$pre_control <- NA

## load all other trait results ###

## attendance
load(file="results/modeloutput/effort/out_attend_with_pre.RData")
cpg_attend_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, intercept, parameter_estimate, parameter, parameter_qval))
cpg_attend_pre$pre_control <- "pre"

load(file="results/modeloutput/effort/out_attend_no_pre.RData")
cpg_attend_no_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos,intercept, parameter_estimate,  parameter, parameter_qval))
cpg_attend_no_pre$pre_control <- "no_pre"

## fight 
load(file="results/modeloutput/effort/out_fight_with_pre.RData")
cpg_fight_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, intercept, parameter_estimate, parameter, parameter_qval))
cpg_fight_pre$pre_control <- "pre"

load(file="results/modeloutput/effort/out_fight_no_pre.RData")
cpg_fight_no_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, intercept, parameter_estimate, parameter, parameter_qval))
cpg_fight_no_pre$pre_control <- "no_pre"

## dist 
load(file="results/modeloutput/effort/out_dist_with_pre.RData")
cpg_dist_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, intercept, parameter_estimate, parameter, parameter_qval))
cpg_dist_pre$pre_control <- "pre"

load(file="results/modeloutput/effort/out_dist_no_pre.RData")
cpg_dist_no_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos,intercept, parameter_estimate,  parameter, parameter_qval))
cpg_dist_no_pre$pre_control <- "no_pre"

## mass 
load(file="results/modeloutput/physio/out_mass_dif_with_pre.RData")
cpg_mass_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos,intercept, parameter_estimate,  parameter, parameter_qval))
cpg_mass_pre$pre_control <- "pre"

load(file="results/modeloutput/physio/out_mass_dif_no_pre.RData")
cpg_mass_no_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos,intercept, parameter_estimate,  parameter, parameter_qval))
cpg_mass_no_pre$pre_control <- "no_pre"

## microf 
load(file="results/modeloutput/physio/out_microf_dif_with_pre.RData")
cpg_microf_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, intercept, parameter_estimate, parameter, parameter_qval))
cpg_microf_pre$pre_control <- "pre"

load(file="results/modeloutput/physio/out_microf_dif_no_pre.RData")
cpg_microf_no_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos,intercept, parameter_estimate,  parameter, parameter_qval))
cpg_microf_no_pre$pre_control <- "no_pre"

## trypa 
load(file="results/modeloutput/physio/out_trypa_dif_with_pre.RData")
cpg_trypa_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, intercept, parameter_estimate, parameter, parameter_qval))
cpg_trypa_pre$pre_control <- "pre"

load(file="results/modeloutput/physio/out_trypa_dif_no_pre.RData")
cpg_trypa_no_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos,intercept, parameter_estimate,  parameter, parameter_qval))
cpg_trypa_no_pre$pre_control <- "no_pre"

## hct 
load(file="results/modeloutput/physio/out_hct_dif_with_pre.RData")
cpg_hct_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos,intercept, parameter_estimate,  parameter, parameter_qval))
cpg_hct_pre$pre_control <- "pre"

load(file="results/modeloutput/physio/out_hct_dif_no_pre.RData")
cpg_hct_no_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, intercept, parameter_estimate, parameter, parameter_qval))
cpg_hct_no_pre$pre_control <- "no_pre"

## ig 
load(file="results/modeloutput/physio/out_ig_dif_with_pre.RData")
cpg_ig_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, intercept, parameter_estimate, parameter, parameter_qval))
cpg_ig_pre$pre_control <- "pre"

load(file="results/modeloutput/physio/out_ig_dif_no_pre.RData")
cpg_ig_no_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos,intercept, parameter_estimate,  parameter, parameter_qval))
cpg_ig_no_pre$pre_control <- "no_pre"

### combine

all_models_sig <- rbind(cpg_attend_no_pre, cpg_attend_pre,
                        cpg_fight_no_pre, cpg_fight_pre,
                        cpg_dist_no_pre, cpg_dist_pre,
                        cpg_mass_no_pre, cpg_mass_pre,
                        cpg_microf_no_pre, cpg_microf_pre,
                        cpg_trypa_no_pre, cpg_trypa_pre,
                        cpg_ig_no_pre, cpg_ig_pre,
                        cpg_hct_no_pre, cpg_hct_pre)



all_models <- rbind(cpg_sig_ams, all_models_sig)

### identify duplicate CpG sites / overlap ####
all_models <- all_models %>% group_by(chr_pos) %>% mutate(n = row_number()) %>% ungroup()
dups <- all_models %>% subset(n > 1) %>% select(chr_pos) %>% unique()

dup_models <- subset(all_models, chr_pos %in% dups$chr_pos) %>% select(-c(n))%>% arrange(chr_pos)

## 9 overlapping!

### Plot raw data ####
#### CpG 1 ####
model_cpg_1 <- subset(dup_models, chr_pos == dups$chr_pos[1]) # attend and AMS
cpg_1 <- subset(delta_meth, chr_pos == dups$chr_pos[1])

ggplot(cpg_1, aes(x = attend_scl, y = delta_meth)) + 
        geom_point(fill=clrs_hunting[1], size=3) + labs(x = "z-transformed attendance", y = expression(Delta*" methylation"),
                     title = paste0("Estimate = ", round(model_cpg_1$parameter_estimate[which(model_cpg_1$parameter=="attend")], 2), 
                                    ", q-value = ", round(model_cpg_1$parameter_qval[which(model_cpg_1$parameter=="attend")], 4)), 
                                    subtitle = cpg_1$chr_pos[1]) +
                                     geom_abline(intercept=model_cpg_1$intercept[which(model_cpg_1$parameter=="attend")], 
                                     slope = model_cpg_1$parameter_estimate[which(model_cpg_1$parameter=="attend")], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_cpg_1_a


ggplot(cpg_1, aes(x = delta_meth, y = MS)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(y = "AMS", x = expression(Delta*" methylation"),
                      title = paste0("Estimate = ", round(model_cpg_1$parameter_estimate[which(model_cpg_1$parameter=="ams")], 2),
                                        ", q-value = ", round(model_cpg_1$parameter_qval[which(model_cpg_1$parameter=="ams")], 2)), 
                      subtitle = cpg_1$chr_pos[1]) +
                                        geom_abline(intercept = model_cpg_1$intercept[which(model_cpg_1$parameter=="ams")],
                                        slope = model_cpg_1$parameter_estimate[which(model_cpg_1$parameter=="ams")], 
                                        color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_cpg_1_b

cowplot::plot_grid(plot_cpg_1_a, plot_cpg_1_b, 
    labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plot_cpg_1

ggsave(plot_cpg_1, file = "plots/model_out/overlap/rawdata_overlap_cpg_1.png", width=14, height=12)                                

#### CpG 2 ####

model_cpg_2 <- subset(dup_models, chr_pos == dups$chr_pos[2]) # dist, ig_dif and ams, both pre and non-pre for both dist and ig_dif. choose 'pre'
cpg_2 <- subset(delta_meth, chr_pos == dups$chr_pos[2])

ggplot(cpg_2, aes(x = dist_scl, y = delta_meth)) + 
        geom_point(fill=clrs_hunting[1], size=3) + labs(x = "z-transformed centrality", y = expression(Delta*" methylation"),
                     title = paste0("Estimate = ", round(model_cpg_2$parameter_estimate[which(model_cpg_2$parameter=="dist" & model_cpg_2$pre_control == "pre")], 2), 
                                    ", q-value = ", round(model_cpg_2$parameter_qval[which(model_cpg_2$parameter=="dist" & model_cpg_2$pre_control == "pre")], 4)), 
                                    subtitle = cpg_2$chr_pos[1]) +
                                     geom_abline(intercept=model_cpg_2$intercept[which(model_cpg_2$parameter=="dist" & model_cpg_2$pre_control == "pre")], 
                                     slope = model_cpg_2$parameter_estimate[which(model_cpg_2$parameter=="dist" & model_cpg_2$pre_control == "pre")], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_cpg_2_a

ggplot(cpg_2, aes(x = ig_dif_scl, y = delta_meth)) + 
        geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression("z-transformed "*Delta*" IgG"), y = expression(Delta*" methylation"),
                     title = paste0("Estimate = ", round(model_cpg_2$parameter_estimate[which(model_cpg_2$parameter=="ig_dif" & model_cpg_2$pre_control == "pre")], 2), 
                                    ", q-value = ", round(model_cpg_2$parameter_qval[which(model_cpg_2$parameter=="ig_dif" & model_cpg_2$pre_control == "pre")], 4)), 
                                    subtitle = cpg_2$chr_pos[1]) +
                                     geom_abline(intercept=model_cpg_2$intercept[which(model_cpg_2$parameter=="ig_dif" & model_cpg_2$pre_control == "pre")], 
                                     slope = model_cpg_2$parameter_estimate[which(model_cpg_2$parameter=="ig_dif" & model_cpg_2$pre_control == "pre")], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_cpg_2_b
                                      

ggplot(cpg_2, aes(x = delta_meth, y = MS)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(y = "AMS", x = expression(Delta*" methylation"),
                      title = paste0("Estimate = ", round(model_cpg_2$parameter_estimate[which(model_cpg_2$parameter=="ams")], 2),
                                        ", q-value = ", round(model_cpg_2$parameter_qval[which(model_cpg_2$parameter=="ams")], 2)), 
                      subtitle = cpg_2$chr_pos[1]) +
                                        geom_abline(intercept = model_cpg_2$intercept[which(model_cpg_2$parameter=="ams")],
                                        slope = model_cpg_2$parameter_estimate[which(model_cpg_2$parameter=="ams")], 
                                        color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_cpg_2_c

cowplot::plot_grid(plot_cpg_2_a, plot_cpg_2_b, plot_cpg_2_c,
    labs="auto", align="hv", axis="lb", ncol=3, label_fontface = "plain", label_size = 22) -> plot_cpg_2

ggsave(plot_cpg_2, file = "plots/model_out/overlap/rawdata_overlap_cpg_2.png", width=18, height=12)                                

#### CpG 3 ####
model_cpg_3 <- subset(dup_models, chr_pos == dups$chr_pos[3]) # mass_dif (no_pre) and AMS
cpg_3 <- subset(delta_meth, chr_pos == dups$chr_pos[3])

ggplot(cpg_3, aes(x = mass_dif_scl, y = delta_meth)) + 
        geom_point(fill=clrs_hunting[1], size=3) + labs(x =  expression("z-transformed "*Delta*" body mass"), y = expression(Delta*" methylation"),
                     title = paste0("Estimate = ", round(model_cpg_3$parameter_estimate[which(model_cpg_3$parameter=="mass_dif")], 2), 
                                    ", q-value = ", round(model_cpg_3$parameter_qval[which(model_cpg_3$parameter=="mass_dif")], 4)), 
                                    subtitle = cpg_3$chr_pos[1]) +
                                     geom_abline(intercept=model_cpg_3$intercept[which(model_cpg_3$parameter=="mass_dif")], 
                                     slope = model_cpg_3$parameter_estimate[which(model_cpg_3$parameter=="mass_dif")], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_cpg_3_a


ggplot(cpg_3, aes(x = delta_meth, y = MS)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(y = "AMS", x = expression(Delta*" methylation"),
                      title = paste0("Estimate = ", round(model_cpg_3$parameter_estimate[which(model_cpg_3$parameter=="ams")], 2),
                                        ", q-value = ", round(model_cpg_3$parameter_qval[which(model_cpg_3$parameter=="ams")], 2)), 
                      subtitle = cpg_3$chr_pos[1]) +
                                        geom_abline(intercept = model_cpg_3$intercept[which(model_cpg_3$parameter=="ams")],
                                        slope = model_cpg_3$parameter_estimate[which(model_cpg_3$parameter=="ams")], 
                                        color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_cpg_3_b

cowplot::plot_grid(plot_cpg_3_a, plot_cpg_3_b, 
    labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plot_cpg_3

ggsave(plot_cpg_3, file = "plots/model_out/overlap/rawdata_overlap_cpg_3.png", width=14, height=12)                                


#### CpG 4 ####
model_cpg_4 <- subset(dup_models, chr_pos == dups$chr_pos[4]) # mass_dif and AMS, both in pre and no_pre, choose 'pre'
cpg_4 <- subset(delta_meth, chr_pos == dups$chr_pos[4])


ggplot(cpg_4, aes(x = mass_dif_scl, y = delta_meth)) + 
        geom_point(fill=clrs_hunting[1], size=3) + labs(x =  expression("z-transformed "*Delta*" body mass"), y = expression(Delta*" methylation"),
                     title = paste0("Estimate = ", round(model_cpg_4$parameter_estimate[which(model_cpg_4$parameter=="mass_dif"& model_cpg_4$pre_control == "pre")], 2), 
                                    ", q-value = ", round(model_cpg_4$parameter_qval[which(model_cpg_4$parameter=="mass_dif"& model_cpg_4$pre_control == "pre")], 4)), 
                                    subtitle = cpg_4$chr_pos[1]) +
                                     geom_abline(intercept=model_cpg_4$intercept[which(model_cpg_4$parameter=="mass_dif"& model_cpg_4$pre_control == "pre")], 
                                     slope = model_cpg_4$parameter_estimate[which(model_cpg_4$parameter=="mass_dif"& model_cpg_4$pre_control == "pre")], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_cpg_4_a


ggplot(cpg_4, aes(x = delta_meth, y = MS)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(y = "AMS", x = expression(Delta*" methylation"),
                      title = paste0("Estimate = ", round(model_cpg_4$parameter_estimate[which(model_cpg_4$parameter=="ams")], 2),
                                        ", q-value = ", round(model_cpg_4$parameter_qval[which(model_cpg_4$parameter=="ams")], 2)), 
                      subtitle = cpg_4$chr_pos[1]) +
                                        geom_abline(intercept = model_cpg_4$intercept[which(model_cpg_4$parameter=="ams")],
                                        slope = model_cpg_4$parameter_estimate[which(model_cpg_4$parameter=="ams")], 
                                        color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_cpg_4_b

cowplot::plot_grid(plot_cpg_4_a, plot_cpg_4_b, 
    labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plot_cpg_4

ggsave(plot_cpg_4, file = "plots/model_out/overlap/rawdata_overlap_cpg_4.png", width=14, height=12)                                                  

#### CpG 5 ####
model_cpg_5 <- subset(dup_models, chr_pos == dups$chr_pos[5]) # attendance and mass_dif
cpg_5 <- subset(delta_meth, chr_pos == dups$chr_pos[5])

ggplot(cpg_5, aes(x = mass_dif_scl, y = delta_meth)) + 
        geom_point(fill=clrs_hunting[1], size=3) + labs(x =  expression("z-transformed "*Delta*" body mass"), y = expression(Delta*" methylation"),
                     title = paste0("Estimate = ", round(model_cpg_5$parameter_estimate[which(model_cpg_5$parameter=="mass_dif")], 2), 
                                    ", q-value = ", round(model_cpg_5$parameter_qval[which(model_cpg_5$parameter=="mass_dif")], 4)), 
                                    subtitle = cpg_5$chr_pos[1]) +
                                     geom_abline(intercept=model_cpg_5$intercept[which(model_cpg_5$parameter=="mass_dif")], 
                                     slope = model_cpg_5$parameter_estimate[which(model_cpg_5$parameter=="mass_dif")], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_cpg_5_a


ggplot(cpg_5, aes(x = attend_scl, y = delta_meth)) + 
        geom_point(fill=clrs_hunting[1], size=3) + labs(x = "z-transformed attendance", y = expression(Delta*" methylation"),
                     title = paste0("Estimate = ", round(model_cpg_5$parameter_estimate[which(model_cpg_5$parameter=="attend")], 2), 
                                    ", q-value = ", round(model_cpg_5$parameter_qval[which(model_cpg_5$parameter=="attend")], 4)), 
                                    subtitle = cpg_5$chr_pos[1]) +
                                     geom_abline(intercept=model_cpg_5$intercept[which(model_cpg_5$parameter=="attend")], 
                                     slope = model_cpg_5$parameter_estimate[which(model_cpg_5$parameter=="attend")], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_cpg_5_b

cowplot::plot_grid(plot_cpg_5_a, plot_cpg_5_b, 
    labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plot_cpg_5

ggsave(plot_cpg_5, file = "plots/model_out/overlap/rawdata_overlap_cpg_5.png", width=14, height=12)                                

#### CpG 6 ####
model_cpg_6 <- subset(dup_models, chr_pos == dups$chr_pos[6]) # microf_dif (pre) and ams
cpg_6 <- subset(delta_meth, chr_pos == dups$chr_pos[6])

ggplot(cpg_6, aes(x = microf_dif_scl, y = delta_meth)) + 
        geom_point(fill=clrs_hunting[1], size=3) + labs(x =  expression("z-transformed "*Delta*" Microfilaria spp."), y = expression(Delta*" methylation"),
                     title = paste0("Estimate = ", round(model_cpg_6$parameter_estimate[which(model_cpg_6$parameter=="microf_dif"& model_cpg_6$pre_control == "pre")], 2), 
                                    ", q-value = ", round(model_cpg_6$parameter_qval[which(model_cpg_6$parameter=="microf_dif"& model_cpg_6$pre_control == "pre")], 4)), 
                                    subtitle = cpg_6$chr_pos[1]) +
                                     geom_abline(intercept=model_cpg_6$intercept[which(model_cpg_6$parameter=="microf_dif"& model_cpg_6$pre_control == "pre")], 
                                     slope = model_cpg_6$parameter_estimate[which(model_cpg_6$parameter=="microf_dif"& model_cpg_6$pre_control == "pre")], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_cpg_6_a


ggplot(cpg_6, aes(x = delta_meth, y = MS)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(y = "AMS", x = expression(Delta*" methylation"),
                      title = paste0("Estimate = ", round(model_cpg_6$parameter_estimate[which(model_cpg_6$parameter=="ams")], 2),
                                        ", q-value = ", round(model_cpg_6$parameter_qval[which(model_cpg_6$parameter=="ams")], 2)), 
                      subtitle = cpg_6$chr_pos[1]) +
                                        geom_abline(intercept = model_cpg_6$intercept[which(model_cpg_6$parameter=="ams")],
                                        slope = model_cpg_6$parameter_estimate[which(model_cpg_6$parameter=="ams")], 
                                        color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_cpg_6_b

cowplot::plot_grid(plot_cpg_6_a, plot_cpg_6_b, 
    labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plot_cpg_6

ggsave(plot_cpg_6, file = "plots/model_out/overlap/rawdata_overlap_cpg_6.png", width=14, height=12)    

#### CpG 7 ####
model_cpg_7 <- subset(dup_models, chr_pos == dups$chr_pos[7]) # microf_dif (pre) and ams
cpg_7 <- subset(delta_meth, chr_pos == dups$chr_pos[7])

ggplot(cpg_7, aes(x = microf_dif_scl, y = delta_meth)) + 
        geom_point(fill=clrs_hunting[1], size=3) + labs(x =  expression("z-transformed "*Delta*" Microfilaria spp."), y = expression(Delta*" methylation"),
                     title = paste0("Estimate = ", round(model_cpg_7$parameter_estimate[which(model_cpg_7$parameter=="microf_dif"& model_cpg_7$pre_control == "pre")], 2), 
                                    ", q-value = ", round(model_cpg_7$parameter_qval[which(model_cpg_7$parameter=="microf_dif"& model_cpg_7$pre_control == "pre")], 4)), 
                                    subtitle = cpg_7$chr_pos[1]) +
                                     geom_abline(intercept=model_cpg_7$intercept[which(model_cpg_7$parameter=="microf_dif"& model_cpg_7$pre_control == "pre")], 
                                     slope = model_cpg_7$parameter_estimate[which(model_cpg_7$parameter=="microf_dif"& model_cpg_7$pre_control == "pre")], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_cpg_7_a


ggplot(cpg_7, aes(x = delta_meth, y = MS)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(y = "AMS", x = expression(Delta*" methylation"),
                      title = paste0("Estimate = ", round(model_cpg_7$parameter_estimate[which(model_cpg_7$parameter=="ams")], 2),
                                        ", q-value = ", round(model_cpg_7$parameter_qval[which(model_cpg_7$parameter=="ams")], 2)), 
                      subtitle = cpg_7$chr_pos[1]) +
                                        geom_abline(intercept = model_cpg_7$intercept[which(model_cpg_7$parameter=="ams")],
                                        slope = model_cpg_7$parameter_estimate[which(model_cpg_7$parameter=="ams")], 
                                        color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_cpg_7_b

cowplot::plot_grid(plot_cpg_7_a, plot_cpg_7_b, 
    labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plot_cpg_7

ggsave(plot_cpg_7, file = "plots/model_out/overlap/rawdata_overlap_cpg_7.png", width=14, height=12)    

#### CpG 8 ####
model_cpg_8 <- subset(dup_models, chr_pos == dups$chr_pos[8]) # igg dif with and without pre, skip

#### CpG 9 ####
model_cpg_9 <- subset(dup_models, chr_pos == dups$chr_pos[9]) # hct_dif (pre) and AMS
cpg_9 <- subset(delta_meth, chr_pos == dups$chr_pos[9])

ggplot(cpg_9, aes(x = hct_dif_scl, y = delta_meth)) + 
        geom_point(fill=clrs_hunting[1], size=3) + labs(x =  expression("z-transformed "*Delta*" HCT"), y = expression(Delta*" methylation"),
                     title = paste0("Estimate = ", round(model_cpg_9$parameter_estimate[which(model_cpg_9$parameter=="hct_dif"& model_cpg_9$pre_control == "pre")], 2), 
                                    ", q-value = ", round(model_cpg_9$parameter_qval[which(model_cpg_9$parameter=="hct_dif"& model_cpg_9$pre_control == "pre")], 4)), 
                                    subtitle = cpg_9$chr_pos[1]) +
                                     geom_abline(intercept=model_cpg_9$intercept[which(model_cpg_9$parameter=="hct_dif"& model_cpg_9$pre_control == "pre")], 
                                     slope = model_cpg_9$parameter_estimate[which(model_cpg_9$parameter=="hct_dif"& model_cpg_9$pre_control == "pre")], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_cpg_9_a

ggplot(cpg_9, aes(x = delta_meth, y = MS)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(y = "AMS", x = expression(Delta*" methylation"),
                      title = paste0("Estimate = ", round(model_cpg_9$parameter_estimate[which(model_cpg_9$parameter=="ams")], 2),
                                        ", q-value = ", round(model_cpg_9$parameter_qval[which(model_cpg_9$parameter=="ams")], 2)), 
                      subtitle = cpg_9$chr_pos[1]) +
                                        geom_abline(intercept = model_cpg_9$intercept[which(model_cpg_9$parameter=="ams")],
                                        slope = model_cpg_9$parameter_estimate[which(model_cpg_9$parameter=="ams")], 
                                        color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_cpg_9_b

cowplot::plot_grid(plot_cpg_9_a, plot_cpg_9_b, 
    labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plot_cpg_9

ggsave(plot_cpg_9, file = "plots/model_out/overlap/rawdata_overlap_cpg_8.png", width=14, height=12)    