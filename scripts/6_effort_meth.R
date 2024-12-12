### load packages
pacman::p_load(tidyverse, data.table, tibble, performance, gaston,
               parallel, performance, lmerTest, tidystats, ggpointdensity)
               
source("scripts/plotting_theme.R")

### load data

load(file="results/modeloutput/changing/changing_sites_glmer.RData")

### load phenotypic data

load(file = "data/phenotypes/fulldata_complete_epi_withdates.RData")

### methylation difference

load(file = "results/modeloutput/all_sites_deltameth.RData")

delta_meth <- subset(delta_meth, chr_pos %in% changing_cpg$chr_pos)

### plot pre_meth effect on delta_meth ####

lmer_pre_delta <- lmerTest::lmer(delta_meth ~ methperc_pre + (1|id), data = delta_meth)

sum_pre_delta <- as.data.frame(summary(lmer_pre_delta)$coef) # negative relationship

ggplot(delta_meth, aes(methperc_pre, delta_meth)) + 
  geom_pointdensity() + 
  scale_color_viridis_c() + 
  geom_abline(intercept=sum_pre_delta$Estimate[1], slope = sum_pre_delta$Estimate[2], 
                                          color="red", linewidth=1)+
  labs(x = "Methylation % pre-lekking", y = expression(Delta*" methylation %")) -> cor_pre_delta

ggplot(delta_meth, aes(methperc_pre, abs(delta_meth))) + 
  geom_pointdensity() + 
  scale_color_viridis_c() + geom_smooth() + 
  labs(x = "Methylation % pre-lekking", y = expression("Absolute "*Delta*" methylation %")) -> cor_pre_delta_abs

cowplot::plot_grid(cor_pre_delta, cor_pre_delta_abs, labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) %>%
  ggsave(file = "plots/explore/pre_vs_delta.png", width=14, height=10)

### combine delta methylation data with site and behaviour info

delta_meth <- left_join(delta_meth, unique(all_pheno_epi[,c("id", "year", "site", "Core")], by = c("id", "year")))

### z-transform the traits before the model that subsets IDs and years where there is
### data for that CpG site

effort <- all_pheno_epi %>% dplyr::select(c("id", "year", "attend", "fight", "dist", "MS")) %>% filter(!is.na(attend)) %>% unique()

effort <- subset(effort, id %in% delta_meth$id)

effort$attend_scl <- scale(effort$attend)
effort$fight_scl <- scale(effort$fight)
effort$dist_scl <- scale(effort$dist)

### combine reproductive effort data with methylation data

delta_meth <- left_join(delta_meth, effort[,c("id", "year", "attend", "fight", "dist", "attend_scl", "fight_scl", "dist_scl")], by = c("id", "year"))
                                           
#### run the model per trait: without pre-lekking ####
source("scripts/function_models.R")


# attendance 

### only select cpg sites with enough data
delta_meth_n_attend <- delta_meth %>% group_by(chr_pos) %>% filter(!is.na(delta_meth)& !is.na(attend)) %>% tally()
delta_meth_n_attend <- subset(delta_meth_n_attend, n > 20)
 
delta_meth_sub_attend <- subset(delta_meth, chr_pos %in% delta_meth_n_attend$chr_pos)
length(unique(delta_meth_sub_attend$chr_pos)) # 808 sites

delta_meth_attend_ls <- delta_meth_sub_attend %>% group_split(chr_pos)

## no pre
m_attend_no_pre <- parallel::mclapply(delta_meth_attend_ls, function_model_delta_pheno, parameter="attend", pre="no_control", mc.cores=4)
m_attend_no_pre_out <- function_process_model(m_attend_no_pre, dir_plots = "plots/model_out/effort", dir_data = "results/modeloutput/effort",
                                            name_file = "attend_no_pre", pretty_name = "Attendance", filter_disp=FALSE) #no sig

nrow(m_attend_no_pre_out$sig) #0

## with pre
m_attend_pre <- parallel::mclapply(delta_meth_attend_ls, function_model_delta_pheno, parameter="attend", pre="control", mc.cores=4)
m_attend_pre_out <- function_process_model(m_attend_pre, dir_plots = "plots/model_out/effort", dir_data = "results/modeloutput/effort",
                                            name_file = "attend_with_pre", pretty_name = "Attendance", filter_disp=FALSE) # n = 3

nrow(m_attend_pre_out$sig) #n=3
#nrow(subset(m_attend_pre_out$sig, chr_pos %in% m_attend_no_pre_out$sig$chr_pos))

#plot raw
system(paste0("rm -i ", getwd(), "/plots/model_out/effort/raw/attend_with_pre_*"))

for (i in 1:nrow(m_attend_pre_out$sig)){
    ggplot(subset(delta_meth, chr_pos == m_attend_pre_out$sig$chr_pos[i]), aes(x = attend_scl, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = "z-transformed attendance", y = expression(Delta*" methylation"),
                      subtitle = m_attend_pre_out$sig$chr_pos[i],
                      title = paste0("Estimate = ", round(m_attend_pre_out$sig$parameter_estimate[i], 2),
                                        ", q-value = ", round(m_attend_pre_out$sig$parameter_qval[i], 2))) +
                                         geom_abline(intercept=m_attend_pre_out$sig$intercept[i], slope = m_attend_pre_out$sig$parameter_estimate[i], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
    ggsave(plot, file = paste0("plots/model_out/effort/raw/attend_with_pre_", i, ".png"), width=8, height=8)
}

# fighting 
### only select cpg sites with enough data
delta_meth_n_fight <- delta_meth %>% group_by(chr_pos) %>% filter(!is.na(delta_meth)& !is.na(fight)) %>% tally()
delta_meth_n_fight <- subset(delta_meth_n_fight, n > 20)
 
delta_meth_sub_fight <- subset(delta_meth, chr_pos %in% delta_meth_n_fight$chr_pos)
length(unique(delta_meth_sub_fight$chr_pos)) # 466 sites

delta_meth_fight_ls <- delta_meth_sub_fight %>% group_split(chr_pos)

## no pre
m_fight_no_pre <- parallel::mclapply(delta_meth_fight_ls, function_model_delta_pheno, parameter="fight", pre="no_control", mc.cores=4)
m_fight_no_pre_out <- function_process_model(m_fight_no_pre, dir_plots = "plots/model_out/effort", dir_data = "results/modeloutput/effort",
                                            name_file = "fight_no_pre", pretty_name = "Fighting", filter_disp=FALSE) 

nrow(m_fight_no_pre_out$sig) #0

system(paste0("rm ", getwd(), "/plots/model_out/effort/raw/fight_no_pre_*"))

#for (i in 1:nrow(m_fight_no_pre_out$sig)){
#    ggplot(subset(delta_meth, chr_pos == m_fight_no_pre_out$sig$chr_pos[i]), aes(x = fight_scl, y = delta_meth)) + 
#    geom_point(fill=clrs_hunting[1], size=3) + labs(x = "z-transformed fighting rate", y = expression(Delta*" methylation"),
#                      subtitle = m_fight_no_pre_out$sig$chr_pos[i],
#                      title = paste0("Estimate = ", round(m_fight_no_pre_out$sig$parameter_estimate[i], 2),
#                                        ", q-value = ", round(m_fight_no_pre_out$sig$parameter_qval[i], 2))) +
#                                         geom_abline(intercept=m_fight_no_pre_out$sig$intercept[i], slope = m_fight_no_pre_out$sig$parameter_estimate[i], 
#                                          color=clrs_hunting[2], linewidth=1)+
#                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
#    ggsave(plot, file = paste0("plots/model_out/effort/raw/fight_no_pre_", i, ".png"), width=8, height=8)
#}

## with pre
m_fight_pre <- parallel::mclapply(delta_meth_fight_ls, function_model_delta_pheno, parameter="fight", pre="control", mc.cores=4)
m_fight_pre_out <- function_process_model(m_fight_pre, dir_plots = "plots/model_out/effort", dir_data = "results/modeloutput/effort",
                                            name_file = "fight_with_pre", pretty_name = "Centrality", filter_disp=FALSE) 

nrow(m_fight_pre_out$sig) #n=3
#nrow(subset(m_fight_pre_out$sig, chr_pos %in% m_fight_no_pre_out$sig$chr_pos))

system(paste0("rm ", getwd(), "/plots/model_out/effort/raw/fight_with_pre_*"))

for (i in 1:nrow(m_fight_pre_out$sig)){
    ggplot(subset(delta_meth, chr_pos == m_fight_pre_out$sig$chr_pos[i]), aes(x = fight_scl, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = "z-transformed fighting rate", y = expression(Delta*" methylation"),
                      subtitle = m_fight_pre_out$sig$chr_pos[i],
                      title = paste0("Estimate = ", round(m_fight_pre_out$sig$parameter_estimate[i], 2),
                                        ", q-value = ", round(m_fight_pre_out$sig$parameter_qval[i], 2))) +
                                         geom_abline(intercept=m_fight_pre_out$sig$intercept[i], slope = m_fight_pre_out$sig$parameter_estimate[i], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
    ggsave(plot, file = paste0("plots/model_out/effort/raw/fight_with_pre_", i, ".png"), width=8, height=8)
}

# centrality 

### only select cpg sites with enough data
delta_meth_n_dist <- delta_meth %>% group_by(chr_pos) %>% filter(!is.na(delta_meth)& !is.na(dist)) %>% tally()
delta_meth_n_dist <- subset(delta_meth_n_dist, n > 20)
 
delta_meth_sub_dist <- subset(delta_meth, chr_pos %in% delta_meth_n_dist$chr_pos)
length(unique(delta_meth_sub_dist$chr_pos)) # 760 sites

delta_meth_dist_ls <- delta_meth_sub_dist %>% group_split(chr_pos)

## no pre
m_dist_no_pre <- parallel::mclapply(delta_meth_dist_ls, function_model_delta_pheno, parameter="dist", pre="no_control", mc.cores=4)
m_dist_no_pre_out <- function_process_model(m_dist_no_pre, dir_plots = "plots/model_out/effort", dir_data = "results/modeloutput/effort",
                                            name_file = "dist_no_pre", pretty_name = "Centrality", filter_disp=FALSE)

nrow(m_dist_no_pre_out$sig) #n=1 sig
nrow(subset(m_dist_no_pre_out$sig, chr_pos %in% m_dist_pre_out$sig$chr_pos)) #same one 1

system(paste0("rm -i ", getwd(), "/plots/model_out/effort/raw/dist_no_pre_*"))

for (i in 1:nrow(m_dist_no_pre_out$sig)){
    ggplot(subset(delta_meth, chr_pos == m_dist_no_pre_out$sig$chr_pos[i]), aes(x = dist_scl, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = "z-transformed centrality", y = expression(Delta*" methylation"),
                      subtitle = m_dist_no_pre_out$sig$chr_pos[i],
                      title = paste0("Estimate = ", round(m_dist_no_pre_out$sig$parameter_estimate[i], 2),
                                        ", q-value = ", round(m_dist_no_pre_out$sig$parameter_qval[i], 2))) +
                                         geom_abline(intercept=m_dist_no_pre_out$sig$intercept[i], slope = m_dist_no_pre_out$sig$parameter_estimate[i], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
    ggsave(plot, file = paste0("plots/model_out/effort/raw/dist_no_pre_", i, ".png"), width=8, height=8)
}

## with pre
m_dist_pre <- parallel::mclapply(delta_meth_dist_ls, function_model_delta_pheno, parameter="dist", pre="control", mc.cores=4)
m_dist_pre_out <- function_process_model(m_dist_pre, dir_plots = "plots/model_out/effort", dir_data = "results/modeloutput/effort",
                                            name_file = "dist_with_pre", pretty_name = "Centrality", filter_disp=FALSE)

nrow(m_dist_pre_out$sig) #n=1 sig
#n = 1, the one is the same as above
nrow(subset(m_dist_pre_out$sig, chr_pos %in% m_dist_no_pre_out$sig$chr_pos))

system(paste0("rm -i ", getwd(), "/plots/model_out/effort/raw/dist_with_pre_*"))

for (i in 1:nrow(m_dist_pre_out$sig)){
    ggplot(subset(delta_meth, chr_pos == m_dist_pre_out$sig$chr_pos[i]), aes(x = dist_scl, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = "z-transformed centrality", y = expression(Delta*" methylation"),
                      subtitle = m_dist_pre_out$sig$chr_pos[i],
                      title = paste0("Estimate = ", round(m_dist_pre_out$sig$parameter_estimate[i], 2),
                                        ", q-value = ", round(m_dist_pre_out$sig$parameter_qval[i], 2))) +
                                         geom_abline(intercept=m_dist_pre_out$sig$intercept[i], slope = m_dist_pre_out$sig$parameter_estimate[i], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
    ggsave(plot, file = paste0("plots/model_out/effort/raw/dist_with_pre_", i, ".png"), width=8, height=8)
}
