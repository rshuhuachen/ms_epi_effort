### load packages
pacman::p_load(tidyverse, data.table, tibble, performance, matrixStats, 
               parallel, performance, lmerTest, tidystats, insight)

source("scripts/plotting_theme.R")

### load data

load(file="results/modeloutput/changing/changing_sites_glmer.RData")

### load phenotypic data

load(file = "data/phenotypes/fulldata_complete_epi_withdates.RData")

load("data/phenotypes/pheno_dif_prepost.RData") ## differences in physiology

### methylation difference

load(file = "results/modeloutput/all_sites_deltameth.RData")

delta_meth <- subset(delta_meth, chr_pos %in% changing_cpg$chr_pos)

### combine delta methylation data with site and behaviour info

delta_meth <- left_join(delta_meth, unique(all_pheno_epi[,c("id", "year", "site", "Core")], by = c("id", "year")))

### z-transform the traits before the model
prepost_dif$mass_dif_scl <- scale(prepost_dif$mass_dif)
prepost_dif$microf_dif_scl <- scale(prepost_dif$microf_dif)
prepost_dif$trypa_dif_scl <- scale(prepost_dif$trypa_dif)
prepost_dif$ig_dif_scl <- scale(prepost_dif$ig_dif)
prepost_dif$hct_dif_scl <- scale(prepost_dif$hct_dif)

### combine data with pre-post delta physio numbers
delta_meth <- left_join(delta_meth, unique(prepost_dif[,c("id", "year", "mass_dif", "microf_dif", "trypa_dif", "ig_dif", "hct_dif",
                        "mass_dif_scl", "microf_dif_scl", "trypa_dif_scl", "ig_dif_scl", "hct_dif_scl")], by = c("id", "year")))


#### run the model per trait####
source("scripts/function_models.R")
source("scripts/function_models_updated.R")

### mass ####

### only select cpg sites with enough data
delta_meth_n_mass <- delta_meth %>% group_by(chr_pos) %>% filter(!is.na(delta_meth)& !is.na(mass_dif)) %>% tally()
delta_meth_n_mass <- subset(delta_meth_n_mass, n > 20)
 
delta_meth_sub_mass <- subset(delta_meth, chr_pos %in% delta_meth_n_mass$chr_pos)
length(unique(delta_meth_sub_mass$chr_pos)) #808 sites

delta_meth_sub_mass_ls <- delta_meth_sub_mass %>% group_split(chr_pos)

### no pre ####
m_mass_no_pre <- parallel::mclapply(delta_meth_sub_mass_ls, function_model_delta_pheno, parameter="mass_dif", pre="no_control", mc.cores=4)
m_mass_no_pre_out <- function_process_model(m_mass_no_pre, dir_plots = "plots/model_out/physio/mass", dir_data = "results/modeloutput/physio",
                                            name_file = "mass_dif_no_pre", pretty_name = "Delta mass", filter_disp=FALSE) #n=5

nrow(m_mass_no_pre_out$sig) # 5
system(paste0("rm -i ", getwd(), "/plots/model_out/physio/mass/mass_no_pre_*"))

for (i in 1:nrow(m_mass_no_pre_out$sig)){
    ggplot(subset(delta_meth, chr_pos == m_mass_no_pre_out$sig$chr_pos[i]), aes(x = mass_dif_scl, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression("z-transformed "*Delta*" body mass"), y = expression(Delta*" methylation"),
                      subtitle = m_mass_no_pre_out$sig$chr_pos[i],
                      title = paste0("Estimate = ", round(m_mass_no_pre_out$sig$parameter_estimate[i], 2),
                                        ", q-value = ", round(m_mass_no_pre_out$sig$parameter_qval[i], 2))) +
                                         geom_abline(intercept=m_mass_no_pre_out$sig$intercept[i], slope = m_mass_no_pre_out$sig$parameter_estimate[i], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
    ggsave(plot, file = paste0("plots/model_out/physio/mass/mass_no_pre_", i, ".png"), width=8, height=8)
}

### with pre ####
m_mass_with_pre <- parallel::mclapply(delta_meth_sub_mass_ls, function_model_delta_pheno, parameter="mass_dif", pre="control", mc.cores=4)
m_mass_with_pre_out <- function_process_model(m_mass_with_pre, dir_plots = "plots/model_out/physio/mass", dir_data = "results/modeloutput/physio",
                                            name_file = "mass_dif_with_pre", pretty_name = "Delta mass", filter_disp=FALSE) 

nrow(m_mass_with_pre_out$sig) # 2
system(paste0("rm -i ", getwd(), "/plots/model_out/physio/mass/mass_with_pre_*"))

nrow(subset(m_mass_no_pre_out$sig, chr_pos %in% m_mass_with_pre_out$sig$chr_pos)) #1 overlap

for (i in 1:nrow(m_mass_with_pre_out$sig)){
    ggplot(subset(delta_meth, chr_pos == m_mass_with_pre_out$sig$chr_pos[i]), aes(x = mass_dif_scl, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression("z-transformed "*Delta*" body mass"), y = expression(Delta*" methylation"),
                      subtitle = m_mass_with_pre_out$sig$chr_pos[i],
                      title = paste0("Estimate = ", round(m_mass_with_pre_out$sig$parameter_estimate[i], 2),
                                        ", q-value = ", round(m_mass_with_pre_out$sig$parameter_qval[i], 2))) +
                                         geom_abline(intercept=m_mass_with_pre_out$sig$intercept[i], slope = m_mass_with_pre_out$sig$parameter_estimate[i], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
    ggsave(plot, file = paste0("plots/model_out/physio/mass/mass_with_pre_", i, ".png"), width=8, height=8)
}

## with pre but remove repeated samples
delta_meth_mass_ls_norepeat <- delta_meth_sub_mass_ls
for (i in 1:length(delta_meth_mass_ls_norepeat)){
  delta_meth_mass_ls_norepeat[[i]] <- delta_meth_mass_ls_norepeat[[i]] %>%
    group_by(id) %>%
    sample_n(1) %>%
    ungroup()
}

m_mass_pre_norepeat <- parallel::mclapply(delta_meth_mass_ls_norepeat, function_model_delta_pheno_norepeat, parameter="mass_dif", pre="control", mc.cores=4)
m_mass_pre_out_norepeat <- function_process_model(m_mass_pre_norepeat, dir_plots = "plots/model_out/physio", dir_data = "results/modeloutput/physio",
                                                  name_file = "mass_with_pre_norepeat", pretty_name = "Delta mass", filter_disp=FALSE) # n = 3

nrow(m_mass_pre_out_norepeat$sig) #n=0

### microf

### only select cpg sites with enough data
delta_meth_n_microf <- delta_meth %>% group_by(chr_pos) %>% filter(!is.na(delta_meth)& !is.na(microf_dif)) %>% tally()
delta_meth_n_microf <- subset(delta_meth_n_microf, n > 20)
 
delta_meth_sub_microf <- subset(delta_meth, chr_pos %in% delta_meth_n_microf$chr_pos)
length(unique(delta_meth_sub_microf$chr_pos)) #149 sites

delta_meth_sub_microf_ls <- delta_meth_sub_microf %>% group_split(chr_pos)

#### no pre ####

m_microf_no_pre <- parallel::mclapply(delta_meth_sub_microf_ls, function_model_delta_pheno, parameter="microf_dif", pre="no_control", mc.cores=4)
m_microf_no_pre_out <- function_process_model(m_microf_no_pre, dir_plots = "plots/model_out/physio/microf", dir_data = "results/modeloutput/physio",
                                            name_file = "microf_dif_no_pre", pretty_name = "Delta Microfilaria spp.", filter_disp=FALSE) #n=1
nrow(m_microf_no_pre_out$sig) #1

system(paste0("rm -i ", getwd(), "/plots/model_out/physio/microf/microf_no_pre_*"))

for (i in 1:nrow(m_microf_no_pre_out$sig)){
    ggplot(subset(delta_meth, chr_pos == m_microf_no_pre_out$sig$chr_pos[i]), aes(x = microf_dif_scl, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression("z-transformed "*Delta*" Microfilaria spp."), y = expression(Delta*" methylation"),
                      subtitle = m_microf_no_pre_out$sig$chr_pos[i],
                      title = paste0("Estimate = ", round(m_microf_no_pre_out$sig$parameter_estimate[i], 2),
                                        ", q-value = ", round(m_microf_no_pre_out$sig$parameter_qval[i], 2))) +
                                         geom_abline(intercept=m_microf_no_pre_out$sig$intercept[i], slope = m_microf_no_pre_out$sig$parameter_estimate[i], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
    ggsave(plot, file = paste0("plots/model_out/physio/microf/microf_no_pre_", i, ".png"), width=8, height=8)
}

#### with pre ####
m_microf_with_pre <- parallel::mclapply(delta_meth_sub_microf_ls, function_model_delta_pheno, parameter="microf_dif", pre="control", mc.cores=4)
m_microf_with_pre_out <- function_process_model(m_microf_with_pre, dir_plots = "plots/model_out/physio/microf", dir_data = "results/modeloutput/physio",
                                            name_file = "microf_dif_with_pre", pretty_name = "Delta Microfilaria spp.", filter_disp=FALSE) #

nrow(m_microf_with_pre_out$sig) #3
nrow(subset(m_microf_with_pre_out$sig, chr_pos %in% m_microf_no_pre_out$sig$chr_pos)) #0 overlap

system(paste0("rm -i ", getwd(), "/plots/model_out/physio/microf/microf_with_pre_*"))

for (i in 1:nrow(m_microf_with_pre_out$sig)){
    ggplot(subset(delta_meth, chr_pos == m_microf_with_pre_out$sig$chr_pos[i]), aes(x = microf_dif_scl, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression("z-transformed "*Delta*" Microfilaria spp."), y = expression(Delta*" methylation"),
                      subtitle = m_microf_with_pre_out$sig$chr_pos[i],
                      title = paste0("Estimate = ", round(m_microf_with_pre_out$sig$parameter_estimate[i], 2),
                                        ", q-value = ", round(m_microf_with_pre_out$sig$parameter_qval[i], 2))) +
                                         geom_abline(intercept=m_microf_with_pre_out$sig$intercept[i], slope = m_microf_with_pre_out$sig$parameter_estimate[i], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
    ggsave(plot, file = paste0("plots/model_out/physio/microf/microf_with_pre_", i, ".png"), width=8, height=8)
}

## with pre but remove repeated samples
delta_meth_microf_ls_norepeat <- delta_meth_sub_microf_ls
for (i in 1:length(delta_meth_microf_ls_norepeat)){
  delta_meth_microf_ls_norepeat[[i]] <- delta_meth_microf_ls_norepeat[[i]] %>%
    group_by(id) %>%
    sample_n(1) %>%
    ungroup()
}

m_microf_pre_norepeat <- parallel::mclapply(delta_meth_microf_ls_norepeat, function_model_delta_pheno_norepeat, parameter="microf_dif", pre="control", mc.cores=4)
m_microf_pre_out_norepeat <- function_process_model(m_microf_pre_norepeat, dir_plots = "plots/model_out/physio", dir_data = "results/modeloutput/physio",
                                                  name_file = "microf_with_pre_norepeat", pretty_name = "Delta Microfilaria", filter_disp=FALSE) # n = 3

nrow(m_microf_pre_out_norepeat$sig) #n=0

### trypa ####

### only select cpg sites with enough data
delta_meth_n_trypa <- delta_meth %>% group_by(chr_pos) %>% filter(!is.na(delta_meth)& !is.na(trypa_dif)) %>% tally()
delta_meth_n_trypa <- subset(delta_meth_n_trypa, n > 20)
 
delta_meth_sub_trypa <- subset(delta_meth, chr_pos %in% delta_meth_n_trypa$chr_pos)
length(unique(delta_meth_sub_trypa$chr_pos)) #164 sites

delta_meth_sub_trypa_ls <- delta_meth_sub_trypa %>% group_split(chr_pos)

#### no pre ####

m_trypa_no_pre <- parallel::mclapply(delta_meth_sub_trypa_ls, function_model_delta_pheno, parameter="trypa_dif", pre="no_control", mc.cores=4)
m_trypa_no_pre_out <- function_process_model(m_trypa_no_pre, dir_plots = "plots/model_out/physio/trypa", dir_data = "results/modeloutput/physio",
                                            name_file = "trypa_dif_no_pre", pretty_name = "Delta Trypanosoma spp.", filter_disp=FALSE) #n=0

nrow(m_trypa_no_pre_out$sig) # 0

#### with pre ####
m_trypa_with_pre <- parallel::mclapply(delta_meth_sub_trypa_ls, function_model_delta_pheno, parameter="trypa_dif", pre="control", mc.cores=4)
m_trypa_with_pre_out <- function_process_model(m_trypa_with_pre, dir_plots = "plots/model_out/physio/trypa", dir_data = "results/modeloutput/physio",
                                            name_file = "trypa_dif_with_pre", pretty_name = "Delta Trypanosoma spp.", filter_disp=FALSE) #n=0


nrow(m_trypa_with_pre_out$sig) # 0


## with pre but remove repeated samples
delta_meth_trypa_ls_norepeat <- delta_meth_sub_trypa_ls
for (i in 1:length(delta_meth_trypa_ls_norepeat)){
  delta_meth_trypa_ls_norepeat[[i]] <- delta_meth_trypa_ls_norepeat[[i]] %>%
    group_by(id) %>%
    sample_n(1) %>%
    ungroup()
}

m_trypa_pre_norepeat <- parallel::mclapply(delta_meth_trypa_ls_norepeat, function_model_delta_pheno_norepeat, parameter="trypa_dif", pre="control", mc.cores=4)
m_trypa_pre_out_norepeat <- function_process_model(m_trypa_pre_norepeat, dir_plots = "plots/model_out/physio", dir_data = "results/modeloutput/physio",
                                                    name_file = "trypa_with_pre_norepeat", pretty_name = "Delta Trypanosoma", filter_disp=FALSE) 

nrow(m_trypa_pre_out_norepeat$sig) #n=0

### ig ####
### only select cpg sites with enough data
delta_meth_n_ig <- delta_meth %>% group_by(chr_pos) %>% filter(!is.na(delta_meth)& !is.na(ig_dif)) %>% tally()
delta_meth_n_ig <- subset(delta_meth_n_ig, n > 20)
 
delta_meth_sub_ig <- subset(delta_meth, chr_pos %in% delta_meth_n_ig$chr_pos)
length(unique(delta_meth_sub_ig$chr_pos)) #775 sites

delta_meth_sub_ig_ls <- delta_meth_sub_ig %>% group_split(chr_pos)

#### no pre ####

m_ig_no_pre <- parallel::mclapply(delta_meth_sub_ig_ls, function_model_delta_pheno, parameter="ig_dif", pre="no_control", mc.cores=4)
m_ig_no_pre_out <- function_process_model(m_ig_no_pre, dir_plots = "plots/model_out/physio/ig", dir_data = "results/modeloutput/physio",
                                            name_file = "ig_dif_no_pre", pretty_name = "Delta IgG", filter_disp=FALSE) 

nrow(m_ig_no_pre_out$sig) #n = 2

system(paste0("rm -i ", getwd(), "/plots/model_out/physio/ig/ig_no_pre_*"))

for (i in 1:nrow(m_ig_no_pre_out$sig)){
    ggplot(subset(delta_meth, chr_pos == m_ig_no_pre_out$sig$chr_pos[i]), aes(x = ig_dif_scl, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression("z-transformed "*Delta*" IgG"), y = expression(Delta*" methylation"),
                      subtitle = m_ig_no_pre_out$sig$chr_pos[i],
                      title = paste0("Estimate = ", round(m_ig_no_pre_out$sig$parameter_estimate[i], 2),
                                        ", q-value = ", round(m_ig_no_pre_out$sig$parameter_qval[i], 2))) +
                                         geom_abline(intercept=m_ig_no_pre_out$sig$intercept[i], slope = m_ig_no_pre_out$sig$parameter_estimate[i], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
    ggsave(plot, file = paste0("plots/model_out/physio/ig/ig_no_pre_", i, ".png"), width=8, height=8)
}

#### with pre ####
m_ig_with_pre <- parallel::mclapply(delta_meth_sub_ig_ls, function_model_delta_pheno, parameter="ig_dif", pre="control", mc.cores=4)
m_ig_with_pre_out <- function_process_model(m_ig_with_pre, dir_plots = "plots/model_out/physio/ig", dir_data = "results/modeloutput/physio",
                                            name_file = "ig_dif_with_pre", pretty_name = "Delta IgG", filter_disp=FALSE) 

nrow(m_ig_with_pre_out$sig) #n=2

nrow(subset(m_ig_with_pre_out$sig, chr_pos %in% m_ig_no_pre_out$sig$chr_pos)) #2 overlap

system(paste0("rm -i ", getwd(), "/plots/model_out/physio/ig/ig_with_pre_*"))

for (i in 1:nrow(m_ig_with_pre_out$sig)){
    ggplot(subset(delta_meth, chr_pos == m_ig_with_pre_out$sig$chr_pos[i]), aes(x = ig_dif_scl, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression("z-transformed "*Delta*" IgG"), y = expression(Delta*" methylation"),
                      subtitle = m_ig_with_pre_out$sig$chr_pos[i],
                      title = paste0("Estimate = ", round(m_ig_with_pre_out$sig$parameter_estimate[i], 2),
                                        ", q-value = ", round(m_ig_with_pre_out$sig$parameter_qval[i], 2))) +
                                         geom_abline(intercept=m_ig_with_pre_out$sig$intercept[i], slope = m_ig_with_pre_out$sig$parameter_estimate[i], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
    ggsave(plot, file = paste0("plots/model_out/physio/ig/ig_with_pre_", i, ".png"), width=8, height=8)
}

## with pre but remove repeated samples
delta_meth_igg_ls_norepeat <- delta_meth_sub_ig_ls
for (i in 1:length(delta_meth_igg_ls_norepeat)){
  delta_meth_igg_ls_norepeat[[i]] <- delta_meth_igg_ls_norepeat[[i]] %>%
    group_by(id) %>%
    sample_n(1) %>%
    ungroup()
}

m_igg_pre_norepeat <- parallel::mclapply(delta_meth_igg_ls_norepeat, function_model_delta_pheno_norepeat, parameter="ig_dif", pre="control", mc.cores=4)
m_igg_pre_out_norepeat <- function_process_model(m_igg_pre_norepeat, dir_plots = "plots/model_out/physio", dir_data = "results/modeloutput/physio",
                                                   name_file = "igg_with_pre_norepeat", pretty_name = "Delta IgG", filter_disp=FALSE) 

nrow(m_igg_pre_out_norepeat$sig) #n=0

### hct
delta_meth_n_hct <- delta_meth %>% group_by(chr_pos) %>% filter(!is.na(delta_meth)& !is.na(hct_dif)) %>% tally()
delta_meth_n_hct <- subset(delta_meth_n_hct, n > 20)
 
delta_meth_sub_hct <- subset(delta_meth, chr_pos %in% delta_meth_n_hct$chr_pos)
length(unique(delta_meth_sub_hct$chr_pos)) #146 sites

delta_meth_sub_hct_ls <- delta_meth_sub_hct %>% group_split(chr_pos)

#### no pre ####

m_hct_no_pre <- parallel::mclapply(delta_meth_sub_hct_ls, function_model_delta_pheno, parameter="hct_dif", pre="no_control", mc.cores=4)
m_hct_no_pre_out <- function_process_model(m_hct_no_pre, dir_plots = "plots/model_out/physio/hct", dir_data = "results/modeloutput/physio",
                                            name_file = "hct_dif_no_pre", pretty_name = "Delta HCT", filter_disp=FALSE) 
nrow(m_hct_no_pre_out$sig) # 0

#### with pre ####
m_hct_with_pre <- parallel::mclapply(delta_meth_sub_hct_ls, function_model_delta_pheno, parameter="hct_dif", pre="control", mc.cores=4)
m_hct_with_pre_out <- function_process_model(m_hct_with_pre, dir_plots = "plots/model_out/physio/hct", dir_data = "results/modeloutput/physio",
                                            name_file = "hct_dif_with_pre", pretty_name = "Delta HCT", filter_disp=FALSE)
nrow(m_hct_with_pre_out$sig) # 1

system(paste0("rm -i ", getwd(), "/plots/model_out/physio/hct/hct_with_pre_*"))

for (i in 1:nrow(m_hct_with_pre_out$sig)){
    ggplot(subset(delta_meth, chr_pos == m_hct_with_pre_out$sig$chr_pos[i]), aes(x = hct_dif_scl, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression("z-transformed "*Delta*" HCT"), y = expression(Delta*" methylation"),
                      subtitle = m_hct_with_pre_out$sig$chr_pos[i],
                      title = paste0("Estimate = ", round(m_hct_with_pre_out$sig$parameter_estimate[i], 2),
                                        ", q-value = ", round(m_hct_with_pre_out$sig$parameter_qval[i], 2))) +
                                         geom_abline(intercept=m_hct_with_pre_out$sig$intercept[i], slope = m_hct_with_pre_out$sig$parameter_estimate[i], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
    ggsave(plot, file = paste0("plots/model_out/physio/hct/hct_with_pre_", i, ".png"), width=8, height=8)
}

## with pre but remove repeated samples
delta_meth_hct_ls_norepeat <- delta_meth_sub_hct_ls
for (i in 1:length(delta_meth_hct_ls_norepeat)){
  delta_meth_hct_ls_norepeat[[i]] <- delta_meth_hct_ls_norepeat[[i]] %>%
    group_by(id) %>%
    sample_n(1) %>%
    ungroup()
}

m_hct_pre_norepeat <- parallel::mclapply(delta_meth_hct_ls_norepeat, function_model_delta_pheno_norepeat, parameter="hct_dif", pre="control", mc.cores=4)
m_hct_pre_out_norepeat <- function_process_model(m_hct_pre_norepeat, dir_plots = "plots/model_out/physio", dir_data = "results/modeloutput/physio",
                                                 name_file = "hct_with_pre_norepeat", pretty_name = "Delta HCT", filter_disp=FALSE) 

nrow(m_hct_pre_out_norepeat$sig) #n=0
