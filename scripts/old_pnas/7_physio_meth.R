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

##### mass ####
## remove repeated samples
delta_meth_mass <- delta_meth %>%
  filter(!is.na(delta_meth)& !is.na(mass_dif))%>%
  group_by(chr_pos, id) %>%
  sample_n(1) %>%
  ungroup()

## only select cpg sites with enough data
delta_meth_n_mass <- delta_meth_mass %>% group_by(chr_pos) %>% tally()
delta_meth_n_mass <- subset(delta_meth_n_mass, n > 20)
 
delta_meth_sub_mass <- subset(delta_meth_mass, chr_pos %in% delta_meth_n_mass$chr_pos)
length(unique(delta_meth_sub_mass$chr_pos)) #607 sites

delta_meth_sub_mass_ls <- delta_meth_sub_mass %>% group_split(chr_pos)
save(delta_meth_sub_mass_ls, file = "data/processed/delta_meth_ls_mass.RData")

## model
m_mass_with_pre <- parallel::mclapply(delta_meth_sub_mass_ls, function_model_delta_pheno_norepeat, parameter="mass_dif", pre="control", mc.cores=4)
m_mass_with_pre_out <- function_process_model(m_mass_with_pre, dir_plots = "plots/model_out/physio/mass", dir_data = "results/modeloutput/physio",
                                            name_file = "mass_dif_with_pre", pretty_name = "Delta mass", filter_disp=FALSE) 

nrow(m_mass_with_pre_out$data) # 607
nrow(m_mass_with_pre_out$sig) # 0

##### microf ####
## remove repeated samples
delta_meth_microf <- delta_meth %>%
  filter(!is.na(delta_meth)& !is.na(microf_dif))%>%
  group_by(chr_pos, id) %>%
  sample_n(1) %>%
  ungroup()

### only select cpg sites with enough data
delta_meth_n_microf <- delta_meth_microf %>% group_by(chr_pos) %>% tally()
delta_meth_n_microf <- subset(delta_meth_n_microf, n > 20)
 
delta_meth_sub_microf <- subset(delta_meth_microf, chr_pos %in% delta_meth_n_microf$chr_pos)
length(unique(delta_meth_sub_microf$chr_pos)) #81 sites

delta_meth_sub_microf_ls <- delta_meth_sub_microf %>% group_split(chr_pos)
save(delta_meth_sub_microf_ls, file = "data/processed/delta_meth_ls_microf.RData")

## model 

m_microf_with_pre <- parallel::mclapply(delta_meth_sub_microf_ls, function_model_delta_pheno_norepeat, parameter="microf_dif", pre="control", mc.cores=4)
m_microf_with_pre_out <- function_process_model(m_microf_no_pre, dir_plots = "plots/model_out/physio/microf", dir_data = "results/modeloutput/physio",
                                            name_file = "microf_dif_with_pre", pretty_name = "Delta Microfilaria spp.", filter_disp=FALSE) #n=1
nrow(m_microf_with_pre_out$data) #81
nrow(m_microf_with_pre_out$sig) #0

##### trypa ####
## remove repeated samples
delta_meth_trypa <- delta_meth %>%
  filter(!is.na(delta_meth)& !is.na(trypa_dif))%>%
  group_by(chr_pos, id) %>%
  sample_n(1) %>%
  ungroup()

### only select cpg sites with enough data
delta_meth_n_trypa <- delta_meth_trypa %>% group_by(chr_pos) %>% tally()
delta_meth_n_trypa <- subset(delta_meth_n_trypa, n > 20)
 
delta_meth_sub_trypa <- subset(delta_meth_trypa, chr_pos %in% delta_meth_n_trypa$chr_pos)
length(unique(delta_meth_sub_trypa$chr_pos)) #108 sites

delta_meth_sub_trypa_ls <- delta_meth_sub_trypa %>% group_split(chr_pos)
save(delta_meth_sub_trypa_ls, file = "data/processed/delta_meth_ls_trypa.RData")

## model
m_trypa_with_pre <- parallel::mclapply(delta_meth_sub_trypa_ls, function_model_delta_pheno_norepeat, parameter="trypa_dif", pre="control", mc.cores=4)
m_trypa_with_pre_out <- function_process_model(m_trypa_with_pre, dir_plots = "plots/model_out/physio/trypa", dir_data = "results/modeloutput/physio",
                                            name_file = "trypa_dif_with_pre", pretty_name = "Delta Trypanosoma spp.", filter_disp=FALSE) #n=0


nrow(m_trypa_with_pre_out$data) # 108
nrow(m_trypa_with_pre_out$sig) # 0

##### ig ####
## remove repeated samples
delta_meth_ig <- delta_meth %>%
  filter(!is.na(delta_meth)& !is.na(ig_dif))%>%
  group_by(chr_pos, id) %>%
  sample_n(1) %>%
  ungroup()

## only select cpg sites with enough data
delta_meth_n_ig <- delta_meth_ig %>% group_by(chr_pos) %>% tally()
delta_meth_n_ig <- subset(delta_meth_n_ig, n > 20)
 
delta_meth_sub_ig <- subset(delta_meth_ig, chr_pos %in% delta_meth_n_ig$chr_pos)
length(unique(delta_meth_sub_ig$chr_pos)) #576 sites

delta_meth_sub_ig_ls <- delta_meth_sub_ig %>% group_split(chr_pos)
save(delta_meth_sub_ig_ls, file = "data/processed/delta_meth_ls_ig.RData")

## model
m_ig_with_pre <- parallel::mclapply(delta_meth_sub_ig_ls, function_model_delta_pheno_norepeat, parameter="ig_dif", pre="control", mc.cores=4)
m_ig_with_pre_out <- function_process_model(m_ig_with_pre, dir_plots = "plots/model_out/physio/ig", dir_data = "results/modeloutput/physio",
                                            name_file = "ig_dif_with_pre", pretty_name = "Delta IgG", filter_disp=FALSE) 

nrow(m_ig_with_pre_out$data) #n=576
nrow(m_ig_with_pre_out$sig) #n=0

### hct
## remove repeated samples
delta_meth_hct <- delta_meth %>%
  filter(!is.na(delta_meth)& !is.na(hct_dif))%>%
  group_by(chr_pos, id) %>%
  sample_n(1) %>%
  ungroup()

## only select cpg sites with enough data
delta_meth_n_hct <- delta_meth_hct %>% group_by(chr_pos) %>% tally()
delta_meth_n_hct <- subset(delta_meth_n_hct, n > 20)
 
delta_meth_sub_hct <- subset(delta_meth_hct, chr_pos %in% delta_meth_n_hct$chr_pos)
length(unique(delta_meth_sub_hct$chr_pos)) #79 sites

delta_meth_sub_hct_ls <- delta_meth_sub_hct %>% group_split(chr_pos)
save(delta_meth_sub_hct_ls, file = "data/processed/delta_meth_ls_hct.RData")

## model 
m_hct_with_pre <- parallel::mclapply(delta_meth_sub_hct_ls, function_model_delta_pheno_norepeat, parameter="hct_dif", pre="control", mc.cores=4)
m_hct_with_pre_out <- function_process_model(m_hct_with_pre, dir_plots = "plots/model_out/physio/hct", dir_data = "results/modeloutput/physio",
                                            name_file = "hct_dif_with_pre", pretty_name = "Delta HCT", filter_disp=FALSE)
nrow(m_hct_with_pre_out$data) # 79
nrow(m_hct_with_pre_out$sig) # 0
