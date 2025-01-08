### load packages
pacman::p_load(tidyverse, data.table, tibble, performance, gaston,
               parallel, performance, lmerTest, tidystats, ggpointdensity)

source("scripts/plotting_theme.R")

### load data

load(file="results/modeloutput/changing/changing_sites_glmer.RData")

### load phenotypic data

load(file = "data/phenotypes/prepost_with_injury.RData")

### methylation difference

load(file = "results/modeloutput/all_sites_deltameth.RData")

delta_meth <- subset(delta_meth, chr_pos %in% changing_cpg$chr_pos) #38618

### combine delta methylation data with site and behaviour info
injury <- subset(prepost, prepost == "post")
injury[which(injury$id == "D229096" & injury$injury_belly == 0),] <- NA #manual correction based on email
injury <- subset(injury, !is.na(epi_nr))
injury <- injury %>% select(c(id, year, site, PC1, PC2))
injury$PC1_scl <- injury$PC1 # don't scale, just to make the function run
injury$PC2_scl <- injury$PC2 # don't scale, just to make the function run
delta_meth <- left_join(delta_meth, injury, by = c("id", "year"))

#### run the model per trait: with pre-lekking ####
source("scripts/function_models_updated.R")

# PC1 

### only select cpg sites with enough data
delta_meth_n_pc1 <- delta_meth %>% group_by(chr_pos) %>% filter(!is.na(delta_meth)& !is.na(PC1)) %>% tally()
delta_meth_n_pc1 <- subset(delta_meth_n_pc1, n > 20)

delta_meth_sub_pc1 <- subset(delta_meth, chr_pos %in% delta_meth_n_pc1$chr_pos)
length(unique(delta_meth_sub_pc1$chr_pos)) # 730 sites

delta_meth_pc1_ls <- delta_meth_sub_pc1 %>% group_split(chr_pos)

# exclude repeated samples
for (i in 1:length(delta_meth_pc1_ls)){
  delta_meth_pc1_ls[[i]] <- delta_meth_pc1_ls[[i]] %>%
    group_by(id) %>%
    sample_n(1) %>%
    ungroup()
}


## with pre
m_pc1_pre <- parallel::mclapply(delta_meth_pc1_ls, function_model_delta_pheno, parameter="PC1", pre="control", mc.cores=4)
m_pc1_pre_out <- function_process_model(m_pc1_pre, dir_plots = "plots/model_out/injury", 
                                        dir_data = "results/modeloutput/injury",
                                           name_file = "pc1_with_pre", pretty_name = "pc1", filter_disp=FALSE)

nrow(m_pc1_pre_out$sig) #n=1

# PC2

### only select cpg sites with enough data
delta_meth_n_pc2 <- delta_meth %>% group_by(chr_pos) %>% filter(!is.na(delta_meth)& !is.na(PC2)) %>% tally()
delta_meth_n_pc2 <- subset(delta_meth_n_pc2, n > 20)

delta_meth_sub_pc2 <- subset(delta_meth, chr_pos %in% delta_meth_n_pc2$chr_pos)
length(unique(delta_meth_sub_pc2$chr_pos)) # 730 sites

delta_meth_pc2_ls <- delta_meth_sub_pc2 %>% group_split(chr_pos)

# exclude repeated samples
for (i in 1:length(delta_meth_pc2_ls)){
  delta_meth_pc2_ls[[i]] <- delta_meth_pc2_ls[[i]] %>%
    group_by(id) %>%
    sample_n(1) %>%
    ungroup()
}

## with pre
m_pc2_pre <- parallel::mclapply(delta_meth_pc2_ls, function_model_delta_pheno, parameter="PC2", pre="control", mc.cores=4)
m_pc2_pre_out <- function_process_model(m_pc2_pre, dir_plots = "plots/model_out/injury", 
                                        dir_data = "results/modeloutput/injury",
                                        name_file = "pc2_with_pre", pretty_name = "pc2", filter_disp=FALSE) 

nrow(m_pc2_pre_out$sig) #n=0

