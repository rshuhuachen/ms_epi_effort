### load packages
pacman::p_load(dplyr, data.table, tibble, performance, gaston,
               parallel, performance, lmerTest, tidystats, ggpointdensity)

source("scripts/plotting_theme.R")

### load data

load(file="results/modeloutput/changing/changing_sites_glmer.RData")

### load phenotypic data

load(file = "data/phenotypes/fulldata_complete_epi_withdates.RData")

### methylation difference

load(file = "results/modeloutput/all_sites_deltameth.RData")

delta_meth <- subset(delta_meth, chr_pos %in% changing_cpg$chr_pos)

### Plot pre_meth effect on delta_meth ####

lmer_pre_delta <- lmerTest::lmer(delta_meth ~ methperc_pre + (1|chr_pos)+ (1|id), data = delta_meth)

sum_pre_delta <- as.data.frame(summary(lmer_pre_delta)$coef) # negative relationship

ggplot(delta_meth, aes(methperc_pre, delta_meth)) + 
  geom_pointdensity() + 
  scale_color_viridis_c() + 
  geom_abline(intercept=sum_pre_delta$Estimate[1], slope = sum_pre_delta$Estimate[2], 
              color="red", linewidth=1)+
  labs(x = "Methylation % pre-lekking", y = expression(Delta*" methylation %")) -> cor_pre_delta

ggplot(delta_meth, aes(methperc_pre, abs(delta_meth))) + 
  geom_pointdensity() + 
  scale_color_viridis_c() + #geom_smooth() + 
  labs(x = "Methylation % pre-lekking", y = expression("Absolute "*Delta*" methylation %")) -> cor_pre_delta_abs

cowplot::plot_grid(cor_pre_delta, cor_pre_delta_abs, labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) %>%
  ggsave(file = "plots/explore/pre_vs_delta.png", width=14, height=10)

### combine delta methylation data with site and behaviour info

delta_meth <- left_join(delta_meth, unique(all_pheno_epi[,c("id", "year", "site", "Core")], by = c("id", "year")))

### combine reproductive effort data with methylation data

effort <- all_pheno_epi %>% dplyr::select(c("id", "year", "attend", "dist", "MS")) %>% unique()

effort <- subset(effort, id %in% delta_meth$id)

delta_meth <- left_join(delta_meth, unique(effort[,c("id", "year", "attend", "dist", "MS")]), by = c("id", "year"))

### Run the model per trait ####
source("scripts/function_models.R")

#### attendance ####
## remove repeated samples
delta_meth_attend <- delta_meth %>%
  filter(!is.na(delta_meth)& !is.na(attend))%>%
  group_by(chr_pos, id) %>%
  sample_n(1) %>%
  ungroup()

## only select cpg sites with enough data
delta_meth_n_attend <- delta_meth_attend %>% group_by(chr_pos) %>% tally()
delta_meth_n_attend <- subset(delta_meth_n_attend, n > 20)

delta_meth_sub_attend <- subset(delta_meth_attend, chr_pos %in% delta_meth_n_attend$chr_pos)
length(unique(delta_meth_sub_attend$chr_pos)) # 607 sites

delta_meth_attend_ls <- delta_meth_sub_attend %>% group_split(chr_pos)

save(delta_meth_sub_attend, file = "data/processed/delta_meth_ls_attend.RData")

## run the model
m_attend_pre <- parallel::mclapply(delta_meth_attend_ls, function_model_delta_pheno_norepeat, parameter="attend", pre="control", mc.cores=4)
m_attend_pre_out <- function_process_model(m_attend_pre, dir_plots = "plots/model_out/effort", dir_data = "results/modeloutput/effort",
                                           name_file = "attend_with_pre", pretty_name = "Attendance", filter_disp=FALSE) 
nrow(m_attend_pre_out$data) #n=602
nrow(m_attend_pre_out$sig) #n=1

almostsig_attend <- subset(m_attend_pre_out$data, parameter_pval < 0.05)
summary(almostsig_attend$parameter_estimate)

#### centrality ####
## remove repeated samples
delta_meth_dist <- delta_meth %>%
  filter(!is.na(delta_meth)& !is.na(dist))%>%
  group_by(chr_pos, id) %>%
  sample_n(1) %>%
  ungroup()

## only select cpg sites with enough data
delta_meth_n_dist <- delta_meth_dist %>% group_by(chr_pos) %>% tally()
delta_meth_n_dist <- subset(delta_meth_n_dist, n > 20)

delta_meth_sub_dist <- subset(delta_meth_dist, chr_pos %in% delta_meth_n_dist$chr_pos)
length(unique(delta_meth_sub_dist$chr_pos)) # 534 sites

delta_meth_dist_ls <- delta_meth_sub_dist %>% group_split(chr_pos)
save(delta_meth_dist_ls, file = "data/processed/delta_meth_ls_dist.RData")

## model
m_dist_pre <- parallel::mclapply(delta_meth_dist_ls, function_model_delta_pheno_norepeat, parameter="dist", pre="control", mc.cores=4)
m_dist_pre_out <- function_process_model(m_dist_pre, dir_plots = "plots/model_out/effort", dir_data = "results/modeloutput/effort",
                                         name_file = "dist_with_pre", pretty_name = "Centrality", filter_disp=FALSE)
nrow(m_dist_pre_out$data) #n=534
nrow(m_dist_pre_out$sig) #n=0

almostsig_dist <- subset(m_dist_pre_out$data, parameter_pval < 0.05)
summary(almostsig_dist$parameter_estimate)

#### mating success ####
## remove repeated samples
delta_meth_MS <- delta_meth %>%
  filter(!is.na(delta_meth)& !is.na(MS))%>%
  group_by(chr_pos, id) %>%
  sample_n(1) %>%
  ungroup()

## only select cpg sites with enough data
delta_meth_n_MS <- delta_meth_MS %>% group_by(chr_pos) %>% tally()
delta_meth_n_MS <- subset(delta_meth_n_MS, n > 20)

delta_meth_sub_MS <- subset(delta_meth_MS, chr_pos %in% delta_meth_n_MS$chr_pos)
length(unique(delta_meth_sub_MS$chr_pos)) # 564 sites
save(delta_meth_sub_MS, file = "data/processed/delta_meth_MS_sub.RData")

delta_meth_MS_ls <-  delta_meth_sub_MS %>% group_split(chr_pos)


## model
m_MS_pre <- parallel::mclapply(delta_meth_MS_ls, function_model_delta_pheno_norepeat, parameter="MS", pre="control", mc.cores=4)
m_MS_pre_out <- function_process_model(m_MS_pre, dir_plots = "plots/model_out/effort", dir_data = "results/modeloutput/effort",
                                         name_file = "MS_with_pre", pretty_name = "Mating success", filter_disp=FALSE)
nrow(m_MS_pre_out$data) #n=563
nrow(m_MS_pre_out$sig) #n=0, was 1 before! due to random sampling of repeats?

almostsig_ms <- subset(m_MS_pre_out$data, parameter_pval < 0.05)
summary(almostsig_ms$parameter_estimate)

### compare almost sig cpg sites ####

almostsig_all <- rbind(almostsig_attend, almostsig_dist, almostsig_ms)
dup <- almostsig_all[duplicated(almostsig_all$chr_pos),]
dup <- left_join(data.frame(chr_pos=dup[,c(1)]), almostsig_all[,c(1,2)],by="chr_pos")

save(dup, file = "results/modeloutput/effort/multiple_cpgs_effort_nofdr.RData")
