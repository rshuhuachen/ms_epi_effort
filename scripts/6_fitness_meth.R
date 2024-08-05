### load packages
pacman::p_load(tidyverse, data.table, tibble, performance, matrixStats, 
               parallel, performance, lmerTest, tidystats, insight, glmmTMB)

### load data

load(file = "results/modeloutput/subset_sites_sig_deltameth.RData")

### load phenotypic data

load("data/phenotypes/fulldata_complete_epi_withdates.RData")
load("data/phenotypes/pheno_dif_prepost.RData")

#combine with site and fitness data
pheno_pre <- subset(all_pheno_epi, prepost=="pre")

delta_meth <- left_join(delta_meth, unique(pheno_pre[,c("id", "year", "MS", "surv")]), by = c("id", "year"))
                                           
delta_meth_ls <- delta_meth %>% group_split(chr_pos)

# function to run the model
function_model_fitness <- function(df){tryCatch({
  chr_pos <- as.character(df[1,1])
  df <- as.data.frame(df)
  df$site <- as.factor(df$site)
  df$id <- as.factor(df$id)
  
  ### AMS
  formula_ams <- formula(paste0("MS ~ scale(delta_meth) + age + scale(methperc_pre) + (1|site/id) "))
  model_ams <- glmmTMB(formula_ams, data=df, family = "poisson", ziformula=~1)
  summary_ams <- summary(model_ams)
  
  #fixed effect
  parameter_estimate <- summary_ams$coefficients$cond["scale(delta_meth)", "Estimate"]
  parameter_se <- summary_ams$coefficients$cond["scale(delta_meth)","Std. Error"]
  parameter_zval <- summary_ams$coefficients$cond["scale(delta_meth)","z value"]
  parameter_pval <- summary_ams$coefficients$cond["scale(delta_meth)", "Pr(>|z|)"]
  
  #age effect
  age_estimate <- summary_ams$coefficients$cond["age", "Estimate"]
  age_se <- summary_ams$coefficients$cond["age", "Std. Error"]
  age_zval <- summary_ams$coefficients$cond["age", "z value"]
  age_pval <- summary_ams$coefficients$cond["age", "Pr(>|z|)"]
  
  #premeth effect
  pre_estimate <- summary_ams$coefficients$cond["scale(methperc_pre)", "Estimate"]
  pre_se <- summary_ams$coefficients$cond["scale(methperc_pre)", "Std. Error"]
  pre_zval <- summary_ams$coefficients$cond["scale(methperc_pre)", "z value"]
  pre_pval <- summary_ams$coefficients$cond["scale(methperc_pre)", "Pr(>|z|)"]
  
  rsqc <- performance::r2(model_ams)$R2_conditional #fixed plus random parameterects relative to overall variance
  rsqm <- performance::r2(model_ams)$R2_marginal #fixed parameterects relative to overall variance
  
  message <- model_ams$fit$message
  
  icc_id_site <-icc(model_ams, by_group = TRUE, tolerance = 0)[1,2]
  icc_site <-icc(model_ams, by_group = TRUE, tolerance = 0)[2,2]
  
 
  ams <- data.frame(chr_pos=as.factor(chr_pos),
                    ams_icc_id_site = as.numeric(icc_id_site),
                    ams_icc_site = as.numeric(icc_site),
                    ams_delta_meth_estimate = as.numeric(parameter_estimate),
                    ams_delta_meth_se = as.numeric(parameter_se),
                    ams_delta_meth_zval = as.numeric(parameter_zval),
                    ams_delta_meth_pval = as.numeric(parameter_pval),
                    ams_age_estimate = as.numeric(age_estimate),
                    ams_age_se = as.numeric(age_se),
                    ams_age_zval = as.numeric(age_zval),
                    ams_age_pval = as.numeric(age_pval),
                    ams_pre_estimate = as.numeric(pre_estimate),
                    ams_pre_se = as.numeric(pre_se),
                    ams_pre_zval = as.numeric(pre_zval),
                    ams_pre_pval = as.numeric(pre_pval),
                    ams_rsqc = as.numeric(rsqc),
                    ams_rsqm = as.numeric(rsqm),
                    ams_message = message              
  ) 
  
  ### surv
  formula_surv <- formula(paste0("surv ~ scale(delta_meth) + age + scale(methperc_pre) + (1|site/id) "))
  model_surv <- glmmTMB(formula_surv, data=df, family = "binomial")
  summary_surv <- summary(model_surv)
  
  #fixed effect
  parameter_estimate <- summary_surv$coefficients$cond["scale(delta_meth)", "Estimate"]
  parameter_se <- summary_surv$coefficients$cond["scale(delta_meth)","Std. Error"]
  parameter_zval <- summary_surv$coefficients$cond["scale(delta_meth)","z value"]
  parameter_pval <- summary_surv$coefficients$cond["scale(delta_meth)", "Pr(>|z|)"]
  
  #age effect
  age_estimate <- summary_surv$coefficients$cond["age", "Estimate"]
  age_se <- summary_surv$coefficients$cond["age", "Std. Error"]
  age_zval <- summary_surv$coefficients$cond["age", "z value"]
  age_pval <- summary_surv$coefficients$cond["age", "Pr(>|z|)"]
  
  #premeth effect
  pre_estimate <- summary_surv$coefficients$cond["scale(methperc_pre)", "Estimate"]
  pre_se <- summary_surv$coefficients$cond["scale(methperc_pre)", "Std. Error"]
  pre_zval <- summary_surv$coefficients$cond["scale(methperc_pre)", "z value"]
  pre_pval <- summary_surv$coefficients$cond["scale(methperc_pre)", "Pr(>|z|)"]
  
  rsqc <- performance::r2(model_surv)$R2_conditional #fixed plus random parameterects relative to overall variance
  rsqm <- performance::r2(model_surv)$R2_marginal #fixed parameterects relative to overall variance
  
  message <- model_surv$fit$message
  
  icc_id_site <-icc(model_surv, by_group = TRUE, tolerance = 0)[1,2]
  icc_site <-icc(model_surv, by_group = TRUE, tolerance = 0)[2,2]
  
  surv <- data.frame(surv_icc_id_site = as.numeric(icc_id_site),
                    surv_icc_site = as.numeric(icc_site),
                    surv_delta_meth_estimate = as.numeric(parameter_estimate),
                    surv_delta_meth_se = as.numeric(parameter_se),
                    surv_delta_meth_zval = as.numeric(parameter_zval),
                    surv_delta_meth_pval = as.numeric(parameter_pval),
                    surv_age_estimate = as.numeric(age_estimate),
                    surv_age_se = as.numeric(age_se),
                    surv_age_zval = as.numeric(age_zval),
                    surv_age_pval = as.numeric(age_pval),
                    surv_pre_estimate = as.numeric(pre_estimate),
                    surv_pre_se = as.numeric(pre_se),
                    surv_pre_zval = as.numeric(pre_zval),
                    surv_pre_pval = as.numeric(pre_pval),
                    surv_rsqc = as.numeric(rsqc),
                    surv_rsqm = as.numeric(rsqm),
                    surv_message = message              
  ) 
  out <- cbind(ams, surv)
  return(out)
  
}, error = function(e){cat("ERROR :", conditionMessage(e), "\n");print(chr_pos)})
}

### run the model 
# run model
delta_out_fitness <- parallel::mclapply(delta_meth_ls, function_model_fitness,mc.cores=12)
delta_out_fitness <- do.call(rbind.data.frame, delta_out_fitness)

# convert to numeric
delta_out_fitness$ams_delta_meth_estimate <- as.numeric(delta_out_fitness$ams_delta_meth_estimate)
delta_out_fitness$surv_delta_meth_estimate <- as.numeric(delta_out_fitness$surv_delta_meth_estimate)

# exclude those with convergence errors
delta_out_ams <- subset(delta_out_fitness, ams_message == "relative convergence (4)")
delta_out_surv <- subset(delta_out_fitness, surv_message == "relative convergence (4)")

# FDR correction
delta_out_ams$ams_delta_meth_qval <- p.adjust(delta_out_ams$ams_delta_meth_pval, method = "fdr", n = nrow(delta_out_ams))
delta_out_surv$surv_delta_meth_qval <- p.adjust(delta_out_surv$surv_delta_meth_pval, method = "fdr", n = nrow(delta_out_surv))

delta_out_ams$chr_pos <- as.factor(delta_out_ams$chr_pos)
delta_out_surv$chr_pos <- as.factor(delta_out_surv$chr_pos)

delta_out_ams <- delta_out_ams %>% select(c(chr_pos, ams_icc_id_site:ams_message, ams_delta_meth_qval))
delta_out_surv <- delta_out_surv %>% select(c(chr_pos, surv_icc_id_site:surv_message, surv_delta_meth_qval))

save(delta_out_ams, file="results/modeloutput/AMS_deltameth_modeloutput_filtered.RData")
save(delta_out_surv, file="results/modeloutput/surv_deltameth_modeloutput_filtered.RData")

nrow(subset(delta_out_ams, ams_delta_meth_qval < 0.05))
nrow(subset(delta_out_surv, surv_delta_meth_qval < 0.05))

### significant ones

cpg_sig_ams <- subset(delta_out_ams, ams_delta_meth_qval < 0.05)

### plotting

source("scripts/plotting_theme.R")

cpg_sig_ams <- cpg_sig_ams %>% arrange(ams_delta_meth_qval)

for (i in 1:10){
    ggplot(subset(delta_meth, chr_pos == cpg_sig_ams$chr_pos[i]), aes(x = delta_meth, y = MS)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression(Delta*" methylation"), y = "Annual mating success",
                      title = paste0("Estimate = ", round(cpg_sig_ams$ams_delta_meth_estimate[i], 2),
                                        ", q-value = ", round(cpg_sig_ams$ams_delta_meth_qval[i], 2))) +
                                        geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
    
    ggsave(plot, file = paste0("plots/model_out/ams/rawdata_plot_AMS_cpg_", i, ".png"), width=8, height=8)}

### is there any overlap with behavioural or physiological traits?
load(file="results/modeloutput/effort_deltameth_modeloutput_filtered.RData")
cpg_effort <- subset(delta_out_all, parameter_qval < 0.05)
  
load(file="results/modeloutput/physio_deltameth_modeloutput_filtered.RData")
cpg_physio <- subset(delta_out_all, parameter_qval < 0.05)
  
subset(cpg_sig_ams, chr_pos %in% cpg_effort$chr_pos)  
subset(cpg_sig_ams, chr_pos %in% cpg_physio$chr_pos)  
subset(cpg_physio, chr_pos %in% cpg_sig_ams$chr_pos)  

### merge the data of those
overlap_ams_physio <- inner_join(cpg_physio, cpg_sig_ams, by="chr_pos")

delta_meth <- left_join(delta_meth, unique(prepost_dif[,c("id", "year",  "trypa_dif")], by = c("id", "year")))
   
ggplot(subset(delta_meth, chr_pos == overlap_ams_physio$chr_pos[1]), aes(x = trypa_dif, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression(Delta*" Trypanosoma spp."), y = expression(Delta*" methylation"),
                      title = paste0("Estimate = ", round(overlap_ams_physio$parameter_estimate[1], 2),
                                        ", q-value = ", round(overlap_ams_physio$parameter_qval[1], 2))) +
                                        geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_trypa_meth_1

ggplot(subset(delta_meth, chr_pos == overlap_ams_physio$chr_pos[1]), aes(x = delta_meth, y = MS)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(y = "Annual mating success", x = expression(Delta*" methylation"),
                      title = paste0("Estimate = ", round(overlap_ams_physio$ams_delta_meth_estimate[1], 2),
                                        ", q-value = ", round(overlap_ams_physio$ams_delta_meth_qval[1], 2))) +
                                        geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_meth_ms_1

ggplot(subset(delta_meth, chr_pos == overlap_ams_physio$chr_pos[2]), aes(x = trypa_dif, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression(Delta*" Trypanosoma spp."), y = expression(Delta*" methylation"),
                      title = paste0("Estimate = ", round(overlap_ams_physio$parameter_estimate[2], 2),
                                        ", q-value = ", round(overlap_ams_physio$parameter_qval[2], 2))) +
                                        geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_trypa_meth_2

ggplot(subset(delta_meth, chr_pos == overlap_ams_physio$chr_pos[2]), aes(x = delta_meth, y = MS)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(y = "Annual mating success", x = expression(Delta*" methylation"),
                      title = paste0("Estimate = ", round(overlap_ams_physio$ams_delta_meth_estimate[2], 2),
                                        ", q-value = ", round(overlap_ams_physio$ams_delta_meth_qval[2], 2))) +
                                        geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1) -> plot_meth_ms_2

cowplot::plot_grid(plot_trypa_meth_1, plot_meth_ms_1, plot_trypa_meth_2, plot_meth_ms_2, 
    labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_physio_ams

ggsave(plots_physio_ams, file = paste0("plots/model_out/rawdata_plot_trypa_delta_AMS.png"), width=12, height=12)
                                
        