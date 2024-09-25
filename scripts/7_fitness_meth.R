### load packages
pacman::p_load(tidyverse, data.table, tibble, performance, gaston, cowplot,
               parallel, performance, lmerTest, tidystats, glmmTMB, DHARMa)

### load data

load(file="results/modeloutput/changing/changing_sites_glmer.RData")

### load phenotypic data

load(file = "data/phenotypes/fulldata_complete_epi_withdates.RData")

load("data/phenotypes/pheno_dif_prepost.RData") ## differences in physiology

### methylation difference

load(file = "results/modeloutput/all_sites_deltameth.RData")

delta_meth_raw <- subset(delta_meth, chr_pos %in% changing_cpg$chr_pos)

### overdisp function
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

### theme
source("scripts/plotting_theme.R")

#combine with site and fitness data
pheno_pre <- subset(all_pheno_epi, prepost=="pre")

delta_meth <- left_join(delta_meth_raw, unique(pheno_pre[,c("id", "year", "MS", "surv", "site", "attend", "fight", "dist")]), by = c("id", "year"))
delta_meth <- left_join(delta_meth, unique(prepost_dif[,c("id", "year", "mass_dif", "trypa_dif", "ig_dif", "hct_dif", "microf_dif")]), by = c("id", "year"))
       
## only select cpg sites with enough data
delta_meth_n <- delta_meth %>% group_by(chr_pos) %>% filter(!is.na(delta_meth) & !is.na(MS) & !is.na(surv)) %>% tally()
delta_meth_n_min20 <- subset(delta_meth_n, n > 20)

delta_meth_sub <- subset(delta_meth, chr_pos %in% delta_meth_n_min20$chr_pos) #777
delta_meth_ls <- delta_meth_sub %>% group_split(chr_pos)

## test for zero inflation
model_zi <- glmmTMB(MS ~ 1 + (1|site/id), family = "poisson", data = prepost_dif)
DHARMa::testZeroInflation(model_zi) #ns

# function to run the model
function_model_fitness <- function(df){tryCatch({
  chr_pos <- as.character(df[1,1])
  df <- as.data.frame(df)
  df$site <- as.factor(df$site)
  df$id <- as.factor(df$id)
  
  ### AMS
  formula_ams <- formula(paste0("MS ~ delta_meth + (1|site/id) "))
  model_ams <- glmmTMB(formula_ams, data=df, family = "poisson")
  summary_ams <- summary(model_ams)
  
  intercept_ams <- summary_ams$coefficients$cond["(Intercept)", "Estimate"]

  #fixed effect
  parameter_estimate <- summary_ams$coefficients$cond["delta_meth", "Estimate"]
  parameter_se <- summary_ams$coefficients$cond["delta_meth","Std. Error"]
  parameter_zval <- summary_ams$coefficients$cond["delta_meth","z value"]
  parameter_pval <- summary_ams$coefficients$cond["delta_meth", "Pr(>|z|)"]
  
  message_ams <- model_ams$fit$message
  dispersion_ams <- overdisp_fun(model_ams)
 
  ams <- data.frame(chr_pos=as.factor(chr_pos),
                    intercept_ams = intercept_ams,
                    ams_delta_meth_estimate = as.numeric(parameter_estimate),
                    ams_delta_meth_se = as.numeric(parameter_se),
                    ams_delta_meth_zval = as.numeric(parameter_zval),
                    ams_delta_meth_pval = as.numeric(parameter_pval),
                    ams_message = message_ams,
                    ams_disp_chi = dispersion_ams[1][[1]],
                    ams_disp_ratio = dispersion_ams[2][[1]],
                    ams_disp_rdf = dispersion_ams[3][[1]],
                    ams_disp_p = dispersion_ams[4][[1]]
  ) 
  
  ### surv
  formula_surv <- formula(paste0("surv ~ delta_meth + (1|site/id) "))
  model_surv <- glmmTMB(formula_surv, data=df, family = "binomial")
  summary_surv <- summary(model_surv)
  
  intercept_surv <- summary_surv$coefficients$cond["(Intercept)", "Estimate"]

  #fixed effect
  parameter_estimate <- summary_surv$coefficients$cond["delta_meth", "Estimate"]
  parameter_se <- summary_surv$coefficients$cond["delta_meth","Std. Error"]
  parameter_zval <- summary_surv$coefficients$cond["delta_meth","z value"]
  parameter_pval <- summary_surv$coefficients$cond["delta_meth", "Pr(>|z|)"]
  
  message <- model_surv$fit$message
  
  dispersion_surv <- overdisp_fun(model_surv)
  
  surv <- data.frame(intercept_surv = intercept_surv,
                  surv_delta_meth_estimate = as.numeric(parameter_estimate),
                    surv_delta_meth_se = as.numeric(parameter_se),
                    surv_delta_meth_zval = as.numeric(parameter_zval),
                    surv_delta_meth_pval = as.numeric(parameter_pval),
                    surv_message = message,
                  surv_disp_chi = dispersion_surv[1][[1]],
                  surv_disp_ratio = dispersion_surv[2][[1]],
                  surv_disp_rdf = dispersion_surv[3][[1]],
                  surv_disp_p = dispersion_surv[4][[1]]
  ) 
  out <- cbind(ams, surv)
  return(out)
  
}, error = function(e){cat("ERROR :", conditionMessage(e), "\n");print(chr_pos)})
}

### run the model 
# run model
delta_out_fitness <- parallel::mclapply(delta_meth_ls, function_model_fitness,mc.cores=4)
delta_out_fitness <- do.call(rbind.data.frame, delta_out_fitness)

save(delta_out_fitness, file="results/modeloutput/fitness/out_fitness_nopre_raw.RData")

#### Subset models and exclude models that did not converge ####

delta_out_ams <- subset(delta_out_fitness, ams_message == "relative convergence (4)")
delta_out_surv <- subset(delta_out_fitness, surv_message == "relative convergence (4)")

#### AMS ####
nrow(delta_out_ams) / nrow(delta_out_fitness) * 100 # retain 96.4%, 749

##### Overdispersion raw data #####

summary(delta_out_ams$ams_disp_ratio)

# histogram dispersion ratio raw 
ggplot(delta_out_ams, aes(x = ams_disp_ratio)) + geom_histogram() + geom_vline(xintercept = 1., col = "red", linetype = "dotted", linewidth=1) +
scale_y_log10() + labs(title = "Histogram dispersion ratio", subtitle= "Raw model output AMS") -> hist_ams_disp_ratio_raw

# histogram p-values raw
ggplot(delta_out_ams, aes(x = ams_delta_meth_pval)) + geom_histogram() + 
  scale_y_continuous(labels = scales::unit_format(unit = "K", scale = 1e-3)) + 
  labs(title = "Histogram p-values", subtitle="Raw model output AMS")-> hist_ams_pvals_raw

plot_grid(hist_ams_disp_ratio_raw, hist_ams_pvals_raw, labs="auto", align="hv", axis="lb", ncol=1, label_fontface = "plain", label_size = 22)-> hists_ams_raw

ggsave(hists_ams_raw, file = "plots/model_out/fitness/ams/hist_raw_ams.png", width = 12, height = 12)

# qqplot raw

png(file = "plots/model_out/fitness/ams/qqplot_raw_ams.png", width = 800, height = 800)
qqplot.pvalues(delta_out_ams$ams_delta_meth_pval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95) 
dev.off()

## filter for 95 percentile
delta_out_ams_clean <- subset(delta_out_ams, ams_disp_ratio < as.vector(quantile(delta_out_ams$ams_disp_ratio, 0.95)))
nrow(delta_out_ams_clean) # 711

# histogram dispersion ratio filtered 
ggplot(delta_out_ams_clean, aes(x = ams_disp_ratio)) + geom_histogram() + 
scale_y_log10() + labs(title = "Histogram dispersion ratio", subtitle= "Filtered model output AMS") -> hist_ams_disp_ratio_filter

# histogram p-values filtered
ggplot(delta_out_ams_clean, aes(x = ams_delta_meth_pval)) + geom_histogram() + 
  scale_y_continuous(labels = scales::unit_format(unit = "K", scale = 1e-3)) + 
  labs(title = "Histogram p-values", subtitle="Filtered model output AMS")-> hist_ams_pvals_filter

plot_grid(hist_ams_disp_ratio_filter, hist_ams_pvals_filter, labs="auto", align="hv", axis="lb", ncol=1, label_fontface = "plain", label_size = 22)-> hists_ams_filter

ggsave(hists_ams_filter, file = "plots/model_out/fitness/ams/hist_filtered_ams.png", width = 12, height = 12)

# qqplot filtered

png(file = "plots/model_out/fitness/ams/qqplot_filtered_ams.png", width = 800, height = 800)
qqplot.pvalues(delta_out_ams_clean$ams_delta_meth_pval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95) 
dev.off()

#### FDR-correction ####

delta_out_ams_clean$ams_delta_meth_qval <- p.adjust(delta_out_ams_clean$ams_delta_meth_pval, method = "fdr", n = nrow(delta_out_ams_clean))

#### Survival #####
nrow(delta_out_surv) / nrow(delta_out_fitness) * 100 # retain 99.7%, 775

##### Overdispersion raw data #####

summary(delta_out_surv$surv_disp_ratio)

# histogram dispersion ratio raw 
ggplot(delta_out_surv, aes(x = surv_disp_ratio)) + geom_histogram() + geom_vline(xintercept = 1., col = "red", linetype = "dotted", linewidth=1) +
scale_y_log10() + labs(title = "Histogram dispersion ratio", subtitle= "Raw model output surv") -> hist_surv_disp_ratio_raw

# histogram p-values raw
ggplot(delta_out_surv, aes(x = surv_delta_meth_pval)) + geom_histogram() + 
  scale_y_continuous(labels = scales::unit_format(unit = "K", scale = 1e-3)) + 
  labs(title = "Histogram p-values", subtitle="Raw model output surv")-> hist_surv_pvals_raw

plot_grid(hist_surv_disp_ratio_raw, hist_surv_pvals_raw, labs="auto", align="hv", axis="lb", ncol=1, label_fontface = "plain", label_size = 22)-> hists_surv_raw

ggsave(hists_surv_raw, file = "plots/model_out/fitness/surv/hist_raw_surv.png", width = 12, height = 12)

# qqplot raw

png(file = "plots/model_out/fitness/surv/qqplot_raw_surv.png", width = 800, height = 800)
qqplot.pvalues(delta_out_surv$surv_delta_meth_pval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95) 
dev.off()

### -> no overdispersion 

#### FDR-correction ####

delta_out_surv$surv_delta_meth_qval <- p.adjust(delta_out_surv$surv_delta_meth_pval, method = "fdr", n = nrow(delta_out_surv))

#### How many significant? ####
nrow(subset(delta_out_ams_clean, ams_delta_meth_qval < 0.05)) #362
nrow(subset(delta_out_surv, surv_delta_meth_qval < 0.05)) #0

#### Save data ####
delta_out_ams <- delta_out_ams_clean %>% select(c(chr_pos, intercept_ams:ams_disp_p, ams_delta_meth_qval))
delta_out_surv <- delta_out_surv %>% select(c(chr_pos, intercept_surv:surv_disp_p, surv_delta_meth_qval))

save(delta_out_ams, file="results/modeloutput/fitness/out_ams_deltameth_filtered.RData")
save(delta_out_surv, file="results/modeloutput/fitness/out_surv_deltameth_filtered.RData")

### Volcano plot ####
source("scripts/plotting_theme.R")

#ams
load(file="results/modeloutput/AMS_deltameth_modeloutput_filtered.RData")

delta_out_ams <- delta_out_ams %>% mutate(sig = case_when(ams_delta_meth_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))

#clrs <- viridisLite::viridis(6)
ggplot(delta_out_ams, aes(x = ams_delta_meth_estimate, y = -log10(ams_delta_meth_qval))) + geom_point(size=4, alpha=0.5, aes(col = sig)) +
    labs(x = expression("Estimate "*Delta*" methylation"), y = "-log10(q-value)") +
    xlim(-21,21)+
    scale_color_manual(values=c(clrs[5], clrs[17])) +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
  #  geom_vline(xintercept = -0.1, col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") -> volcano_ams

ggsave(volcano_ams, file = "plots/model_out/fitness/ams/volcano_ams.png", width=10, height=10)    

# surv
delta_out_surv <- delta_out_surv %>% mutate(sig = case_when(surv_delta_meth_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))

ggplot(delta_out_surv, aes(x = surv_delta_meth_estimate, y = -log10(surv_delta_meth_qval))) + geom_point(size=4, alpha=0.5, aes(col = sig)) +
    labs(x = expression("Estimate "*Delta*" methylation"), y = "-log10(q-value)") +
    xlim(-10,10)+
    scale_color_manual(values=c("grey60", clrs[4])) +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
  #  geom_vline(xintercept = -0.1, col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") -> volcano_surv

ggsave(volcano_surv, file = "plots/model_out/fitness/surv/volcano_surv.png", width=10, height=10)    

### significant ones

cpg_sig_ams <- subset(delta_out_ams, ams_delta_meth_qval < 0.05) #362

### plotting

source("scripts/plotting_theme.R")
cpg_sig_ams$intercept_ams <- as.numeric(cpg_sig_ams$intercept_ams)
cpg_sig_ams <- cpg_sig_ams %>% arrange(ams_delta_meth_qval)

for (i in c(1:21)){
    ggplot(subset(delta_meth, chr_pos == cpg_sig_ams$chr_pos[i]), aes(x = delta_meth, y = MS)) + 
      geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression(Delta*" methylation"), 
      y = "Annual mating success",
      title = paste0("Estimate = ", round(cpg_sig_ams$ams_delta_meth_estimate[i], 2), ", q-value = ", 
      round(cpg_sig_ams$ams_delta_meth_qval[i], 4))) +
       geom_abline(intercept=cpg_sig_ams$intercept_ams[i], slope = cpg_sig_ams$ams_delta_meth_estimate[i], 
       color=clrs_hunting[2], linewidth=1)+
       geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
    ggsave(plot, file = paste0("plots/model_out/fitness/ams/rawdata_cpg_", i, ".png"), width=10, height=10)     
   }

### with traits ####

# function to run the model
function_model_surv_attend <- function(df){tryCatch({
  chr_pos <- as.character(df[1,1])
  df <- as.data.frame(df)
  df$site <- as.factor(df$site)
  df$id <- as.factor(df$id)
  
  ### surv
  formula_surv <- formula(paste0("surv ~ delta_meth*(attend+mass_dif) + (1|site/id) "))
  model_surv <- glmmTMB(formula_surv, data=df, family = "binomial")
  summary_surv <- summary(model_surv)
  
  intercept_surv <- summary_surv$coefficients$cond["(Intercept)", "Estimate"]
  
  #delta effect
  delta_estimate <- summary_surv$coefficients$cond["delta_meth", "Estimate"]
  delta_se <- summary_surv$coefficients$cond["delta_meth","Std. Error"]
  delta_zval <- summary_surv$coefficients$cond["delta_meth","z value"]
  delta_pval <- summary_surv$coefficients$cond["delta_meth", "Pr(>|z|)"]
  
  #attend effect
  attend_estimate <- summary_surv$coefficients$cond["attend", "Estimate"]
  attend_se <- summary_surv$coefficients$cond["attend","Std. Error"]
  attend_zval <- summary_surv$coefficients$cond["attend","z value"]
  attend_pval <- summary_surv$coefficients$cond["attend", "Pr(>|z|)"]
  
  #mass effect 
  mass_estimate <- summary_surv$coefficients$cond["mass_dif", "Estimate"]
  mass_se <- summary_surv$coefficients$cond["mass_dif","Std. Error"]
  mass_zval <- summary_surv$coefficients$cond["mass_dif","z value"]
  mass_pval <- summary_surv$coefficients$cond["mass_dif", "Pr(>|z|)"]
  
  #delta:attend
  inter_delta_attend_estimate <- summary_surv$coefficients$cond["delta_meth:attend", "Estimate"]
  inter_delta_attend_se <- summary_surv$coefficients$cond["delta_meth:attend","Std. Error"]
  inter_delta_attend_zval <- summary_surv$coefficients$cond["delta_meth:attend","z value"]
  inter_delta_attend_pval <- summary_surv$coefficients$cond["delta_meth:attend", "Pr(>|z|)"]
  
  #delta:mass
  inter_delta_mass_estimate <- summary_surv$coefficients$cond["delta_meth:mass_dif", "Estimate"]
  inter_delta_mass_se <- summary_surv$coefficients$cond["delta_meth:mass_dif","Std. Error"]
  inter_delta_mass_zval <- summary_surv$coefficients$cond["delta_meth:mass_dif","z value"]
  inter_delta_mass_pval <- summary_surv$coefficients$cond["delta_meth:mass_dif", "Pr(>|z|)"]
  
  message <- model_surv$fit$message
  
  dispersion_surv <- overdisp_fun(model_surv)
  
  surv <- data.frame(chr_pos = chr_pos,
                     intercept_surv = intercept_surv,
                     delta_estimate = as.numeric(delta_estimate),
                     delta_se = as.numeric(delta_se),
                     delta_zval = as.numeric(delta_zval),
                     delta_pval = as.numeric(delta_pval),
                     attend_estimate = attend_estimate,
                     attend_se = attend_se,
                     attend_zval = attend_zval,
                     attend_pval = attend_pval,
                     mass_estimate = mass_estimate,
                     mass_se = mass_se,
                     mass_zval = mass_zval,
                     mass_pval = mass_pval,
                     inter_delta_attend_estimate = inter_delta_attend_estimate,
                     inter_delta_attend_se = inter_delta_attend_se,
                     inter_delta_attend_zval = inter_delta_attend_zval,
                     inter_delta_attend_pval = inter_delta_attend_pval,
                     inter_delta_mass_estimate = inter_delta_mass_estimate,
                     inter_delta_mass_se = inter_delta_mass_se,
                     inter_delta_mass_zval = inter_delta_mass_zval,
                     inter_delta_mass_pval = inter_delta_mass_pval,
                     surv_message = message,
                     surv_disp_chi = dispersion_surv[1][[1]],
                     surv_disp_ratio = dispersion_surv[2][[1]],
                     surv_disp_rdf = dispersion_surv[3][[1]],
                     surv_disp_p = dispersion_surv[4][[1]]
  ) 
  
  return(surv)
  
}, error = function(e){cat("ERROR :", conditionMessage(e), "\n");print(chr_pos)})
}


surv_attend_delta <- parallel::mclapply(delta_meth_ls, function_model_surv_attend,mc.cores=4)
surv_attend_delta_raw <- do.call(rbind.data.frame, surv_attend_delta) #777

## only those that converged
surv_attend_delta_conv <- subset(surv_attend_delta_raw, surv_message == "relative convergence (4)")

#hist p-vals
ggplot(surv_attend_delta_conv, aes((delta_pval))) + geom_histogram() -> hist_p_delta
ggplot(surv_attend_delta_conv, aes((attend_pval))) + geom_histogram() -> hist_p_attend
ggplot(surv_attend_delta_conv, aes((mass_pval))) + geom_histogram() -> hist_p_mass
ggplot(surv_attend_delta_conv, aes((inter_delta_attend_pval))) + geom_histogram() -> hist_p_i_delta_attend
ggplot(surv_attend_delta_conv, aes((inter_delta_mass_pval))) + geom_histogram() -> hist_p_i_delta_mass

#disp
ggplot(surv_attend_delta_conv, aes(as.numeric(surv_disp_ratio))) + geom_histogram()-> hist_disp

cowplot::plot_grid(hist_p_delta, hist_p_attend,hist_p_mass, hist_p_i_delta_attend,  hist_p_i_delta_mass , hist_disp,
  labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> summary_model_hist

ggsave(summary_model_hist, file = "plots/model_out/fitness/surv/hist_raw_surv_with_pheno.png", width=10, height=18)

#qqplot
png(file = "plots/model_out/fitness/surv/qqplot_filtered_surv_with_pheno_delta.png", width = 800, height = 800)
qqplot.pvalues(surv_attend_delta_conv$delta_pval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95)
dev.off()

png(file = "plots/model_out/fitness/surv/qqplot_filtered_surv_with_pheno_attend.png", width = 800, height = 800)
qqplot.pvalues(surv_attend_delta_conv$attend_pval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95)
dev.off()

png(file = "plots/model_out/fitness/surv/qqplot_filtered_surv_with_pheno_mass.png", width = 800, height = 800)
qqplot.pvalues(surv_attend_delta_conv$mass_pval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95)
dev.off()

png(file = "plots/model_out/fitness/surv/qqplot_filtered_surv_with_pheno_i_delta_attend.png", width = 800, height = 800)
qqplot.pvalues(surv_attend_delta_conv$inter_delta_attend_pval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95)
dev.off()

png(file = "plots/model_out/fitness/surv/qqplot_filtered_surv_with_pheno_i_delta_mass.png", width = 800, height = 800)
qqplot.pvalues(surv_attend_delta_conv$inter_delta_mass_pval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95)
dev.off()


# FDR correction
surv_attend_delta_conv$delta_qval <- p.adjust(surv_attend_delta_conv$delta_pval, method = "fdr", n = nrow(surv_attend_delta_conv))
surv_attend_delta_conv$attend_qval <- p.adjust(surv_attend_delta_conv$attend_pval, method = "fdr", n = nrow(surv_attend_delta_conv))
surv_attend_delta_conv$mass_qval <- p.adjust(surv_attend_delta_conv$mass_pval, method = "fdr", n = nrow(surv_attend_delta_conv))
surv_attend_delta_conv$inter_delta_attend_qval <- p.adjust(surv_attend_delta_conv$inter_delta_attend_pval, method = "fdr", n = nrow(surv_attend_delta_conv))
surv_attend_delta_conv$inter_delta_mass_qval <- p.adjust(surv_attend_delta_conv$inter_delta_mass_pval, method = "fdr", n = nrow(surv_attend_delta_conv))

nrow(subset(surv_attend_delta_conv, delta_qval < 0.05)) # n = 0
nrow(subset(surv_attend_delta_conv, attend_qval < 0.05)) # n = 0
nrow(subset(surv_attend_delta_conv, mass_qval < 0.05)) # n = 0
nrow(subset(surv_attend_delta_conv, inter_delta_attend_qval < 0.05)) # n = 0
nrow(subset(surv_attend_delta_conv, inter_delta_mass_qval < 0.05)) # n = 0

#### Annotate AMS CpG sites ####

### Packages ####
pacman::p_load(genomation, GenomicFeatures, rtracklayer, 
               GenomicRanges)


### Combine all sites vs changing sites
cpg_all <- delta_out_ams %>% dplyr::select(c(chr_pos, ams_delta_meth_qval))
names(cpg_all)[2] <- "parameter_qval"
cpg_all$parameter <- "all"

cpg_ams_select <- cpg_sig_ams %>% dplyr::select(c(chr_pos, ams_delta_meth_qval))
names(cpg_ams_select)[2] <- "parameter_qval"
cpg_ams_select$parameter <- "ams"

all_models_sig <- rbind(cpg_all, cpg_ams_select)

### Rename chr_pos and divide ###
all_models_sig$chr_pos <- gsub("__", ";", all_models_sig$chr_pos)
all_models_sig$chr_pos <- gsub("HRSCAF_", "HRSCAF=", all_models_sig$chr_pos, )

# Extract the numbers following HRSCAF=XXX_number
# Split the chr_pos column into two columns based on the first "_"
split_chr_pos <- strsplit(all_models_sig$chr_pos, "_", fixed = TRUE)

all_models_sig$chr <- paste0(sapply(split_chr_pos, "[", 1), "_",
                             sapply(split_chr_pos, "[", 2))

all_models_sig$pos <- sapply(split_chr_pos, "[", 3)

all_models_sig <- all_models_sig %>% 
  relocate(chr, .after = chr_pos) %>%
  relocate(pos, .after = chr_pos)

#revert scafnames
all_models_sig$chr_pos <- gsub(";","__", all_models_sig$chr_pos)
all_models_sig$chr_pos <- gsub("HRSCAF=", "HRSCAF_", all_models_sig$chr_pos)

all_models_sig$chr <- gsub(";","__", all_models_sig$chr)
all_models_sig$chr <- gsub("HRSCAF=", "HRSCAF_", all_models_sig$chr)

### Load annotation data
annotation_dir <- "~/PhD_grouse/grouse-annotation/output"

promoter=unique(gffToGRanges(paste0(annotation_dir, "/promoters.gff3")))
genes=unique(gffToGRanges(paste0(annotation_dir, "/genes.gff3")))
TSS=unique(gffToGRanges(paste0(annotation_dir, "/TSS.gff3")))
exons_gene=unique(gffToGRanges(paste0(annotation_dir, "/exons_gene.gff3")))
introns=unique(gffToGRanges(paste0(annotation_dir, "/introns_transcripts.gff3")))
downstream=unique(gffToGRanges(paste0(annotation_dir, "/downstream.gff3")))
upstream=unique(gffToGRanges(paste0(annotation_dir, "/upstream.gff3")))
threeUTR =unique(gffToGRanges(paste0(annotation_dir, "/threeUTRs.gff3")))
fiveUTR=unique(gffToGRanges(paste0(annotation_dir, "/fiveUTRs.gff3")))

#### Annotate ####
all_models_sig$end <- all_models_sig$pos
all_models_sig$start <- all_models_sig$pos
sig_gr <- as(all_models_sig, "GRanges")

sig_promoter <- subsetByOverlaps(sig_gr, promoter) %>% as.data.frame() %>%
  add_column("region" = "promoter", .after="parameter") %>% 
  dplyr::select(-c(seqnames:strand)) 

sig_gene <- as.data.frame(subsetByOverlaps(sig_gr, genes)) %>% as.data.frame() %>%
  add_column("region" = "gene", .after="parameter") %>% 
  dplyr::select(-c(seqnames:strand)) 

sig_tss <- as.data.frame(subsetByOverlaps(sig_gr, TSS)) %>% as.data.frame() %>%
  add_column("region" = "TSS", .after="parameter") %>%
  dplyr::select(-c(seqnames:strand)) 

sig_exon <- as.data.frame(subsetByOverlaps(sig_gr, exons_gene)) %>% as.data.frame() %>%
  add_column("region" = "exon", .after="parameter") %>% 
  dplyr::select(-c(seqnames:strand)) 

sig_intron <- as.data.frame(subsetByOverlaps(sig_gr, introns))  %>% as.data.frame() %>%
  add_column("region" = "intron", .after="parameter") %>% 
  dplyr::select(-c(seqnames:strand)) 

sig_down <- as.data.frame(subsetByOverlaps(sig_gr, downstream)) %>% as.data.frame() %>%
  add_column("region" = "downstream", .after="parameter") %>% 
  dplyr::select(-c(seqnames:strand)) 

sig_up <- as.data.frame(subsetByOverlaps(sig_gr, upstream))  %>% as.data.frame() %>%
  add_column("region" = "upstream", .after="parameter") %>% 
  dplyr::select(-c(seqnames:strand)) 

sig_threeUTR <- as.data.frame(subsetByOverlaps(sig_gr, threeUTR))  %>% as.data.frame() %>%
  add_column("region" = "threeUTR", .after="parameter") %>% 
  dplyr::select(-c(seqnames:strand)) 

sig_fiveUTR <- as.data.frame(subsetByOverlaps(sig_gr, fiveUTR))  %>% as.data.frame() %>%
  add_column("region" = "fiveUTR", .after="parameter") %>% 
  dplyr::select(-c(seqnames:strand)) 

all_models_sig_annotated <- rbind(sig_promoter, sig_gene,
                                  sig_tss, sig_exon, sig_intron, sig_down,
                                  sig_up, sig_threeUTR,  sig_fiveUTR)

summary(as.factor(all_models_sig_annotated$region))

save(all_models_sig_annotated, file="results/modeloutput/fitness/annotated_sig_cpg_ams.RData")

#### Summarise number of sites per region ####
sum_annotated <- as.data.frame(table(as.factor(all_models_sig_annotated$region), all_models_sig_annotated$parameter))
names(sum_annotated) <- c("region", "model", "n")

sum_annotated$model <- gsub("all", "All", sum_annotated$model)
sum_annotated$model <- gsub("ams", "Annual mating success", sum_annotated$model)

sum_annotated$region <- gsub("downstream", "Downstream", sum_annotated$region)
sum_annotated$region <- gsub("upstream", "Upstream", sum_annotated$region)
sum_annotated$region <- gsub("exon", "Exon", sum_annotated$region)
sum_annotated$region <- gsub("fiveUTR", "5' UTR", sum_annotated$region)
sum_annotated$region <- gsub("gene", "Gene body", sum_annotated$region)
sum_annotated$region <- gsub("intron", "Intron", sum_annotated$region)
sum_annotated$region <- gsub("promoter", "Promoter", sum_annotated$region)
sum_annotated$region <- gsub("threeUTR", "3' UTR", sum_annotated$region)

sum_annotated$region <- factor(sum_annotated$region, levels = c("3' UTR", "5' UTR", "Downstream", "Upstream", "Gene body", "Exon", "Intron", "Promoter", "TSS"))

# add total sig CpGs
sum_annotated <- sum_annotated %>% mutate(n_total = case_when(
  model == "All" ~ nrow(delta_out_ams),
  model == "Annual mating success" ~ nrow(cpg_sig_ams)))

sum_annotated <- sum_annotated %>% mutate(perc = n / n_total * 100)

write.csv(sum_annotated, file="results/modeloutput/fitness/summary_regions_sig_cpgs_ams.csv", row.names=F, quote=F)

#### Plot number of sites per region ####
clrs <- viridisLite::viridis(6)
ggplot(subset(sum_annotated, region != "5' UTR" & region != "3' UTR"), aes(x = region, y = perc, fill = model)) + geom_bar(stat="identity", position="dodge") + 
  labs(y="Percentage of CpG sites", x="Region", fill = "Subset")+ coord_flip() + 
  scale_fill_manual(values=c(clrs[5], clrs[17])) -> num_perc_region
ggsave(num_perc_region, file="plots/model_out/fitness/perc_sig_per_region_ams.png", width=10, height=12)


#### Based on the 'similar to' column ####

sig_promoter <- mergeByOverlaps(sig_gr, promoter) 
sig_promoter_df <- unique(data.frame(chr_pos = sig_promoter$chr_pos,
                                     # chr = sig_promoter$sig_gr@seqnames,
                                     # pos = sig_promoter$sig_gr$pos,
                                      parameter = sig_promoter$parameter,
                                      parameter_qval = sig_promoter$parameter_qval,
                                      gene_id = sig_promoter$gene_id,
                                      region = "promoter"))

sig_gene <- mergeByOverlaps(sig_gr, genes) 
sig_gene_df <- unique(data.frame(chr_pos = sig_gene$chr_pos,
                                          #   chr = sig_gene$sig_gr@seqnames,
                                          #   pos = sig_gene$sig_gr$pos,
                                             parameter = sig_gene$parameter,
                                      parameter_qval = sig_gene$parameter_qval,
                                           gene_id = sig_gene$gene_id,
                                             region = "gene"))

sig_TSS <- mergeByOverlaps(sig_gr, TSS) 
sig_TSS_df <- unique(data.frame(chr_pos = sig_TSS$chr_pos,
                                        # chr = sig_TSS$sig_gr@seqnames,
                                        # pos = sig_TSS$sig_gr$pos,
                                         parameter = sig_TSS$parameter,
                                      parameter_qval = sig_TSS$parameter_qval,
                                       gene_id = unlist(sig_TSS$gene_id),
                                         region = "TSS"))

sig_exon <- mergeByOverlaps(sig_gr, exons_gene) 
sig_exon_df <- unique(data.frame(chr_pos = sig_exon$chr_pos,
                                       # chr = sig_exon$sig_gr@seqnames,
                                       # pos = sig_exon$sig_gr$pos,
                                        parameter = sig_exon$parameter,
                                      parameter_qval = sig_exon$parameter_qval,
                                        gene_id = sig_exon$exon_name,
                                        region = "exon"))

sig_intron <- mergeByOverlaps(sig_gr, introns) 
sig_intron_df <- unique(data.frame(chr_pos = sig_intron$chr_pos,
                                      #   chr = sig_intron$sig_gr@seqnames,
                                      #   pos = sig_intron$sig_gr$pos,
                                         parameter = sig_intron$parameter,
                                      parameter_qval = sig_intron$parameter_qval,
                                        gene_id = NA,
                                         region = "intron"))

sig_down <- mergeByOverlaps(sig_gr, downstream) 
sig_down_df <- unique(data.frame(chr_pos = sig_down$chr_pos,
                                    #    chr = sig_down$sig_gr@seqnames,
                                     #   pos = sig_down$sig_gr$pos,
                                        parameter = sig_down$parameter,
                                      parameter_qval = sig_down$parameter_qval,
                                         gene_id = sig_down$gene_id,
                                        region = "downstream"))

sig_up <- mergeByOverlaps(sig_gr, upstream) 
sig_up_df <- unique(data.frame(chr_pos = sig_up$chr_pos,
                                    #     chr = sig_up$sig_gr@seqnames,
                                      #   pos = sig_up$sig_gr$pos,
                                        parameter = sig_up$parameter,
                                      parameter_qval = sig_up$parameter_qval,
                                       gene_id = sig_up$gene_id,
                                         region = "upstream"))


all_models_sig_annotated <- rbind(sig_promoter_df, sig_gene_df,
                                          sig_TSS_df, sig_exon_df, sig_down_df,
                                  sig_up_df) # left out intron due to error

### merge with 'similar to' column
lookup <- fread("/home/nioo/rebeccash/PhD_grouse/grouse-annotation/data/lookup_ANN_gene_id.txt")
names(lookup) <- c("original", "similar")

all_models_sig_annotated <- left_join(all_models_sig_annotated, lookup, by = c("gene_id" = "original"))

all_models_sig_annotated <- subset(all_models_sig_annotated,
                                           !is.na(similar) &
                                             !grepl("LOC", similar))

all_models_sig_annotated$similar <- toupper(all_models_sig_annotated$similar)
subset(all_models_sig_annotated, parameter == "all") %>% arrange(parameter_qval) %>%
  dplyr::select(similar) %>% unique() %>% write.csv("results/modeloutput/fitness/all_gene_ids_changing_similar.csv", quote=F, row.names=F, col.names = F)

subset(all_models_sig_annotated, parameter == "ams") %>% arrange(parameter_qval) %>%
  dplyr::select(similar) %>% unique() %>% write.csv("results/modeloutput/fitness/gene_ids_sig_ams_similar.csv", quote=F, row.names=F, col.names = F)
