### load packages
pacman::p_load(tidyverse, data.table, tibble, performance, matrixStats, 
               parallel, performance, lmerTest, tidystats, insight)

### load data

load(file = "results/modeloutput/subset_sites_sig_deltameth.RData")

### load phenotypic data

load("data/phenotypes/fulldata_complete_epi_withdates.RData")
load("data/phenotypes/pheno_dif_prepost.RData")

#combine with site info

delta_meth <- left_join(delta_meth, unique(all_pheno_epi[,c("id", "year", "site", "Core")], by = c("id", "year")))

#combine with delta physiology
prepost_dif$mass_dif_scl <- scale(prepost_dif$mass_dif)
prepost_dif$microf_dif_scl <- scale(prepost_dif$microf_dif)
prepost_dif$trypa_dif_scl <- scale(prepost_dif$trypa_dif)
prepost_dif$ig_dif_scl <- scale(prepost_dif$ig_dif)
prepost_dif$hct_dif_scl <- scale(prepost_dif$hct_dif)

delta_meth <- left_join(delta_meth, unique(prepost_dif[,c("id", "year", "mass_dif", "microf_dif", "trypa_dif", "ig_dif", "hct_dif",
                        "mass_dif_scl", "microf_dif_scl", "trypa_dif_scl", "ig_dif_scl", "hct_dif_scl")], by = c("id", "year")))

delta_meth_ls <- delta_meth %>% group_split(chr_pos)

# function to run the model
function_model_delta <- function(df, parameter){tryCatch({
  chr_pos <- as.character(df[1,1])
  df <- as.data.frame(df)
  df$methperc_pre_scl <- scale(df$methperc_pre)

  formula <- formula(paste0("scale(delta_meth) ~ ", parameter, "_dif_scl + age + scale(methperc_pre) + (1|site/id) "))
  
  model <- lmerTest::lmer(formula, data=df)
  summary <- summary(model)
  
  overdisp.lmer_fun <- function(model) {
    vpars <- function(m) {
      nrow(m)*(nrow(m)+1)/2
    }
    model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
    rdf <- nrow(model.frame(model))-model.df
    rp <- residuals(model,type="pearson")
    Pearson.chisq <- sum(rp^2)
    prat <- Pearson.chisq/rdf
    pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
    data.frame(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
  }
  
  #fixed effect
  parameter_estimate <- summary$coefficients[2,1]
  parameter_se <- summary$coefficients[2,2]
  parameter_df <- summary$coefficients[2,3]
  parameter_tval <- summary$coefficients[2,4]
  parameter_pval <- summary$coefficients[2,5]
  
  #age effect
  age_estimate <- summary$coefficients["age", "Estimate"]
  age_se <- summary$coefficients["age", "Std. Error"]
  age_df <- summary$coefficients["age", "df"]
  age_tval <- summary$coefficients["age", "t value"]
  age_pval <- summary$coefficients["age", "Pr(>|t|)"]
  
  #premeth effect
  pre_estimate <- summary$coefficients["scale(methperc_pre)", "Estimate"]
  pre_se <- summary$coefficients["scale(methperc_pre)", "Std. Error"]
  pre_df <- summary$coefficients["scale(methperc_pre)", "df"]
  pre_tval <- summary$coefficients["scale(methperc_pre)", "t value"]
  pre_pval <- summary$coefficients["scale(methperc_pre)", "Pr(>|t|)"]
  
  rsqc <- performance::r2(model)$R2_conditional #fixed plus random parameterects relative to overall variance
  rsqm <- performance::r2(model)$R2_marginal #fixed parameterects relative to overall variance
  
  dispersion.chisq <- overdisp.lmer_fun(model)[1,1]
  dispersion.ratio <- overdisp.lmer_fun(model)[1,2]
  dispersion.rdf <- overdisp.lmer_fun(model)[1,3]
  dispersion.pval <- overdisp.lmer_fun(model)[1,4]
  
  isSingular <- isSingular(model)
  
  icc_id_site <-icc(model, by_group = TRUE, tolerance = 0)[1,2]
  icc_site <-icc(model, by_group = TRUE, tolerance = 0)[2,2]
  
  return(data.frame(chr_pos=as.factor(chr_pos),
                    parameter = as.factor(parameter),
                    icc_id_site = as.numeric(icc_id_site),
                    icc_site = as.numeric(icc_site),
                    parameter_estimate = as.numeric(parameter_estimate),
                    parameter_se = as.numeric(parameter_se),
                    parameter_df = as.numeric(parameter_df),
                    parameter_tval = as.numeric(parameter_tval),
                    parameter_pval = as.numeric(parameter_pval),
                    age_estimate = as.numeric(age_estimate),
                    age_se = as.numeric(age_se),
                    age_df = as.numeric(age_df),
                    age_tval = as.numeric(age_tval),
                    age_pval = as.numeric(age_pval),
                    pre_estimate = as.numeric(pre_estimate),
                    pre_se = as.numeric(pre_se),
                    pre_df = as.numeric(pre_df),
                    pre_tval = as.numeric(pre_tval),
                    pre_pval = as.numeric(pre_pval),
                    rsqc = as.numeric(rsqc),
                    rsqm = as.numeric(rsqm),
                    dispersion.chisq = as.numeric(dispersion.chisq),
                    dispersion.ratio = as.numeric(dispersion.ratio),
                    dispersion.rdf = as.numeric(dispersion.rdf),
                    dispersion.pval = as.numeric(dispersion.pval),
                    isSingular = as.logical(isSingular)
                    
  ))
}, error = function(e){cat("ERROR :", conditionMessage(e), "\n");print(chr_pos)})
}

### run the model per trait

### mass
# run model
delta_out_mass <- parallel::mclapply(delta_meth_ls, function_model_delta, parameter="mass",mc.cores=12)
delta_out_mass <- do.call(rbind.data.frame, delta_out_mass)

# convert to numeric
delta_out_mass$parameter_pval <- as.numeric(delta_out_mass$parameter_pval)
delta_out_mass$age_pval <- as.numeric(delta_out_mass$age_pval)

# exclude those with overdispersion
delta_out_mass <- subset(delta_out_mass, dispersion.ratio < 1.1 & dispersion.pval > 0.05)

# FDR correction
delta_out_mass$parameter_qval <- p.adjust(delta_out_mass$parameter_pval, method = "fdr", n = nrow(delta_out_mass))

delta_out_mass$age_qval <- p.adjust(delta_out_mass$age_pval, method = "fdr", n = nrow(delta_out_mass))

### microf
# run model
delta_out_microf <- parallel::mclapply(delta_meth_ls, function_model_delta, parameter="microf",mc.cores=12)
delta_out_microf <- do.call(rbind.data.frame, delta_out_microf)

# to numeric
delta_out_microf$parameter_pval <- as.numeric(delta_out_microf$parameter_pval)
delta_out_microf$age_pval <- as.numeric(delta_out_microf$age_pval)

# exclude overdispersion
delta_out_microf <- subset(delta_out_microf, isSingular == FALSE & dispersion.ratio < 1.1 & dispersion.pval > 0.05)

# FDR correction
delta_out_microf$parameter_qval <- p.adjust(delta_out_microf$parameter_pval, method = "fdr", n = nrow(delta_out_microf))
delta_out_microf$age_qval <- p.adjust(delta_out_microf$age_pval, method = "fdr", n = nrow(delta_out_microf))

### trypa
# run model
delta_out_trypa <- parallel::mclapply(delta_meth_ls, function_model_delta, parameter="trypa",mc.cores=12)
delta_out_trypa <- do.call(rbind.data.frame, delta_out_trypa)

# as numeric
delta_out_trypa$parameter_pval <- as.numeric(delta_out_trypa$parameter_pval)
delta_out_trypa$age_pval <- as.numeric(delta_out_trypa$age_pval)

# exclude overdispersion
delta_out_trypa <- subset(delta_out_trypa, dispersion.ratio < 1.1 & dispersion.pval > 0.05)

# FDR correction
delta_out_trypa$parameter_qval <- p.adjust(delta_out_trypa$parameter_pval, method = "fdr", n = nrow(delta_out_trypa))
delta_out_trypa$age_qval <- p.adjust(delta_out_trypa$age_pval, method = "fdr", n = nrow(delta_out_trypa))

### ig
# run model
delta_out_ig <- parallel::mclapply(delta_meth_ls, function_model_delta, parameter="ig",mc.cores=12)
delta_out_ig <- do.call(rbind.data.frame, delta_out_ig)

# as numeric
delta_out_ig$parameter_pval <- as.numeric(delta_out_ig$parameter_pval)
delta_out_ig$age_pval <- as.numeric(delta_out_ig$age_pval)

# exclude overdispersion
delta_out_ig <- subset(delta_out_ig, dispersion.ratio < 1.1 & dispersion.pval > 0.05)

# FDR correction
delta_out_ig$parameter_qval <- p.adjust(delta_out_ig$parameter_pval, method = "fdr", n = nrow(delta_out_ig))
delta_out_ig$age_qval <- p.adjust(delta_out_ig$age_pval, method = "fdr", n = nrow(delta_out_ig))

### hct
# run model
delta_out_hct <- parallel::mclapply(delta_meth_ls, function_model_delta, parameter="hct",mc.cores=12)
delta_out_hct <- do.call(rbind.data.frame, delta_out_hct)

# as numeric
delta_out_hct$parameter_pval <- as.numeric(delta_out_hct$parameter_pval)
delta_out_hct$age_pval <- as.numeric(delta_out_hct$age_pval)

# exclude overdispersion
delta_out_hct <- subset(delta_out_hct, dispersion.ratio < 1.1 & dispersion.pval > 0.05)

# FDR correction
delta_out_hct$parameter_qval <- p.adjust(delta_out_hct$parameter_pval, method = "fdr", n = nrow(delta_out_hct))
delta_out_hct$age_qval <- p.adjust(delta_out_hct$age_pval, method = "fdr", n = nrow(delta_out_hct))

### combine
delta_out_all <- rbind(delta_out_mass, delta_out_microf, delta_out_trypa, delta_out_ig, delta_out_hct)
delta_out_all$chr_pos <- as.factor(delta_out_all$chr_pos)
delta_out_all$parameter <- as.factor(delta_out_all$parameter)
delta_out_all$isSingular <- as.logical(delta_out_all$isSingular)
delta_out_all[c(3:25, 27,28)] <- lapply(delta_out_all[c(3:25, 27:28)], as.numeric)

save(delta_out_all, file="results/modeloutput/physio_deltameth_modeloutput_filtered.RData")

### volcano plot
source("scripts/plotting_theme.R")

load(file="results/modeloutput/physio_deltameth_modeloutput_filtered.RData")

delta_out_all$parameter <- gsub("mass", "Delta body mass",delta_out_all$parameter)
delta_out_all$parameter <- gsub("microf", "Delta Microfilaria spp.", delta_out_all$parameter)
delta_out_all$parameter <- gsub("trypa", "Delta Trypanosoma spp.", delta_out_all$parameter)
delta_out_all$parameter <- gsub("ig", "Delta IgG", delta_out_all$parameter)
delta_out_all$parameter <- gsub("hct", "Delta HCT", delta_out_all$parameter)

delta_out_all <- delta_out_all %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))

clrs <- viridisLite::viridis(6)
ggplot(delta_out_all, aes(x = parameter_estimate, y = -log10(parameter_qval))) + geom_point(size=4, alpha=0.5, aes(col = sig)) +
    facet_wrap(~parameter, ncol=1) +
    xlim(-1,1)+
    ylim(0,3)+
    labs(x = "Estimate", y = "-log10(q-value)") +
    scale_color_manual(values=c("grey60", clrs[4])) +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = -0.1, col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0.1, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") -> volcano_physio

ggsave(volcano_physio, file = "plots/model_out/volcano_physio.png", width=8, height=18)    

### significant ones

cpg_sig_microf <- subset(delta_out_all, parameter_qval < 0.05 & parameter == "Delta Microfilaria spp." & abs(parameter_estimate) > 0.1)
cpg_sig_trypa <- subset(delta_out_all, parameter_qval < 0.05 & parameter == "Delta Trypanosoma spp."& abs(parameter_estimate) > 0.1)
cpg_sig_ig <- subset(delta_out_all, parameter_qval < 0.05 & parameter == "Delta IgG"& abs(parameter_estimate) > 0.1)
cpg_sig_hct <- subset(delta_out_all, parameter_qval < 0.05 & parameter == "Delta HCT"& abs(parameter_estimate) > 0.1)

### plotting

source("scripts/plotting_theme.R")

### microf

ggplot(subset(delta_meth, chr_pos == cpg_sig_microf$chr_pos[1]), aes(x = microf_dif, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression(Delta*" Microfilaria spp."), y = expression(Delta*" methylation"),
                      title = paste0("Estimate = ", round(cpg_sig_microf$parameter_estimate[1], 2),
                                        ", q-value = ", round(cpg_sig_microf$parameter_qval[1], 2))) +
                                        geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot_microf
    
ggsave(plot_microf, file = paste0("plots/model_out/rawdata_plot_microf_dif_cpg_1.png"), width=8, height=8)

### hct plot all 3 in one
ggplot(subset(delta_meth, chr_pos == cpg_sig_hct$chr_pos[1]), aes(x = hct_dif, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression(Delta*" Haematocrit"), y = expression(Delta*" methylation"),
                      title = paste0("Estimate = ", round(cpg_sig_hct$parameter_estimate[1], 2),
                                        ", q-value = ", round(cpg_sig_hct$parameter_qval[1], 2))) +
                                        geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot_hct_1

ggplot(subset(delta_meth, chr_pos == cpg_sig_hct$chr_pos[2]), aes(x = hct_dif, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression(Delta*" Haematocrit"), y = expression(Delta*" methylation"),
                      title = paste0("Estimate = ", round(cpg_sig_hct$parameter_estimate[2], 2),
                                        ", q-value = ", round(cpg_sig_hct$parameter_qval[2], 2))) +
                                        geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot_hct_2

ggplot(subset(delta_meth, chr_pos == cpg_sig_hct$chr_pos[3]), aes(x = hct_dif, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression(Delta*" Haematocrit"), y = expression(Delta*" methylation"),
                      title = paste0("Estimate = ", round(cpg_sig_hct$parameter_estimate[3], 2),
                                        ", q-value = ", round(cpg_sig_hct$parameter_qval[3], 2))) +
                                        geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot_hct_3

cowplot::plot_grid(plot_hct_1, plot_hct_2, plot_hct_3, labs="auto", align="hv", axis="lb", ncol=1, label_fontface = "plain", label_size = 22) -> plots_hct

ggsave(plots_hct, file = paste0("plots/model_out/rawdata_plot_hct_dif.png"), width=8, height=14)

## trypa
# sort by significance
cpg_sig_trypa <- cpg_sig_trypa %>% arrange(parameter_qval)

# loop over all
for (i in 1:nrow(cpg_sig_trypa)){
    ggplot(subset(delta_meth, chr_pos == cpg_sig_trypa$chr_pos[i]), aes(x = trypa_dif, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression(Delta*" Trypanosoma spp."), y = expression(Delta*" methylation"),
                      title = paste0("Estimate = ", round(cpg_sig_trypa$parameter_estimate[i], 2),
                                        ", q-value = ", round(cpg_sig_trypa$parameter_qval[i], 2))) +
                                        geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
    
    ggsave(plot, file = paste0("plots/model_out/trypa/rawdata_plot_trypa_dif_cpg_", i, ".png"), width=8, height=8)}

# top 4
ggplot(subset(delta_meth, chr_pos == cpg_sig_trypa$chr_pos[1]), aes(x = trypa_dif, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression(Delta*" Trypanosoma spp."), y = expression(Delta*" methylation"),
                      title = paste0("Estimate = ", round(cpg_sig_trypa$parameter_estimate[1], 2),
                                        ", q-value = ", round(cpg_sig_trypa$parameter_qval[1], 2))) +
                                        geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot_trypa_1

ggplot(subset(delta_meth, chr_pos == cpg_sig_trypa$chr_pos[2]), aes(x = trypa_dif, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression(Delta*" Trypanosoma spp."), y = expression(Delta*" methylation"),
                      title = paste0("Estimate = ", round(cpg_sig_trypa$parameter_estimate[2], 2),
                                        ", q-value = ", round(cpg_sig_trypa$parameter_qval[2], 2))) +
                                        geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot_trypa_2

ggplot(subset(delta_meth, chr_pos == cpg_sig_trypa$chr_pos[3]), aes(x = trypa_dif, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression(Delta*" Trypanosoma spp."), y = expression(Delta*" methylation"),
                      title = paste0("Estimate = ", round(cpg_sig_trypa$parameter_estimate[3], 2),
                                        ", q-value = ", round(cpg_sig_trypa$parameter_qval[3], 2))) +
                                        geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot_trypa_3

ggplot(subset(delta_meth, chr_pos == cpg_sig_trypa$chr_pos[4]), aes(x = trypa_dif, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression(Delta*" Trypanosoma spp."), y = expression(Delta*" methylation"),
                      title = paste0("Estimate = ", round(cpg_sig_trypa$parameter_estimate[4], 2),
                                        ", q-value = ", round(cpg_sig_trypa$parameter_qval[4], 2))) +
                                        geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot_trypa_4

cowplot::plot_grid(plot_trypa_1, plot_trypa_2, plot_trypa_3, plot_trypa_4, labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_trypa

ggsave(plots_trypa, file = paste0("plots/model_out/rawdata_plot_trypa_dif_cpg_top.png"), width=14, height=14)

### IgG
# sort by significance
cpg_sig_ig <- cpg_sig_ig %>% arrange(parameter_qval)

# loop over all
for (i in 1:nrow(cpg_sig_ig)){
    ggplot(subset(delta_meth, chr_pos == cpg_sig_ig$chr_pos[i]), aes(x = ig_dif, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression(Delta*" IgG"), y = expression(Delta*" methylation"),
                      title = paste0("Estimate = ", round(cpg_sig_ig$parameter_estimate[i], 2),
                                        ", q-value = ", round(cpg_sig_ig$parameter_qval[i], 2))) +
                                        geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
    
    ggsave(plot, file = paste0("plots/model_out/igg/rawdata_plot_igg_dif_cpg_", i, ".png"), width=8, height=8)
}

# top 4
ggplot(subset(delta_meth, chr_pos == cpg_sig_ig$chr_pos[1]), aes(x = ig_dif, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression(Delta*" IgG"), y = expression(Delta*" methylation"),
                      title = paste0("Estimate = ", round(cpg_sig_ig$parameter_estimate[1], 2),
                                        ", q-value = ", round(cpg_sig_ig$parameter_qval[1], 2))) +
                                        geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot_igg_1

ggplot(subset(delta_meth, chr_pos == cpg_sig_ig$chr_pos[2]), aes(x = ig_dif, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression(Delta*" IgG"), y = expression(Delta*" methylation"),
                      title = paste0("Estimate = ", round(cpg_sig_ig$parameter_estimate[2], 2),
                                        ", q-value = ", round(cpg_sig_ig$parameter_qval[2], 2))) +
                                        geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot_igg_2

ggplot(subset(delta_meth, chr_pos == cpg_sig_ig$chr_pos[3]), aes(x = ig_dif, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression(Delta*" IgG"), y = expression(Delta*" methylation"),
                      title = paste0("Estimate = ", round(cpg_sig_ig$parameter_estimate[3], 2),
                                        ", q-value = ", round(cpg_sig_ig$parameter_qval[3], 2))) +
                                        geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot_igg_3

ggplot(subset(delta_meth, chr_pos == cpg_sig_ig$chr_pos[4]), aes(x = ig_dif, y = delta_meth)) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression(Delta*" IgG"), y = expression(Delta*" methylation"),
                      title = paste0("Estimate = ", round(cpg_sig_ig$parameter_estimate[4], 2),
                                        ", q-value = ", round(cpg_sig_ig$parameter_qval[4], 2))) +
                                        geom_smooth(method="lm", color=clrs_hunting[2], linewidth=1) +
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot_igg_4

cowplot::plot_grid(plot_igg_1, plot_igg_2, plot_igg_3, plot_igg_4, labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_igg

ggsave(plots_igg, file = paste0("plots/model_out/rawdata_plot_igg_dif_cpg_top.png"), width=10, height=14)
