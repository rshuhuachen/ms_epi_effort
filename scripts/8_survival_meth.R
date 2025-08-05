### load packages
pacman::p_load(dplyr, data.table, tibble, performance, gaston, cowplot,
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


function_model_surv <- function(df){tryCatch({
  chr_pos <- as.character(df[1,1])
  df <- as.data.frame(df)
  df$site <- as.factor(df$site)
  df$id <- as.factor(df$id)
  
  ### surv
  formula_surv <- formula(paste0("surv ~ delta_meth + age + (1|site) "))
  model_surv <- glmmTMB(formula_surv, data=df, family = "binomial", REML=FALSE)
  summary_surv <- summary(model_surv)
  
  intercept_surv <- summary_surv$coefficients$cond["(Intercept)", "Estimate"]
  
  #fixed effect
  parameter_estimate <- summary_surv$coefficients$cond["delta_meth", "Estimate"]
  parameter_se <- summary_surv$coefficients$cond["delta_meth","Std. Error"]
  parameter_zval <- summary_surv$coefficients$cond["delta_meth","z value"]
  parameter_pval <- summary_surv$coefficients$cond["delta_meth", "Pr(>|z|)"]
  
  age_estimate <- summary_surv$coefficients$cond["ageYearling", "Estimate"]
  age_se <- summary_surv$coefficients$cond["ageYearling","Std. Error"]
  age_zval <- summary_surv$coefficients$cond["ageYearling","z value"]
  age_pval <- summary_surv$coefficients$cond["ageYearling", "Pr(>|z|)"]
  
  message <- model_surv$fit$message
  
  dispersion_surv <- overdisp_fun(model_surv)
  
  out <- data.frame(chr_pos=as.factor(chr_pos),
                    intercept_surv = intercept_surv,
                    surv_delta_meth_estimate = as.numeric(parameter_estimate),
                    surv_delta_meth_se = as.numeric(parameter_se),
                    surv_delta_meth_zval = as.numeric(parameter_zval),
                    surv_delta_meth_pval = as.numeric(parameter_pval),
                    surv_age_estimate = as.numeric(age_estimate),
                    surv_age_se = as.numeric(age_se),
                    surv_age_zval = as.numeric(age_zval),
                    surv_age_pval = as.numeric(age_pval),
                    surv_message = message,
                    surv_disp_chi = dispersion_surv[1][[1]],
                    surv_disp_ratio = dispersion_surv[2][[1]],
                    surv_disp_rdf = dispersion_surv[3][[1]],
                    surv_disp_p = dispersion_surv[4][[1]]
  ) 
  
  return(out)
  
}, error = function(e){cat("ERROR :", conditionMessage(e), "\n");print(chr_pos)})
}

### run the model 

### survival
# exclude repeated samples
delta_meth_surv <- delta_meth %>%
  filter(!is.na(delta_meth) & !is.na(surv)) %>% 
  group_by(chr_pos, id) %>%
  sample_n(1) %>%
  ungroup()

## only select cpg sites with enough data
delta_meth_n_surv <- delta_meth_surv %>% group_by(chr_pos) %>% tally()
delta_meth_n_min_surv <- subset(delta_meth_n_surv, n > 20)

delta_meth_sub_surv <- subset(delta_meth_surv, chr_pos %in% delta_meth_n_min_surv$chr_pos) #399
delta_meth_ls_surv <- delta_meth_sub_surv %>% group_split(chr_pos)

# save this data set for later due to randomisation
save(delta_meth_ls_surv, file = "data/processed/delta_meth_ls_surv.RData")

delta_out_surv_raw <- parallel::mclapply(delta_meth_ls_surv, function_model_surv,mc.cores=4)
delta_out_surv <- do.call(rbind.data.frame, delta_out_surv_raw)

save(delta_out_surv, file="results/modeloutput/fitness/out_surv_nopre_raw.RData")

#### Subset models and exclude models that did not converge ####

delta_out_surv <- subset(delta_out_surv, surv_message == "relative convergence (4)")

#### Survival #####
nrow(delta_out_surv) # 584

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
nrow(subset(delta_out_surv, surv_delta_meth_qval < 0.05)) #0

subset(delta_out_surv, surv_delta_meth_pval < 0.05) -> almost_sig_surv
summary(as.numeric(almost_sig_surv$surv_delta_meth_estimate))

load(dup, file = "results/modeloutput/effort/multiple_cpgs_effort_nofdr.RData")

#### Save data ####
delta_out_surv <- delta_out_surv %>% select(c(chr_pos, intercept_surv:surv_disp_p, surv_delta_meth_qval))

save(delta_out_surv, file="results/modeloutput/fitness/out_surv_deltameth_filtered.RData")

### Volcano plot ####
source("scripts/plotting_theme.R")

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
