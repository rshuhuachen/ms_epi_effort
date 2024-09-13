### load packages
pacman::p_load(tidyverse, data.table, tibble, performance, matrixStats, gaston,
               parallel, performance, lmerTest, tidystats, insight, glmmTMB, DHARMa)

### load data

load(file = "results/modeloutput/subset_sites_sig_prepost.RData")

### load phenotypic data

load("data/phenotypes/fulldata_complete_epi_withdates.RData")
load("data/phenotypes/pheno_dif_prepost.RData")

### load delta meth
load(file = "results/modeloutput/all_sites_deltameth.RData")
delta_meth <- subset(delta_meth, chr_pos %in% changing_cpg$chr_pos)

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

delta_meth <- left_join(delta_meth, unique(pheno_pre[,c("id", "year", "MS", "surv", "site", "attend", "fight", "dist")]), by = c("id", "year"))
delta_meth <- left_join(delta_meth, unique(prepost_dif[,c("id", "year", "mass_dif", "trypa_dif", "ig_dif", "hct_dif", "microf_dif")]), by = c("id", "year"))
       
## first only select cpg sites with enough data
delta_meth_n <- delta_meth %>% group_by(chr_pos) %>% filter(!is.na(delta_meth)) %>% tally()
delta_meth_n_min25 <- subset(delta_meth_n, n > 25)

delta_meth_sub <- subset(delta_meth, chr_pos %in% delta_meth_n_min25$chr_pos)
delta_meth_ls <- delta_meth_sub %>% group_split(chr_pos)

## test for zero inflation
model_zi <- glmmTMB(MS ~ 1 + (1|id), family = "poisson", data = prepost_dif)
DHARMa::testZeroInflation(model_zi) #ns

# function to run the model
function_model_fitness <- function(df){tryCatch({
  chr_pos <- as.character(df[1,1])
  df <- as.data.frame(df)
  df$site <- as.factor(df$site)
  df$id <- as.factor(df$id)
  
  ### AMS
  formula_ams <- formula(paste0("MS ~ delta_meth + (1|id) "))
  model_ams <- glmmTMB(formula_ams, data=df, family = "poisson")
  summary_ams <- summary(model_ams)
  
  intercept_ams <- summary_ams$coefficients$cond["(Intercept)", "Estimate"]

  #fixed effect
  parameter_estimate <- summary_ams$coefficients$cond["delta_meth", "Estimate"]
  parameter_se <- summary_ams$coefficients$cond["delta_meth","Std. Error"]
  parameter_zval <- summary_ams$coefficients$cond["delta_meth","z value"]
  parameter_pval <- summary_ams$coefficients$cond["delta_meth", "Pr(>|z|)"]
  
  #age effect
 # age_estimate <- summary_ams$coefficients$cond["ageYearling", "Estimate"]
#  age_se <- summary_ams$coefficients$cond["ageYearling", "Std. Error"]
#  age_zval <- summary_ams$coefficients$cond["ageYearling", "z value"]
#  age_pval <- summary_ams$coefficients$cond["ageYearling", "Pr(>|z|)"]
  
  #premeth effect
 # pre_estimate <- summary_ams$coefficients$cond["scale(methperc_pre)", "Estimate"]
 # pre_se <- summary_ams$coefficients$cond["scale(methperc_pre)", "Std. Error"]
 # pre_zval <- summary_ams$coefficients$cond["scale(methperc_pre)", "z value"]
 # pre_pval <- summary_ams$coefficients$cond["scale(methperc_pre)", "Pr(>|z|)"]
  
  #rsqc <- performance::r2(model_ams)$R2_conditional #fixed plus random parameterects relative to overall variance
  #rsqm <- performance::r2(model_ams)$R2_marginal #fixed parameterects relative to overall variance
  
  message <- model_ams$fit$message
  
 # icc_id_site <-icc(model_ams, by_group = TRUE, tolerance = 0)[1,2]
#  icc_site <-icc(model_ams, by_group = TRUE, tolerance = 0)[2,2]
  overdisp_fun <- function(model) {
    rdf <- df.residual(model)
    rp <- residuals(model,type="pearson")
    Pearson.chisq <- sum(rp^2)
    prat <- Pearson.chisq/rdf
    pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
    c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
  }
  dispersion_ams <- overdisp_fun(model_ams)
 
  ams <- data.frame(chr_pos=as.factor(chr_pos),
                    intercept_ams = intercept_ams,
                    #ams_icc_id_site = as.numeric(icc_id_site),
                   # ams_icc_site = as.numeric(icc_site),
                    ams_delta_meth_estimate = as.numeric(parameter_estimate),
                    ams_delta_meth_se = as.numeric(parameter_se),
                    ams_delta_meth_zval = as.numeric(parameter_zval),
                    ams_delta_meth_pval = as.numeric(parameter_pval),
                  #  ams_age_yearling_estimate = as.numeric(age_estimate),
                  #  ams_age_yearling_se = as.numeric(age_se),
                  #  ams_age_yearling_zval = as.numeric(age_zval),
                  #  ams_age_yearling_pval = as.numeric(age_pval),
                    #ams_pre_estimate = as.numeric(pre_estimate),
                   # ams_pre_se = as.numeric(pre_se),
                   # ams_pre_zval = as.numeric(pre_zval),
                   # ams_pre_pval = as.numeric(pre_pval),
                   # ams_rsqc = as.numeric(rsqc),
                   # ams_rsqm = as.numeric(rsqm),
                    ams_message = message,
                   ams_disp_chi = dispersion_ams[1][[1]],
                   ams_disp_ratio = dispersion_ams[2][[1]],
                   ams_disp_rdf = dispersion_ams[3][[1]],
                   ams_disp_p = dispersion_ams[4][[1]]
  ) 
  
  ### surv
  formula_surv <- formula(paste0("surv ~ delta_meth + (1|id) "))
  model_surv <- glmmTMB(formula_surv, data=df, family = "binomial")
  summary_surv <- summary(model_surv)
  
  intercept_surv <- summary_surv$coefficients$cond["(Intercept)", "Estimate"]

  #fixed effect
  parameter_estimate <- summary_surv$coefficients$cond["delta_meth", "Estimate"]
  parameter_se <- summary_surv$coefficients$cond["delta_meth","Std. Error"]
  parameter_zval <- summary_surv$coefficients$cond["delta_meth","z value"]
  parameter_pval <- summary_surv$coefficients$cond["delta_meth", "Pr(>|z|)"]
  
  #age effect
 # age_estimate <- summary_surv$coefficients$cond["ageYearling", "Estimate"]
#  age_se <- summary_surv$coefficients$cond["ageYearling", "Std. Error"]
#  age_zval <- summary_surv$coefficients$cond["ageYearling", "z value"]
#  age_pval <- summary_surv$coefficients$cond["ageYearling", "Pr(>|z|)"]
  
  #premeth effect
  #pre_estimate <- summary_surv$coefficients$cond["scale(methperc_pre)", "Estimate"]
  #pre_se <- summary_surv$coefficients$cond["scale(methperc_pre)", "Std. Error"]
  #pre_zval <- summary_surv$coefficients$cond["scale(methperc_pre)", "z value"]
  #pre_pval <- summary_surv$coefficients$cond["scale(methperc_pre)", "Pr(>|z|)"]
  
  #rsqc <- performance::r2(model_surv)$R2_conditional #fixed plus random parameterects relative to overall variance
  #rsqm <- performance::r2(model_surv)$R2_marginal #fixed parameterects relative to overall variance
  
  message <- model_surv$fit$message
  
  dispersion_surv <- overdisp_fun(model_surv)
  
  ##icc_id_site <-icc(model_surv, by_group = TRUE, tolerance = 0)[1,2]
  #icc_site <-icc(model_surv, by_group = TRUE, tolerance = 0)[2,2]
  
  surv <- data.frame(intercept_surv = intercept_surv,
                   # surv_icc_id_site = as.numeric(icc_id_site),
                    #surv_icc_site = as.numeric(icc_site),
                    surv_delta_meth_estimate = as.numeric(parameter_estimate),
                    surv_delta_meth_se = as.numeric(parameter_se),
                    surv_delta_meth_zval = as.numeric(parameter_zval),
                    surv_delta_meth_pval = as.numeric(parameter_pval),
                  #  surv_age_yearling_estimate = as.numeric(age_estimate),
                  #  surv_age_yearling_se = as.numeric(age_se),
                  #  surv_age_yearling_zval = as.numeric(age_zval),
                  #  surv_age_yearling_pval = as.numeric(age_pval),
                  #  surv_pre_estimate = as.numeric(pre_estimate),
                   # surv_pre_se = as.numeric(pre_se),
                   # surv_pre_zval = as.numeric(pre_zval),
                  #  surv_pre_pval = as.numeric(pre_pval),
                   # surv_rsqc = as.numeric(rsqc),
                  #  surv_rsqm = as.numeric(rsqm),
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

ggplot(delta_out_fitness, aes(as.numeric(surv_delta_meth_pval))) + geom_histogram()
ggplot(delta_out_fitness, aes(as.numeric(ams_delta_meth_pval))) + geom_histogram()

ggplot(delta_out_fitness, aes(as.numeric(surv_disp_ratio))) + geom_histogram()
ggplot(delta_out_fitness, aes(as.numeric(ams_disp_ratio))) + geom_histogram()

qqplot.pvalues(delta_out_fitness$surv_delta_meth_pval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95)
qqplot.pvalues(delta_out_fitness$ams_delta_meth_pval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95)

# convert to numeric
delta_out_fitness$ams_delta_meth_estimate <- as.numeric(delta_out_fitness$ams_delta_meth_estimate)
delta_out_fitness$surv_delta_meth_estimate <- as.numeric(delta_out_fitness$surv_delta_meth_estimate)

# exclude those with convergence errors
delta_out_ams <- subset(delta_out_fitness, ams_message == "relative convergence (4)")
delta_out_surv <- subset(delta_out_fitness, surv_message == "relative convergence (4)")

qqplot.pvalues(delta_out_surv$surv_delta_meth_pval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95)
qqplot.pvalues(delta_out_ams$ams_delta_meth_pval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95)

## exclude those with dispersion problems

# how ??
delta_out_ams_dis <- subset(delta_out_ams, ams_disp_ratio < 1.1 & ams_disp_p > 0.05)

qqplot.pvalues(delta_out_ams_dis$ams_delta_meth_pval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95)

# FDR correction
delta_out_ams$ams_delta_meth_qval <- p.adjust(delta_out_ams$ams_delta_meth_pval, method = "fdr", n = nrow(delta_out_ams))
delta_out_surv$surv_delta_meth_qval <- p.adjust(delta_out_surv$surv_delta_meth_pval, method = "fdr", n = nrow(delta_out_surv))

# delta_out_ams$ams_pre_qval <- p.adjust(delta_out_ams$ams_pre_pval, method = "fdr", n = nrow(delta_out_ams))
# delta_out_surv$surv_pre_qval <- p.adjust(delta_out_surv$surv_pre_pval, method = "fdr", n = nrow(delta_out_surv))

delta_out_ams$chr_pos <- as.factor(delta_out_ams$chr_pos)
delta_out_surv$chr_pos <- as.factor(delta_out_surv$chr_pos)

delta_out_ams <- delta_out_ams %>% select(c(chr_pos, intercept_ams:ams_message, ams_delta_meth_qval))
delta_out_surv <- delta_out_surv %>% select(c(chr_pos, intercept_surv:surv_message, surv_delta_meth_qval))

save(delta_out_ams, file="results/modeloutput/AMS_deltameth_modeloutput_filtered.RData")
save(delta_out_surv, file="results/modeloutput/surv_deltameth_modeloutput_filtered.RData")

nrow(subset(delta_out_ams, ams_delta_meth_qval < 0.05)) # n = 45
nrow(subset(delta_out_surv, surv_delta_meth_qval < 0.05)) # n = 0

### volcano plot
source("scripts/plotting_theme.R")

## mean delta methylation %

mean_delta_meth <- delta_meth_sub %>% group_by(chr_pos) %>% summarise_at(vars(delta_meth), funs(mean(., na.rm=TRUE)))
names(mean_delta_meth)[2] <- "mean_delta_meth"

#ams
load(file="results/modeloutput/AMS_deltameth_modeloutput_filtered.RData")

delta_out_ams_data <- left_join(delta_out_ams, mean_delta_meth, by = "chr_pos")
delta_out_ams_data <- delta_out_ams_data %>% mutate(sig = case_when(ams_delta_meth_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))

clrs <- viridisLite::viridis(6)
ggplot(delta_out_ams_data, aes(x = mean_delta_meth, y = -log10(ams_delta_meth_qval))) + geom_point(size=4, alpha=0.5, aes(col = sig)) +
    labs(x = "Estimate delta methylation", y = "-log10(q-value)") +
    #xlim(-1,1)+
    scale_color_manual(values=c("grey60", clrs[4])) +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = -0.1, col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0.1, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") -> volcano_ams

# surv
load(file="results/modeloutput/surv_deltameth_modeloutput_filtered.RData")

delta_out_surv <- delta_out_surv %>% mutate(sig = case_when(surv_delta_meth_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))

ggplot(delta_out_surv, aes(x = surv_delta_meth_estimate, y = -log10(surv_delta_meth_qval))) + geom_point(size=4, alpha=0.5, aes(col = sig)) +
    labs(x = "Estimate delta methylation", y = "-log10(q-value)") +
    xlim(-1,1)+
    scale_color_manual(values=c("grey60", clrs[4])) +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = -0.1, col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0.1, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") -> volcano_surv

cowplot::plot_grid(volcano_ams, volcano_surv, labs=c("a) Annual mating succes", "b) Survival"), 
        align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_fitness

ggsave(plots_fitness, file = "plots/model_out/volcano_fitness.png", width=14, height=10)    

### significant ones

cpg_sig_ams <- subset(delta_out_ams, ams_delta_meth_qval < 0.05) #481
cpg_sig_surv <- subset(delta_out_surv, surv_delta_meth_qval < 0.05) #0

### plotting

source("scripts/plotting_theme.R")
cpg_sig_ams$intercept_ams <- as.numeric(cpg_sig_ams$intercept_ams)
cpg_sig_ams <- cpg_sig_ams %>% arrange(ams_delta_meth_qval)

list_plot_ams <- list()
for (i in 1:21){
    ggplot(subset(delta_meth, chr_pos == cpg_sig_ams$chr_pos[i]), aes(x = scale(delta_meth), y = MS)) + 
      geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression("z-transformed "*Delta*" methylation"), 
      y = "Annual mating success",
      title = paste0("Estimate = ", round(cpg_sig_ams$ams_delta_meth_estimate[i], 2), ", q-value = ", 
      round(cpg_sig_ams$ams_delta_meth_qval[i], 4))) +
       geom_abline(intercept=cpg_sig_ams$intercept_ams[i], slope = cpg_sig_ams$ams_delta_meth_estimate[i], 
       color=clrs_hunting[2], linewidth=1)+
       geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
    list_plot_ams[[i]] <- plot   
   }

cowplot::plot_grid(list_plot_ams[[1]], list_plot_ams[[2]], list_plot_ams[[3]], list_plot_ams[[4]], list_plot_ams[[5]], list_plot_ams[[6]], 
                    list_plot_ams[[7]], list_plot_ams[[8]], list_plot_ams[[9]], list_plot_ams[[10]],
        labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_ams_a

cowplot::plot_grid(list_plot_ams[[11]], list_plot_ams[[12]], list_plot_ams[[13]], list_plot_ams[[14]], list_plot_ams[[15]], list_plot_ams[[16]], 
                    list_plot_ams[[17]], list_plot_ams[[18]], list_plot_ams[[19]], list_plot_ams[[20]],list_plot_ams[[21]],
        labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_ams_b

ggsave(plots_ams_a, file = paste0("plots/model_out/rawdata_plot_AMS_cpg_top10.png"), width=14, height=22)
ggsave(plots_ams_b, file = paste0("plots/model_out/rawdata_plot_AMS_cpg_top20.png"), width=14, height=22)

### post meth on surv ####

# function to run the model
function_model_surv_post <- function(df){tryCatch({
  chr_pos <- as.character(df[1,1])
  df <- as.data.frame(df)
  df$site <- as.factor(df$site)
  df$id <- as.factor(df$id)
  
  ### surv
  formula_surv <- formula(paste0("surv ~ methperc_post + (1|id) "))
  model_surv <- glmmTMB(formula_surv, data=df, family = "binomial")
  summary_surv <- summary(model_surv)
  
  intercept_surv <- summary_surv$coefficients$cond["(Intercept)", "Estimate"]
  
  #fixed effect
  parameter_estimate <- summary_surv$coefficients$cond["methperc_post", "Estimate"]
  parameter_se <- summary_surv$coefficients$cond["methperc_post","Std. Error"]
  parameter_zval <- summary_surv$coefficients$cond["methperc_post","z value"]
  parameter_pval <- summary_surv$coefficients$cond["methperc_post", "Pr(>|z|)"]
  
  message <- model_surv$fit$message
  
  dispersion_surv <- overdisp_fun(model_surv)
  
  surv <- data.frame(chr_pos = chr_pos,
                     intercept_surv = intercept_surv,
                     surv_post_meth_estimate = as.numeric(parameter_estimate),
                     surv_post_meth_se = as.numeric(parameter_se),
                     surv_post_meth_zval = as.numeric(parameter_zval),
                     surv_post_meth_pval = as.numeric(parameter_pval),
                     surv_message = message,
                     surv_disp_chi = dispersion_surv[1][[1]],
                     surv_disp_ratio = dispersion_surv[2][[1]],
                     surv_disp_rdf = dispersion_surv[3][[1]],
                     surv_disp_p = dispersion_surv[4][[1]]
  ) 
  
  return(surv)
  
}, error = function(e){cat("ERROR :", conditionMessage(e), "\n");print(chr_pos)})
}

delta_out_surv_post <- parallel::mclapply(delta_meth_ls, function_model_surv_post,mc.cores=4)
delta_out_surv_post <- do.call(rbind.data.frame, delta_out_surv_post)

ggplot(delta_out_surv_post, aes(as.numeric(surv_post_meth_pval))) + geom_histogram()

ggplot(delta_out_surv_post, aes(as.numeric(surv_disp_ratio))) + geom_histogram()

qqplot.pvalues(delta_out_surv_post$surv_post_meth_pval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95)

# convert to numeric
delta_out_surv_post$surv_post_meth_estimate <- as.numeric(delta_out_surv_post$surv_post_meth_estimate)

# FDR correction
delta_out_surv_post$surv_post_meth_qval <- p.adjust(delta_out_surv_post$surv_post_meth_pval, method = "fdr", n = nrow(delta_out_surv_post))

nrow(subset(delta_out_surv_post, surv_post_meth_pval < 0.05)) # n = 7
nrow(subset(delta_out_surv_post, surv_post_meth_qval < 0.05)) # n = 0

## mean delta methylation %

mean_delta_meth <- delta_meth_sub %>% group_by(chr_pos) %>% summarise_at(vars(delta_meth), funs(mean(., na.rm=TRUE)))
names(mean_delta_meth)[2] <- "mean_delta_meth"

#ams
load(file="results/modeloutput/AMS_deltameth_modeloutput_filtered.RData")

delta_out_surv_post <- left_join(delta_out_surv_post, mean_delta_meth, by = "chr_pos")
delta_out_surv_post <- delta_out_surv_post %>% mutate(sig = case_when(surv_post_meth_pval < 0.05 ~ "sig", TRUE ~ "nonsig"))

sig_surv_post <- subset(delta_out_surv_post, surv_post_meth_pval < 0.05)

ggplot(subset(delta_meth_sub, chr_pos == sig_surv_post$chr_pos[1]), aes(x = methperc_post, y = surv)) + geom_boxplot() + geom_point()
ggplot(subset(delta_meth_sub, chr_pos == sig_surv_post$chr_pos[2]), aes(x = methperc_post, y = surv)) + geom_boxplot() + geom_point()
ggplot(subset(delta_meth_sub, chr_pos == sig_surv_post$chr_pos[3]), aes(x = methperc_post, y = surv)) + geom_boxplot() + geom_point()
ggplot(subset(delta_meth_sub, chr_pos == sig_surv_post$chr_pos[4]), aes(x = methperc_post, y = surv)) + geom_boxplot() + geom_point()
ggplot(subset(delta_meth_sub, chr_pos == sig_surv_post$chr_pos[5]), aes(x = methperc_post, y = surv)) + geom_boxplot() + geom_point()
ggplot(subset(delta_meth_sub, chr_pos == sig_surv_post$chr_pos[6]), aes(x = methperc_post, y = surv)) + geom_boxplot() + geom_point()
ggplot(subset(delta_meth_sub, chr_pos == sig_surv_post$chr_pos[7]), aes(x = methperc_post, y = surv)) + geom_boxplot() + geom_point()

### with traits ####

# function to run the model

# function to run the model
function_model_surv_attend <- function(df){tryCatch({
  chr_pos <- as.character(df[1,1])
  df <- as.data.frame(df)
  df$site <- as.factor(df$site)
  df$id <- as.factor(df$id)
  
  ### surv
  formula_surv <- formula(paste0("surv ~ delta_meth*attend + (1|id) "))
  model_surv <- glmmTMB(formula_surv, data=df, family = "binomial")
  summary_surv <- summary(model_surv)
  
  intercept_surv <- summary_surv$coefficients$cond["(Intercept)", "Estimate"]
  
  #fixed effect
  delta_estimate <- summary_surv$coefficients$cond["delta_meth", "Estimate"]
  delta_se <- summary_surv$coefficients$cond["delta_meth","Std. Error"]
  delta_zval <- summary_surv$coefficients$cond["delta_meth","z value"]
  delta_pval <- summary_surv$coefficients$cond["delta_meth", "Pr(>|z|)"]
  
  attend_estimate <- summary_surv$coefficients$cond["attend", "Estimate"]
  attend_se <- summary_surv$coefficients$cond["attend","Std. Error"]
  attend_zval <- summary_surv$coefficients$cond["attend","z value"]
  attend_pval <- summary_surv$coefficients$cond["attend", "Pr(>|z|)"]
  
  interact_estimate <- summary_surv$coefficients$cond["delta_meth:attend", "Estimate"]
  interact_se <- summary_surv$coefficients$cond["delta_meth:attend","Std. Error"]
  interact_zval <- summary_surv$coefficients$cond["delta_meth:attend","z value"]
  interact_pval <- summary_surv$coefficients$cond["delta_meth:attend", "Pr(>|z|)"]
  
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
                     interact_estimate = interact_estimate,
                     interact_se = interact_se,
                     interact_zval = interact_zval,
                     interact_pval = interact_pval,
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
surv_attend_delta <- do.call(rbind.data.frame, surv_attend_delta)

#hist p-vals
ggplot(surv_attend_delta, aes(as.numeric(delta_pval))) + geom_histogram()
ggplot(surv_attend_delta, aes(as.numeric(attend_pval))) + geom_histogram()
ggplot(surv_attend_delta, aes(as.numeric(interact_pval))) + geom_histogram()

#disp
ggplot(surv_attend_delta, aes(as.numeric(surv_disp_ratio))) + geom_histogram()

#qqplot
qqplot.pvalues(surv_attend_delta$delta_pval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95)
qqplot.pvalues(surv_attend_delta$attend_pval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95)
qqplot.pvalues(surv_attend_delta$interact_pval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95)

# look at warnings
summary(as.factor(surv_attend_delta$surv_message))

# FDR correction
delta_out_surv_post$surv_post_meth_qval <- p.adjust(delta_out_surv_post$surv_post_meth_pval, method = "fdr", n = nrow(delta_out_surv_post))

nrow(subset(delta_out_surv_post, surv_post_meth_pval < 0.05)) # n = 7
nrow(subset(delta_out_surv_post, surv_post_meth_qval < 0.05)) # n = 0

## mean delta methylation %

mean_delta_meth <- delta_meth_sub %>% group_by(chr_pos) %>% summarise_at(vars(delta_meth), funs(mean(., na.rm=TRUE)))
names(mean_delta_meth)[2] <- "mean_delta_meth"

#ams
load(file="results/modeloutput/AMS_deltameth_modeloutput_filtered.RData")

delta_out_surv_post <- left_join(delta_out_surv_post, mean_delta_meth, by = "chr_pos")
delta_out_surv_post <- delta_out_surv_post %>% mutate(sig = case_when(surv_post_meth_pval < 0.05 ~ "sig", TRUE ~ "nonsig"))

sig_surv_post <- subset(delta_out_surv_post, surv_post_meth_pval < 0.05)

ggplot(subset(delta_meth_sub, chr_pos == sig_surv_post$chr_pos[1]), aes(x = methperc_post, y = surv)) + geom_boxplot() + geom_point()
