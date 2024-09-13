### load packages
pacman::p_load(tidyverse, data.table, tibble, performance, matrixStats, 
               parallel, performance, lmerTest, tidystats, insight)

### load data

load(file = "results/modeloutput/subset_sites_sig_deltameth.RData")
load(file = "results/modeloutput/subset_sites_sig_prepost.RData")
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

  formula <- formula(paste0("scale(delta_meth) ~ ", parameter, "_dif_scl + age + scale(methperc_pre) + (1|id) + (1|site) "))
  
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
  
  intercept = summary$coefficients["(Intercept)","Estimate"]
  #fixed effect
  parameter_estimate <- summary$coefficients[2,1]
  parameter_se <- summary$coefficients[2,2]
  parameter_df <- summary$coefficients[2,3]
  parameter_tval <- summary$coefficients[2,4]
  parameter_pval <- summary$coefficients[2,5]
  
  #age effect
  age_yearling_estimate <- summary$coefficients["ageYearling", "Estimate"]
  age_yearling_se <- summary$coefficients["ageYearling", "Std. Error"]
  age_yearling_df <- summary$coefficients["ageYearling", "df"]
  age_yearling_tval <- summary$coefficients["ageYearling", "t value"]
  age_yearling_pval <- summary$coefficients["ageYearling", "Pr(>|t|)"]
  
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
  

  if(is.null(summary(model)$optinfo$conv$lme4$messages )){
    convergence <- NA
  }

  if(!is.null(summary(model)$optinfo$conv$lme4$messages )){
    convergence <- summary(model)$optinfo$conv$lme4$messages
  }

  icc_id <-icc(model, by_group = TRUE, tolerance = 0)[1,2]
  icc_site <-icc(model, by_group = TRUE, tolerance = 0)[2,2]
  
  return(data.frame(chr_pos=as.factor(chr_pos),
                    parameter = as.factor(parameter),
                    intercept= as.numeric(intercept),
                    icc_id = as.numeric(icc_id),
                    icc_site = as.numeric(icc_site),
                    parameter_estimate = as.numeric(parameter_estimate),
                    parameter_se = as.numeric(parameter_se),
                    parameter_df = as.numeric(parameter_df),
                    parameter_tval = as.numeric(parameter_tval),
                    parameter_pval = as.numeric(parameter_pval),
                    age_yearling_estimate = as.numeric(age_yearling_estimate),
                    age_yearling_se = as.numeric(age_yearling_se),
                    age_yearling_df = as.numeric(age_yearling_df),
                    age_yearling_tval = as.numeric(age_yearling_tval),
                    age_yearling_pval = as.numeric(age_yearling_pval),
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
                    isSingular = as.logical(isSingular),
                    convergence = convergence
                    
  ))
}, error = function(e){cat("ERROR :", conditionMessage(e), "\n");print(chr_pos)})
}

### run the model per trait

### mass
# run model
delta_out_mass <- parallel::mclapply(delta_meth_ls, function_model_delta, parameter="mass",mc.cores=4)
# some have multiple convergence warnings, exclude them, and some do not have enough data for site:id, exclude
errors <- NULL
for (i in 1:length(delta_out_mass)){
  length <- length(delta_out_mass[[i]])
  if(length != 28){
    errors <- c(errors, i)
  }
}

# some have wrong col names
errors_cols <- NULL
names <- names(delta_out_mass[[1]])
for (i in 1:length(delta_out_mass)){
  wrongnames <- names == names(delta_out_mass[[i]])
  if((FALSE %in% wrongnames) == TRUE){
    errors_cols <- c(errors_cols, i)
  }
}

delta_out_mass <- delta_out_mass[-errors]
delta_out_mass <- delta_out_mass[-errors_cols]

# some have wrong col names
errors_cols <- NULL
names <- names(delta_out_mass[[1]])
for (i in 1:length(delta_out_mass)){
  wrongnames <- names == names(delta_out_mass[[i]])
  if((FALSE %in% wrongnames) == TRUE){
    errors_cols <- c(errors_cols, i)
  }
}
delta_out_mass <- delta_out_mass[-errors_cols]

delta_out_mass <- do.call(rbind.data.frame, delta_out_mass)

# convert to numeric
delta_out_mass$parameter_pval <- as.numeric(delta_out_mass$parameter_pval)
delta_out_mass$age_yearling_pval <- as.numeric(delta_out_mass$age_yearling_pval)

# exclude those with overdispersion
delta_out_mass <- subset(delta_out_mass, dispersion.ratio < 1.1 & dispersion.pval > 0.05)

# FDR correction
delta_out_mass$parameter_qval <- p.adjust(delta_out_mass$parameter_pval, method = "fdr", n = nrow(delta_out_mass))
delta_out_mass$age_yearling_qval <- p.adjust(delta_out_mass$age_yearling_pval, method = "fdr", n = nrow(delta_out_mass))

### microf
# run model
delta_out_microf <- parallel::mclapply(delta_meth_ls, function_model_delta, parameter="microf",mc.cores=4)

# some have multiple convergence warnings, exclude them, and some do not have enough data for site:id, exclude
errors <- NULL
for (i in 1:length(delta_out_microf)){
  length <- length(delta_out_microf[[i]])
  if(length != 28){
    errors <- c(errors, i)
  }
}

# some have wrong col names
errors_cols <- NULL
names <- names(delta_out_microf[[4]])
for (i in 1:length(delta_out_microf)){
  wrongnames <- names == names(delta_out_microf[[i]])
  if((FALSE %in% wrongnames) == TRUE){
    errors_cols <- c(errors_cols, i)
  }
}

delta_out_microf <- delta_out_microf[-errors]
delta_out_microf <- delta_out_microf[-errors_cols]

# some have wrong col names
errors_cols <- NULL
names <- names(delta_out_microf[[1]])
for (i in 1:length(delta_out_microf)){
  wrongnames <- names == names(delta_out_microf[[i]])
  if((FALSE %in% wrongnames) == TRUE){
    errors_cols <- c(errors_cols, i)
  }
}
delta_out_microf <- delta_out_microf[-errors_cols]
delta_out_microf <- do.call(rbind.data.frame, delta_out_microf) #1632

# to numeric
delta_out_microf$parameter_pval <- as.numeric(delta_out_microf$parameter_pval)
delta_out_microf$age_yearling_pval <- as.numeric(delta_out_microf$age_yearling_pval)

# exclude overdispersion
delta_out_microf <- subset(delta_out_microf, isSingular == FALSE & dispersion.ratio < 1.1 & dispersion.pval > 0.05)
#445

# FDR correction
delta_out_microf$parameter_qval <- p.adjust(delta_out_microf$parameter_pval, method = "fdr", n = nrow(delta_out_microf))
delta_out_microf$age_yearling_qval <- p.adjust(delta_out_microf$age_yearling_pval, method = "fdr", n = nrow(delta_out_microf))

### trypa
# run model
delta_out_trypa <- parallel::mclapply(delta_meth_ls, function_model_delta, parameter="trypa",mc.cores=4)

# some have multiple convergence warnings, exclude them, and some do not have enough data for site:id, exclude
errors <- NULL
for (i in 1:length(delta_out_trypa)){
  length <- length(delta_out_trypa[[i]])
  if(length != 28){
    errors <- c(errors, i)
  }
}

# some have wrong col names
errors_cols <- NULL
names <- names(delta_out_trypa[[4]])
for (i in 1:length(delta_out_trypa)){
  wrongnames <- names == names(delta_out_trypa[[i]])
  if((FALSE %in% wrongnames) == TRUE){
    errors_cols <- c(errors_cols, i)
  }
}

delta_out_trypa <- delta_out_trypa[-errors]
delta_out_trypa <- delta_out_trypa[-errors_cols]

# some have wrong col names
errors_cols <- NULL
names <- names(delta_out_trypa[[1]])
for (i in 1:length(delta_out_trypa)){
  wrongnames <- names == names(delta_out_trypa[[i]])
  if((FALSE %in% wrongnames) == TRUE){
    errors_cols <- c(errors_cols, i)
  }
}
delta_out_trypa <- delta_out_trypa[-errors_cols]

delta_out_trypa <- do.call(rbind.data.frame, delta_out_trypa)

# as numeric
delta_out_trypa$parameter_pval <- as.numeric(delta_out_trypa$parameter_pval)
delta_out_trypa$age_yearling_pval <- as.numeric(delta_out_trypa$age_yearling_pval)

# exclude overdispersion
delta_out_trypa <- subset(delta_out_trypa, dispersion.ratio < 1.1 & dispersion.pval > 0.05)

# FDR correction
delta_out_trypa$parameter_qval <- p.adjust(delta_out_trypa$parameter_pval, method = "fdr", n = nrow(delta_out_trypa))
delta_out_trypa$age_yearling_qval <- p.adjust(delta_out_trypa$age_yearling_pval, method = "fdr", n = nrow(delta_out_trypa))

### ig
# run model
delta_out_ig <- parallel::mclapply(delta_meth_ls, function_model_delta, parameter="ig",mc.cores=4)

# some have multiple convergence warnings, exclude them, and some do not have enough data for site:id, exclude
errors <- NULL
for (i in 1:length(delta_out_ig)){
  length <- length(delta_out_ig[[i]])
  if(length != 28){
    errors <- c(errors, i)
  }
}

# some have wrong col names
errors_cols <- NULL
names <- names(delta_out_ig[[4]])
for (i in 1:length(delta_out_ig)){
  wrongnames <- names == names(delta_out_ig[[i]])
  if((FALSE %in% wrongnames) == TRUE){
    errors_cols <- c(errors_cols, i)
  }
}

delta_out_ig <- delta_out_ig[-errors]
delta_out_ig <- delta_out_ig[-errors_cols]

# some have wrong col names
errors_cols <- NULL
names <- names(delta_out_ig[[1]])
for (i in 1:length(delta_out_ig)){
  wrongnames <- names == names(delta_out_ig[[i]])
  if((FALSE %in% wrongnames) == TRUE){
    errors_cols <- c(errors_cols, i)
  }
}
delta_out_ig <- delta_out_ig[-errors_cols]

delta_out_ig <- do.call(rbind.data.frame, delta_out_ig)#3422

# as numeric
delta_out_ig$parameter_pval <- as.numeric(delta_out_ig$parameter_pval)
delta_out_ig$age_yearling_pval <- as.numeric(delta_out_ig$age_yearling_pval)

# exclude overdispersion
delta_out_ig <- subset(delta_out_ig, dispersion.ratio < 1.1 & dispersion.pval > 0.05)#3364

# FDR correction
delta_out_ig$parameter_qval <- p.adjust(delta_out_ig$parameter_pval, method = "fdr", n = nrow(delta_out_ig))
delta_out_ig$age_yearling_qval <- p.adjust(delta_out_ig$age_yearling_pval, method = "fdr", n = nrow(delta_out_ig))

### hct
# run model
delta_out_hct <- parallel::mclapply(delta_meth_ls, function_model_delta, parameter="hct",mc.cores=4)
# some have multiple convergence warnings, exclude them, and some do not have enough data for site:id, exclude
errors <- NULL
for (i in 1:length(delta_out_hct)){
  length <- length(delta_out_hct[[i]])
  if(length != 28){
    errors <- c(errors, i)
  }
}

# some have wrong col names
errors_cols <- NULL
names <- names(delta_out_hct[[4]])
for (i in 1:length(delta_out_hct)){
  wrongnames <- names == names(delta_out_hct[[i]])
  if((FALSE %in% wrongnames) == TRUE){
    errors_cols <- c(errors_cols, i)
  }
}

delta_out_hct <- delta_out_hct[-errors]
delta_out_hct <- delta_out_hct[-errors_cols]

# some have wrong col names
errors_cols <- NULL
names <- names(delta_out_hct[[1]])
for (i in 1:length(delta_out_hct)){
  wrongnames <- names == names(delta_out_hct[[i]])
  if((FALSE %in% wrongnames) == TRUE){
    errors_cols <- c(errors_cols, i)
  }
}
delta_out_hct <- delta_out_hct[-errors_cols]

delta_out_hct <- do.call(rbind.data.frame, delta_out_hct) #2797

# as numeric
delta_out_hct$parameter_pval <- as.numeric(delta_out_hct$parameter_pval)
delta_out_hct$age_yearling_pval <- as.numeric(delta_out_hct$age_yearling_pval)

# exclude overdispersion
delta_out_hct <- subset(delta_out_hct, dispersion.ratio < 1.1 & dispersion.pval > 0.05) #2651

# FDR correction
delta_out_hct$parameter_qval <- p.adjust(delta_out_hct$parameter_pval, method = "fdr", n = nrow(delta_out_hct))
delta_out_hct$age_yearling_qval <- p.adjust(delta_out_hct$age_yearling_pval, method = "fdr", n = nrow(delta_out_hct))

### combine
delta_out_all <- rbind(delta_out_mass, delta_out_microf, delta_out_trypa, delta_out_ig, delta_out_hct)
delta_out_all$chr_pos <- as.factor(delta_out_all$chr_pos)
delta_out_all$parameter <- as.factor(delta_out_all$parameter)
delta_out_all$isSingular <- as.logical(delta_out_all$isSingular)

delta_out_all <- subset(delta_out_all, convergence == "boundary (singular) fit: see help('isSingular')" | is.na(convergence))

save(delta_out_all, file="results/modeloutput/physio_deltameth_modeloutput_filtered.RData")

### volcano plot
source("scripts/plotting_theme.R")

load(file="results/modeloutput/physio_deltameth_modeloutput_filtered.RData")

delta_out_all$parameter <- gsub("body mass", "Delta body mass",delta_out_all$parameter)
delta_out_all$parameter <- gsub("microf", "Delta Microfilaria spp.", delta_out_all$parameter)
delta_out_all$parameter <- gsub("trypa", "Delta Trypanosoma spp.", delta_out_all$parameter)
delta_out_all$parameter <- gsub("ig", "Delta IgG", delta_out_all$parameter)
delta_out_all$parameter <- gsub("hct", "Delta HCT", delta_out_all$parameter)

delta_out_all <- delta_out_all %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))

clrs <- viridisLite::viridis(6)
ggplot(delta_out_all, aes(x = parameter_estimate, y = -log10(parameter_qval))) + geom_point(size=4, alpha=0.5, aes(col = sig)) +
    facet_wrap(~parameter, ncol=1, scales="free") +
   # xlim(-1,1)+
   # ylim(0,3)+
    labs(x = "Estimate", y = "-log10(q-value)") +
    scale_color_manual(values=c("grey60", clrs[4])) +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = -0.1, col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0.1, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") -> volcano_physio

ggsave(volcano_physio, file = "plots/model_out/volcano_physio.png", width=8, height=18)    

### significant ones
cpg_sig_mass <- subset(delta_out_all, parameter_qval < 0.05 & parameter == "Delta body mass" & abs(parameter_estimate) > 0.1) #8
cpg_sig_microf <- subset(delta_out_all, parameter_qval < 0.05 & parameter == "Delta Microfilaria spp." & abs(parameter_estimate) > 0.1) #12
cpg_sig_trypa <- subset(delta_out_all, parameter_qval < 0.05 & parameter == "Delta Trypanosoma spp."& abs(parameter_estimate) > 0.1) #21
cpg_sig_ig <- subset(delta_out_all, parameter_qval < 0.05 & parameter == "Delta IgG"& abs(parameter_estimate) > 0.1) #5
cpg_sig_hct <- subset(delta_out_all, parameter_qval < 0.05 & parameter == "Delta HCT"& abs(parameter_estimate) > 0.1) #23

### plotting

source("scripts/plotting_theme.R")

### microf
list_plot_mass <- list()
for (i in 1:nrow(cpg_sig_mass)){
    ggplot(subset(delta_meth, chr_pos == cpg_sig_mass$chr_pos[i]), aes(x = mass_dif_scl, y = scale(delta_meth))) + 
      geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression("z-transformed "*Delta*" body mass"), y = expression("z-transformed "*Delta*" methylation"),
                      title = paste0("Estimate = ", round(cpg_sig_mass$parameter_estimate[i], 2), ", q-value = ", round(cpg_sig_mass$parameter_qval[i], 4))) +
                                        geom_abline(intercept=cpg_sig_mass$intercept[i], slope = cpg_sig_mass$parameter_estimate[i], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
    list_plot_mass[[i]] <- plot   
   }

cowplot::plot_grid(list_plot_mass[[1]], list_plot_mass[[2]], list_plot_mass[[3]], list_plot_mass[[4]], list_plot_mass[[5]], 
list_plot_mass[[6]], list_plot_mass[[7]], list_plot_mass[[8]], 
        labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_mass


ggsave(plots_mass, file = paste0("plots/model_out/rawdata_plot_mass.png"), width=14, height=22)

### microf
list_plot_microf <- list()
for (i in 1:nrow(cpg_sig_microf)){
    ggplot(subset(delta_meth, chr_pos == cpg_sig_microf$chr_pos[i]), aes(x = microf_dif_scl, y = scale(delta_meth))) + 
      geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression("z-transformed "*Delta*" Microfilaria spp."), y = expression("z-transformed "*Delta*" methylation"),
                      title = paste0("Estimate = ", round(cpg_sig_microf$parameter_estimate[i], 2), ", q-value = ", round(cpg_sig_microf$parameter_qval[i], 4))) +
                                        geom_abline(intercept=cpg_sig_microf$intercept[i], slope = cpg_sig_microf$parameter_estimate[i], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
    list_plot_microf[[i]] <- plot   
   }

cowplot::plot_grid(list_plot_microf[[1]], list_plot_microf[[2]], list_plot_microf[[3]], list_plot_microf[[4]], list_plot_microf[[5]], list_plot_microf[[6]], 
        labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_microf_a

cowplot::plot_grid(list_plot_microf[[7]], list_plot_microf[[8]], list_plot_microf[[9]], list_plot_microf[[10]], list_plot_microf[[11]], list_plot_microf[[12]], 
        labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_microf_b

ggsave(plots_microf_a, file = paste0("plots/model_out/rawdata_plot_microf_a.png"), width=14, height=20)
ggsave(plots_microf_b, file = paste0("plots/model_out/rawdata_plot_microf_b.png"), width=14, height=20)


### trypa
cpg_sig_trypa <- cpg_sig_trypa %>% arrange(parameter_qval)

list_plot_trypa <- list()
for (i in 1:nrow(cpg_sig_trypa)){
    ggplot(subset(delta_meth, chr_pos == cpg_sig_trypa$chr_pos[i]), aes(x = trypa_dif_scl, y = scale(delta_meth))) + 
      geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression("z-transformed "*Delta*" Trypanosoma spp."), y = expression("z-transformed "*Delta*" methylation"),
                      title = paste0("Estimate = ", round(cpg_sig_trypa$parameter_estimate[i], 2), ", q-value = ", round(cpg_sig_trypa$parameter_qval[i], 4))) +
                                        geom_abline(intercept=cpg_sig_trypa$intercept[i], slope = cpg_sig_trypa$parameter_estimate[i], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
    list_plot_trypa[[i]] <- plot   
   }

cowplot::plot_grid(list_plot_trypa[[1]], list_plot_trypa[[2]], list_plot_trypa[[3]], list_plot_trypa[[4]], list_plot_trypa[[5]], list_plot_trypa[[6]], 
                    list_plot_trypa[[7]], list_plot_trypa[[8]], list_plot_trypa[[9]], list_plot_trypa[[10]],
        labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_trypa_a

cowplot::plot_grid(list_plot_trypa[[11]], list_plot_trypa[[12]], list_plot_trypa[[13]], list_plot_trypa[[14]], list_plot_trypa[[15]], list_plot_trypa[[16]], 
                    list_plot_trypa[[17]], list_plot_trypa[[18]], list_plot_trypa[[19]], list_plot_trypa[[20]],list_plot_trypa[[21]],
        labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_trypa_b

ggsave(plots_trypa_a, file = paste0("plots/model_out/rawdata_plot_trypa_a.png"), width=14, height=22)
ggsave(plots_trypa_b, file = paste0("plots/model_out/rawdata_plot_trypa_b.png"), width=14, height=22)

### igg
list_plot_igg <- list()
for (i in 1:nrow(cpg_sig_ig)){
    ggplot(subset(delta_meth, chr_pos == cpg_sig_ig$chr_pos[i]), aes(x = ig_dif_scl, y = scale(delta_meth))) + 
      geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression("z-transformed "*Delta*" IgG"), y = expression("z-transformed "*Delta*" methylation"),
                      title = paste0("Estimate = ", round(cpg_sig_ig$parameter_estimate[i], 2), ", q-value = ", round(cpg_sig_ig$parameter_qval[i], 4))) +
                                        geom_abline(intercept=cpg_sig_ig$intercept[i], slope = cpg_sig_ig$parameter_estimate[i], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
    list_plot_igg[[i]] <- plot   
   }

cowplot::plot_grid(list_plot_igg[[1]], list_plot_igg[[2]], list_plot_igg[[3]], list_plot_igg[[4]], list_plot_igg[[5]], 
        labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_igg

ggsave(plots_igg, file = paste0("plots/model_out/rawdata_plot_igg.png"), width=14, height=20)

### hct plot all 3 in one
cpg_sig_hct <- cpg_sig_hct %>% arrange(parameter_qval)

list_plot_hct <- list()
for (i in 1:nrow(cpg_sig_hct)){
    ggplot(subset(delta_meth, chr_pos == cpg_sig_hct$chr_pos[i]), aes(x = trypa_dif_scl, y = scale(delta_meth))) + 
      geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression("z-transformed "*Delta*" HCT"), y = expression("z-transformed "*Delta*" methylation"),
                      title = paste0("Estimate = ", round(cpg_sig_hct$parameter_estimate[i], 2), ", q-value = ", round(cpg_sig_hct$parameter_qval[i], 4))) +
                                        geom_abline(intercept=cpg_sig_hct$intercept[i], slope = cpg_sig_hct$parameter_estimate[i], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
    list_plot_hct[[i]] <- plot   
   }

cowplot::plot_grid(list_plot_hct[[1]], list_plot_hct[[2]], list_plot_hct[[3]], list_plot_hct[[4]], list_plot_hct[[5]], list_plot_hct[[6]], 
                    list_plot_hct[[7]], list_plot_hct[[8]], list_plot_hct[[9]], list_plot_hct[[10]],
        labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_hct_a

cowplot::plot_grid(list_plot_hct[[11]], list_plot_hct[[12]], list_plot_hct[[13]], list_plot_hct[[14]], list_plot_hct[[15]], list_plot_hct[[16]], 
                    list_plot_hct[[17]], list_plot_hct[[18]], list_plot_hct[[19]], list_plot_hct[[20]],list_plot_hct[[21]],
        labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_hct_b

ggsave(plots_hct_a, file = paste0("plots/model_out/rawdata_plot_hct_a.png"), width=14, height=22)
ggsave(plots_hct_b, file = paste0("plots/model_out/rawdata_plot_hct_b.png"), width=14, height=24)
