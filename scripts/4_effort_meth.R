### load packages
pacman::p_load(tidyverse, data.table, tibble, performance, matrixStats, 
               parallel, performance, lmerTest, tidystats, insight, effects)

### load data

load(file="results/modeloutput/subset_sites_sig_prepost.RData")

### load phenotypic data

load("data/phenotypes/fulldata_complete_epi_withdates.RData")

### Calculate delta methylation by matching up pre-post ####

delta_meth <- left_join(subset(changing_cpg, prepost == "pre"),
                            subset(changing_cpg, prepost == "post")[,c("chr_pos", "lib_id", "epi_nr", "lib", "methperc", "cov", "id", "year", "fulldate")],
                            by = c("chr_pos", "id", "year"), suffix = c("_pre", "_post"))

delta_meth <- delta_meth %>% dplyr::select(-c(numC, numT, n_sample, prepost))
delta_meth <- delta_meth %>% relocate(c(id, year, born:age), .before=lib_id_pre)

delta_meth <- delta_meth %>% mutate(delta_meth = methperc_post - methperc_pre, .after =born)
delta_meth <- delta_meth %>% mutate(diff_date = fulldate_post - fulldate_pre)
delta_meth$diff_date <- as.numeric(delta_meth$diff_date)

save(delta_meth, file = "results/modeloutput/subset_sites_sig_deltameth.RData")

#Now, we can use this calculation of delta methylation as a predictor of reproductive effort, and loop over CpG sites with this model.

#combine with site and behaviour info

delta_meth <- left_join(delta_meth, unique(all_pheno_epi[,c("id", "year", "site", "Core")], by = c("id", "year")))

#z-transform the traits before the model that subsets IDs and years where there is
#data for that CpG site

effort <- all_pheno_epi %>% dplyr::select(c("id", "year", "attend", "fight", "dist", "MS")) %>% filter(!is.na(attend)) %>% unique()

effort <- subset(effort, id %in% delta_meth$id)

effort$attend_scl <- scale(effort$attend)
effort$fight_scl <- scale(effort$fight)
effort$dist_scl <- scale(effort$dist)

#combine effort data with methylation data

delta_meth <- left_join(delta_meth, effort[,c("id", "year", "attend", "fight", "dist", "attend_scl", "fight_scl", "dist_scl")], by = c("id", "year"))
                                           
delta_meth_ls <- delta_meth %>% group_split(chr_pos)

# function to run the model
function_model_delta <- function(df, parameter){tryCatch({
  chr_pos <- as.character(df[1,1])
  df <- as.data.frame(df)
  df$methperc_pre_scl <- scale(df$methperc_pre)

  formula <- formula(paste0("scale(delta_meth) ~ ", parameter, "_scl + age + scale(methperc_pre) + (1|site/id) "))
  
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
  
  intercept = summary$coefficients["(Intercept)", "Estimate"]

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

  icc_id_site <-icc(model, by_group = TRUE, tolerance = 0)[1,2]
  icc_site <-icc(model, by_group = TRUE, tolerance = 0)[2,2]
  
  return(data.frame(chr_pos=as.factor(chr_pos),
                    parameter = as.factor(parameter),
                    intercept = as.numeric(intercept),
                    icc_id_site = as.numeric(icc_id_site),
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

### centrality
# run model
delta_out_dist <- parallel::mclapply(delta_meth_ls, function_model_delta, parameter="dist",mc.cores=12)

# some have multiple convergence warnings, exclude them, and some do not have enough data for site:id, exclude
errors <- NULL
for (i in 1:length(delta_out_dist)){
  length <- length(delta_out_dist[[i]])
  if(length != 28){
    errors <- c(errors, i)
  }
}

# some have wrong col names
errors_cols <- NULL
names <- names(delta_out_dist[[1]])
for (i in 1:length(delta_out_dist)){
  wrongnames <- names == names(delta_out_dist[[i]])
  if((FALSE %in% wrongnames) == TRUE){
    errors_cols <- c(errors_cols, i)
  }
}

delta_out_dist <- delta_out_dist[-errors]
delta_out_dist <- delta_out_dist[-errors_cols]

# some have wrong col names
errors_cols <- NULL
names <- names(delta_out_dist[[1]])
for (i in 1:length(delta_out_dist)){
  wrongnames <- names == names(delta_out_dist[[i]])
  if((FALSE %in% wrongnames) == TRUE){
    errors_cols <- c(errors_cols, i)
  }
}
delta_out_dist <- delta_out_dist[-errors_cols]

delta_out_dist <- do.call(rbind.data.frame, delta_out_dist) #n = 3399

# convert to numeric
delta_out_dist$parameter_pval <- as.numeric(delta_out_dist$parameter_pval)
delta_out_dist$age_yearling_pval <- as.numeric(delta_out_dist$age_yearling_pval)

# exclude those with overdispersion
delta_out_dist <- subset(delta_out_dist, dispersion.ratio < 1.1 & dispersion.pval > 0.05) # 3355

# FDR correction
delta_out_dist$parameter_qval <- p.adjust(delta_out_dist$parameter_pval, method = "fdr", n = nrow(delta_out_dist))
delta_out_dist$age_yearling_qval <- p.adjust(delta_out_dist$age_yearling_pval, method = "fdr", n = nrow(delta_out_dist))

### attendance
# run model
delta_out_attend <- parallel::mclapply(delta_meth_ls, function_model_delta, parameter="attend",mc.cores=12)

# some have multiple convergence warnings, exclude them, and some do not have enough data for site:id, exclude
errors <- NULL
for (i in 1:length(delta_out_attend)){
  length <- length(delta_out_attend[[i]])
  if(length != 28){
    errors <- c(errors, i)
  }
}

# some have wrong col names
errors_cols <- NULL
names <- names(delta_out_attend[[1]])
for (i in 1:length(delta_out_attend)){
  wrongnames <- names == names(delta_out_attend[[i]])
  if((FALSE %in% wrongnames) == TRUE){
    errors_cols <- c(errors_cols, i)
  }
}

delta_out_attend <- delta_out_attend[-errors]
delta_out_attend <- delta_out_attend[-errors_cols]

# run again?
errors_cols <- NULL
names <- names(delta_out_attend[[1]])
for (i in 1:length(delta_out_attend)){
  wrongnames <- names == names(delta_out_attend[[i]])
  if((FALSE %in% wrongnames) == TRUE){
    errors_cols <- c(errors_cols, i)
  }
}
delta_out_attend <- delta_out_attend[-errors_cols]

delta_out_attend <- do.call(rbind.data.frame, delta_out_attend)

# to numeric
delta_out_attend$parameter_pval <- as.numeric(delta_out_attend$parameter_pval)
delta_out_attend$age_yearling_pval <- as.numeric(delta_out_attend$age_yearling_pval)

# exclude overdispersion
delta_out_attend <- subset(delta_out_attend, isSingular == FALSE & dispersion.ratio < 1.1 & dispersion.pval > 0.05)

# FDR correction
delta_out_attend$parameter_qval <- p.adjust(delta_out_attend$parameter_pval, method = "fdr", n = nrow(delta_out_attend))
delta_out_attend$age_yearling_qval <- p.adjust(delta_out_attend$age_yearling_pval, method = "fdr", n = nrow(delta_out_attend))

### fighting
# run model
delta_out_fight <- parallel::mclapply(delta_meth_ls, function_model_delta, parameter="fight",mc.cores=12)

# some have multiple convergence warnings, exclude them, and some do not have enough data for site:id, exclude
errors <- NULL
for (i in 1:length(delta_out_fight)){
  length <- length(delta_out_fight[[i]])
  if(length != 28){
    errors <- c(errors, i)
  }
}

# some have wrong col names
errors_cols <- NULL
names <- names(delta_out_fight[[1]])
for (i in 1:length(delta_out_fight)){
  wrongnames <- names == names(delta_out_fight[[i]])
  if((FALSE %in% wrongnames) == TRUE){
    errors_cols <- c(errors_cols, i)
  }
}

delta_out_fight <- delta_out_fight[-errors]
delta_out_fight <- delta_out_fight[-errors_cols]

# some have wrong col names
errors_cols <- NULL
names <- names(delta_out_fight[[1]])
for (i in 1:length(delta_out_fight)){
  wrongnames <- names == names(delta_out_fight[[i]])
  if((FALSE %in% wrongnames) == TRUE){
    errors_cols <- c(errors_cols, i)
  }
}
delta_out_fight <- delta_out_fight[-errors_cols]
delta_out_fight <- do.call(rbind.data.frame, delta_out_fight)

# as numeric
delta_out_fight$parameter_pval <- as.numeric(delta_out_fight$parameter_pval)
delta_out_fight$age_yearling_pval <- as.numeric(delta_out_fight$age_yearling_pval)

# exclude overdispersion
delta_out_fight <- subset(delta_out_fight, dispersion.ratio < 1.1 & dispersion.pval > 0.05)

# FDR correction
delta_out_fight$parameter_qval <- p.adjust(delta_out_fight$parameter_pval, method = "fdr", n = nrow(delta_out_fight))
delta_out_fight$age_yearling_qval <- p.adjust(delta_out_fight$age_yearling_pval, method = "fdr", n = nrow(delta_out_fight))

delta_out_all <- rbind(delta_out_dist, delta_out_attend, delta_out_fight)
delta_out_all$chr_pos <- as.factor(delta_out_all$chr_pos)
delta_out_all$parameter <- as.factor(delta_out_all$parameter)
delta_out_all$isSingular <- as.logical(delta_out_all$isSingular)
delta_out_all[c(3:26, 29,30)] <- lapply(delta_out_all[c(3:26, 29,30)], as.numeric)

delta_out_all <- subset(delta_out_all, convergence == "boundary (singular) fit: see help('isSingular')" | is.na(convergence))
save(delta_out_all, file="results/modeloutput/effort_deltameth_modeloutput_filtered.RData")


### volcano plot
source("scripts/plotting_theme.R")

load(file="results/modeloutput/effort_deltameth_modeloutput_filtered.RData")

delta_out_all$parameter <- gsub("dist", "Centrality",delta_out_all$parameter)
delta_out_all$parameter <- gsub("attend", "Attendance", delta_out_all$parameter)
delta_out_all$parameter <- gsub("fight", "Fighting", delta_out_all$parameter)

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
    theme(legend.position="none") -> volcano_effort

ggsave(volcano_effort, file = "plots/model_out/volcano_effort.png", width=8, height=12)    

### significant ones

cpg_sig <- subset(delta_out_all, parameter_qval < 0.05)
cpg_sig_dist <- subset(cpg_sig, parameter == "Centrality") # 7
cpg_sig_attend <- subset(cpg_sig, parameter == "Attendance") #3
cpg_sig_fight <- subset(cpg_sig, parameter == "Fighting") # 35

### plotting

# as separate plots and all together

# dist
list_plot_dist <- list()
for (i in 1:nrow(cpg_sig_dist)){
    ggplot(subset(delta_meth, chr_pos == cpg_sig_dist$chr_pos[i]), aes(x = dist_scl, y = scale(delta_meth))) + 
      geom_point(fill=clrs_hunting[1], size=3) + labs(x = "z-transformed centrality", y = expression("z-tranformed"*delta*" methylation"),
                      title = paste0("Estimate = ", round(cpg_sig_dist$parameter_estimate[i], 2),", 
SE = ", round(cpg_sig_dist$parameter_se[i],2), ", q-value = ", round(cpg_sig_dist$parameter_qval[i], 2))) +
                                        geom_abline(intercept=cpg_sig_dist$intercept[i], slope = cpg_sig_dist$parameter_estimate[i], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
    list_plot_dist[[i]] <- plot   
   # ggsave(plot, file = paste0("plots/model_out/rawdata_plot_dist_cpg_", i, ".png"), width=8, height=8)
    }

cowplot::plot_grid(list_plot_dist[[1]], list_plot_dist[[2]], list_plot_dist[[3]], list_plot_dist[[4]], list_plot_dist[[5]], list_plot_dist[[6]], list_plot_dist[[7]], 
        labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_dist

ggsave(plots_dist, file = paste0("plots/model_out/rawdata_plot_dist.png"), width=14, height=20)

# attend
list_plot_attend <- list()
for (i in 1:nrow(cpg_sig_attend)){
    ggplot(subset(delta_meth, chr_pos == cpg_sig_attend$chr_pos[i]), aes(x = attend_scl, y = scale(delta_meth))) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = "z-transfomred attendance", y = expression("z-transformed"*delta*" methylation"),
                      title = paste0("Estimate = ", round(cpg_sig_attend$parameter_estimate[i], 2),
                                        ", q-value = ", round(cpg_sig_attend$parameter_qval[i], 2))) +
                                         geom_abline(intercept=cpg_sig_attend$intercept[i], slope = cpg_sig_attend$parameter_estimate[i], 
                                          color=clrs_hunting[2], linewidth=1)+ geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
    list_plot_attend[[i]] <- plot   
   # ggsave(plot, file = paste0("plots/model_out/rawdata_plot_attend_cpg_", i, ".png"), width=8, height=8)
}

cowplot::plot_grid(list_plot_attend[[1]], list_plot_attend[[2]], list_plot_attend[[3]], 
        labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_attend

ggsave(plots_attend, file = paste0("plots/model_out/rawdata_plot_attend.png"), width=14, height=14)

#fight
list_plot_fight <- list()
for (i in 1:nrow(cpg_sig_fight)){
    ggplot(subset(delta_meth, chr_pos == cpg_sig_fight$chr_pos[i]), aes(x = fight_scl, y = scale(delta_meth))) + 
    geom_point(fill=clrs_hunting[1], size=3) + labs(x = "z-transformed fighting rate", y = expression("z-transformed"*delta*" methylation"),
                      title = paste0("Estimate = ", round(cpg_sig_fight$parameter_estimate[i], 2),
                                        ", q-value = ", round(cpg_sig_fight$parameter_qval[i], 2))) +
                                         geom_abline(intercept=cpg_sig_fight$intercept[i], slope = cpg_sig_fight$parameter_estimate[i], 
                                          color=clrs_hunting[2], linewidth=1)+
                                        geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
    list_plot_fight[[i]] <- plot   
    #ggsave(plot, file = paste0("plots/model_out/rawdata_plot_fight_cpg_", i, ".png"), width=8, height=8)
}

cowplot::plot_grid(list_plot_fight[[1]], list_plot_fight[[2]], list_plot_fight[[3]], list_plot_fight[[4]], 
                    list_plot_fight[[5]], list_plot_fight[[6]], list_plot_fight[[7]], 
                    list_plot_fight[[8]], list_plot_fight[[9]], list_plot_fight[[10]], 
        labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_fight_a

cowplot::plot_grid(list_plot_fight[[11]], list_plot_fight[[12]], list_plot_fight[[13]], list_plot_fight[[14]], 
                    list_plot_fight[[15]], list_plot_fight[[16]], list_plot_fight[[17]], 
                    list_plot_fight[[18]], list_plot_fight[[19]], list_plot_fight[[20]], 
        labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_fight_b

cowplot::plot_grid(list_plot_fight[[21]], list_plot_fight[[22]], list_plot_fight[[23]], list_plot_fight[[24]], 
                    list_plot_fight[[25]], list_plot_fight[[26]], list_plot_fight[[27]], 
                    list_plot_fight[[28]], list_plot_fight[[29]], list_plot_fight[[30]], 
        labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_fight_c

cowplot::plot_grid(list_plot_fight[[31]], list_plot_fight[[32]], list_plot_fight[[33]], list_plot_fight[[34]], 
                    list_plot_fight[[35]], 
        labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_fight_d

ggsave(plots_fight_a, file = paste0("plots/model_out/rawdata_plot_fight_a.png"), width=14, height=22)
ggsave(plots_fight_b, file = paste0("plots/model_out/rawdata_plot_fight_b.png"), width=14, height=22)
ggsave(plots_fight_c, file = paste0("plots/model_out/rawdata_plot_fight_c.png"), width=14, height=22)
ggsave(plots_fight_d, file = paste0("plots/model_out/rawdata_plot_fight_d.png"), width=14, height=14)
