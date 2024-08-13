
##### This script identifies which CpG sites change from pre- to post-lekking

### Packages ####
#devtools::install_github("mastoffel/rptR", build_vignettes = TRUE)
pacman::p_load(tidyverse, data.table, tibble, performance, matrixStats, 
               parallel, performance, lmerTest, tidystats, insight, rptR)

### Plotting ###
source("scripts/plotting_theme.R")

### Data ####

#### Load data ####
### converted prepost meth data
load(file = "data/processed/methylkit_prepost_long_onlyvar_thres0.3_min_0.5_group.RData") #prepost_long only variation
prepost_long <- prepost_long_clean
rm(prepost_long_clean)

## phenotype data ##
load("data/phenotypes/fulldata_complete_epi_withdates.RData")
prepost <- subset(all_pheno_epi, !is.na(prepost))

rm(all_pheno_epi)

### merge with some metadata

prepost_long <- left_join(prepost_long, prepost[,c("id", "year", "Core", "born", "fulldate")], 
                          by = c("id", "year", "fulldate"))

prepost_long <- prepost_long %>% mutate(age_year = as.factor(case_when(Core == "Core" ~ year - born,
                                                        Core == "No core" ~ NA)),
                                        age = as.factor(case_when(Core == "Core" & (year - born > 1) ~ "Adult",
                                                        Core == "Core" & (year - born == 1) ~ "Yearling",
                                                        Core == "No core" ~ "Adult")))

### convert data to a list, one per CpG site
data <- prepost_long %>% group_split(chr_pos)

### define function to collect overdispersion statistics
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

### build function to run the model

## lmer version
function_model_lmer <- function(df){tryCatch({
  chr_pos <- as.character(df[1,1])
  df <- as.data.frame(df)
  df$prepost <- as.factor(df$prepost)
  df$id <- as.factor(df$id)
  
  # model
  model <- lmerTest::lmer(methperc ~ prepost + (1|id), df)
  
  #fixed effects
  prepost_estimate <- summary(model)$coefficients[2,1]
  prepost_se <- summary(model)$coefficients[2,2]
  prepost_tval <- summary(model)$coefficients[2,4]
  prepost_pval <-  summary(model)$coefficients[2,5]
  
  #random effects 
  id_sd <- attributes(VarCorr(model)$id)$stddev
  id_variance <- data.frame(VarCorr(model), comp="Variance")[1,4]
  
  rsqc <- as.vector(performance::r2(model)$R2_conditional) #fixed plus random effects relative to overall variance
  rsqm <- as.vector(performance::r2(model)$R2_marginal) #fixed effects relative to overall variance
  
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
  icc_id <- icc(model, by_group = TRUE, tolerance = 0)[1,2]
  
  return(data.frame(chr_pos=chr_pos, 
                    prepost_estimate = prepost_estimate,
                    prepost_se = prepost_se,
                    prepost_tval = prepost_tval,
                    prepost_pval = prepost_pval,
                    id_sd = id_sd,
                    id_variance = id_variance,
                    rsqc = rsqc,
                    rsqm = rsqm,
                    dispersion.chisq = dispersion.chisq,
                    dispersion.ratio = dispersion.ratio,
                    dispersion.rdf = dispersion.rdf,
                    dispersion.pval = dispersion.pval,
                    isSingular = isSingular,
                    convergence = convergence,
                    icc_id = icc_id
                    ))
}, error = function(e){cat("ERROR :", conditionMessage(e), "\n");print(chr_pos)})
}

### build function to run the model

## glmer version
function_model_glmer <- function(df){tryCatch({
  chr_pos <- as.character(df[1,1])
  df <- as.data.frame(df)
  df$prepost <- as.factor(df$prepost)
  df$id <- as.factor(df$id)
  
  # model
  model <- lme4::glmer(cbind(numC, numT) ~ prepost + (1|id), family = "binomial", df)
  
  #fixed effects
  prepost_estimate <- summary(model)$coefficients[2,1]
  prepost_se <- summary(model)$coefficients[2,2]
  prepost_zval <- summary(model)$coefficients[2,3]
  prepost_pval <-  summary(model)$coefficients[2,4]
  
  #random effects 
  id_sd <- attributes(VarCorr(model)$id)$stddev
  id_variance <- data.frame(VarCorr(model), comp="Variance")[1,4]
  
  rsqc <- performance::r2(model)$R2_conditional #fixed plus random effects relative to overall variance
  rsqm <- performance::r2(model)$R2_marginal #fixed effects relative to overall variance
  
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
  
  icc_id <- icc(model, by_group = TRUE, tolerance = 0)[1,2]
  
  return(data.frame(chr_pos=chr_pos, 
                    prepost_estimate = prepost_estimate,
                    prepost_se = prepost_se,
                    prepost_zval = prepost_zval,
                    prepost_pval = prepost_pval,
                    id_sd = id_sd,
                    id_variance = id_variance,
                    rsqc = rsqc,
                    rsqm = rsqm,
                    dispersion.chisq = dispersion.chisq,
                    dispersion.ratio = dispersion.ratio,
                    dispersion.rdf = dispersion.rdf,
                    dispersion.pval = dispersion.pval,
                    isSingular = isSingular,
                    convergence = convergence,
                    icc_id = icc_id
                    ))
}, error = function(e){cat("ERROR :", conditionMessage(e), "\n");print(chr_pos)})
}

### run the model in parallel per CpG site (list item)
out_glmer <- parallel::mclapply(data, function_model_glmer, mc.cores=4) #274188

# some have multiple convergence warnings, exclude them
errors <- NULL
for (i in 1:length(out_glmer)){
  length <- length(out_glmer[[i]])
  if(length != 16){
    errors <- c(errors, i)
  }
}

out_glmer <- out_glmer[-errors]

# some have wrong col names
errors_cols <- NULL
names <- names(out_glmer[[1]])
for (i in 1:length(out_glmer)){
  wrongnames <- names == names(out_glmer[[i]])
  if((FALSE %in% wrongnames) == TRUE){
    errors_cols <- c(errors_cols, i)
  }
}
# all are due to large eigenvalues, unindentifiable model

out_glmer <- out_glmer[-errors_cols]

out_glmer <- do.call(rbind.data.frame, out_glmer)
save(out_glmer, file="results/modeloutput/prepost_modeloutput_glmer_min0.75_raw.RData")

### apply a FDR multiple-testing correction
out_glmer <- subset(out_glmer, dispersion.ratio < 1.1 & dispersion.pval > 0.05 & is.na(convergence)) # 168720, 179311 with cov not nT
out_glmer$prepost_qval <- p.adjust(out_glmer$prepost_pval, method = "fdr", n = nrow(out_glmer))

save(out_glmer, file="results/modeloutput/prepost_modeloutput_glmer_min0.75.RData")

## lmer
out_lmer <- parallel::mclapply(data, function_model_lmer, mc.cores=4) #274188
# some have multiple convergence warnings, exclude them
errors <- NULL
for (i in 1:length(out_lmer)){
  length <- length(out_lmer[[i]])
  if(length != 16){
    errors <- c(errors, i)
  }
}

# some have wrong col names
errors_cols <- NULL
names <- names(out_lmer[[1]])
for (i in 1:length(out_lmer)){
  wrongnames <- names == names(out_lmer[[i]])
  if((FALSE %in% wrongnames) == TRUE){
    errors_cols <- c(errors_cols, i)
  }
}
# both are due to large eigenvalues, unindentifiable model

out_lmer <- out_lmer[-errors]
out_lmer <- out_lmer[-errors_cols]
out_lmer <- do.call(rbind.data.frame, out_lmer)
save(out_lmer, file="results/modeloutput/prepost_modeloutput_lmer_min0.75_raw.RData")

### apply a FDR multiple-testing correction
out_lmer <- subset(out_lmer, dispersion.ratio < 1.1 & dispersion.pval > 0.05 & is.na(convergence)) # 151449
out_lmer$prepost_qval <- p.adjust(out_lmer$prepost_pval, method = "fdr", n = nrow(out_lmer))

save(out_lmer, file="results/modeloutput/prepost_modeloutput_lmer_min0.75.RData")

## save the epi data only from cpg's that are sig
sub_glmer_prepost <- subset(out_glmer, prepost_qval < 0.05 & abs(prepost_estimate) >= 0.1)
# N = 3563

sub_lmer_prepost <- subset(out_lmer, prepost_qval < 0.05 & abs(prepost_estimate) >= 0.1)
# N = 199 of which 103 in glmer

changing_cpg <- subset(prepost_long, chr_pos %in% sub_glmer_prepost$chr_pos)
save(changing_cpg, file="results/modeloutput/subset_sites_sig_prepost.RData")

changing_cpg_lmer <- subset(prepost_long, chr_pos %in% sub_lmer_prepost$chr_pos)
save(changing_cpg_lmer, file="results/modeloutput/subset_sites_sig_prepost_lmer.RData")

### volcano plot
source("scripts/plotting_theme.R")

load(file="results/modeloutput/prepost_modeloutput_glmer_min0.75.RData")
out_glmer <- out_glmer %>% mutate(sig = as.factor(case_when(abs(as.numeric(prepost_estimate)) > 0.1 & prepost_qval < 0.05 ~ "sig", TRUE ~ "nonsig")))

clrs <- viridisLite::viridis(6)
ggplot(out_glmer, aes(x = as.numeric(prepost_estimate), y = -log10(as.numeric(prepost_qval)))) + geom_point(size=4, alpha=0.5, aes(col = as.factor(sig))) +
    labs(x = "Estimate time period", y = "-log10(q-value)") +
    #xlim(-1, 1)+
    scale_color_manual(values=c("grey60", clrs[4])) +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = -0.1, col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0.1, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") -> volcano_change

ggsave(volcano_change, file = "plots/model_out/volcano_change.png", width=10, height=10)    


# lmer
out_lmer <- out_lmer %>% mutate(sig = as.factor(case_when(abs(as.numeric(prepost_estimate)) > 0.1 & prepost_qval < 0.05 ~ "sig", TRUE ~ "nonsig")))

ggplot(out_lmer, aes(x = as.numeric(prepost_estimate), y = -log10(as.numeric(prepost_qval)))) + geom_point(size=4, alpha=0.5, aes(col = as.factor(sig))) +
    labs(x = "Estimate time period", y = "-log10(q-value)") +
    #xlim(-1, 1)+
    scale_color_manual(values=c("grey60", clrs[4])) +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = -0.1, col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0.1, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") -> volcano_change_lmer

ggsave(volcano_change_lmer, file = "plots/model_out/volcano_change_lmer.png", width=10, height=10)    


### plot the raw data of the five most sig cpg sites
out_glmer <- out_glmer %>% arrange(prepost_qval)

changing_cpg$prepost <- factor(changing_cpg$prepost, levels = c("pre", "post"))
changing_cpg$prepost <- factor(changing_cpg$prepost, levels = c("pre", "post"), labels = c("Pre-lekking", "Post-lekking"))
changing_cpg$id_year <- as.factor(paste0(changing_cpg$id, "_", changing_cpg$year))

subset(changing_cpg, chr_pos == out_glmer$chr_pos[1]) %>%
  arrange(id, year) %>%
  ggplot(., aes(x = prepost, y = methperc))+
  geom_boxplot(linewidth=1, outlier.shape=NA) + 
  geom_path(aes(group = id_year), alpha = 0.8, col = "grey60", position = position_jitter(width = 0.1, seed = 3922)) +
  geom_point(aes(alpha = 0.8, size=cov), col = clrs[4], position = position_jitter(width = 0.1, seed = 3922)) + 
  labs(x = "Time period", y = "Methylation percentage", subtitle = out_glmer$chr_pos[1]) +
  theme(legend.position="none") -> plot_top_cpg_1


subset(changing_cpg, chr_pos == out_glmer$chr_pos[2]) %>%
  arrange(id, year) %>%
  ggplot(., aes(x = prepost, y = methperc))+
  geom_boxplot(linewidth=1, outlier.shape=NA) + 
  geom_path(aes(group = id_year), alpha = 0.8, col = "grey60", position = position_jitter(width = 0.1, seed = 3922)) +
  geom_point(aes(alpha = 0.8, size=cov), col = clrs[4], position = position_jitter(width = 0.1, seed = 3922)) + 
  labs(x = "Time period", y = "Methylation percentage", subtitle = out_glmer$chr_pos[2]) +
  theme(legend.position="none") -> plot_top_cpg_2

subset(changing_cpg, chr_pos == out_glmer$chr_pos[3]) %>%
  arrange(id, year) %>%
  ggplot(., aes(x = prepost, y = methperc))+
  geom_boxplot(linewidth=1, outlier.shape=NA) + 
  geom_path(aes(group = id_year), alpha = 0.8, col = "grey60", position = position_jitter(width = 0.1, seed = 3922)) +
  geom_point(aes(alpha = 0.8, size=cov), col = clrs[4], position = position_jitter(width = 0.1, seed = 3922)) + 
  labs(x = "Time period", y = "Methylation percentage", subtitle = out_glmer$chr_pos[3]) +
  theme(legend.position="none") -> plot_top_cpg_3

subset(changing_cpg, chr_pos == out_glmer$chr_pos[4]) %>%
  arrange(id, year) %>%
  ggplot(., aes(x = prepost, y = methperc))+
  geom_boxplot(linewidth=1, outlier.shape=NA) + 
  geom_path(aes(group = id_year), alpha = 0.8, col = "grey60", position = position_jitter(width = 0.1, seed = 3922)) +
  geom_point(aes(alpha = 0.8, size=cov), col = clrs[4], position = position_jitter(width = 0.1, seed = 3922)) + 
  labs(x = "Time period", y = "Methylation percentage", subtitle = out_glmer$chr_pos[4]) +
  theme(legend.position="none") -> plot_top_cpg_4

  subset(changing_cpg, chr_pos == out_glmer$chr_pos[5]) %>%
  arrange(id, year) %>%
  ggplot(., aes(x = prepost, y = methperc))+
  geom_boxplot(linewidth=1, outlier.shape=NA) + 
  geom_path(aes(group = id_year), alpha = 0.8, col = "grey60", position = position_jitter(width = 0.1, seed = 3922)) +
  geom_point(aes(alpha = 0.8, size=cov), col = clrs[4], position = position_jitter(width = 0.1, seed = 3922)) + 
  labs(x = "Time period", y = "Methylation percentage", subtitle = out_glmer$chr_pos[5]) +
  theme(legend.position="none") -> plot_top_cpg_5

  subset(changing_cpg, chr_pos == out_glmer$chr_pos[6]) %>%
  arrange(id, year) %>%
  ggplot(., aes(x = prepost, y = methperc))+
  geom_boxplot(linewidth=1, outlier.shape=NA) + 
  geom_path(aes(group = id_year), alpha = 0.8, col = "grey60", position = position_jitter(width = 0.1, seed = 3922)) +
  geom_point(aes(alpha = 0.8, size=cov), col = clrs[4], position = position_jitter(width = 0.1, seed = 3922)) + 
  labs(x = "Time period", y = "Methylation percentage", subtitle = out_glmer$chr_pos[6]) +
  theme(legend.position="none") -> plot_top_cpg_6

cowplot::plot_grid(plot_top_cpg_1, plot_top_cpg_2, plot_top_cpg_3, plot_top_cpg_4, 
                    plot_top_cpg_5, plot_top_cpg_6, labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_change

ggsave(plots_change, file = "plots/model_out/plot_top_change_cpg.png", width=16, height=20)    
