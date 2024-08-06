
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
load(file = "data/processed/methylkit_prepost_long_onlyvar_min0.75.RData") #prepost_long only variation

## phenotype data ##
load("data/phenotypes/fulldata_complete_epi_withdates.RData")
prepost <- subset(all_pheno_epi, !is.na(prepost))

rm(all_pheno_epi)

### merge with some metadata

prepost_long <- left_join(prepost_long, prepost[,c("id", "prepost", "year", "born", "fulldate")], 
                          by = c("id", "year", "fulldate"))

prepost_long$age <- prepost_long$year - prepost_long$born

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
                    icc_id = icc_id
                    ))
}, error = function(e){cat("ERROR :", conditionMessage(e), "\n");print(chr_pos)})
}

### run the model in parallel per CpG site (list item)
out_glmer <- parallel::mclapply(data, function_model_glmer, mc.cores=4)
out_glmer <- do.call(rbind.data.frame, out_glmer)

### apply a FDR multiple-testing correction
out_glmer$prepost_qval <- p.adjust(out_glmer$prepost_pval, method = "fdr", n = nrow(out_glmer))
out_glmer <- subset(out_glmer, dispersion.ratio < 1.1 & dispersion.pval > 0.05)

save(out_glmer, file="results/modeloutput/prepost_modeloutput_glmer_min0.75.RData")

## save the epi data only from cpg's that are sig
sub_glmer_prepost <- subset(out_glmer, prepost_qval < 0.05 & dispersion.ratio < 1.1 & dispersion.pval > 0.05 & abs(prepost_estimate >= 0.1))
# N = 3,418

changing_cpg <- subset(prepost_long, chr_pos %in% sub_glmer_prepost$chr_pos)
save(changing_cpg, file="results/modeloutput/subset_sites_sig_prepost.RData")

### volcano plot
source("scripts/plotting_theme.R")

load(file="results/modeloutput/prepost_modeloutput_glmer_min0.75.RData")
out_glmer <- out_glmer %>% mutate(sig = as.factor(case_when(abs(as.numeric(prepost_estimate)) > 0.1 & prepost_qval < 0.05 ~ "sig", TRUE ~ "nonsig")))

clrs <- viridisLite::viridis(6)
ggplot(out_glmer, aes(x = as.numeric(prepost_estimate), y = -log10(as.numeric(prepost_qval)))) + geom_point(size=4, alpha=0.5, aes(col = as.factor(sig))) +
    labs(x = "Estimate time period", y = "-log10(q-value)") +
    xlim(-1, 1)+
    scale_color_manual(values=c("grey60", clrs[4])) +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = -0.1, col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0.1, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") -> volcano_change

ggsave(volcano_change, file = "plots/model_out/volcano_change.png", width=10, height=10)    
