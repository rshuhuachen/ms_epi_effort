##### Extra ######
#### lmer #####

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
out_lmer_raw <- do.call(rbind.data.frame, out_lmer)
save(out_lmer_raw, file="results/modeloutput/prepost_modeloutput_lmer_min0.75_raw.RData")

# qq plot raw data
png(file = "plots/model_out/qqplot_changing_qqplot_lmer_raw.png", width = 1000, height = 1000)
qqplot.pvalues(out_lmer_raw$prepost_pval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95) 
dev.off()

# quite some overdispersion but not as much

### apply a FDR multiple-testing correction

## two options: dispersion filter by ratio < 1.1 (threshold)
out_lmer_threshold <- subset(out_lmer_raw, dispersion.ratio < 1.1 & dispersion.pval > 0.05 & is.na(convergence)) # 168720, 179311 with cov not nT
out_lmer_threshold$prepost_qval <- p.adjust(out_lmer_threshold$prepost_pval, method = "fdr", n = nrow(out_lmer_threshold))

nrow(out_lmer_threshold)/nrow(out_lmer_raw) # = 0.60
nrow(subset(out_lmer_threshold, prepost_qval < 0.05)) # N = 347

# qq plot
png(file = "plots/model_out/qqplot_changing_qqplot_lmer_threshold.png", width = 1000, height = 1000)
qqplot.pvalues(out_lmer_threshold$prepost_qval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95) 
dev.off()

## second option: within the 90% quantiles
out_lmer_perc <- subset(out_lmer_raw, dispersion.ratio < as.vector(quantile(out_lmer_raw$dispersion.ratio, 0.975)) & dispersion.ratio > as.vector(quantile(out_lmer_raw$dispersion.ratio, 0.025)))
out_lmer_perc$prepost_qval <- p.adjust(out_lmer_perc$prepost_pval, method = "fdr", n = nrow(out_lmer_perc))
nrow(out_lmer_perc)/nrow(out_lmer_raw) # = 0.90 (obvs)
nrow(subset(out_lmer_perc, prepost_qval < 0.05)) # N = 337

png(file = "plots/model_out/qqplot_changing_qqplot_lmer_90percentile.png", width = 1000, height = 1000)
qqplot.pvalues(out_lmer_perc$prepost_qval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95) 
dev.off()



out_lmer <- left_join(out_lmer_threshold, mean_delta_meth, by = "chr_pos")
sub_lmer_prepost <- subset(out_lmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.1)
# N = 204 
nrow(subset(sub_glmer_prepost, chr_pos %in% sub_lmer_prepost$chr_pos)) #105
nrow(subset(sub_lmer_prepost, chr_pos %in% sub_glmer_prepost$chr_pos)) #105
changing_cpg_lmer <- subset(prepost_long, chr_pos %in% sub_lmer_prepost$chr_pos)
save(changing_cpg_lmer, file="results/modeloutput/subset_sites_sig_prepost_lmer.RData")
