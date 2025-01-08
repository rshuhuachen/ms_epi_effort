#### function to run the models that predict delta_meth based on phenotype ####
function_model_delta_pheno <- function(df, parameter, pre){tryCatch({
  chr_pos <- as.character(df[1,1])
  df <- as.data.frame(df)
  df$methperc_pre_scl <- scale(df$methperc_pre)
  
  if (pre == "control"){
    formula <- formula(paste0("delta_meth ~ ", parameter, "_scl + methperc_pre + (1|site/id) "))}
  
  if (pre == "no_control"){
    formula <- formula(paste0("delta_meth ~ ", parameter, "_scl + (1|site/id) "))}
  
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
  
  if (pre == "control"){
    #premeth effect
    pre_estimate <- summary$coefficients["methperc_pre", "Estimate"]
    pre_se <- summary$coefficients["methperc_pre", "Std. Error"]
    pre_df <- summary$coefficients["methperc_pre", "df"]
    pre_tval <- summary$coefficients["methperc_pre", "t value"]
    pre_pval <- summary$coefficients["methperc_pre", "Pr(>|t|)"]
  }
  
  if (pre == "no_control"){
    pre_estimate <- NA
    pre_se <- NA
    pre_df <- NA
    pre_tval <- NA
    pre_pval <- NA
  }
  
  rsqc <- performance::r2(model)$R2_conditional #fixed plus random parameterects relative to overall variance
  rsqm <- performance::r2(model)$R2_marginal #fixed parameterects relative to overall variance
  
  dispersion.chisq <- overdisp.lmer_fun(model)[1,1]
  dispersion.ratio <- overdisp.lmer_fun(model)[1,2]
  dispersion.rdf <- overdisp.lmer_fun(model)[1,3]
  dispersion.pval <- overdisp.lmer_fun(model)[1,4]
  
  isSingular <- isSingular(model)
  
  if(is.null(summary(model)$optinfo$conv$lme4$messages == TRUE)){
    convergence <- NA
  }
  
  if(!is.null(summary(model)$optinfo$conv$lme4$messages == TRUE)){
    convergence <- summary(model)$optinfo$conv$lme4$messages
    if (length(convergence) == 1){
      convergence <- convergence
    }
    else if (length(convergence) == 2){
      convergence <- paste0(convergence[1], "_", convergence[2])
    }
    else if (length(convergence) == 3){
      convergence <- paste0(convergence[1], "_", convergence[2], "_", convergence[3])
    }
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
}, error = function(e){cat("ERROR :", conditionMessage(e), "\n");print(paste0(chr_pos, " ", conditionMessage(e)))})
}

#### function to process the model output

function_process_model <- function(list, dir_plots, dir_data, name_file, pretty_name, filter_disp){
  pacman::p_load(tidyverse, cowplot)
  
  # some have multiple convergence warnings, exclude them, and some do not have enough data for site:id, exclude
  errors <- NULL
  for (i in 1:length(list)){
    length <- length(list[[i]])
    if(length != 23){
      errors <- c(errors, i)
    }
  }
  
  # some have wrong col names
  errors_cols <- NULL
  names <- names(list[[2]])
  for (i in 1:length(list)){
    wrongnames <- names == names(list[[i]])
    if((FALSE %in% wrongnames) == TRUE){
      errors_cols <- c(errors_cols, i)
    }
  }
  
  if (length(errors) != 0){list <- list[-errors]}
  if (length(errors_cols) != 0){list <- list[-errors_cols]}
  
  # some still have wrong col names
  errors_cols <- NULL
  names <- names(list[[1]])
  for (i in 1:length(list)){
    wrongnames <- names == names(list[[i]])
    if((FALSE %in% wrongnames) == TRUE){
      errors_cols <- c(errors_cols, i)
    }
  }
  if (length(errors_cols) != 0){list <- list[-errors_cols]}
  
  ## unlist the list
  
  data <- do.call(rbind.data.frame, list) 
  
  # convert all to numeric
  nums <- c(3:21)
  data[nums] <- lapply(data[nums], as.numeric)
  
  data$chr_pos <- as.factor(data$chr_pos)
  data$parameter <- as.factor(data$parameter)
  data$isSingular <- as.logical(data$isSingular)
  
  # exclude those with convergence error
  data <- subset(data, !grepl("converge", convergence))
  
  # plot overdispersion raw data
  
  ggplot(data, aes(x = dispersion.ratio)) + geom_histogram() + labs(title = "Histogram dispersion ratio", 
                                                                    subtitle= paste0("Raw model output ", data$parameter[1], "; n = ", nrow(data))) -> hist_disp
  
  ggplot(data, aes(x = parameter_pval))+ geom_histogram() + 
    scale_y_continuous(labels = scales::unit_format(unit = "K", scale = 1e-3)) + 
    labs(title = "Histogram p-values", subtitle= paste0("Raw model output ", data$parameter[1])) -> hist_pval
  
  plot_grid(hist_disp, hist_pval, labs="auto", align="hv", axis="lb", ncol=1, label_fontface = "plain", label_size = 22)-> hists_glmer_raw
  ggsave(hists_glmer_raw, file = paste0(dir_plots, "/hist_raw_disp_pval_", name_file, ".png"), width = 12, height = 12)
  
  if(filter_disp == TRUE){
    # exclude those with overdispersion
    data <- subset(data, dispersion.ratio > as.vector(quantile(data$dispersion.ratio, 0.025, na.rm=T)))
    
    ggplot(data, aes(x = dispersion.ratio)) + geom_histogram() + labs(title = "Histogram dispersion ratio", 
                                                                      subtitle= paste0("Raw model output ", data$parameter[1], "; n = ", nrow(data))) -> hist_disp_filter
    
    ggplot(data, aes(x = parameter_qval))+ geom_histogram() + 
      scale_y_continuous(labels = scales::unit_format(unit = "K", scale = 1e-3)) + 
      labs(title = "Histogram p-values", subtitle= paste0("Raw model output ", data$parameter[1])) -> hist_pval_filter
    
    plot_grid(hist_disp_filter, hist_pval_filter, labs="auto", align="hv", axis="lb", ncol=1, label_fontface = "plain", label_size = 22)-> hists_glmer_filter
    ggsave(hists_glmer_filter, file = paste0(dir_plots, "/hist_filter_disp_pval_", name_file, ".png"), width = 12, height = 12)
  }
  
  # FDR correction
  data$parameter_qval <- p.adjust(data$parameter_pval, method = "fdr", n = nrow(data))
  save(data, file = paste0(dir_data, "/out_", name_file, ".RData"))
  
  # subset significant ones 
  
  sig <- subset(data, parameter_qval < 0.05)
  
  # volcano plot
  
  source("scripts/plotting_theme.R")
  
  data$parameter <- gsub(as.character(data$parameter[1]), pretty_name, data$parameter)
  data <- data %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))
  clrs <- viridisLite::viridis(6)
  
  ggplot(data, aes(x = parameter_estimate, y = -log10(parameter_qval))) + geom_point(size=4, alpha=0.5, aes(col = sig)) +
    facet_wrap(~parameter, ncol=1, scales="free") +
    labs(x = "Estimate", y = "-log10(q-value)", subtitle = paste0("Number of significant sites: ", nrow(sig))) +
    scale_color_manual(values=c("grey60", clrs[4])) +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = -0.1, col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0.1, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") -> volcano
  
  ggsave(volcano, file = paste0(dir_plots, "/volcano_", name_file, ".png"), width=10, height=10)    
  
  out <- list(data = data, volcano = volcano, sig = sig, hists_glmer_raw = hists_glmer_raw)
  
  return(out)
  
}