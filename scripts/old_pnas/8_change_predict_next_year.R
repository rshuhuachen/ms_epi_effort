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

### merge data with 'next year' phenotypes
load(file="data/phenotypes/phenotypes_long_formodels_withsurv.RData")

pheno_next <- pheno_long_models_ly %>% select(c(id, site, lastyear, lifespan, eyec, blue, lyre, attend, fight, dist, surv))
pheno_next <- left_join(pheno_next, unique(delta_meth_raw[,c("id", "born")]), by = "id")
pheno_next$dead <- pheno_next$born + pheno_next$lifespan
pheno_next <- pheno_next %>% select(c(id, lastyear, site, born, lifespan, dead, lastyear, eyec:surv))
pheno_next$year <- pheno_next$lastyear + 1

# merge metadata
delta_meth_nextyear <- left_join(delta_meth_raw[,-6], pheno_next[,c("id", "site", "lastyear", "born", "lifespan", "dead", "year", "surv")], by = (c("id", "year")))
delta_meth_nextyear <- left_join(delta_meth_nextyear, pheno_next[,c("id", "lastyear", "eyec", "blue", "lyre", "attend", "fight", "dist")], by = (c("id", "year" = "lastyear")))

#check how much is missing
delta_meth_nextyear <- delta_meth_nextyear %>% rename(eyec_nextyear = eyec,
                                                      blue_nextyear = blue,
                                                      lyre_nextyear = lyre,
                                                      attend_nextyear = attend,
                                                      fight_nextyear = fight,
                                                      dist_nextyear = dist)

delta_meth_nextyear %>% subset(year!=dead) %>% select(c(id, year, eyec_nextyear:dist_nextyear)) %>% unique() -> pheno_data

write.csv(pheno_data, "data/phenotypes/next_year.csv", quote=F, row.names = F)

delta_meth_nextyear <- delta_meth_nextyear %>%
  group_by(chr_pos, id) %>%
  sample_n(1) %>%
  ungroup()


## only select cpg sites with enough data
delta_meth_nextyear_n <- delta_meth_nextyear %>% group_by(chr_pos) %>% tally()
delta_meth_n_min10 <- subset(delta_meth_nextyear_n, n > 10)

delta_meth_nextyear_sub <- subset(delta_meth_nextyear, chr_pos %in% delta_meth_n_min10$chr_pos) #564

delta_meth_nextyear_ls <- delta_meth_nextyear %>% group_split(chr_pos)

### select one for those repeated samples ####
function_model_nextyear <- function(df, parameter, pre){tryCatch({
  chr_pos <- as.character(df[1,1])
  df <- as.data.frame(df)
  df$methperc_pre_scl <- scale(df$methperc_pre)
  
  if (pre == "control"){
    formula <- formula(paste0(parameter, "_nextyear ~ delta_meth + methperc_pre + (1|site) "))}
  
  if (pre == "no_control"){
    formula <- formula(paste0(parameter, "_nextyear ~ delta_meth + (1|site) "))}
  
  model <- tryCatch(lmerTest::lmer(formula, data=df, REML=FALSE),
                    warning = function(w) {
                      warning_message <<- conditionMessage(w)  # capture the warning message
                    },
                    error = function(w) {
                      error_message <<- errorCondition(w)  # capture the error message
                    })
  
  if (exists("warning_message") == FALSE) {
    warning_message = NA
  }
  
  if (exists("error_message") == TRUE) {
    
    if(grepl("problems: id:site",error_message)){
      if (pre == "control"){
        formula <- formula(paste0(parameter, "_nextyear ~ delta_meth + methperc_pre + (1|site) "))}
      
      if (pre == "no_control"){
        formula <- formula(paste0(parameter, "_nextyear ~ delta_meth + (1|site) "))
        
        
      }
      error_message <-  error_message$message 
      model <- tryCatch(lmerTest::lmer(formula, data=df),
                        warning = function(w) {
                          warning_message <<- conditionMessage(w)  # capture the warning message
                        },
                        error = function(w) {
                          error_message <<- errorCondition(w)  # capture the error message
                        })
    }}
  
  if (grepl("Model failed|unable to evaluate scaled gradient|Model may not have converged*", warning_message) == TRUE) {
    
    if (pre == "control"){
      formula <- formula(paste0(parameter, "_nextyear ~ delta_meth + methperc_pre + (1|site) "))}
    
    if (pre == "no_control"){
      formula <- formula(paste0(parameter, "_nextyear ~ delta_meth + (1|site) "))
      
    }
    
    model <- tryCatch(lmerTest::lmer(formula, data=df),
                      warning = function(w) {
                        warning_message <<- conditionMessage(w)  # capture the warning message
                      },
                      error = function(w) {
                        error_message <<- errorCondition(w)  # capture the error message
                      })
  }
  
  if (exists("error_message") == FALSE) {
    error_message = NA
  }
  
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
  
  if(is.null(summary(model)$optinfo$conv$lme4$messages) == TRUE){
    convergence <- NA
  }
  
  if(!is.null(summary(model)$optinfo$conv$lme4$messages) == TRUE){
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
  
  icc_site <-icc(model, by_group = TRUE, tolerance = 0)[1,2]
  
  return(data.frame(chr_pos=as.factor(chr_pos),
                    parameter = as.factor(parameter),
                    intercept = as.numeric(intercept),
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
                    convergence = convergence,
                    warning = warning_message,
                    error = as.character(error_message)
                    
  ))
}, error = function(e){cat("ERROR :", conditionMessage(e), "\n");print(paste0(chr_pos, " ", conditionMessage(e)))})
}


### lyre
delta_out_next_lyre_raw <- parallel::mclapply(delta_meth_nextyear_ls, function_model_nextyear, pre = "no_control", parameter = "lyre",mc.cores=4)
delta_out_next_lyre <- do.call(rbind.data.frame, delta_out_next_lyre_raw)
delta_out_next_lyre$parameter_qval <- p.adjust(delta_out_next_lyre$parameter_pval, method = "fdr", n = nrow(delta_out_next_lyre))

nrow(subset(delta_out_next_lyre, parameter_pval < 0.05))
nrow(subset(delta_out_next_lyre, parameter_qval < 0.05))

### blue
delta_out_next_blue_raw <- parallel::mclapply(delta_meth_nextyear_ls, function_model_nextyear, pre = "no_control", parameter = "blue",mc.cores=4)
delta_out_next_blue <- do.call(rbind.data.frame, delta_out_next_blue_raw)
delta_out_next_blue$parameter_qval <- p.adjust(delta_out_next_blue$parameter_pval, method = "fdr", n = nrow(delta_out_next_blue))

nrow(subset(delta_out_next_blue, parameter_pval < 0.05))
nrow(subset(delta_out_next_blue, parameter_qval < 0.05))

### eyec
delta_out_next_eyec_raw <- parallel::mclapply(delta_meth_nextyear_ls, function_model_nextyear, pre = "no_control", parameter = "eyec",mc.cores=4)
delta_out_next_eyec <- do.call(rbind.data.frame, delta_out_next_eyec_raw)
delta_out_next_eyec$parameter_qval <- p.adjust(delta_out_next_eyec$parameter_pval, method = "fdr", n = nrow(delta_out_next_eyec))

nrow(subset(delta_out_next_eyec, parameter_pval < 0.05))
nrow(subset(delta_out_next_eyec, parameter_qval < 0.05))

### attend
delta_out_next_attend_raw <- parallel::mclapply(delta_meth_nextyear_ls, function_model_nextyear, pre = "no_control", parameter = "attend",mc.cores=4)
delta_out_next_attend <- do.call(rbind.data.frame, delta_out_next_attend_raw)
delta_out_next_attend$parameter_qval <- p.adjust(delta_out_next_attend$parameter_pval, method = "fdr", n = nrow(delta_out_next_attend))

nrow(subset(delta_out_next_attend, parameter_pval < 0.05))
nrow(subset(delta_out_next_attend, parameter_qval < 0.05))

### fight
delta_out_next_fight_raw <- parallel::mclapply(delta_meth_nextyear_ls, function_model_nextyear, pre = "no_control", parameter = "fight",mc.cores=4)
delta_out_next_fight <- do.call(rbind.data.frame, delta_out_next_fight_raw)
delta_out_next_fight$parameter_qval <- p.adjust(delta_out_next_fight$parameter_pval, method = "fdr", n = nrow(delta_out_next_fight))

nrow(subset(delta_out_next_fight, parameter_pval < 0.05))
nrow(subset(delta_out_next_fight, parameter_qval < 0.05))

### dist
delta_out_next_dist_raw <- parallel::mclapply(delta_meth_nextyear_ls, function_model_nextyear, pre = "no_control", parameter = "dist",mc.cores=4)
delta_out_next_dist <- do.call(rbind.data.frame, delta_out_next_dist_raw)
delta_out_next_dist$parameter_qval <- p.adjust(delta_out_next_dist$parameter_pval, method = "fdr", n = nrow(delta_out_next_dist))

nrow(subset(delta_out_next_dist, parameter_pval < 0.05))
nrow(subset(delta_out_next_dist, parameter_qval < 0.05))

sig_dist <- subset(delta_out_next_dist, parameter_qval < 0.05)

ggplot(subset(delta_meth_nextyear_sub, chr_pos == as.character(sig_dist$chr_pos[1])), aes(x=delta_meth, y = dist_nextyear)) + 
  geom_point() + geom_smooth(method='lm')

ggplot(subset(delta_meth_nextyear_sub, chr_pos == as.character(sig_dist$chr_pos[2])), aes(x=delta_meth, y = dist_nextyear)) + 
  geom_point() + geom_smooth(method='lm')

ggplot(subset(delta_meth_nextyear_sub, chr_pos == as.character(sig_dist$chr_pos[3])), aes(x=delta_meth, y = dist_nextyear)) + 
  geom_point() + geom_smooth(method='lm')

ggplot(subset(delta_meth_nextyear_sub, chr_pos == as.character(sig_dist$chr_pos[4])), aes(x=delta_meth, y = dist_nextyear)) + 
  geom_point() + geom_smooth(method='lm')

ggplot(subset(delta_meth_nextyear_sub, chr_pos == as.character(sig_dist$chr_pos[5])), aes(x=delta_meth, y = dist_nextyear)) + 
  geom_point() + geom_smooth(method='lm')

### need to subset data according to phenotypic data available