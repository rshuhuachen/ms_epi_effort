### load packages ###

pacman::p_load(tidyverse, brms, bayesplot, cowplot)

### load posteriors ####
load(file="results/modeloutput/brms_fit_attend.RData") #fit_attend
load(file="results/modeloutput/brms_fit_fight.RData") #fit_fight
load(file="results/modeloutput/brms_fit_dist.RData") #fit_dist
load(file="results/modeloutput/brms_fit_mass.RData") #fit_mass

### model performance / diagnosis ###

list_models <- list(fit_attend = fit_attend,
                    fit_fight = fit_fight,
                    fit_dist = fit_dist,
                    fit_mass = fit_mass)


for (i in 1:length(list_models)){
  pdf(file = paste0("plots/pheno/diagnose_", names(list_models[i]), ".pdf"))
  # extract beta's for some plots
  betas <- subset(variables(list_models[[i]]), grepl("b_", variables(list_models[[i]])) & !grepl("_Intercept", variables(list_models[[i]])))
  #autocor
  print(mcmc_acf(list_models[[i]], lags = 10), pars = c(betas))
  #trace betas
  print(mcmc_trace(list_models[[i]], pars = c(betas)))
    
  #rhat
  print(mcmc_rhat(brms::rhat(list_models[[i]])))
  #neff
  print(mcmc_neff(neff_ratio(list_models[[i]])))
  
  # areas
  print(mcmc_areas(list_models[[i]], pars = c(betas)))
  
  dev.off()
  i = i+1
  }

### plot posteriors ####
### load theme ####
source("scripts/plotting_theme.R")

### load function to plot and save posteriors
source('scripts/plots/function_plot_posterior.R')

### attendance ####
post_attend_mass <- plot_posterior(posteriors = fit_attend, response = "attend", predictor = "scalemass", name = "attend_mass", label = "bottom", minx = -0.3, maxx = 0.3)

post_ms_attend <- plot_posterior(posteriors = fit_attend, response = "MS", predictor = "scaleattend", name = "MS_attend", label = "bottom", minx = -3, maxx = 3)
post_surv_attend <- plot_posterior(posteriors = fit_attend, response = "surv", predictor = "scaleattend", name = "surv_attend", label = "bottom", minx = -1.5, maxx=1.5)

post_ms_mass_a <- plot_posterior(posteriors = fit_attend, response = "MS", predictor = "scalemass", name = "MS_mass_a", label = "bottom", min = -1.5, maxx = 1.5)
post_surv_mass_a <- plot_posterior(posteriors = fit_attend, response = "surv", predictor = "scalemass", name = "surv_mass_a", label = "bottom")

### fight ####
post_fight_mass <- plot_posterior(posteriors = fit_fight, response = "fight", predictor = "scalemass", name = "fight_mass", label = "bottom", min = -0.3, maxx = 0.3)

post_ms_fight <- plot_posterior(posteriors = fit_fight, response = "MS", predictor = "scalefight", name = "MS_fight", label = "bottom", minx = -3, maxx=3)
post_surv_fight <- plot_posterior(posteriors = fit_fight, response = "surv", predictor = "scalefight", name = "surv_fight", label = "bottom", minx = -1.5, maxx=1.5)

post_ms_mass_f <- plot_posterior(posteriors = fit_fight, response = "MS", predictor = "scalemass", name = "MS_mass_f", label = "bottom", min = -1.5, maxx = 1.5)
post_surv_mass_f <- plot_posterior(posteriors = fit_fight, response = "surv", predictor = "scalemass", name = "surv_mass_f", label = "bottom")

### dist ####
post_dist_mass <- plot_posterior(posteriors = fit_dist, response = "dist", predictor = "scalemass", name = "dist_mass", label = "bottom", min = -6, maxx = 1)

post_ms_dist <- plot_posterior(posteriors = fit_dist, response = "MS", predictor = "scaledist", name = "MS_dist", label = "bottom", minx = -3, maxx=3)
post_surv_dist <- plot_posterior(posteriors = fit_dist, response = "surv", predictor = "scaledist", name = "surv_dist", label = "bottom", minx = -1.5, maxx=1.5)

post_ms_mass_d <- plot_posterior(posteriors = fit_dist, response = "MS", predictor = "scalemass", name = "MS_mass_d", label = "bottom", minx = -1.5, maxx=1.5)
post_surv_mass_d <- plot_posterior(posteriors = fit_dist, response = "surv", predictor = "scalemass", name = "surv_mass_d", label = "bottom")

### mass ####
post_mass_dif_attend <- plot_posterior(posteriors = fit_mass, response = "massdif", predictor = "scaleattend", name = "mass_dif_attend", label = "bottom", minx = -90, maxx=30)
post_mass_dif_fight <- plot_posterior(posteriors = fit_mass, response = "massdif", predictor = "scalefight", name = "mass_dif_fight", label = "bottom", minx=-20, maxx=30)
#post_mass_dif_dist <- plot_posterior(posteriors = fit_mass, response = "massdif", predictor = "scaledist", name = "mass_dist", label = "bottom")

post_ms_mass_dif <- plot_posterior(posteriors = fit_mass, response = "MS", predictor = "scalemass_dif", name = "MS_mass_dif", label = "bottom")
post_surv_mass_dif <- plot_posterior(posteriors = fit_mass, response = "surv", predictor = "scalemass_dif", name = "surv_mass_dif", label = "bottom", minx = -1, maxx= 5)

