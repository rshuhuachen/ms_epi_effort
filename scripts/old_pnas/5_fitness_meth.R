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
       
## test for zero inflation
model_zi <- glmmTMB(MS ~ 1 + (1|site/id), family = "poisson", data = prepost_dif)
DHARMa::testZeroInflation(model_zi) #ns

# function to run the model
function_model_ams <- function(df){tryCatch({
  chr_pos <- as.character(df[1,1])
  df <- as.data.frame(df)
  df$site <- as.factor(df$site)
  df$id <- as.factor(df$id)
  
  ### AMS
  formula_ams <- formula(paste0("MS ~ delta_meth + (1|site) "))
  model_ams <- glmmTMB(formula_ams, data=df, family = "poisson", REML=FALSE)
  summary_ams <- summary(model_ams)
  
  intercept_ams <- summary_ams$coefficients$cond["(Intercept)", "Estimate"]

  #fixed effect
  parameter_estimate <- summary_ams$coefficients$cond["delta_meth", "Estimate"]
  parameter_se <- summary_ams$coefficients$cond["delta_meth","Std. Error"]
  parameter_zval <- summary_ams$coefficients$cond["delta_meth","z value"]
  parameter_pval <- summary_ams$coefficients$cond["delta_meth", "Pr(>|z|)"]
  
  message_ams <- model_ams$fit$message
  dispersion_ams <- overdisp_fun(model_ams)
 
  out <- data.frame(chr_pos=as.factor(chr_pos),
                    intercept_ams = intercept_ams,
                    ams_delta_meth_estimate = as.numeric(parameter_estimate),
                    ams_delta_meth_se = as.numeric(parameter_se),
                    ams_delta_meth_zval = as.numeric(parameter_zval),
                    ams_delta_meth_pval = as.numeric(parameter_pval),
                    ams_message = message_ams,
                    ams_disp_chi = dispersion_ams[1][[1]],
                    ams_disp_ratio = dispersion_ams[2][[1]],
                    ams_disp_rdf = dispersion_ams[3][[1]],
                    ams_disp_p = dispersion_ams[4][[1]]
  ) 
  
  return(out)
  
}, error = function(e){cat("ERROR :", conditionMessage(e), "\n");print(chr_pos)})
}

function_model_surv <- function(df){tryCatch({
  chr_pos <- as.character(df[1,1])
  df <- as.data.frame(df)
  df$site <- as.factor(df$site)
  df$id <- as.factor(df$id)
  
  ### surv
  formula_surv <- formula(paste0("surv ~ delta_meth + (1|site) "))
  model_surv <- glmmTMB(formula_surv, data=df, family = "binomial", REML=FALSE)
  summary_surv <- summary(model_surv)
  
  intercept_surv <- summary_surv$coefficients$cond["(Intercept)", "Estimate"]

  #fixed effect
  parameter_estimate <- summary_surv$coefficients$cond["delta_meth", "Estimate"]
  parameter_se <- summary_surv$coefficients$cond["delta_meth","Std. Error"]
  parameter_zval <- summary_surv$coefficients$cond["delta_meth","z value"]
  parameter_pval <- summary_surv$coefficients$cond["delta_meth", "Pr(>|z|)"]
  
  message <- model_surv$fit$message
  
  dispersion_surv <- overdisp_fun(model_surv)
  
  out <- data.frame(chr_pos=as.factor(chr_pos),
                    intercept_surv = intercept_surv,
                    surv_delta_meth_estimate = as.numeric(parameter_estimate),
                    surv_delta_meth_se = as.numeric(parameter_se),
                    surv_delta_meth_zval = as.numeric(parameter_zval),
                    surv_delta_meth_pval = as.numeric(parameter_pval),
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

### ams
# exclude repeated samples
set.seed(1908)
delta_meth_ams <- delta_meth %>%
    filter(!is.na(delta_meth) & !is.na(MS)) %>% 
    group_by(chr_pos, id) %>%
    sample_n(1) %>%
    ungroup()

## only select cpg sites with enough data
delta_meth_n_ams <- delta_meth_ams %>% group_by(chr_pos) %>% tally()
delta_meth_n_min20_ams <- subset(delta_meth_n_ams, n > 20)

delta_meth_sub_ams <- subset(delta_meth_ams, chr_pos %in% delta_meth_n_min20_ams$chr_pos) #564
delta_meth_ls_ams <- delta_meth_sub_ams %>% group_split(chr_pos)

# save this data set for later due to randomisation
save(delta_meth_ls_ams, file = "data/processed/delta_meth_ls_ams.RData")

delta_out_ams_raw <- parallel::mclapply(delta_meth_ls_ams, function_model_ams,mc.cores=4)
delta_out_ams <- do.call(rbind.data.frame, delta_out_ams_raw)

save(delta_out_ams, file="results/modeloutput/fitness/out_ams_nopre_raw.RData")

### survival
# exclude repeated samples
delta_meth_surv <- delta_meth %>%
  filter(!is.na(delta_meth) & !is.na(surv)) %>% 
  group_by(chr_pos, id) %>%
  sample_n(1) %>%
  ungroup()

## only select cpg sites with enough data
delta_meth_n_surv <- delta_meth_surv %>% group_by(chr_pos) %>% tally()
delta_meth_n_min20_surv <- subset(delta_meth_n_surv, n > 20)

delta_meth_sub_surv <- subset(delta_meth_surv, chr_pos %in% delta_meth_n_min20_surv$chr_pos) #607
delta_meth_ls_surv <- delta_meth_sub_surv %>% group_split(chr_pos)

# save this data set for later due to randomisation
save(delta_meth_ls_surv, file = "data/processed/delta_meth_ls_surv.RData")

delta_out_surv_raw <- parallel::mclapply(delta_meth_ls_surv, function_model_surv,mc.cores=4)
delta_out_surv <- do.call(rbind.data.frame, delta_out_surv_raw)

save(delta_out_surv, file="results/modeloutput/fitness/out_surv_nopre_raw.RData")

#### Subset models and exclude models that did not converge ####

delta_out_ams <- subset(delta_out_ams, ams_message == "relative convergence (4)")
delta_out_surv <- subset(delta_out_surv, surv_message == "relative convergence (4)")

#### AMS ####
nrow(delta_out_ams)  # retain 511

##### Overdispersion raw data #####

summary(delta_out_ams$ams_disp_ratio)

# histogram dispersion ratio raw 
ggplot(delta_out_ams, aes(x = ams_disp_ratio)) + geom_histogram() + geom_vline(xintercept = 1., col = "red", linetype = "dotted", linewidth=1) +
scale_y_log10() + labs(title = "Histogram dispersion ratio", subtitle= "Raw model output AMS") -> hist_ams_disp_ratio_raw

# histogram p-values raw
ggplot(delta_out_ams, aes(x = ams_delta_meth_pval)) + geom_histogram() + 
  scale_y_continuous(labels = scales::unit_format(unit = "K", scale = 1e-3)) + 
  labs(title = "Histogram p-values", subtitle="Raw model output AMS")-> hist_ams_pvals_raw

plot_grid(hist_ams_disp_ratio_raw, hist_ams_pvals_raw, labs="auto", align="hv", axis="lb", ncol=1, label_fontface = "plain", label_size = 22)-> hists_ams_raw

ggsave(hists_ams_raw, file = "plots/model_out/fitness/ams/hist_raw_ams.png", width = 12, height = 12)

# qqplot raw

png(file = "plots/model_out/fitness/ams/qqplot_raw_ams.png", width = 800, height = 800)
qqplot.pvalues(delta_out_ams$ams_delta_meth_pval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95) 
dev.off()

## filter for 95 percentile
delta_out_ams_clean <- subset(delta_out_ams, ams_disp_ratio < as.vector(quantile(delta_out_ams$ams_disp_ratio, 0.95)))
nrow(delta_out_ams_clean) # 485

# histogram dispersion ratio filtered 
ggplot(delta_out_ams_clean, aes(x = ams_disp_ratio)) + geom_histogram() + 
scale_y_log10() + labs(title = "Histogram dispersion ratio", subtitle= "Filtered model output AMS") -> hist_ams_disp_ratio_filter

# histogram p-values filtered
ggplot(delta_out_ams_clean, aes(x = ams_delta_meth_pval)) + geom_histogram() + 
  scale_y_continuous(labels = scales::unit_format(unit = "K", scale = 1e-3)) + 
  labs(title = "Histogram p-values", subtitle="Filtered model output AMS")-> hist_ams_pvals_filter

plot_grid(hist_ams_disp_ratio_filter, hist_ams_pvals_filter, labs="auto", align="hv", axis="lb", ncol=1, label_fontface = "plain", label_size = 22)-> hists_ams_filter

ggsave(hists_ams_filter, file = "plots/model_out/fitness/ams/hist_filtered_ams.png", width = 12, height = 12)

# qqplot filtered

png(file = "plots/model_out/fitness/ams/qqplot_filtered_ams.png", width = 800, height = 800)
qqplot.pvalues(delta_out_ams_clean$ams_delta_meth_pval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95) 
dev.off()

#### FDR-correction ####

delta_out_ams_clean$ams_delta_meth_qval <- p.adjust(delta_out_ams_clean$ams_delta_meth_pval, method = "fdr", n = nrow(delta_out_ams_clean))

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
nrow(subset(delta_out_ams_clean, ams_delta_meth_qval < 0.05)) #203
nrow(subset(delta_out_surv, surv_delta_meth_qval < 0.05)) #0

#### Save data ####
delta_out_ams <- delta_out_ams_clean %>% select(c(chr_pos, intercept_ams:ams_disp_p, ams_delta_meth_qval))
delta_out_surv <- delta_out_surv %>% select(c(chr_pos, intercept_surv:surv_disp_p, surv_delta_meth_qval))

save(delta_out_ams, file="results/modeloutput/fitness/out_ams_deltameth_filtered.RData")
save(delta_out_surv, file="results/modeloutput/fitness/out_surv_deltameth_filtered.RData")

### Volcano plot ####
source("scripts/plotting_theme.R")

#ams
load(file="results/modeloutput/fitness/out_ams_deltameth_filtered.RData")

delta_out_ams <- delta_out_ams %>% mutate(sig = case_when(ams_delta_meth_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))

#clrs <- viridisLite::viridis(6)
ggplot(delta_out_ams, aes(x = ams_delta_meth_estimate, y = -log10(ams_delta_meth_qval))) + geom_point(size=4, alpha=0.5, aes(col = sig)) +
    labs(x = expression("Estimate "*Delta*" methylation"), y = "-log10(q-value)") +
    xlim(-21,21)+
    scale_color_manual(values=c(clrs[5], clrs[17])) +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
  #  geom_vline(xintercept = -0.1, col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") -> volcano_ams

ggsave(volcano_ams, file = "plots/model_out/fitness/ams/volcano_ams.png", width=10, height=10)    

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

### significant ones

cpg_sig_ams <- subset(delta_out_ams, ams_delta_meth_qval < 0.05) #203

### plotting

source("scripts/plotting_theme.R")
cpg_sig_ams$intercept_ams <- as.numeric(cpg_sig_ams$intercept_ams)
cpg_sig_ams <- cpg_sig_ams %>% arrange(ams_delta_meth_qval)

for (i in c(1:21)){
    ggplot(subset(delta_meth, chr_pos == cpg_sig_ams$chr_pos[i]), aes(x = delta_meth, y = MS)) + 
      geom_point(fill=clrs_hunting[1], size=3) + labs(x = expression(Delta*" methylation"), 
      y = "Annual mating success",
      title = paste0("Estimate = ", round(cpg_sig_ams$ams_delta_meth_estimate[i], 2), ", q-value = ", 
      round(cpg_sig_ams$ams_delta_meth_qval[i], 4))) +
       geom_abline(intercept=cpg_sig_ams$intercept_ams[i], slope = cpg_sig_ams$ams_delta_meth_estimate[i], 
       color=clrs_hunting[2], linewidth=1)+
       geom_hline(yintercept=0, color=clrs_hunting[3], linetype="dotted", linewidth =1)-> plot
    ggsave(plot, file = paste0("plots/model_out/fitness/ams/rawdata_cpg_", i, ".png"), width=10, height=10)     
   }

#### Summarise if it's more optimal to change or remain stable, and whether to increase or decrease ####
delta_meth_ls_ams_unlist <- do.call(rbind.data.frame, delta_meth_ls_ams)

categories <- data.frame()
plots_raw <- list()
plots_raw_abs <- list()

for (i in 1:nrow(cpg_sig_ams)){
  subset <- subset(delta_meth_ls_ams_unlist, chr_pos == cpg_sig_ams$chr_pos[i])
  subset$delta_meth_abs <- abs(subset$delta_meth)
  
  model <- glmmTMB(MS ~ delta_meth + (1|site), family = "poisson", REML=FALSE, data=subset)
  
  subset$predict <- predict(model)
  
  plot <- ggplot(subset, aes(x = delta_meth, y = MS)) + geom_point() + 
    geom_abline(intercept=cpg_sig_ams$intercept_ams[i], slope = cpg_sig_ams$ams_delta_meth_estimate[i]) + 
    geom_hline(yintercept=0, col = "red", linetype="dotted")
  
  plot_abs <- ggplot(subset, aes(x = delta_meth_abs, y = MS)) + geom_point() + 
    geom_smooth(method='lm') + geom_hline(yintercept=0, col = "red", linetype="dotted")
  
  predicted_nochange <- subset %>% arrange(delta_meth_abs) %>% head(1) %>% select(predict)
  predicted_change <- subset %>% arrange(delta_meth_abs) %>% tail(1) %>% select(predict)
  
  change <- case_when(predicted_change > predicted_nochange ~ "change",
                      predicted_nochange > predicted_change ~ "don't change")
  
  if(change == "change"){
    predicted_decrease <- subset %>% arrange(delta_meth) %>% head(1) %>% select(predict)
    predicted_increase <- subset %>% arrange(delta_meth) %>% tail(1) %>% select(predict)
    
    direction <- case_when(predicted_increase > predicted_decrease ~ "increase",
                           predicted_decrease > predicted_increase ~ "decrease")
  }
  if(change == "don't change"){
    direction <- NA
  }
  
  ### collect all outputs for a list
  
  # category
  category <- data.frame(chr_pos = cpg_sig_ams$chr_pos[i],
                         category = change,
                         direction = direction)
  categories <- rbind(categories, category)
  
  # plot
  plots_raw[[i]] <- plot
  plots_raw_abs[[i]] <- plot_abs
  
  category_optimal_list <- list(categories = categories, 
                                plots_raw = plots_raw, 
                                plots_raw_abs = plots_raw_abs)
  
}

# count different categories of CpG sites
n_changing <- length(unique(delta_meth$chr_pos)) # 1026, changing sites
n_data <- length(delta_meth_ls_ams) # 564, sites with enough data
nodata = n_changing - n_data # changing sites - sites that have enough data = no data
nomodel <- n_data - nrow(delta_out_ams_clean) # with data - converged = not converged
nopredict <- nrow(delta_out_ams_clean) - nrow(cpg_sig_ams) # converged - sig = not sig
dontchange <- nrow(subset(category_optimal_list$categories, category == "don't change"))
inc <- nrow(subset(category_optimal_list$categories, direction == "increase"))
dec <- nrow(subset(category_optimal_list$categories, direction == "decrease"))

summary_cpgs <- data.frame(what = c("Not enough data", "Convergence problems", 
                                    "Not predictive of MS", "No change predicts higher MS",
                                    "Increase predicts higher MS", "Decrease predicts higher MS"),
                           n = c(nodata, nomodel, nopredict, dontchange, inc, dec))

summary_cpgs$what <- factor(summary_cpgs$what, 
                            levels = c("Not enough data", "Convergence problems", 
                                       "Not predictive of MS", "No change predicts higher MS",
                                       "Increase predicts higher MS", "Decrease predicts higher MS"))

# prepare data for a pie chart

summary_cpgs <- summary_cpgs %>% 
  arrange(desc(what)) %>%
  mutate(prop = n / sum(summary_cpgs$n) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

library(ggrepel)

save(summary_cpgs, file = "results/modeloutput/fitness/cpg_categories_for_ams.RData")

# plot
ggplot(summary_cpgs, aes(x="", y=prop, fill=what)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values=c("grey90", "grey70", "#536B74", clrs[3], clrs[6], clrs[2])) +
  geom_label_repel(aes(y = ypos, label = paste0("N = ", n), size = 10), nudge_x = 0.7, show.legend=F)+
  theme_void() + 
  labs(fill = "CpG site category")+
  theme(title = element_text(size=20),
        plot.title = element_text(hjust = 0.5, margin=margin(0,0,15,0)),
        plot.subtitle = element_text(size=16, family = "Arial"),
        text=element_text(size=18, family = "Arial"),
        legend.text =  element_text(size = 18, family = "Arial"),
        legend.title = element_text(size = 18, family = "Arial"),
        strip.text = element_text(size = 18, family = "Arial"),
        plot.margin = margin(1,1,1,1, "cm"),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA)) -> pie_cpg_cat
  
ggsave(pie_cpg_cat, file = "plots/model_out/fitness/ams/pie_cpg_cat.png", width=10, height=10)    

#### Annotate AMS CpG sites ####

### Packages ####
pacman::p_load(genomation, GenomicFeatures, rtracklayer, 
               GenomicRanges)


### Combine all sites vs changing sites
cpg_all <- delta_out_ams %>% dplyr::select(c(chr_pos, ams_delta_meth_qval))
names(cpg_all)[2] <- "parameter_qval"
cpg_all$parameter <- "all"

cpg_ams_select <- cpg_sig_ams %>% dplyr::select(c(chr_pos, ams_delta_meth_qval))
names(cpg_ams_select)[2] <- "parameter_qval"
cpg_ams_select$parameter <- "ams"

all_models_sig <- rbind(cpg_all, cpg_ams_select)

### Rename chr_pos and divide ###
all_models_sig$chr_pos <- gsub("__", ";", all_models_sig$chr_pos)
all_models_sig$chr_pos <- gsub("HRSCAF_", "HRSCAF=", all_models_sig$chr_pos, )

# Extract the numbers following HRSCAF=XXX_number
# Split the chr_pos column into two columns based on the first "_"
split_chr_pos <- strsplit(all_models_sig$chr_pos, "_", fixed = TRUE)

all_models_sig$chr <- paste0(sapply(split_chr_pos, "[", 1), "_",
                             sapply(split_chr_pos, "[", 2))

all_models_sig$pos <- sapply(split_chr_pos, "[", 3)

all_models_sig <- all_models_sig %>% 
  relocate(chr, .after = chr_pos) %>%
  relocate(pos, .after = chr_pos)

#revert scafnames
all_models_sig$chr_pos <- gsub(";","__", all_models_sig$chr_pos)
all_models_sig$chr_pos <- gsub("HRSCAF=", "HRSCAF_", all_models_sig$chr_pos)

all_models_sig$chr <- gsub(";","__", all_models_sig$chr)
all_models_sig$chr <- gsub("HRSCAF=", "HRSCAF_", all_models_sig$chr)

### Load annotation data
annotation_dir <- "~/PhD_grouse/grouse-annotation/output"

promoter=unique(gffToGRanges(paste0(annotation_dir, "/promoters.gff3")))
genes=unique(gffToGRanges(paste0(annotation_dir, "/genes.gff3")))
TSS=unique(gffToGRanges(paste0(annotation_dir, "/TSS.gff3")))
exons_gene=unique(gffToGRanges(paste0(annotation_dir, "/exons_gene.gff3")))
introns=unique(gffToGRanges(paste0(annotation_dir, "/introns_transcripts.gff3")))
downstream=unique(gffToGRanges(paste0(annotation_dir, "/downstream.gff3")))
upstream=unique(gffToGRanges(paste0(annotation_dir, "/upstream.gff3")))
threeUTR =unique(gffToGRanges(paste0(annotation_dir, "/threeUTRs.gff3")))
fiveUTR=unique(gffToGRanges(paste0(annotation_dir, "/fiveUTRs.gff3")))

#### Annotate ####
all_models_sig$end <- all_models_sig$pos
all_models_sig$start <- all_models_sig$pos
sig_gr <- as(all_models_sig, "GRanges")

sig_promoter <- mergeByOverlaps(promoter, sig_gr)
sig_promoter <- as.data.frame(sig_promoter@listData)
sig_promoter <- sig_promoter %>% add_column("region" = "promoter", .after="parameter") %>% 
  mutate(distance = as.numeric(promoter.end) - as.numeric(pos) ) %>% 
  dplyr::select(c(chr_pos, pos, parameter, region, parameter_qval,  ID, distance)) 

sig_gene <- mergeByOverlaps(sig_gr, genes)
sig_gene <- as.data.frame(sig_gene@listData)
sig_gene <- sig_gene %>% add_column("region" = "gene", .after="parameter") %>% 
  mutate(distance = as.numeric(pos) - as.numeric(genes.start)) %>% 
  dplyr::select(c(chr_pos, pos, parameter, region, parameter_qval,  ID, distance)) 

sig_tss <- mergeByOverlaps(sig_gr, TSS)
sig_tss <- as.data.frame(sig_tss@listData)
sig_tss <- sig_tss %>% add_column("region" = "TSS", .after="parameter") %>% 
  mutate(distance = as.numeric(TSS.end) - as.numeric(pos) ) %>% 
  dplyr::select(c(chr_pos, pos, parameter, region, parameter_qval,  ID, distance)) 

sig_exon <- mergeByOverlaps(sig_gr, exons_gene)
sig_exon <- as.data.frame(sig_exon@listData)
sig_exon <- sig_exon %>% add_column("region" = "exon", .after="parameter") %>% 
  mutate(distance = as.numeric(pos) - as.numeric(exons_gene.start)) %>% 
  dplyr::select(c(chr_pos, pos, parameter, region, parameter_qval,  ID, distance)) 

sig_intron <- mergeByOverlaps(sig_gr, introns)
sig_intron <- as.data.frame(sig_intron@listData)
sig_intron <- sig_intron %>% add_column("region" = "intron", .after="parameter") %>% 
  mutate(distance = as.numeric(pos) - as.numeric(introns.start)) %>% 
  dplyr::select(c(chr_pos, pos, parameter, region, parameter_qval,  ID, distance)) 

sig_down <- mergeByOverlaps(downstream, sig_gr)
sig_down <- as.data.frame(sig_down@listData)
sig_down <- sig_down %>% add_column("region" = "downstream", .after="parameter") %>% 
  mutate(distance = as.numeric(pos) - as.numeric(downstream.start)) %>% 
  dplyr::select(c(chr_pos, pos, parameter, region, parameter_qval, ID, distance)) 

sig_up <- mergeByOverlaps(upstream, sig_gr)
sig_up <- as.data.frame(sig_up@listData)
sig_up <- sig_up %>% add_column("region" = "upstream", .after="parameter") %>% 
  mutate(distance =  as.numeric(upstream.end) - as.numeric(pos)) %>% 
  dplyr::select(c(chr_pos, pos, parameter, region, parameter_qval,  ID, distance)) 
# 
# sig_threeUTR <- mergeByOverlaps(sig_gr, threeUTR)
# sig_threeUTR <- as.data.frame(sig_threeUTR@listData)
# sig_threeUTR <- sig_threeUTR %>% add_column("region" = "threeUTR", .after="parameter") %>% 
#   mutate(distance =  as.numeric(threeUTR.start) - as.numeric(pos)) %>% 
#   dplyr::select(c(chr_pos, pos, parameter, region, parameter_qval,  ID, distance)) 
# 
# sig_fiveUTR <- mergeByOverlaps(sig_gr, fiveUTR)
# sig_fiveUTR <- as.data.frame(sig_fiveUTR@listData)
# sig_fiveUTR <- sig_fiveUTR %>% add_column("region" = "fiveUTR", .after="parameter") %>% 
#   mutate(distance =  as.numeric(pos) - as.numeric(fiveUTR.start)) %>%
#   dplyr::select(c(chr_pos, pos, parameter, region, parameter_qval,  ID, distance)) 

all_models_sig_annotated_raw <- rbind(sig_promoter, sig_gene,
                                      sig_tss, sig_exon, sig_intron, sig_down,
                                      sig_up)


summary(as.factor(all_models_sig_annotated_raw$region))

save(all_models_sig_annotated_raw, file="results/modeloutput/fitness/annotated_sig_cpg_ams_raw.RData")

#### Priority workflow: TSS > promoter > exon/intron > down/up ####

prioritize_region <- function(region) {
  # Create a lookup table for region priorities
  priority_table <- c(TSS = 1, promoter = 2, exon = 3, intron = 3, downstream = 4, upstream = 4, gene = 5)
  
  # Get the priority for the given region
  priority <- priority_table[region]
  
  # Handle missing regions (assign a low priority)
  priority[is.na(priority)] <- 7
  
  priority
}


all_models_sig_annotated <- all_models_sig_annotated_raw %>%
  group_by(chr_pos, parameter) %>% #doesn't really matter to group it by gene id, same priority applies
  mutate(region_priority = prioritize_region(region)) %>%
  # Select the row with the highest priority (lowest numeric value)
  slice_min(region_priority, with_ties=T) %>%
  # If multiple rows have the same priority (the case for up/downstream), select the one with the lowest distance
  slice_min(distance, with_ties=T) %>%
  ungroup() %>%
  dplyr::select(-region_priority)

save(all_models_sig_annotated, file="results/modeloutput/fitness/annotated_sig_cpg_ams.RData")

#### Summarise number of sites per region ####
sum_annotated <- as.data.frame(table(as.factor(all_models_sig_annotated$region), all_models_sig_annotated$parameter))
names(sum_annotated) <- c("region", "model", "n")

sum_annotated$model <- gsub("all", "All", sum_annotated$model)
sum_annotated$model <- gsub("ams", "Annual mating success", sum_annotated$model)

sum_annotated$region <- gsub("downstream", "Downstream", sum_annotated$region)
sum_annotated$region <- gsub("upstream", "Upstream", sum_annotated$region)
sum_annotated$region <- gsub("exon", "Exon", sum_annotated$region)
sum_annotated$region <- gsub("fiveUTR", "5' UTR", sum_annotated$region)
sum_annotated$region <- gsub("gene", "Gene body", sum_annotated$region)
sum_annotated$region <- gsub("intron", "Intron", sum_annotated$region)
sum_annotated$region <- gsub("promoter", "Promoter", sum_annotated$region)
sum_annotated$region <- gsub("threeUTR", "3' UTR", sum_annotated$region)

sum_annotated$region <- factor(sum_annotated$region, levels = c("3' UTR", "5' UTR", "Downstream", "Upstream", "Gene body", "Exon", "Intron", "Promoter", "TSS"))

# add total sig CpGs
sum_annotated <- sum_annotated %>% mutate(n_total = case_when(
  model == "All" ~ nrow(delta_out_ams),
  model == "Annual mating success" ~ nrow(cpg_sig_ams)))

sum_annotated <- sum_annotated %>% mutate(perc = n / n_total * 100)

write.csv(sum_annotated, file="results/modeloutput/fitness/summary_regions_sig_cpgs_ams.csv", row.names=F, quote=F)

#### Plot number of sites per region ####
clrs <- viridisLite::viridis(6)
ggplot(subset(sum_annotated, region != "5' UTR" & region != "3' UTR"), aes(x = region, y = perc, fill = model)) + geom_bar(stat="identity", position="dodge") + 
  labs(y="Percentage of CpG sites", x="Region", fill = "Subset")+ coord_flip() + 
  scale_fill_manual(values=c(clrs[5], clrs[17])) -> num_perc_region
ggsave(num_perc_region, file="plots/model_out/fitness/perc_sig_per_region_ams.png", width=10, height=12)


### merge with 'similar to' column to get gene IDs
lookup <- fread("/home/nioo/rebeccash/PhD_grouse/grouse-annotation/data/lookup_ANN_gene_id.txt")
names(lookup) <- c("original", "similar")

all_models_sig_annotated_id <- left_join(all_models_sig_annotated, lookup, by = c("ID" = "original"))

all_models_sig_annotated_id <- subset(all_models_sig_annotated_id,
                                      !is.na(similar) &
                                        !grepl("LOC", similar))

all_models_sig_annotated_id$similar <- toupper(all_models_sig_annotated_id$similar)

subset(all_models_sig_annotated_id, parameter == "all") %>% arrange(parameter_qval) %>%
  dplyr::select(similar) %>% unique() %>% write.csv("results/modeloutput/fitness/all_gene_ids_changing_similar.csv", quote=F, row.names=F, col.names = F)

subset(all_models_sig_annotated_id, parameter == "ams") %>% arrange(parameter_qval) %>%
  dplyr::select(similar) %>% unique() %>% write.csv("results/modeloutput/fitness/gene_ids_sig_ams_similar.csv", quote=F, row.names=F, col.names = F)
