
##### This script identifies which CpG sites change from pre- to post-lekking

### Packages ####
pacman::p_load(dplyr, data.table, tibble, performance,  
               parallel, performance, lmerTest, tidystats, cowplot, gaston)

### Plotting ###
source("scripts/plotting_theme.R")

### Data ####

#### Load data ####
### converted prepost meth data
load(file = "data/processed/methylkit_prepost_long_onlyvar_thres0.3_min_0.5_group_15X.RData") #prepost_long only variation
prepost_long <- prepost_long_clean_15X
rm(prepost_long_clean_15X)

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
# Core referes to whether the male belongs to the core dataset. 
# This core dataset includes males that were first captured as yearlings, and therefore we know their exact age
# Non-core males were caught as adullts, and therefore their exact age can only be estimated from plumage but is not known
# with certainty

# swap levels
prepost_long$prepost <- factor(prepost_long$prepost, levels = c("pre", "post"))

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

#### glmer version ####
function_model_glmer <- function(df){tryCatch({
  chr_pos <- as.character(df[1,1])
  df <- as.data.frame(df)
  df$prepost <- as.factor(df$prepost)
  df$id <- as.factor(df$id)
  
  # model
  model <- lme4::glmer(cbind(numC, numT) ~ prepost + (1|id), family = "binomial", df)
  
  #fixed effects
  prepost_estimate <- summary(model)$coefficients["prepostpost","Estimate"]
  prepost_se <- summary(model)$coefficients["prepostpost","Std. Error"]
  prepost_zval <- summary(model)$coefficients["prepostpost","z value"]
  prepost_pval <-  summary(model)$coefficients["prepostpost","Pr(>|z|)"]
  
  #random effects 
  id_sd <- attributes(VarCorr(model)$id)$stddev
  id_variance <- data.frame(VarCorr(model), comp="Variance")[1,"vcov"]
  
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
out_glmer <- parallel::mclapply(data, function_model_glmer, mc.cores=4) 

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

out_glmer_raw <- do.call(rbind.data.frame, out_glmer)
save(out_glmer_raw, file="results/modeloutput/prepost_modeloutput_glmer_min0.75_15X_raw.RData")

#### Exclude models that did not converge ####
out_glmer_raw_conv <- subset(out_glmer_raw, convergence == "boundary (singular) fit: see help('isSingular')" | is.na(convergence))
nrow(out_glmer_raw_conv) / nrow(out_glmer_raw) * 100 # retain 97.5%, 345937 out of 354765

##### Overdispersion raw data #####

## explore dispersion ratio and make plots

summary(out_glmer_raw_conv$dispersion.ratio)

# histogram dispersion ratio raw 
ggplot(out_glmer_raw_conv, aes(x = dispersion.ratio)) + geom_histogram() + geom_vline(xintercept = 1., col = "red", linetype = "dotted", linewidth=1) +
  scale_y_log10() + labs(title = "Histogram dispersion ratio", subtitle= "Raw model output changing CpGs") -> hist_glmer_disp_ratio_raw

# histogram p-values raw
ggplot(out_glmer_raw_conv, aes(x = prepost_pval)) + geom_histogram() + 
  scale_y_continuous(labels = scales::unit_format(unit = "K", scale = 1e-3)) + 
  labs(title = "Histogram p-values", subtitle="Raw model output changing CpGs")-> hist_glmer_pvals_raw

plot_grid(hist_glmer_disp_ratio_raw, hist_glmer_pvals_raw, labs="auto", align="hv", axis="lb", ncol=1, label_fontface = "plain", label_size = 22)-> hists_glmer_raw

ggsave(hists_glmer_raw, file = "plots/model_out/changing/hist_raw_glmer_15X.png", width = 12, height = 12)

# qqplot raw

png(file = "plots/model_out/changing/qqplot_raw_glmer_15X.png", width = 800, height = 800)
qqplot.pvalues(out_glmer_raw_conv$prepost_pval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95) 
dev.off()

## filter for 95 percentile
out_glmer<- subset(out_glmer_raw_conv, dispersion.ratio < as.vector(quantile(out_glmer_raw_conv$dispersion.ratio, 0.95)))
nrow(out_glmer) # 235714

# histogram dispersion ratio filtered 
ggplot(out_glmer, aes(x = dispersion.ratio)) + geom_histogram() + 
  scale_y_log10() + labs(title = "Histogram dispersion ratio", subtitle= "Filtered model output changing CpGs") -> hist_glmer_disp_ratio_filter

# histogram p-values filtered
ggplot(out_glmer, aes(x = prepost_pval)) + geom_histogram() + 
  scale_y_continuous(labels = scales::unit_format(unit = "K", scale = 1e-3)) + 
  labs(title = "Histogram p-values", subtitle="Filtered model output changing CpGs")-> hist_glmer_pvals_filter

plot_grid(hist_glmer_disp_ratio_filter, hist_glmer_pvals_filter, labs="auto", align="hv", axis="lb", ncol=1, label_fontface = "plain", label_size = 22)-> hists_glmer_filter

ggsave(hists_glmer_filter, file = "plots/model_out/changing/hist_filtered_glmer_15X.png", width = 12, height = 12)

# qqplot filtered

png(file = "plots/model_out/changing/qqplot_filtered_glmer_15X.png", width = 800, height = 800)
qqplot.pvalues(out_glmer$prepost_pval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95) 
dev.off()

#### FDR-correction ####

out_glmer$prepost_qval <- p.adjust(out_glmer$prepost_pval, method = "fdr", n = nrow(out_glmer))

#### Filter for mean delta methylation ####

load(file = "results/modeloutput/all_sites_deltameth.RData")

### Calculate average delta_meth per CpG site

mean_delta_meth <- delta_meth %>% group_by(chr_pos) %>% summarise_at(vars(delta_meth), funs(mean(., na.rm=TRUE)))
names(mean_delta_meth)[2] <- "mean_delta_meth"

out_glmer <- left_join(out_glmer, mean_delta_meth, by = "chr_pos")

### Filter min absolute mean methylation of 10%

sub_glmer_prepost <- subset(out_glmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.1)
nrow(sub_glmer_prepost) #326

sensitivity_threshold <- data.frame(threshold = seq(0.05, 0.25, by = 0.01),
                                    n_sig = c(nrow(subset(out_glmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.05)),
                                              nrow(subset(out_glmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.06)),
                                              nrow(subset(out_glmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.07)),
                                              nrow(subset(out_glmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.08)),
                                              nrow(subset(out_glmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.09)),
                                              nrow(subset(out_glmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.10)),
                                              nrow(subset(out_glmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.11)),
                                              nrow(subset(out_glmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.12)),
                                              nrow(subset(out_glmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.13)),
                                              nrow(subset(out_glmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.14)),
                                              nrow(subset(out_glmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.15)),
                                              nrow(subset(out_glmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.16)),
                                              nrow(subset(out_glmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.17)),
                                              nrow(subset(out_glmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.18)),
                                              nrow(subset(out_glmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.19)),
                                              nrow(subset(out_glmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.20)),
                                              nrow(subset(out_glmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.21)),
                                              nrow(subset(out_glmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.22)),
                                              nrow(subset(out_glmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.23)),
                                              nrow(subset(out_glmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.24)),
                                              nrow(subset(out_glmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.25))))

source("scripts/plotting_theme.R")
ggplot(sensitivity_threshold, aes(x = threshold, y = n_sig)) + geom_point() + geom_line() + 
  labs(x = "Threshold", y = "Number of dynamic CpG sites") -> plot_thres

ggsave(plot_thres, file = "plots/final/supp/sfig_sensitivity_15X.png", width=10,height=6)

### Save original data (per CpG site per individual) for models but only subset significant CpG sites
changing_cpg <- subset(prepost_long, chr_pos %in% sub_glmer_prepost$chr_pos)
save(changing_cpg, file="results/modeloutput/changing/changing_sites_glmer_15X.RData")

### Save the model output
save(out_glmer, file="results/modeloutput/changing/modeloutput_glmer_15X.RData")

#### Plotting results ####

### Volcano plot

### Add column indicationg it's significant or not

out_glmer <- out_glmer %>% mutate(sig = as.factor(case_when(abs(mean_delta_meth) >= 0.1 & prepost_qval < 0.05 ~ "sig", TRUE ~ "nonsig")))

ggplot(out_glmer, aes(x = mean_delta_meth, y = -log10(as.numeric(prepost_qval)))) + 
  geom_point(size=4, alpha=0.5, aes(col = as.factor(sig))) +
  labs(x = expression("Mean "*Delta*" methylation %"), y = "-log10(q-value)") +
  scale_color_manual(values=c(clrs[5], clrs[17])) +
  geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = -0.1, col = "darkred", linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = 0.1, col = "darkred", linetype = "dotted", linewidth = 1) +
  theme(legend.position="none") -> volcano_change

ggsave(volcano_change, file = "plots/model_out/changing/volcano_change_15X.png", width=10, height=10)    

### How many get up- and how many get downregulated?
sub_glmer_prepost <- sub_glmer_prepost %>% mutate(posneg = as.factor(case_when(mean_delta_meth > 0 ~ "pos", mean_delta_meth < 0 ~ "neg")))
summary(sub_glmer_prepost$posneg) #29 neg, 297 pos

### Manhattan plot

# load scaffold numbers
load("data/scaffold_names_dovetail.RData")

# Split the chr_pos column into two columns based on the first "_"
split_chr_pos <- strsplit(out_glmer$chr_pos, "_", fixed = TRUE)

# Extract the numbers following HRSCAF=XXX_number
out_glmer$chr <- paste0(sapply(split_chr_pos, "[", 1), "_",
                        sapply(split_chr_pos, "[", 2), ";", 
                        sapply(split_chr_pos, "[", 4), "=",
                        sapply(split_chr_pos, "[", 5))

out_glmer$pos <- as.numeric(sapply(split_chr_pos, "[", 6))

# join
out_glmer <- left_join(out_glmer, genome[,c("contig", "scaf_nr")], by = c("chr" = "contig"))

# plot 
# lmer
test <- sample_n(out_glmer, 100)
out_glmer <- out_glmer %>% mutate(col = case_when(scaf_nr %% 2 == 0 ~ "even",
                                                  TRUE ~ "odd"))

shade <- out_glmer %>%
  subset(scaf_nr <= 30) %>%
  group_by(scaf_nr) %>%
  summarise(min = min(pos), max = max(pos)) %>%
  mutate(min = case_when(scaf_nr == 2 | scaf_nr == 4 | scaf_nr == 6 | scaf_nr == 8 | scaf_nr == 10 |
                           scaf_nr == 12 | scaf_nr == 14 | scaf_nr == 16 | scaf_nr == 18 | scaf_nr == 20 |
                           scaf_nr == 22 | scaf_nr == 24 | scaf_nr == 26 | scaf_nr == 28 | scaf_nr == 30 ~ 0,
                         TRUE ~ min)) %>%
  mutate(max = case_when(scaf_nr == 2 | scaf_nr == 4 | scaf_nr == 6 | scaf_nr == 8 | scaf_nr == 10 |
                           scaf_nr == 12 | scaf_nr == 14 | scaf_nr == 16 | scaf_nr == 18 | scaf_nr == 20 |
                           scaf_nr == 22 | scaf_nr == 24 | scaf_nr == 26 | scaf_nr == 28 | scaf_nr == 30  ~ 0,
                         TRUE ~ max))

#clrs <- viridisLite::viridis(6)
out_glmer %>% subset(scaf_nr <= 30) %>% 
  ggplot(aes(x = pos, y = -log10(as.numeric(prepost_pval)))) + 
  geom_point(size=5, alpha=0.5, aes(col = as.factor(col), fill = as.factor(col))) +
  facet_grid(~scaf_nr,scales = 'free_x', space = 'free_x', switch = 'x') +
  labs(x = "Scaffold number", y = expression(-log[10]*"(p-value)")) +
  #geom_rect(data=shade, aes(xmin=min, xmax=max, ymin=0, ymax=-log10(as.numeric(test$prepost_pval))), 
  #        alpha=0.5, fill = "#eceff4") + # "#f7f7f7" "#eceff4"
  #xlim(-1, 1)+
  scale_color_manual(values=c(clrs[5], clrs[17])) +
  scale_fill_manual(values=alpha(c(clrs[5], clrs[17]), 0.5)) +
  geom_hline(yintercept = -log10(0.05/nrow(out_glmer)), col = "darkred", linetype = "dotted", linewidth = 1) +
  theme(axis.text.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        plot.margin = margin(r = 0.5, l = 0.1, b = 0.1, t = 0.1, unit = "cm"),
        axis.line.x = element_blank(),
        legend.position="none",
        axis.title.x = element_text(margin=margin(t=10)),
        axis.title.y = element_text(margin=margin(r=5)),
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank()) -> manhattan_change

ggsave(manhattan_change, file = "plots/model_out/changing/manhattan_glmer_15X.png", width=26, height=10)    

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
  geom_point(aes(alpha = 0.8, size=cov), col = clrs[17], position = position_jitter(width = 0.1, seed = 3922)) + 
  labs(x = "Time period", y = "Methylation percentage", subtitle = out_glmer$chr_pos[1]) +
  theme(legend.position="none") +
  ylim(0,1)-> plot_top_cpg_1

subset(changing_cpg, chr_pos == out_glmer$chr_pos[2]) %>%
  arrange(id, year) %>%
  ggplot(., aes(x = prepost, y = methperc))+
  geom_boxplot(linewidth=1, outlier.shape=NA) + 
  geom_path(aes(group = id_year), alpha = 0.8, col = "grey60", position = position_jitter(width = 0.1, seed = 3922)) +
  geom_point(aes(alpha = 0.8, size=cov), col = clrs[17], position = position_jitter(width = 0.1, seed = 3922)) + 
  labs(x = "Time period", y = "Methylation percentage", subtitle = out_glmer$chr_pos[2]) +
  theme(legend.position="none") +
  ylim(0,1)-> plot_top_cpg_2

subset(changing_cpg, chr_pos == out_glmer$chr_pos[3]) %>%
  arrange(id, year) %>%
  ggplot(., aes(x = prepost, y = methperc))+
  geom_boxplot(linewidth=1, outlier.shape=NA) + 
  geom_path(aes(group = id_year), alpha = 0.8, col = "grey60", position = position_jitter(width = 0.1, seed = 3922)) +
  geom_point(aes(alpha = 0.8, size=cov), col = clrs[17], position = position_jitter(width = 0.1, seed = 3922)) + 
  labs(x = "Time period", y = "Methylation percentage", subtitle = out_glmer$chr_pos[3]) +
  theme(legend.position="none") +
  ylim(0,1)-> plot_top_cpg_3

subset(changing_cpg, chr_pos == out_glmer$chr_pos[4]) %>%
  arrange(id, year) %>%
  ggplot(., aes(x = prepost, y = methperc))+
  geom_boxplot(linewidth=1, outlier.shape=NA) + 
  geom_path(aes(group = id_year), alpha = 0.8, col = "grey60", position = position_jitter(width = 0.1, seed = 3922)) +
  geom_point(aes(alpha = 0.8, size=cov), col = clrs[17], position = position_jitter(width = 0.1, seed = 3922)) + 
  labs(x = "Time period", y = "Methylation percentage", subtitle = out_glmer$chr_pos[4]) +
  theme(legend.position="none") +
  ylim(0,1)-> plot_top_cpg_4

subset(changing_cpg, chr_pos == out_glmer$chr_pos[5]) %>%
  arrange(id, year) %>%
  ggplot(., aes(x = prepost, y = methperc))+
  geom_boxplot(linewidth=1, outlier.shape=NA) + 
  geom_path(aes(group = id_year), alpha = 0.8, col = "grey60", position = position_jitter(width = 0.1, seed = 3922)) +
  geom_point(aes(alpha = 0.8, size=cov), col = clrs[17], position = position_jitter(width = 0.1, seed = 3922)) + 
  labs(x = "Time period", y = "Methylation percentage", subtitle = out_glmer$chr_pos[5]) +
  theme(legend.position="none") +
  ylim(0,1)-> plot_top_cpg_5

subset(changing_cpg, chr_pos == out_glmer$chr_pos[6]) %>%
  arrange(id, year) %>%
  ggplot(., aes(x = prepost, y = methperc))+
  geom_boxplot(linewidth=1, outlier.shape=NA) + 
  geom_path(aes(group = id_year), alpha = 0.8, col = "grey60", position = position_jitter(width = 0.1, seed = 3922)) +
  geom_point(aes(alpha = 0.8, size=cov), col = clrs[17], position = position_jitter(width = 0.1, seed = 3922)) + 
  labs(x = "Time period", y = "Methylation percentage", subtitle = out_glmer$chr_pos[6]) +
  theme(legend.position="none") +
  ylim(0,1)-> plot_top_cpg_6

cowplot::plot_grid(plot_top_cpg_1, plot_top_cpg_2, plot_top_cpg_3, plot_top_cpg_4, 
                   plot_top_cpg_5, plot_top_cpg_6, labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_change

ggsave(plots_change, file = "plots/model_out/changing/rawdata_top_sig_15X.png", width=16, height=20)    


