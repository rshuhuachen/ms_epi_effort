
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

out_glmer_raw <- do.call(rbind.data.frame, out_glmer)
save(out_glmer_raw, file="results/modeloutput/prepost_modeloutput_glmer_min0.75_raw.RData")

## explore dispersion ratio
summary(out_glmer_raw$dispersion.ratio)
ggplot(out_glmer_raw, aes(x = dispersion.ratio)) + geom_histogram() + geom_vline(xintercept = 1., col = "red", linetype = "dotted") +
scale_y_log10()-> hist_glmer_disp_ratio
ggsave(hist_glmer_disp_ratio, file = "plots/model_out/changing_histogram_dispersion_glmer_raw.png", width = 8, height = 8)

## qqplot without filtering for overdispersion
pacman::p_load(gaston)

png(file = "plots/model_out/qqplot_changing_qqplot_glmer_raw.png", width = 800, height = 800)
qqplot.pvalues(out_glmer_raw$prepost_pval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95) 
dev.off()

### apply a FDR multiple-testing correction

## two options: dispersion filter by ratio < 1.1 (threshold)
out_glmer_threshold <- subset(out_glmer_raw, dispersion.ratio < 1.1 & dispersion.pval > 0.05 & is.na(convergence)) # 168720, 179311 with cov not nT
out_glmer_threshold$prepost_qval <- p.adjust(out_glmer_threshold$prepost_pval, method = "fdr", n = nrow(out_glmer_threshold))

nrow(out_glmer_threshold)/nrow(out_glmer_raw) # = 0.53
nrow(subset(out_glmer_threshold, prepost_qval < 0.05)) # N = 3563

# qq plot
png(file = "plots/model_out/qqplot_changing_qqplot_glmer_threshold.png", width = 800, height = 800)
qqplot.pvalues(out_glmer_threshold$prepost_qval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95) 
dev.off()

## second option: within the 90% quantiles
out_glmer_perc <- subset(out_glmer_raw, dispersion.ratio < as.vector(quantile(out_glmer_raw$dispersion.ratio, 0.95)) & dispersion.ratio > as.vector(quantile(out_glmer_raw$dispersion.ratio, 0.05)))
out_glmer_perc$prepost_qval <- p.adjust(out_glmer_perc$prepost_pval, method = "fdr", n = nrow(out_glmer_perc))
nrow(out_glmer_perc)/nrow(out_glmer_raw) # = 0.90 (obvs)
nrow(subset(out_glmer_perc, prepost_qval < 0.05)) # N = 19110

png(file = "plots/model_out/qqplot_changing_qqplot_glmer_90percentile.png", width = 800, height = 800)
qqplot.pvalues(out_glmer_perc$prepost_qval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95) 
dev.off()

## save the one we choose: threshold
save(out_glmer_threshold, file="results/modeloutput/prepost_modeloutput_glmer_min0.75.RData")

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

## after filtering, quite some underdispersion?

#### Filter for mean delta methylation ####

load(file = "results/modeloutput/all_sites_deltameth.RData")

### Calculate average delta_meth per CpG site

mean_delta_meth <- delta_meth %>% group_by(chr_pos) %>% summarise_at(vars(delta_meth), funs(mean(., na.rm=TRUE)))
names(mean_delta_meth)[2] <- "mean_delta_meth"

out_glmer <- left_join(out_glmer_threshold, mean_delta_meth, by = "chr_pos")
out_lmer <- left_join(out_lmer_threshold, mean_delta_meth, by = "chr_pos")

### Filter min absolute mean methylation of 10%

## save the epi data only from cpg's that are sig
sub_glmer_prepost <- subset(out_glmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.1)
# N = 3563 but with 10% mean difference only 434

sub_lmer_prepost <- subset(out_lmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.1)
# N = 204 
nrow(subset(sub_glmer_prepost, chr_pos %in% sub_lmer_prepost$chr_pos)) #105
nrow(subset(sub_lmer_prepost, chr_pos %in% sub_glmer_prepost$chr_pos)) #105

changing_cpg <- subset(prepost_long, chr_pos %in% sub_glmer_prepost$chr_pos)
save(changing_cpg, file="results/modeloutput/subset_sites_sig_prepost.RData")

changing_cpg_lmer <- subset(prepost_long, chr_pos %in% sub_lmer_prepost$chr_pos)
save(changing_cpg_lmer, file="results/modeloutput/subset_sites_sig_prepost_lmer.RData")

### volcano plot
source("scripts/plotting_theme.R")

out_glmer <- out_glmer %>% mutate(sig = as.factor(case_when(abs(mean_delta_meth) >= 0.1 & prepost_qval < 0.05 ~ "sig", TRUE ~ "nonsig")))

clrs <- viridisLite::viridis(6)
ggplot(out_glmer, aes(x = mean_delta_meth, y = -log10(as.numeric(prepost_qval)))) + 
    geom_point(size=4, alpha=0.5, aes(col = as.factor(sig))) +
    labs(x = expression("Mean "*Delta*" methylation %"), y = "-log10(q-value)") +
    #xlim(-1, 1)+
    scale_color_manual(values=c("grey60", clrs[4])) +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = -0.1, col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0.1, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") -> volcano_change

ggsave(volcano_change, file = "plots/model_out/volcano_change.png", width=10, height=10)    

### manhattan plot


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
                                        
clrs <- viridisLite::viridis(6)
out_glmer %>% subset(scaf_nr <= 30) %>% 
  ggplot(aes(x = pos, y = -log10(as.numeric(prepost_pval)))) + 
    geom_point(size=5, alpha=0.5, aes(col = as.factor(col), fill = as.factor(col))) +
    facet_grid(~scaf_nr,scales = 'free_x', space = 'free_x', switch = 'x') +
    labs(x = "Scaffold number", y = expression(-log[10]*"(p-value)")) +
    #geom_rect(data=shade, aes(xmin=min, xmax=max, ymin=0, ymax=-log10(as.numeric(test$prepost_pval))), 
    #        alpha=0.5, fill = "#eceff4") + # "#f7f7f7" "#eceff4"
    #xlim(-1, 1)+
    scale_color_manual(values=c(clrs[2], clrs[4])) +
    scale_fill_manual(values=alpha(c(clrs[2], clrs[4]), 0.5)) +
    geom_hline(yintercept = -log10(0.05/nrow(out_glmer)), col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(axis.text.x = element_blank(),
    panel.spacing = unit(0, "lines"),
    plot.margin = margin(r = 0.5, l = 0.1, b = 0.1, t = 0.1, unit = "cm"),
    axis.line.x = element_blank(),
    legend.position="none",
    axis.title.x = element_text(margin=margin(t=10)),
    axis.title.y = element_text(margin=margin(r=5)),
     axis.ticks.x = element_blank()
    axis.line.y = element_blank()) -> manhattan_change

ggsave(manhattan_change, file = "plots/model_out/manhattan_change_glmer.png", width=26, height=10)    


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
  theme(legend.position="none") +
  ylim(0,1)-> plot_top_cpg_1


subset(changing_cpg, chr_pos == out_glmer$chr_pos[2]) %>%
  arrange(id, year) %>%
  ggplot(., aes(x = prepost, y = methperc))+
  geom_boxplot(linewidth=1, outlier.shape=NA) + 
  geom_path(aes(group = id_year), alpha = 0.8, col = "grey60", position = position_jitter(width = 0.1, seed = 3922)) +
  geom_point(aes(alpha = 0.8, size=cov), col = clrs[4], position = position_jitter(width = 0.1, seed = 3922)) + 
  labs(x = "Time period", y = "Methylation percentage", subtitle = out_glmer$chr_pos[2]) +
  theme(legend.position="none") +
  ylim(0,1)-> plot_top_cpg_2

subset(changing_cpg, chr_pos == out_glmer$chr_pos[3]) %>%
  arrange(id, year) %>%
  ggplot(., aes(x = prepost, y = methperc))+
  geom_boxplot(linewidth=1, outlier.shape=NA) + 
  geom_path(aes(group = id_year), alpha = 0.8, col = "grey60", position = position_jitter(width = 0.1, seed = 3922)) +
  geom_point(aes(alpha = 0.8, size=cov), col = clrs[4], position = position_jitter(width = 0.1, seed = 3922)) + 
  labs(x = "Time period", y = "Methylation percentage", subtitle = out_glmer$chr_pos[3]) +
  theme(legend.position="none") +
  ylim(0,1)-> plot_top_cpg_3

subset(changing_cpg, chr_pos == out_glmer$chr_pos[4]) %>%
  arrange(id, year) %>%
  ggplot(., aes(x = prepost, y = methperc))+
  geom_boxplot(linewidth=1, outlier.shape=NA) + 
  geom_path(aes(group = id_year), alpha = 0.8, col = "grey60", position = position_jitter(width = 0.1, seed = 3922)) +
  geom_point(aes(alpha = 0.8, size=cov), col = clrs[4], position = position_jitter(width = 0.1, seed = 3922)) + 
  labs(x = "Time period", y = "Methylation percentage", subtitle = out_glmer$chr_pos[4]) +
  theme(legend.position="none") +
  ylim(0,1)-> plot_top_cpg_4

  subset(changing_cpg, chr_pos == out_glmer$chr_pos[5]) %>%
  arrange(id, year) %>%
  ggplot(., aes(x = prepost, y = methperc))+
  geom_boxplot(linewidth=1, outlier.shape=NA) + 
  geom_path(aes(group = id_year), alpha = 0.8, col = "grey60", position = position_jitter(width = 0.1, seed = 3922)) +
  geom_point(aes(alpha = 0.8, size=cov), col = clrs[4], position = position_jitter(width = 0.1, seed = 3922)) + 
  labs(x = "Time period", y = "Methylation percentage", subtitle = out_glmer$chr_pos[5]) +
  theme(legend.position="none") +
  ylim(0,1)-> plot_top_cpg_5

subset(changing_cpg, chr_pos == out_glmer$chr_pos[6]) %>%
  arrange(id, year) %>%
  ggplot(., aes(x = prepost, y = methperc))+
  geom_boxplot(linewidth=1, outlier.shape=NA) + 
  geom_path(aes(group = id_year), alpha = 0.8, col = "grey60", position = position_jitter(width = 0.1, seed = 3922)) +
  geom_point(aes(alpha = 0.8, size=cov), col = clrs[4], position = position_jitter(width = 0.1, seed = 3922)) + 
  labs(x = "Time period", y = "Methylation percentage", subtitle = out_glmer$chr_pos[6]) +
  theme(legend.position="none") +
  ylim(0,1)-> plot_top_cpg_6

cowplot::plot_grid(plot_top_cpg_1, plot_top_cpg_2, plot_top_cpg_3, plot_top_cpg_4, 
                    plot_top_cpg_5, plot_top_cpg_6, labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> plots_change

ggsave(plots_change, file = "plots/model_out/plot_top_change_cpg.png", width=16, height=20)    
