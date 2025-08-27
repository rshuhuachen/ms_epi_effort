
##### This script identifies which CpG sites change from pre- to post-lekking

### Packages ####
#devtools::install_github("mastoffel/rptR", build_vignettes = TRUE)
pacman::p_load(dplyr, data.table, tibble, performance,  
               parallel, performance, lmerTest, tidystats, cowplot, gaston)
#matrixStats, insight, rptR, 

### Plotting ###
source("scripts/plotting_theme.R")
#clrs <- viridisLite::viridis(6)

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

## for sup table: export relevant info ###
head(prepost)

table <- prepost %>% dplyr::select(c(id, year, site, fulldate, prepost, attend, dist, MS, surv))
ny <- read.csv("data/phenotypes/data_for_nextyear_corrected.csv")
table <- left_join(table, ny[,c("id", "year", "age_cat", "blue_nextyear", "lyre_nextyear")], by = c("id", "year"))
table <- table %>% dplyr::select(c(id, site, age_cat, prepost, fulldate, attend, dist, MS, surv, blue_nextyear, lyre_nextyear))
write.csv(table, file = "data/metadata_samples.csv", quote=F, row.names = F)

### merge with some metadata

prepost_long <- left_join(prepost_long, prepost[,c("id", "year", "Core", "born", "fulldate")], 
                          by = c("id", "year", "fulldate"))

prepost_long <- prepost_long %>% mutate(age_year = as.factor(case_when(Core == "Core" ~ year - born,
                                                        Core == "No core" ~ NA)),
                                        age = as.factor(case_when(Core == "Core" & (year - born > 1) ~ "Adult",
                                                        Core == "Core" & (year - born == 1) ~ "Yearling",
                                                        Core == "No core" ~ "Adult")))

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
save(out_glmer_raw, file="results/modeloutput/prepost_modeloutput_glmer_min0.75_raw.RData")

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

ggsave(hists_glmer_raw, file = "plots/model_out/changing/hist_raw_glmer.png", width = 12, height = 12)

# qqplot raw

png(file = "plots/model_out/changing/qqplot_raw_glmer.png", width = 800, height = 800)
qqplot.pvalues(out_glmer_raw_conv$prepost_pval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95) 
dev.off()

## filter for 95 percentile
out_glmer<- subset(out_glmer_raw_conv, dispersion.ratio < as.vector(quantile(out_glmer_raw_conv$dispersion.ratio, 0.95)))
nrow(out_glmer) # 328640

# histogram dispersion ratio filtered 
ggplot(out_glmer, aes(x = dispersion.ratio)) + geom_histogram() + 
scale_y_log10() + labs(title = "Histogram dispersion ratio", subtitle= "Filtered model output changing CpGs") -> hist_glmer_disp_ratio_filter

# histogram p-values filtered
ggplot(out_glmer, aes(x = prepost_pval)) + geom_histogram() + 
  scale_y_continuous(labels = scales::unit_format(unit = "K", scale = 1e-3)) + 
  labs(title = "Histogram p-values", subtitle="Filtered model output changing CpGs")-> hist_glmer_pvals_filter

plot_grid(hist_glmer_disp_ratio_filter, hist_glmer_pvals_filter, labs="auto", align="hv", axis="lb", ncol=1, label_fontface = "plain", label_size = 22)-> hists_glmer_filter

ggsave(hists_glmer_filter, file = "plots/model_out/changing/hist_filtered_glmer.png", width = 12, height = 12)

# qqplot filtered

png(file = "plots/model_out/changing/qqplot_filtered_glmer.png", width = 800, height = 800)
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
nrow(sub_glmer_prepost)
# 1026

### Save original data (per CpG site per individual) for models but only subset significant CpG sites
changing_cpg <- subset(prepost_long, chr_pos %in% sub_glmer_prepost$chr_pos)
save(changing_cpg, file="results/modeloutput/changing/changing_sites_glmer.RData")

### Save the model output
save(out_glmer, file="results/modeloutput/changing/modeloutput_glmer.RData")

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

ggsave(volcano_change, file = "plots/model_out/changing/volcano_change.png", width=10, height=10)    

### How many get up- and how many get downregulated?
sub_glmer_prepost <- sub_glmer_prepost %>% mutate(posneg = as.factor(case_when(mean_delta_meth > 0 ~ "pos", mean_delta_meth < 0 ~ "neg")))
summary(sub_glmer_prepost$posneg)

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

ggsave(manhattan_change, file = "plots/model_out/changing/manhattan_glmer.png", width=26, height=10)    

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

ggsave(plots_change, file = "plots/model_out/changing/rawdata_top_sig.png", width=16, height=20)    

#### Annotation of changing CpG sites: gene regions #####

### Packages ####
pacman::p_load(genomation, GenomicFeatures, rtracklayer, 
               GenomicRanges)


### Combine all sites vs changing sites
cpg_all <- out_glmer %>% dplyr::select(c(chr_pos, prepost_qval))
names(cpg_all)[2] <- "parameter_qval"
cpg_all$parameter <- "all"

cpg_changing_select <- sub_glmer_prepost %>% dplyr::select(c(chr_pos, prepost_qval))
names(cpg_changing_select)[2] <- "parameter_qval"
cpg_changing_select$parameter <- "time_period"

all_models_sig <- rbind(cpg_all, cpg_changing_select)

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

save(all_models_sig_annotated_raw, file="results/modeloutput/changing/annotated_sig_cpg_raw.RData")

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

save(all_models_sig_annotated, file="results/modeloutput/changing/annotated_sig_cpg.RData")

### merge with 'similar to' column to get gene IDs
lookup <- fread("/home/nioo/rebeccash/PhD_grouse/grouse-annotation/data/lookup_ANN_gene_id.txt")
names(lookup) <- c("original", "similar")

all_models_sig_annotated_id <- left_join(all_models_sig_annotated, lookup, by = c("ID" = "original"))

all_models_sig_annotated_id <- subset(all_models_sig_annotated_id,
                                   !is.na(similar) &
                                     !grepl("LOC", similar))

all_models_sig_annotated_id$similar <- toupper(all_models_sig_annotated_id$similar)

subset(all_models_sig_annotated_id, parameter == "all") %>% arrange(parameter_qval) %>%
  dplyr::select(similar) %>% unique() %>% write.csv("results/modeloutput/all_gene_ids_similar.csv", quote=F, row.names=F, col.names = F)

subset(all_models_sig_annotated_id, parameter == "time_period") %>% arrange(parameter_qval) %>%
  dplyr::select(similar) %>% unique() %>% write.csv("results/modeloutput/changing/gene_ids_sig_changing_similar.csv", quote=F, row.names=F, col.names = F)

subset(all_models_sig_annotated_id, parameter == "time_period") %>% arrange(parameter_qval) -> annotated_changing
save(annotated_changing, file = "results/modeloutput/changing/gene_ids_sig_changing_similar.RData")

#### Summarise number of sites per region ####
sum_annotated <- as.data.frame(table(as.factor(all_models_sig_annotated$region), all_models_sig_annotated$parameter))
names(sum_annotated) <- c("region", "model", "n")

sum_annotated$model <- gsub("all", "All", sum_annotated$model)
sum_annotated$model <- gsub("time_period", "Time period", sum_annotated$model)

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
  model == "All" ~ nrow(out_glmer),
  model == "Time period" ~ nrow(sub_glmer_prepost)))

sum_annotated <- sum_annotated %>% mutate(perc = n / n_total * 100)

write.csv(sum_annotated, file="results/modeloutput/changing/summary_regions_sig_cpgs.csv", row.names=F, quote=F)

#### Plot number of sites per region ####

ggplot(subset(sum_annotated, region != "5' UTR" & region != "3' UTR"), aes(x = region, y = perc, fill = model)) + geom_bar(stat="identity", position="dodge") + 
  labs(y="Percentage of CpG sites", x="Region", fill = "Subset")+ coord_flip() + 
  scale_fill_manual(values=c("grey60", clrs[4])) -> num_perc_region
ggsave(num_perc_region, file="plots/model_out/changing/perc_sig_per_region.png", width=10, height=12)

ggplot(sum_annotated, aes(x = region, y = n)) + geom_bar(stat="identity") + 
  labs(y="Number of CpG sites", x="Region")+
  facet_wrap(~model, ncol=2, scales="free") + coord_flip() -> num_region_cpg

ggsave(num_region_cpg, file="plots/model_out/changing/n_sig_cpg_per_region.png", width=10, height=18)

ggplot(sum_annotated, aes(x = region, y = perc)) + geom_bar(stat="identity") + ylim(0,100)+
  labs(y="Percentage of CpG sites", x="Region")+
  facet_wrap(~model, ncol=2,) + coord_flip() -> perc_region_cpg

ggsave(perc_region_cpg, file="plots/model_out/changing/perc_sig_cpg_per_region.png", width=20, height=12)
