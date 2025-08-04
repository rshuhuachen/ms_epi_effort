
## load packages
pacman::p_load(tidyverse, data.table, methylKit, tibble, matrixStats, ggpointdensity)

source("scripts/plotting_theme.R")

## load data
load("data/processed/methylkit_prepost_raw.RData")
data <- getData(ltet_meth_unite)

#### Get summary statistics ####

# select columns for specific data: number of sites that are methylated, coverage
# number of columns in data
ncol <- ncol(data) 

# first column with data respectively
first_cov <- 5
first_meth_n <- 6

# select columns 
cols_cov <- seq(first_cov, ncol-(2-(first_cov-5)), 3)
cols_meth_n <- seq(first_meth_n, ncol-(2-(first_meth_n-5)), 3)

# make summary across all cpg sites
summary <- data.frame(site = c(paste(data$chr, data$start, sep="_")))
summary <- summary %>% mutate(mean_cov = rowMeans(data[cols_cov], na.rm=T),
                              sd_cov = rowSds(as.matrix(data[cols_cov]), na.rm=T),
                              mean_n_meth = rowMeans(data[cols_meth_n], na.rm=T),
                              sd_n_meth = rowSds(as.matrix(data[cols_meth_n]), na.rm=T),
                              na =  rowSums(is.na(data[cols_cov])),
                              n =  rowSums(!is.na(data[cols_cov])))

# calculate percentage methylation
sum_meth_prop <- data.frame(site = c(paste(data$chr, data$start, sep="_")))
for (i in 1:length(cols_cov)){
  sum_meth_prop <- sum_meth_prop %>% mutate(methProp = data[cols_meth_n[i]]/data[cols_cov[i]])
  names(sum_meth_prop)[i+1] <- ltet_meth_unite@sample.ids[i] 
}                              

# add summary percentage methylation
summary <- summary %>% mutate(mean_perc_meth = rowMeans(sum_meth_prop[,-1], na.rm=T),
                              sd_perc_meth = rowSds(as.matrix(sum_meth_prop[,-1]), na.rm=T), .before = mean_cov)


# site to rownames
summary <- summary %>% remove_rownames %>% column_to_rownames(var = "site")
write.csv(summary, file = "results/qc/summary_coverage_meth_unfiltered.csv", row.names=T, quote=F)
write.csv(summary[c(1:10),], file = "results/tables/summary_coverage_meth_unfiltered_n10.csv", row.names=T, quote=F)

### Plot distributions of mean methylation etc. ####

ggplot(summary, aes(mean_cov)) + geom_histogram() + labs(x = "Mean coverage") -> hist_mean_cov
ggplot(summary, aes(sd_cov)) + geom_histogram() + labs(x = "Mean coverage SD") -> hist_sd_cov
ggplot(summary, aes(mean_n_meth)) + geom_histogram()+ labs(x = "Mean number of methylated C's") -> hist_mean_n_meth
ggplot(summary, aes(sd_n_meth)) + geom_histogram() + labs(x = "Mean SD of number of methylated C's") -> hist_sd_n_meth
ggplot(summary, aes(n)) + geom_histogram() + labs(x = "Mean sample size") -> hist_n
ggplot(summary, aes(mean_perc_meth)) + geom_histogram() + labs(x = "Mean methylation %") -> hist_mean_perc_meth

cowplot::plot_grid(hist_mean_cov, hist_sd_cov, hist_mean_n_meth, hist_sd_n_meth, hist_n, hist_mean_perc_meth,
                  labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> hist_plots

ggsave(hist_plots, file="plots/explore/hist_summary_stats_cov_meth.png", width=14, height=12)

#### Plot relationshp coverage and dna methylation #####

# geom pointdensity to get an idea of which points are where since many are overlapping
# on a subset of random CpGs
random_n <- 100000
ggplot(data=summary[sample(nrow(summary),random_n),], aes(x=mean_cov, y=mean_perc_meth)) + labs(x = "Mean coverage", y = "Mean methylation percentage") +
  scale_x_continuous(trans = scales::log_trans(),
                     breaks = scales::log_breaks()) +
  geom_pointdensity()  +
  scale_color_viridis_c() -> plot_cov_meth 

ggsave(plot_cov_meth, file = "plots/explore/plot_coverage_mean_methylation_random.png", width = 10, height = 10)

#### PCA ####
# create dataset for PCA with only complete data
data_pca <- data.frame(matrix(unlist(sum_meth_prop), nrow=nrow(sum_meth_prop)),stringsAsFactors=FALSE)
data_pca <- data_pca[complete.cases(data_pca),] 
nrow(data_pca) 
data_pca <- data_pca[,-1]
data_pca <- lapply(data_pca, as.numeric)

# conduct pca and save plots
PCA <- prcomp(t(as.data.frame(data_pca)), center=F, scale=F) # t() transposes the matrix meth_PCA to get one coordinate for each id

save(PCA, file = "results/pca/pca.RData")

png(file = "plots/explore/pca_biplot.png", width=1000, height = 1000)
biplot(PCA,scale = T,center=F,cex=0.25)
dev.off()

ggplot() + geom_point(aes(x=PCA$x[,1],y=PCA$x[,2]), size = 3) + labs(x = "PC1", y = "PC2") -> pca_plot
ggsave(pca_plot, file = "plots/explore/pca.png", width = 10, height=10)

# get eigenvalues and percentage explained
fviz_eig(PCA)

eigs <- PCA$sdev^2
var <- eigs/sum(eigs)
explained<-100*eigs/sum(eigs)

head(explained) #PC1 explains 97% of data

# PCA coloured by library, site

# first collect data required to plot the PCs

# data on lek site, sampling 
load("data/phenotypes/fulldata_complete_epi_withdates.RData")
meta_ltet <- all_pheno_epi %>% dplyr::select(c(epi_nr, id, site, prepost, year)) %>% filter(!is.na(prepost))

# data on QC
load(file="/home/nioo/rebeccash/PhD_grouse/methylation_grouse/data/genomics/qc_epi.RData")

qc_ltet <- qc %>% dplyr::select(c(sample_id, Sample, lib, n_in_lib, conc_std, batch, extraction_batch))
qc_ltet <- subset(qc_ltet, Sample %in% meta_ltet$epi_nr)

# exclude repeats
qc_ltet <- subset(qc_ltet, sample_id != "lib99_1")
qc_ltet <- subset(qc_ltet, sample_id != "lib20_119")
qc_ltet <- subset(qc_ltet, sample_id != "lib20_191")
qc_ltet <- subset(qc_ltet, sample_id != "lib7_250")

# combine data with PC loadings
merge_pca <- data.frame(sample_id = ltet_meth_unite@sample.ids)
merge_pca <- left_join(merge_pca, qc_ltet, by = "sample_id")
merge_pca <- left_join(merge_pca, meta_ltet, by = c("Sample" = "epi_nr"))
merge_pca$pc1 <- PCA$x[,1]
merge_pca$pc2 <- PCA$x[,2]

save(merge_pca, file = "results/pca/quality_check_pca_meanmeth.RData")

# and with mean methylation level
summary_per_sample <- data.frame(sample_id = ltet_meth_unite@sample.ids)
summary_per_sample <- summary_per_sample %>% mutate(mean_perc_meth = colMeans(sum_meth_prop[,-1], na.rm=T))

merge_pca <- left_join(merge_pca, summary_per_sample, by = "sample_id")

# plot PCs
ggplot(merge_pca, aes(x = pc1, y = pc2)) + geom_point(size=3, aes(col = site)) + 
    labs(x = "PC 1", y = "PC 2", col = "Lek site") +
    scale_color_viridis_d() -> pca_site

ggplot(merge_pca, aes(x = pc1, y = pc2)) + geom_point(size=3, aes(col = prepost)) + 
    labs(x = "PC 1", y = "PC 2", col = "Time period") +
    scale_color_viridis_d() -> pca_prepost

ggplot(merge_pca, aes(x = pc1, y = pc2)) + geom_point(size=3, aes(col = as.factor(year))) + 
    labs(x = "PC 1", y = "PC 2", col = "Year") +
    scale_color_viridis_d() -> pca_year
    
ggplot(merge_pca, aes(x = pc1, y = pc2)) + geom_point(size=3, aes(col = lib)) + 
    labs(x = "PC 1", y = "PC 2", col = "Library") +
    scale_color_viridis_d() -> pca_lib

ggplot(merge_pca, aes(x = pc1, y = pc2)) + geom_point(size=3, aes(col = as.factor(n_in_lib))) + 
    labs(x = "PC 1", y = "PC 2", subtitle = "Sample within library") +
    theme(legend.position="none")+
    scale_color_viridis_d() -> pca_lib_n

merge_pca$conc_std[which(merge_pca$conc_std > 400)] <- NA
ggplot(merge_pca, aes(x = pc1, y = pc2)) + geom_point(size=3, aes(col = conc_std)) + 
    labs(x = "PC 1", y = "PC 2", col = "Concentration") +
    scale_color_viridis_c() -> pca_conc

ggplot(merge_pca, aes(x = pc1, y = pc2)) + geom_point(size=3, aes(col = as.factor(batch))) + 
    labs(x = "PC 1", y = "PC 2", col = "Batch") +
    scale_color_viridis_d() -> pca_batch

ggplot(merge_pca, aes(x = pc1, y = pc2)) + geom_point(size=3, aes(col = as.factor(extraction_batch))) + 
    labs(x = "PC 1", y = "PC 2", col = "Batch") +
    scale_color_viridis_d() -> pca_batch_ext    

cowplot::plot_grid(pca_site, pca_prepost, pca_year, pca_lib, pca_lib_n, pca_conc, pca_batch, pca_batch_ext,
                  labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) -> pca_plots

ggsave(pca_plots, file = "plots/explore/pca_plots_factors.png", width=12, height=20)


#### LM ####

lmer_lib <- lmerTest::lmer(mean_perc_meth ~ lib + (1|id), data = merge_pca)
lmer_null  <- lmerTest::lmer(mean_perc_meth ~ 1 + (1|id), data = merge_pca)
anova(lmer_lib, lmer_null)
