#### Packages #####
pacman::p_load(tidyverse, cowplot, data.table)

#### Theme ####
source("scripts/plotting_theme.R")

##### Fig 1 - changing CpG sites ####

##### Load data #####
load(file="results/modeloutput/changing/changing_sites_glmer.RData") # changing_cpg

changing_cpg$prepost <- factor(changing_cpg$prepost, levels = c("pre", "post"))
changing_cpg$prepost <- factor(changing_cpg$prepost, levels = c("pre", "post"), labels = c("Pre-lekking", "Post-lekking"))
changing_cpg$id_year <- as.factor(paste0(changing_cpg$id, "_", changing_cpg$year))

load(file="results/modeloutput/changing/modeloutput_glmer.RData") # out_glmer

sub_glmer_prepost <- subset(out_glmer, prepost_qval < 0.05 & abs(mean_delta_meth) >= 0.1)

#### volcano ####

out_glmer <- out_glmer %>% mutate(sig = as.factor(case_when(abs(mean_delta_meth) >= 0.1 & prepost_qval < 0.05 ~ "sig", TRUE ~ "nonsig")))

p_cutoff <- max(out_glmer$prepost_pval[out_glmer$prepost_qval < 0.05])

ggplot(out_glmer, aes(x = mean_delta_meth*100, y = -log10(as.numeric(prepost_pval)))) + 
    geom_point(size=4, alpha=0.5, aes(col = as.factor(sig), fill = as.factor(sig))) +
    labs(x = expression("Mean "*Delta*" methylation %"), y = expression(-log[10]*"("*italic(p*")"))) +
    scale_color_manual(values=c(clrs[5], clr_sig)) +
    scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.5)) +
    geom_hline(yintercept = -log10(p_cutoff), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = -10, col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 10, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") -> fig1_volcano


fig1_volcano

#### manhattan #####

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

out_glmer <- out_glmer %>% mutate(col_sig = case_when(col == "even" & sig == "sig"~ "even_sig",
                                                      col == "odd"  & sig == "sig" ~ "odd_sig",
                                                      col == "even"  & sig == "nonsig" ~ "even_nonsig",
                                                      col == "odd"  & sig == "nonsig" ~ "odd_nonsig")) 

out_glmer$col_sig <- factor(out_glmer$col_sig, levels=c("odd_sig","even_sig", "odd_nonsig", "even_nonsig"))
                                
out_glmer %>% subset(scaf_nr <= 10) %>% 
  ggplot(aes(x = pos, y = -log10(as.numeric(prepost_pval)))) + 
  geom_point(data = subset(out_glmer, col_sig == "odd_nonsig" & scaf_nr <= 10),
             colour = "#8E8E99", fill = alpha("#8E8E99", 0.8), size=3, shape=21) +
  geom_point(data = subset(out_glmer, col_sig == "even_nonsig" & scaf_nr <= 10),
             colour = "#C9C9CF", fill = alpha("#C9C9CF", 0.8), size=3, shape=21) +
  geom_point(data = subset(out_glmer, col_sig == "odd_sig" & scaf_nr <= 10), 
             colour = "#D34E38", fill = alpha("#D34E38", 0.8), size=3, shape=21) +
  geom_point(data = subset(out_glmer, col_sig == "even_sig" & scaf_nr <= 10),
             colour = "#E5988A", fill = alpha("#E5988A", 0.8), size=3, shape=21) +
  facet_grid(~scaf_nr,scales = 'free_x', space = 'free_x', switch = 'x') +
    labs(x = "Scaffold number", y = expression(-log[10]*"("*italic(p*")")))+
    geom_hline(yintercept = -log10(p_cutoff), col = clr_high, linewidth = 1,linetype="dotted") +
    theme(axis.text.x = element_blank(),
    panel.spacing = unit(0, "lines"),
   axis.line.x = element_blank(),
    legend.position="none",
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank()) -> fig1_manhattan

fig1_manhattan

ggsave(fig1_manhattan, file="plots/test.png", height=10, width=10)

#### Raw data #### 
out_glmer <- out_glmer %>% arrange(prepost_qval)
changing_cpg$methperc <- changing_cpg$methperc*100
cpg_sig_increase <- subset(out_glmer, mean_delta_meth > 0 & sig == "sig") %>% head(n=1)
cpg_sig_increase2 <- subset(out_glmer, mean_delta_meth > 0 & sig == "sig") %>% tail(n=1)
cpg_sig_decrease <- subset(out_glmer, mean_delta_meth < 0 & sig == "sig") %>% head(n=1)
#cpg_sig_decrease2 <- subset(out_glmer, mean_delta_meth < 0 & sig == "sig") %>% slice_sample(n=1) #ScEsiA3_19552__HRSCAF_24415_11239844
cpg_sig_decrease2 <- subset(out_glmer, chr_pos == "ScEsiA3_19552__HRSCAF_24415_11239844")

subset(changing_cpg, chr_pos == cpg_sig_increase$chr_pos) %>%
  arrange(id, year) %>%
  ggplot(., aes(x = prepost, y = methperc))+
  geom_boxplot(linewidth=1, outlier.shape=NA) + 
  geom_path(aes(group = id_year), alpha = 0.8, linewidth=0.8,col = "grey60", position = position_jitter(width = 0.1, seed = 3922)) +
  geom_point(aes(alpha = 0.8, col = prepost), size=4, position = position_jitter(width = 0.1, seed = 3922)) + 
  scale_color_manual(values = clr_prepost)+
  labs(x = "Time period", y = "Methylation %") +
  theme(legend.position="none") +
  ylim(0,100)-> plot_top_cpg_1

subset(changing_cpg, chr_pos == cpg_sig_increase2$chr_pos) %>%
  arrange(id, year) %>%
  ggplot(., aes(x = prepost, y = methperc))+
  geom_boxplot(linewidth=1, outlier.shape=NA) + 
  geom_path(aes(group = id_year), alpha = 0.8, linewidth=0.8,col = "grey60", position = position_jitter(width = 0.1, seed = 3922)) +
  geom_point(aes(alpha = 0.8,  col = prepost), size=4, position = position_jitter(width = 0.1, seed = 3922)) + 
  scale_color_manual(values = clr_prepost)+
  labs(x = "Time period", y = "Methylation %") +
  theme(legend.position="none") +
  ylim(0,100)-> plot_top_cpg_2

subset(changing_cpg, chr_pos == cpg_sig_decrease$chr_pos) %>%
  arrange(id, year) %>%
  ggplot(., aes(x = prepost, y = methperc))+
  geom_boxplot(linewidth=1, outlier.shape=NA) + 
  geom_path(aes(group = id_year), alpha = 0.8, linewidth=0.8, col = "grey60", position = position_jitter(width = 0.1, seed = 3922)) +
  geom_point(aes(alpha = 0.8,col = prepost), size=4, position = position_jitter(width = 0.1, seed = 3922)) + 
  scale_color_manual(values = clr_prepost)+
  labs(x = "Time period", y = "Methylation %") +
  theme(legend.position="none") +
  ylim(0,100)-> plot_top_cpg_3

subset(changing_cpg, chr_pos == cpg_sig_decrease2$chr_pos) %>%
  arrange(id, year) %>%
  ggplot(., aes(x = prepost, y = methperc))+
  geom_boxplot(linewidth=1, outlier.shape=NA) + 
  geom_path(aes(group = id_year), alpha = 0.8, linewidth=0.8,col = "grey60", position = position_jitter(width = 0.1, seed = 3922)) +
  geom_point(aes(alpha = 0.8, col = prepost), size=4, position = position_jitter(width = 0.1, seed = 3922)) + 
  scale_color_manual(values = clr_prepost)+
  labs(x = "Time period", y = "Methylation %") +
  theme(legend.position="none") +
  ylim(0,100)-> plot_top_cpg_4


cowplot::plot_grid(plot_top_cpg_1, plot_top_cpg_2, plot_top_cpg_3, plot_top_cpg_4, 
                    align="hv", axis="lb", ncol=2) -> fig1_raw

ggsave(fig1_raw, file="plots/test.png", height=10, width=10)

### sites per region ###

sum_annotated <- read.csv(file="results/modeloutput/changing/summary_regions_sig_cpgs.csv")
sum_annotated$region <- factor(sum_annotated$region, levels=c("Upstream", "Downstream", "Intron", "Exon", "Promoter", "TSS"))

ggplot(subset(sum_annotated, region != "5' UTR" & region != "3' UTR"), aes(x = region, y = perc)) + 
  geom_bar(stat="identity", position="dodge", aes(fill = model, col = model)) + 
  labs(y="Percentage of CpG sites", x="Region", fill = "Subset")+ coord_flip() + 
  scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.8), label = c("All", "Significantly changing")) +
  scale_color_manual(values=c(clrs[5], clr_sig)) + 
  guides(color = "none")+
  ylim(0, 28) +
  geom_text(aes(label = paste0(round(perc, 1), " %"), x = region, y = perc, group=model), 
            hjust=-0.2, size = 6, position=position_dodge(width=0.9)) +
  theme(legend.position = "top") +
  geom_text(label = "*", x = 1, y = 14.2, hjust=-9, vjust = 0.8, size = 6) + #upstream
  geom_text(label = "*", x = 6, y = 6.9, hjust=-9, size = 6, vjust = 0.8) + #tss
  geom_text(label = "*", x = 5, y = 18.6, hjust=-9, size = 6, vjust = 0.8) + #promo
  geom_text(label = "*", x = 3, y = 15, hjust=-9, size = 6, vjust = 0.8) -> fig1_region #intron
  

ggsave(fig1_region, file="plots/test.png", height=10, width=10)

# sig dif from zero with binomial test?
annotated_wide <- spread(dplyr::select(sum_annotated,-c("n_total", "n")), model, perc)
annotated_wide <- left_join(annotated_wide, subset(sum_annotated, model == "Time period"), by = c("region", "Time period" = "perc"))

binom.test(x = annotated_wide$n[which(annotated_wide$region == "Downstream")], 
            n = annotated_wide$n_total[which(annotated_wide$region == "Downstream")], 
            p = annotated_wide$All[which(annotated_wide$region == "Downstream")]/100) 

binom.test(x = annotated_wide$n[which(annotated_wide$region == "Exon")], 
            n = annotated_wide$n_total[which(annotated_wide$region == "Exon")], 
            p = annotated_wide$All[which(annotated_wide$region == "Exon")]/100)

binom.test(x = annotated_wide$n[which(annotated_wide$region == "Intron")], 
            n = annotated_wide$n_total[which(annotated_wide$region == "Intron")], 
            p = annotated_wide$All[which(annotated_wide$region == "Intron")]/100) # sig

binom.test(x = annotated_wide$n[which(annotated_wide$region == "Promoter")], 
            n = annotated_wide$n_total[which(annotated_wide$region == "Promoter")], 
            p = annotated_wide$All[which(annotated_wide$region == "Promoter")]/100) # sig

binom.test(x = annotated_wide$n[which(annotated_wide$region == "TSS")], 
            n = annotated_wide$n_total[which(annotated_wide$region == "TSS")], 
            p = annotated_wide$All[which(annotated_wide$region == "TSS")]/100) # sig

binom.test(x = annotated_wide$n[which(annotated_wide$region == "Upstream")], 
            n = annotated_wide$n_total[which(annotated_wide$region == "Upstream")], 
            p = annotated_wide$All[which(annotated_wide$region == "Upstream")]/100) # sig

#### assembly fig ####

plot_grid(fig1_volcano, fig1_region, ncol=2, labels=c("a", "b"), align="hv", axis="lb", label_fontface = "plain", label_size = 22) -> fig1_top
plot_grid(fig1_manhattan, fig1_raw, ncol=1, labels=c("c", "d"), label_fontface = "plain", label_size = 22, rel_heights=c(1,2)) -> fig1_bottom

plot_grid(fig1_top, fig1_bottom, ncol=1, align="hv", axis="lb", rel_heights=c(1,2)) -> fig1

ggsave(fig1, file="plots/final/main/fig_1_changing.png", width=15, height=20)
ggsave(fig1, file="plots/final/main/fig_1_changing.pdf", width=15, height=20, device = cairo_pdf)
#ggsave(fig1, file="plots/test.png", height=16, width=18)


