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

ggplot(out_glmer, aes(x = mean_delta_meth, y = -log10(as.numeric(prepost_qval)))) + 
    geom_point(size=4, alpha=0.5, aes(col = as.factor(sig), fill = as.factor(sig))) +
    labs(x = expression("Mean "*Delta*" methylation %"), y = "-log10(q-value)") +
    scale_color_manual(values=c(clrs[5], clrs[17])) +
    scale_fill_manual(values=alpha(c(clrs[5], clrs[17]), 0.5)) +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = -0.1, col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0.1, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") -> fig1_volcano

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

out_glmer <- out_glmer %>% mutate(col_sig = case_when(col == "even" & prepost_qval <0.05 ~ "sig_even",
                                                      col == "odd" & prepost_qval <0.05 ~ "sig_odd",  
                                                      TRUE ~ "nonsig"))                                        

                                
out_glmer %>% subset(scaf_nr <= 15) %>% 
  ggplot(aes(x = pos, y = -log10(as.numeric(prepost_qval)))) + 
    geom_point(size=5, alpha=0.5, aes(col = as.factor(col_sig), fill = as.factor(col_sig))) +
    facet_grid(~scaf_nr,scales = 'free_x', space = 'free_x', switch = 'x') +
    labs(x = "Scaffold number", y = expression(-log[10]*"(p-value)")) +
    scale_color_manual(values=c(clrs[5], clrs[17], clrs[1])) +
    scale_fill_manual(values=alpha(c(clrs[5], clrs[17], clrs[1]), 0.5)) +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(axis.text.x = element_blank(),
    panel.spacing = unit(0, "lines"),
   axis.line.x = element_blank(),
    legend.position="none",
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank()) -> fig1_manhattan

ggsave(fig1_manhattan, file="plots/test.png", height=10, width=10)

#### Raw data #### 
out_glmer <- out_glmer %>% arrange(prepost_qval)

subset(changing_cpg, chr_pos == out_glmer$chr_pos[1]) %>%
  arrange(id, year) %>%
  ggplot(., aes(x = prepost, y = methperc))+
  geom_boxplot(linewidth=1, outlier.shape=NA) + 
  geom_path(aes(group = id_year), alpha = 0.8, linewidth=0.8,col = "grey60", position = position_jitter(width = 0.1, seed = 3922)) +
  geom_point(aes(alpha = 0.8, col = prepost), size=4, position = position_jitter(width = 0.1, seed = 3922)) + 
  scale_color_manual(values = c(clrs[1], clrs[17]))+
  labs(x = "Time period", y = "Methylation percentage") +
  theme(legend.position="none") +
  ylim(0,1)-> plot_top_cpg_1

subset(changing_cpg, chr_pos == out_glmer$chr_pos[2]) %>%
  arrange(id, year) %>%
  ggplot(., aes(x = prepost, y = methperc))+
  geom_boxplot(linewidth=1, outlier.shape=NA) + 
  geom_path(aes(group = id_year), alpha = 0.8, linewidth=0.8,col = "grey60", position = position_jitter(width = 0.1, seed = 3922)) +
  geom_point(aes(alpha = 0.8,  col = prepost), size=4, position = position_jitter(width = 0.1, seed = 3922)) + 
  scale_color_manual(values = c(clrs[1], clrs[17]))+
  labs(x = "Time period", y = "Methylation percentage") +
  theme(legend.position="none") +
  ylim(0,1)-> plot_top_cpg_2

subset(changing_cpg, chr_pos == out_glmer$chr_pos[3]) %>%
  arrange(id, year) %>%
  ggplot(., aes(x = prepost, y = methperc))+
  geom_boxplot(linewidth=1, outlier.shape=NA) + 
  geom_path(aes(group = id_year), alpha = 0.8, linewidth=0.8, col = "grey60", position = position_jitter(width = 0.1, seed = 3922)) +
  geom_point(aes(alpha = 0.8,col = prepost), size=4, position = position_jitter(width = 0.1, seed = 3922)) + 
  scale_color_manual(values = c(clrs[1], clrs[17]))+
  labs(x = "Time period", y = "Methylation percentage") +
  theme(legend.position="none") +
  ylim(0,1)-> plot_top_cpg_3

subset(changing_cpg, chr_pos == out_glmer$chr_pos[4]) %>%
  arrange(id, year) %>%
  ggplot(., aes(x = prepost, y = methperc))+
  geom_boxplot(linewidth=1, outlier.shape=NA) + 
  geom_path(aes(group = id_year), alpha = 0.8, linewidth=0.8,col = "grey60", position = position_jitter(width = 0.1, seed = 3922)) +
  geom_point(aes(alpha = 0.8, col = prepost), size=4, position = position_jitter(width = 0.1, seed = 3922)) + 
  scale_color_manual(values = c(clrs[1], clrs[17]))+
  labs(x = "Time period", y = "Methylation percentage") +
  theme(legend.position="none") +
  ylim(0,1)-> plot_top_cpg_4


cowplot::plot_grid(plot_top_cpg_1, plot_top_cpg_2, plot_top_cpg_3, plot_top_cpg_4, 
                    align="hv", axis="lb", ncol=2) -> fig1_raw

ggsave(fig1_raw, file="plots/test.png", height=10, width=10)

### sites per region ###

sum_annotated <- read.csv(file="results/modeloutput/changing/summary_regions_sig_cpgs.csv")

ggplot(subset(sum_annotated, region != "5' UTR" & region != "3' UTR"), aes(x = region, y = perc)) + 
    geom_bar(stat="identity", position="dodge", aes(fill = model, col = model)) + 
  labs(y="Percentage of CpG sites", x="Region", fill = "Subset")+ coord_flip() + 
  scale_fill_manual(values=alpha(c(clrs[5], clrs[17]), 0.8)) +
  scale_color_manual(values=c(clrs[5], clrs[17])) + 
  guides(color = "none")+
  ylim(0, 28) +
  geom_text(aes(label = paste0(round(perc, 1), " %"), x = region, y = perc, group=model), 
              hjust=-0.2, size = 6, position=position_dodge(width=0.9)) -> fig1_region

ggsave(fig1_region, file="plots/test.png", height=10, width=10)

#### assembly fig ####

plot_grid(fig1_volcano, fig1_region, ncol=2, labels=c("a", "b"), align="hv", axis="lb", label_fontface = "plain", label_size = 22) -> fig1_top
plot_grid(fig1_manhattan, fig1_raw, ncol=1, labels=c("c", "d"), label_fontface = "plain", label_size = 22, rel_heights=c(1,2)) -> fig1_bottom

plot_grid(fig1_top, fig1_bottom, ncol=1, align="hv", axis="lb", rel_heights=c(1,2)) -> fig1

ggsave(fig1, file="plots/final/main/fig_changing.png", width=16, height=18)


