#### Packages #####
pacman::p_load(tidyverse, cowplot, data.table, ggrepel, emmeans, lme4)

#### Theme ####
source("scripts/plotting_theme.R")

#### Raw data (big) for rerunning models #####

load(file="results/modeloutput/changing/changing_sites_glmer.RData")

### load phenotypic data

load(file = "data/phenotypes/fulldata_complete_epi_withdates.RData")

### methylation difference

load(file = "results/modeloutput/all_sites_deltameth.RData")

delta_meth <- subset(delta_meth, chr_pos %in% changing_cpg$chr_pos)

### combine delta methylation data with site and behaviour info

delta_meth <- left_join(delta_meth, unique(all_pheno_epi[,c("id", "year", "site", "Core", "surv", "MS")], by = c("id", "year")))
  
### load model output fitness ####
## AMS
load(file="results/modeloutput/fitness/out_ams_deltameth_filtered.RData")
## surv 
load(file="results/modeloutput/fitness/out_surv_deltameth_filtered.RData")

### Volcano plot ####
source("scripts/plotting_theme.R")

#ams
delta_out_ams <- delta_out_ams %>% mutate(sig = case_when(ams_delta_meth_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))

ggplot(delta_out_ams, aes(x = ams_delta_meth_estimate, y = -log10(ams_delta_meth_qval))) + 
    geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
    scale_color_manual(values=c(clrs[5], clr_sig)) +
    scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
    xlim(-21, 21)+
    labs(x = expression("Estimate "*Delta*" methylation"), y = "-log10(q-value)", title = "Annual mating success") +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") -> volcano_ams

ggsave(plot, file = "plots/test.png", width=10, height=10)    

# raw data top 5

cpg_sig_ams <- subset(delta_out_ams, ams_delta_meth_qval < 0.05) #362

### plotting
cpg_sig_ams$intercept_ams <- as.numeric(cpg_sig_ams$intercept_ams)
cpg_sig_ams <- cpg_sig_ams %>% arrange(ams_delta_meth_qval)

list_plot_ams <- list()
for (i in 1:2){
        sub <- subset(delta_meth, chr_pos == cpg_sig_ams$chr_pos[[i]] & !is.na(delta_meth) & !is.na(MS))
        
        model <- glmmTMB::glmmTMB(MS ~ delta_meth + (1|site/id), family = "poisson",
                            data = sub)
                            
        gr <- ref_grid(model, cov.keep= c('delta_meth'))
        predict <- as.data.frame(emmeans(gr, spec="delta_meth", level=0.95))
        
        ggplot(sub, aes(x = delta_meth, y = MS)) + 
                geom_ribbon(data= predict, aes(ymin = asymp.LCL, ymax = asymp.UCL, y= NULL), fill= clrs[5], alpha = 0.6) +
                geom_line(data= predict, aes(y = emmean), col = "black", linewidth=1.5) + 
                geom_point(size = 7, fill = clr_high, alpha = 0.6, col = clr_high) + 
                labs(x = "Annual mating success", x = expression(Delta*" methylation %")) -> plot

        ggsave(plot, file= paste0("plots/model_out/fitness/raw_AMS_", i, ".png"), width=8, height=8)

        list_plot_ams[[i]] <- plot
}
   
# surv
delta_out_surv <- delta_out_surv %>% mutate(sig = case_when(surv_delta_meth_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))

ggplot(delta_out_surv, aes(x = surv_delta_meth_estimate, y = -log10(surv_delta_meth_qval))) + 
    geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
    scale_color_manual(values=c(clrs[5], clr_sig)) +
    scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
    labs(x = expression("Estimate "*Delta*" methylation"), y = "-log10(q-value)", 
         title = "Survival") +
    xlim(-42, 42)+
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") -> volcano_surv

ggsave(volcano_surv, file = "plots/test.png", width=10, height=10)    

### significant ones
delta_out_surv <- delta_out_surv %>% arrange(surv_delta_meth_qval)

list_plot_surv <- list()
for (i in 1:2){
        sub <- subset(delta_meth, chr_pos == delta_out_surv$chr_pos[[i]] & !is.na(delta_meth) & !is.na(surv))
        
        ggplot(sub, aes(x = delta_meth, y = surv)) + 
        geom_boxplot(lwd=1) + 
        geom_point(size = 7, fill = clr_high, alpha = 0.6, col = clr_high, position=position_jitter(height=0.1)) + 
        labs(y = "Survival", x = expression(Delta*" methylation %")) -> plot

        ggsave(plot, file= paste0("plots/model_out/fitness/raw_surv_", i, ".png"), width=8, height=8)

        list_plot_surv[[i]] <- plot
}

ggsave(plot, file = "plots/test.png", width=10, height=10)    

### for AMS: manhattan ###

# load scaffold numbers
load("data/scaffold_names_dovetail.RData")

# Split the chr_pos column into two columns based on the first "_"
delta_out_ams$chr_pos <- as.character(delta_out_ams$chr_pos)
split_chr_pos <- strsplit(delta_out_ams$chr_pos, "_", fixed = TRUE)

# Extract the numbers following HRSCAF=XXX_number
delta_out_ams$chr <- paste0(sapply(split_chr_pos, "[", 1), "_",
                        sapply(split_chr_pos, "[", 2), ";", 
                        sapply(split_chr_pos, "[", 4), "=",
                        sapply(split_chr_pos, "[", 5))

delta_out_ams$pos <- as.numeric(sapply(split_chr_pos, "[", 6))

# join
delta_out_ams <- left_join(delta_out_ams, genome[,c("contig", "scaf_nr")], by = c("chr" = "contig"))

# plot 
delta_out_ams <- delta_out_ams %>% mutate(col = case_when(scaf_nr %% 2 == 0 ~ "even",
                                                  TRUE ~ "odd"))

delta_out_ams <- delta_out_ams %>% mutate(col_sig = case_when(col == "even" & ams_delta_meth_qval <0.05 ~ "sig_even",
                                                      col == "odd" & ams_delta_meth_qval <0.05 ~ "sig_odd",  
                                                      TRUE ~ "nonsig"))                                        


delta_out_ams %>% subset(scaf_nr <= 10) %>% 
  ggplot(aes(x = pos, y = -log10(as.numeric(ams_delta_meth_qval)))) + 
  geom_point(size=5, alpha=0.5, aes(col = as.factor(col_sig), fill = as.factor(col_sig))) +
  facet_grid(~scaf_nr,scales = 'free_x', space = 'free_x', switch = 'x') +
  labs(x = "Scaffold number", y = expression(-log[10]*"(p-value)")) +
  scale_color_manual(values=c(clrs[5], "#E28979", clr_sig)) +
  scale_fill_manual(values=alpha(c(clrs[5], "#E28979", clr_sig), 0.5)) +
  geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1.5) +
  theme(axis.text.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        axis.line.x = element_blank(),
        legend.position="none",
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank()) -> fig_ams_manhattan

fig_ams_manhattan

### assemble #####

plot_grid(volcano_ams, list_plot_ams[[1]], 
          volcano_surv, list_plot_surv[[1]], align="vh", axis="lb",
          ncol=2, labels="auto", label_fontface = "plain", label_size = 22) -> fig

plot_grid(volcano_ams, 
          volcano_surv, align="vh", axis="lb",
          ncol=2, labels="auto", label_fontface = "plain", label_size = 22) -> fig

ggsave(fig, file="plots/final/main/fig_fitness.png", width=14, height=8)
