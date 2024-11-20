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

delta_meth <- left_join(delta_meth, unique(all_pheno_epi[,c("id", "year", "site", "Core")], by = c("id", "year")))

### z-transform the traits before the model that subsets IDs and years where there is
### data for that CpG site

effort <- all_pheno_epi %>% dplyr::select(c("id", "year", "attend", "fight", "dist", "MS")) %>% filter(!is.na(attend)) %>% unique()

effort <- subset(effort, id %in% delta_meth$id)

effort$attend_scl <- scale(effort$attend)
effort$fight_scl <- scale(effort$fight)
effort$dist_scl <- scale(effort$dist)

### combine reproductive effort data with methylation data

delta_meth <- left_join(delta_meth, effort[,c("id", "year", "attend", "fight", "dist", "attend_scl", "fight_scl", "dist_scl")], by = c("id", "year"))
  
#### Fig 2 - reproductive effort ####
##### With pre-lekking ######

## attendance
load(file="results/modeloutput/effort/out_attend_with_pre.RData")
attend <- data

## fight 
load(file="results/modeloutput/effort/out_fight_with_pre.RData")
fight <- data

## dist 
load(file="results/modeloutput/effort/out_dist_with_pre.RData")
dist <- data

rm(data)

### Make 3 volcano plots ####

## attend 
labels_attend <- data.frame(chr_pos = c("ScEsiA3_17655__HRSCAF_20552_371333",
                                        "ScEsiA3_18278__HRSCAF_21663_133779007",
                                        "ScEsiA3_18752__HRSCAF_22883_2574288"),
                            lab = c("CpG A: upstream from GPR87",
                                    "CpG B: exon of ADCK2",
                                    "CpG C: downstream from UNK"))

attend <- attend %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))
attend <- left_join(attend, labels_attend)

ggplot(attend, aes(x = parameter_estimate, y = -log10(parameter_qval))) + 
    geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
    scale_color_manual(values=c(clrs[5], clr_sig)) +
    scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
    labs(x = expression(paste(beta, " estimate")), y = "-log10(q-value)", title = "Lek attendance") +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") +
    xlim(c(-0.3,0.3))+
    geom_label_repel(aes(label = lab, x = parameter_estimate, y = -log10(parameter_qval)), 
              nudge_x = .05, nudge_y = 2, size = 6) -> volcano_attend

ggsave(volcano_attend, file="plots/test.png", width=10, height=10)

#### fitting model raw data ####

list_plot_attend <- list()
for (i in 1:nrow(labels_attend)){
        sub <- subset(delta_meth, chr_pos == labels_attend$chr_pos[[i]] & !is.na(delta_meth) & !is.na(attend_scl))
        
        model <- lmerTest::lmer(delta_meth ~ attend_scl + methperc_pre + (1|site), 
                            data = sub)
                            
        gr <- ref_grid(model, cov.keep= c('attend_scl'))
        predict <- as.data.frame(emmeans(gr, spec="attend_scl", level=0.95))
        
        ggplot(sub, aes(x = attend_scl, y = delta_meth)) + 
                geom_ribbon(data= predict, aes(ymin = lower.CL, ymax = upper.CL, y= NULL), fill= clrs[5], alpha = 0.6) +
                geom_line(data= predict, aes(y = emmean), col = "black", linewidth=1.5) + 
                geom_point(size = 7, fill = clr_high, alpha = 0.6, col = clr_high) + 
                labs(x = "z-transformed attendance", y = expression(Delta*" methylation %"), 
                        title = labels_attend$lab[[i]]) -> plot

        ggsave(plot, file= paste0("plots/pheno/raw_attend_", i, ".png"), width=8, height=8)

        list_plot_attend[[i]] <- plot
}
   
## fight 
labels_fight <- data.frame(chr_pos = c("ScEsiA3_15486__HRSCAF_17393_6262304",
                                        "ScEsiA3_18290__HRSCAF_21698_5473558",
                                        "ScEsiA3_21978__HRSCAF_26928_5532732"),
                            lab = c("CpG D: exon of RHOF",
                                    "CpG E: upstream from SLC38A10",
                                    "CpG F: intron of TRAK1"))

fight <- fight %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))
fight <- left_join(fight, labels_fight)

ggplot(fight, aes(x = parameter_estimate, y = -log10(parameter_qval))) + 
    geom_point(size=6, alpha=0.5, aes(col = sig, fill = sig)) +
    labs(x = expression(paste(beta, " estimate")), y = "-log10(q-value)", title = "Fighting rate") +
    scale_color_manual(values=c(clrs[5], clr_sig)) +
    scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") +
    xlim(c(-0.3,0.3))+
    geom_label_repel(aes(label = lab, x = parameter_estimate, y = -log10(parameter_qval)), 
              nudge_x = .05, nudge_y = c(40,14,40), size = 6) -> volcano_fight

ggsave(volcano_fight, file="plots/test.png", width=10, height=10)

list_plot_fight <- list()
for (i in 1:nrow(labels_fight)){
        sub <- subset(delta_meth, chr_pos == labels_fight$chr_pos[[i]] & !is.na(delta_meth) & !is.na(fight_scl))
        
        model <- lmerTest::lmer(delta_meth ~ fight_scl + methperc_pre + (1|site), 
                            data = sub)
                            
        gr <- ref_grid(model, cov.keep= c('fight_scl'))
        predict <- as.data.frame(emmeans(gr, spec="fight_scl", level=0.95))
        
        ggplot(sub, aes(x = fight_scl, y = delta_meth)) + 
                geom_ribbon(data= predict, aes(ymin = lower.CL, ymax = upper.CL, y= NULL), fill= clrs[5], alpha = 0.6) +
                geom_line(data= predict, aes(y = emmean), col = "black", linewidth=1.5) + 
                geom_point(size = 7, fill = clr_high, alpha = 0.6, col = clr_high) + 
                labs(x = "z-transformed fighting rate", y = expression(Delta*" methylation %"), 
                        title = labels_fight$lab[[i]]) -> plot

        ggsave(plot, file= paste0("plots/pheno/raw_fight_", i, ".png"), width=8, height=8)

        list_plot_fight[[i]] <- plot
}

   
## dist 
labels_dist <- data.frame(chr_pos = c("ScEsiA3_21979__HRSCAF_26929_1282292"),
                            lab = c("CpG G: upstream from NFIC"))

dist <- dist %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))
dist <- left_join(dist, labels_dist)

ggplot(dist, aes(x = parameter_estimate, y = -log10(parameter_qval))) + 
    geom_point(size=6, alpha=0.5, aes(col = sig, fill = sig)) +
    labs(x = expression(paste(beta, " estimate")), y = "-log10(q-value)", title = "Lek centrality") +
    scale_color_manual(values=c(clrs[5], clr_sig)) +
    scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") +
    xlim(c(-0.3,0.3))+
    geom_label_repel(aes(label = lab, x = parameter_estimate, y = -log10(parameter_qval)), 
              nudge_x = .07, nudge_y = -1, size = 6) -> volcano_dist

ggsave(volcano_dist, file="plots/test.png", width=10, height=10)

list_plot_dist <- list()
for (i in 1:nrow(labels_dist)){
        sub <- subset(delta_meth, chr_pos == labels_dist$chr_pos[[i]] & !is.na(delta_meth) & !is.na(dist_scl))
        
        model <- lmerTest::lmer(delta_meth ~ dist_scl + methperc_pre + (1|site), 
                            data = sub)
                            
        gr <- ref_grid(model, cov.keep= c('dist_scl'))
        predict <- as.data.frame(emmeans(gr, spec="dist_scl", level=0.95))
        
        ggplot(sub, aes(x = dist_scl, y = delta_meth)) + 
                geom_ribbon(data= predict, aes(ymin = lower.CL, ymax = upper.CL, y= NULL), fill= clrs[5], alpha = 0.6) +
                geom_line(data= predict, aes(y = emmean), col = "black", linewidth=1.5) + 
                geom_point(size = 7, fill = clr_high, alpha = 0.6, col = clr_high) + 
                labs(x = "z-transformed centrality", y = expression(Delta*" methylation %"), 
                        title = labels_dist$lab[[i]]) -> plot

        ggsave(plot, file= paste0("plots/pheno/raw_dist_", i, ".png"), width=8, height=8)

        list_plot_dist[[i]] <- plot
}


#### assembly fig ####

volcano_attend
volcano_fight
volcano_dist
list_plot_attend[[1]]
list_plot_attend[[2]]
list_plot_attend[[3]]
list_plot_fight[[1]]
list_plot_fight[[2]]
list_plot_fight[[3]]
list_plot_dist[[1]]

all_raw <- plot_grid(list_plot_attend[[1]], list_plot_attend[[2]], list_plot_attend[[3]], 
                        list_plot_fight[[1]], list_plot_fight[[2]], list_plot_fight[[3]],
                        list_plot_dist[[1]], ncol=3, align="hv", axis="lb") 
ggsave(all_raw, file="plots/final/supp/raw_effort_with_pre.png", width=26, height=24)

plot_grid(volcano_attend, list_plot_attend[[1]], 
          volcano_fight, list_plot_fight[[2]],
          volcano_dist, list_plot_dist[[1]],
          ncol=2, labels="auto", label_fontface = "plain", label_size = 22) -> fig2

ggsave(fig2, file="plots/final/main/fig_effort_with_pre.png", width=16, height=20)

##### WithOUT pre-lekking ######

## attendance
load(file="results/modeloutput/effort/out_attend_no_pre.RData")
attend <- data

## fight 
load(file="results/modeloutput/effort/out_fight_no_pre.RData")
fight <- data

## dist 
load(file="results/modeloutput/effort/out_dist_no_pre.RData")
dist <- data

rm(data)

### Make 3 volcano plots ####

## attend 
attend <- attend %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))

ggplot(attend, aes(x = parameter_estimate, y = -log10(parameter_qval))) + 
    geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
    scale_color_manual(values=c(clrs[5], clr_sig)) +
    scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
    labs(x = expression(paste(beta, " estimate")), y = "-log10(q-value)", title = "Lek attendance") +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") +
    xlim(c(-0.3,0.3)) -> volcano_attend

ggsave(volcano_attend, file="plots/test.png", width=10, height=10)

## no sig so no raw data

## fight 

fight <- fight %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))


ggplot(fight, aes(x = parameter_estimate, y = -log10(parameter_qval))) + 
    geom_point(size=6, alpha=0.5, aes(col = sig, fill = sig)) +
    labs(x = expression(paste(beta, " estimate")), y = "-log10(q-value)", title = "Fighting rate") +
    scale_color_manual(values=c(clrs[5], clr_sig)) +
    scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") +
    xlim(c(-0.3,0.3)) -> volcano_fight

ggsave(volcano_fight, file="plots/test.png", width=10, height=10)

## no sig so no raw data
   
## dist 
labels_dist <- data.frame(chr_pos = c("ScEsiA3_21979__HRSCAF_26929_1282292"),
                            lab = c("CpG G: upstream of NFIC"))

dist <- dist %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))
dist <- left_join(dist, labels_dist)

ggplot(dist, aes(x = parameter_estimate, y = -log10(parameter_qval))) + 
    geom_point(size=6, alpha=0.5, aes(col = sig, fill = sig)) +
    labs(x = expression(paste(beta, " estimate")), y = "-log10(q-value)", title = "Lek centrality") +
    scale_color_manual(values=c(clrs[5], clr_sig)) +
    scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") +
    xlim(c(-0.3,0.3))+
    geom_label_repel(aes(label = lab, x = parameter_estimate, y = -log10(parameter_qval)), 
              nudge_x = .07, nudge_y = -1, size = 6) -> volcano_dist

ggsave(volcano_dist, file="plots/test.png", width=10, height=10)

## same sig as the one with pre-lekking so don't plot here

#### assembly fig ####

plot_grid(volcano_attend, 
          volcano_fight, 
          volcano_dist, 
          ncol=1, labels="auto", label_fontface = "plain", label_size = 22) -> fig2_no_pre

ggsave(fig2_no_pre, file="plots/final/main/fig_effort_no_pre.png", width=12, height=20)
