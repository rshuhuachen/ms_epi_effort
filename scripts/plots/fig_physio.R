#### Packages #####
pacman::p_load(tidyverse, cowplot, data.table, ggrepel, emmeans, lme4)

#### Theme ####
source("scripts/plotting_theme.R")

#### Raw data (big) for rerunning models #####

load(file="results/modeloutput/changing/changing_sites_glmer.RData")

### load phenotypic data

load(file = "data/phenotypes/fulldata_complete_epi_withdates.RData")

load("data/phenotypes/pheno_dif_prepost.RData") ## differences in physiology

### methylation difference

load(file = "results/modeloutput/all_sites_deltameth.RData")

delta_meth <- subset(delta_meth, chr_pos %in% changing_cpg$chr_pos)

### combine delta methylation data with site and behaviour info

delta_meth <- left_join(delta_meth, unique(all_pheno_epi[,c("id", "year", "site", "Core")], by = c("id", "year")))

### z-transform the traits before the model
prepost_dif$mass_dif_scl <- scale(prepost_dif$mass_dif)
prepost_dif$microf_dif_scl <- scale(prepost_dif$microf_dif)
prepost_dif$trypa_dif_scl <- scale(prepost_dif$trypa_dif)
prepost_dif$ig_dif_scl <- scale(prepost_dif$ig_dif)
prepost_dif$hct_dif_scl <- scale(prepost_dif$hct_dif)

### combine data with pre-post delta physio numbers
delta_meth <- left_join(delta_meth, unique(prepost_dif[,c("id", "year", "mass_dif", "microf_dif", "trypa_dif", "ig_dif", "hct_dif",
                        "mass_dif_scl", "microf_dif_scl", "trypa_dif_scl", "ig_dif_scl", "hct_dif_scl")], by = c("id", "year")))

#### Fig - physiological changes ####
##### With pre-lekking ######

## mass
load(file="results/modeloutput/physio/out_mass_dif_with_pre.RData")
mass_dif <- data

## microf 
load(file="results/modeloutput/physio/out_microf_dif_with_pre.RData")
microf_dif <- data

## trypa 
load(file="results/modeloutput/physio/out_trypa_dif_with_pre.RData")
trypa_dif <- data
    
## hct 
load(file="results/modeloutput/physio/out_hct_dif_with_pre.RData")
hct_dif <- data

## ig 
load(file="results/modeloutput/physio/out_ig_dif_with_pre.RData")
igg_dif <- data

rm(data)

### Make 5 volcano plots ####

## mass_dif 
labels_mass_dif <- data.frame(chr_pos = c("ScEsiA3_17616__HRSCAF_20466_3956821",
                                        "ScEsiA3_18278__HRSCAF_21663_133779007"),
                            lab = c("CpG H*: downstream from ID1",
                                    "CpG B: exon of ADCK2"))

mass_dif <- mass_dif %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))
mass_dif <- left_join(mass_dif, labels_mass_dif)

ggplot(mass_dif, aes(x = parameter_estimate, y = -log10(parameter_qval))) + 
    geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
    scale_color_manual(values=c(clrs[5], clr_sig)) +
    scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
    labs(x = expression(beta*" estimate "*Delta*" body mass"), y = "-log10(q-value)") +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") +
    xlim(c(-0.3,0.3))+
    geom_label_repel(aes(label = lab, x = parameter_estimate, y = -log10(parameter_qval)), 
              nudge_x = .05, nudge_y = 1, size = 6) -> volcano_mass_dif

ggsave(volcano_mass_dif, file="plots/test.png", width=10, height=10)

#### fitting model raw data ####

list_plot_mass_dif <- list()
for (i in 1:nrow(labels_mass_dif)){
        sub <- subset(delta_meth, chr_pos == labels_mass_dif$chr_pos[[i]] & !is.na(delta_meth) & !is.na(mass_dif_scl))
        
        model <- lmerTest::lmer(delta_meth ~ mass_dif_scl + methperc_pre + (1|site), 
                            data = sub)
                            
        gr <- ref_grid(model, cov.keep= c('mass_dif_scl'))
        predict <- as.data.frame(emmeans(gr, spec="mass_dif_scl", level=0.95))
        
        ggplot(sub, aes(x = mass_dif_scl, y = delta_meth*100)) + 
                geom_ribbon(data= predict, aes(ymin = lower.CL*100, ymax = upper.CL*100, y= NULL), fill= clrs[5], alpha = 0.6) +
                geom_line(data= predict, aes(y = emmean*100), col = "black", linewidth=1.5) + 
                geom_point(size = 7, fill = clr_high, alpha = 0.6, col = clr_high) + 
                geom_hline(yintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
                scale_x_reverse()+
                labs(x = expression("z-transformed "*Delta*" body mass"), y = expression(Delta*" methylation %"), 
                        title = labels_mass_dif$lab[[i]]) -> plot

        ggsave(plot, file= paste0("plots/model_out/physio/raw_mass_dif_", i, ".png"), width=8, height=8)

        list_plot_mass_dif[[i]] <- plot
}

### microf ####

## microf_dif 
labels_microf_dif <- data.frame(chr_pos = c("ScEsiA3_16766__HRSCAF_19082_28300178",
                                        "ScEsiA3_16771__HRSCAF_19097_2512626",
                                        "ScEsiA3_21978__HRSCAF_26928_6268951"),
                            lab = c("CpG I*",
                                    "CpG J: exon of JAM3",
                                    "CpG K*: TSS of PPP1R1B"))

microf_dif <- microf_dif %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))
microf_dif <- left_join(microf_dif, labels_microf_dif)

ggplot(microf_dif, aes(x = parameter_estimate, y = -log(parameter_qval, base=exp(20)))) + 
    geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
    scale_color_manual(values=c(clrs[5], clr_sig)) +
    scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
    labs(x = expression(beta*" estimate "*Delta*italic(" Microfilaria spp.")), y = "-log20(q-value)") +
    geom_hline(yintercept = -log(0.05,base=exp(20)), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") +
    xlim(c(-0.3,0.3))+
    geom_label_repel(aes(label = lab, x = parameter_estimate, y = -log(parameter_qval, base=exp(20))), 
              nudge_x = c(-0.05, 0, 0.05), nudge_y = c(10,10,10), size = 6) -> volcano_microf_dif
#I, J
#volcano_microf_dif

ggsave(volcano_microf_dif, file="plots/test.png", width=10, height=10)

#### fitting model raw data ####

list_plot_microf_dif <- list()
for (i in 1:nrow(labels_microf_dif)){
        sub <- subset(delta_meth, chr_pos == labels_microf_dif$chr_pos[[i]] & !is.na(delta_meth) & !is.na(microf_dif_scl))
        
        model <- lmerTest::lmer(delta_meth ~ microf_dif_scl + methperc_pre + (1|site), 
                            data = sub)
                            
        gr <- ref_grid(model, cov.keep= c('microf_dif_scl'))
        predict <- as.data.frame(emmeans(gr, spec="microf_dif_scl", level=0.95))
        
        ggplot(sub, aes(x = microf_dif_scl, y = delta_meth*100)) + 
                geom_ribbon(data= predict, aes(ymin = lower.CL*100, ymax = upper.CL*100, y= NULL), fill= clrs[5], alpha = 0.6) +
                geom_line(data= predict, aes(y = emmean*100), col = "black", linewidth=1.5) + 
                geom_point(size = 7, fill = clr_high, alpha = 0.6, col = clr_high) + 
                geom_hline(yintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
                labs(x = expression("z-transformed "*Delta*italic(" Microfilaria spp.")), y = expression(Delta*" methylation %"), 
                        title = labels_microf_dif$lab[[i]]) -> plot

        ggsave(plot, file= paste0("plots/model_out/physio/raw_microf_dif_", i, ".png"), width=8, height=8)

        list_plot_microf_dif[[i]] <- plot
}

#### trypa ####
## trypa_dif 

trypa_dif <- trypa_dif %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))

ggplot(trypa_dif, aes(x = parameter_estimate, y = -log10(parameter_qval))) + 
    geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
    scale_color_manual(values=c(clrs[5], clr_sig)) +
    scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
    labs(x = expression(beta*" estimate "*Delta*italic("Trypanosoma spp.")), y = "-log10(q-value)") +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") +
    xlim(c(-0.3,0.3))-> volcano_trypa_dif

ggsave(volcano_trypa_dif, file="plots/final/supp/volcano_trypa_with_pre.png", width=10, height=10)

### HCT ####

## hct_dif 
hct_dif <- hct_dif %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))
labels_hct_dif <- data.frame(chr_pos = c("ScEsiA3_15486__HRSCAF_17393_15878140"),
                            lab = c("CpG M*: promoter of GAL"))

hct_dif <- left_join(hct_dif, labels_hct_dif)

ggplot(hct_dif, aes(x = parameter_estimate, y = -log10(parameter_qval))) + 
    geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
    scale_color_manual(values=c(clrs[5], clr_sig)) +
    scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
    labs(x = expression(beta*" estimate "*Delta*" HCT"), y = "-log10(q-value)") +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") +
    xlim(c(-0.3,0.3))+
    geom_label_repel(aes(label = lab, x = parameter_estimate, y = -log10(parameter_qval)), 
              nudge_x = .05, nudge_y = -5, size = 6) -> volcano_hct_dif

ggsave(volcano_hct_dif, file="plots/test.png", width=10, height=10)

#### fitting model raw data ####

list_plot_hct_dif <- list()
for (i in 1:nrow(labels_hct_dif)){
        sub <- subset(delta_meth, chr_pos == labels_hct_dif$chr_pos[[i]] & !is.na(delta_meth) & !is.na(hct_dif_scl))
        
        model <- lmerTest::lmer(delta_meth ~ hct_dif_scl + methperc_pre + (1|site), 
                            data = sub)
                            
        gr <- ref_grid(model, cov.keep= c('hct_dif_scl'))
        predict <- as.data.frame(emmeans(gr, spec="hct_dif_scl", level=0.95))
        
        ggplot(sub, aes(x = hct_dif_scl, y = delta_meth*100)) + 
                geom_ribbon(data= predict, aes(ymin = lower.CL*100, ymax = upper.CL*100, y= NULL), fill= clrs[5], alpha = 0.6) +
                geom_line(data= predict, aes(y = emmean*100), col = "black", linewidth=1.5) + 
                geom_point(size = 7, fill = clr_high, alpha = 0.6, col = clr_high) + 
                geom_hline(yintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
                labs(x = expression("z-transformed "*Delta*" HCT"), y = expression(Delta*" methylation %"), 
                        title = labels_hct_dif$lab[[i]]) -> plot

        ggsave(plot, file= paste0("plots/model_out/physio/raw_hct_dif_", i, ".png"), width=8, height=8)

        list_plot_hct_dif[[i]] <- plot
}

#### IgG ####

## igg_dif 
igg_dif <- igg_dif %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))
labels_igg_dif <- data.frame(chr_pos = c("ScEsiA3_21979__HRSCAF_26929_1282292",
                                         "ScEsiA3_21979__HRSCAF_26929_1282301"),
                            lab = c("CpG G*: upstream from NFIC",
                                    "CpG L: upstream from NFIC"))

igg_dif <- left_join(igg_dif, labels_igg_dif)

ggplot(igg_dif, aes(x = parameter_estimate, y = -log10(parameter_qval))) + 
    geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
    scale_color_manual(values=c(clrs[5], clr_sig)) +
    scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
    labs(x = expression(beta*" estimate "*Delta*" IgG"), y = "-log10(q-value)") +
    geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
    theme(legend.position="none") +
    xlim(c(-0.3,0.3))+
    geom_label_repel(aes(label = lab, x = parameter_estimate, y = -log10(parameter_qval)), 
              nudge_x = .1, nudge_y = c(1, -1), size = 6) -> volcano_igg_dif

ggsave(volcano_igg_dif, file="plots/test.png", width=10, height=10)

#### fitting model raw data ####

list_plot_igg_dif <- list()
for (i in 1:nrow(labels_igg_dif)){
        sub <- subset(delta_meth, chr_pos == labels_igg_dif$chr_pos[[i]] & !is.na(delta_meth) & !is.na(ig_dif_scl))
        
        model <- lmerTest::lmer(delta_meth ~ ig_dif_scl + methperc_pre + (1|site), 
                            data = sub)
                            
        gr <- ref_grid(model, cov.keep= c('ig_dif_scl'))
        predict <- as.data.frame(emmeans(gr, spec="ig_dif_scl", level=0.95))
        
        ggplot(sub, aes(x = ig_dif_scl, y = delta_meth*100)) + 
                geom_ribbon(data= predict, aes(ymin = lower.CL*100, ymax = upper.CL*100, y= NULL), fill= clrs[5], alpha = 0.6) +
                geom_line(data= predict, aes(y = emmean*100), col = "black", linewidth=1.5) + 
                geom_point(size = 7, fill = clr_high, alpha = 0.6, col = clr_high) + 
                geom_hline(yintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
                labs(x = expression("z-transformed "*Delta*" IgG"), y = expression(Delta*" methylation %"), 
                        title = labels_igg_dif$lab[[i]]) -> plot

        ggsave(plot, file= paste0("plots/model_out/physio/raw_igg_dif_", i, ".png"), width=8, height=8)

        list_plot_igg_dif[[i]] <- plot
}


#### assemble figure

all_raw <- plot_grid(list_plot_mass_dif[[1]], list_plot_mass_dif[[2]], 
                     list_plot_microf_dif[[1]], list_plot_microf_dif[[2]], list_plot_microf_dif[[3]],
                     list_plot_igg_dif[[1]], list_plot_igg_dif[[2]], 
                     list_plot_hct_dif[[1]], ncol=2, align="hv", axis="lb", labels = "auto", label_fontface = "plain", label_size = 22) 
ggsave(all_raw, file="plots/final/supp/raw_physio_with_pre.png", width=16, height=22)

plot_grid(volcano_mass_dif, list_plot_mass_dif[[1]], 
          volcano_microf_dif, list_plot_microf_dif[[3]],
          volcano_igg_dif, list_plot_igg_dif[[2]],
          volcano_hct_dif,list_plot_hct_dif[[1]], 
          ncol=2, labels= c("a)", " ", "b)", " ", "c)", " ", "d)", " "), label_fontface = "plain", label_size = 22) -> fig3

ggsave(fig3, file="plots/final/main/fig_physio_with_pre.png", width=14, height=20)

# ##### WithOUT pre-lekking ######
# 
# ## mass
# load(file="results/modeloutput/physio/out_mass_dif_no_pre.RData")
# mass_dif <- data
# 
# ## microf 
# load(file="results/modeloutput/physio/out_microf_dif_no_pre.RData")
# microf_dif <- data
# 
# ## trypa 
# load(file="results/modeloutput/physio/out_trypa_dif_no_pre.RData")
# trypa_dif <- data
#     
# ## hct 
# load(file="results/modeloutput/physio/out_hct_dif_no_pre.RData")
# hct_dif <- data
# 
# ## ig 
# load(file="results/modeloutput/physio/out_ig_dif_no_pre.RData")
# igg_dif <- data
# 
# rm(data)
# 
# ### Make 5 volcano plots ####
# 
# ## mass_dif 
# mass_dif <- mass_dif %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))
# labels_mass_dif <- data.frame(chr_pos = c("ScEsiA3_15486__HRSCAF_17393_53786697",
#                                         "ScEsiA3_16858__HRSCAF_19404_1346",
#                                         "ScEsiA3_17616__HRSCAF_20466_3956821",
#                                         "ScEsiA3_18752__HRSCAF_22883_4031831",
#                                         "ScEsiA3_21976__HRSCAF_26926_4403331"),
#                             lab = c("CpG N: intron of Six4",
#                                     "CpG O: unannotated", 
#                                     "CpG H*: upstream of ID1",
#                                     "CpG P: intron of TET3",
#                                     "CpG Q: downstream from KIAA0319L"))
# 
# mass_dif <- left_join(mass_dif, labels_mass_dif)
# 
# ggplot(mass_dif, aes(x = parameter_estimate, y = -log10(parameter_qval))) + 
#     geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
#     scale_color_manual(values=c(clrs[5], clr_sig)) +
#     scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
#     labs(x = expression(paste(beta, " estimate")), y = "-log10(q-value)", title = expression(Delta~" body mass")) +
#     geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
#     geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
#     theme(legend.position="none") +
#     xlim(c(-0.3,0.3))+
#     geom_label_repel(aes(label = lab, x = parameter_estimate, y = -log10(parameter_qval)), 
#               nudge_x = .09, nudge_y = 2, size = 6) -> volcano_mass_dif
# 
# ggsave(volcano_mass_dif, file="plots/test.png", width=10, height=10)
# 
# #### fitting model raw data ####
# 
# list_plot_mass_dif <- list()
# for (i in 1:nrow(labels_mass_dif)){
#         sub <- subset(delta_meth, chr_pos == labels_mass_dif$chr_pos[[i]] & !is.na(delta_meth) & !is.na(mass_dif_scl))
#         
#         model <- lmerTest::lmer(delta_meth ~ mass_dif_scl + methperc_pre + (1|site), 
#                             data = sub)
#                             
#         gr <- ref_grid(model, cov.keep= c('mass_dif_scl'))
#         predict <- as.data.frame(emmeans(gr, spec="mass_dif_scl", level=0.95))
#         
#         ggplot(sub, aes(x = mass_dif_scl, y = delta_meth)) + 
#                 geom_ribbon(data= predict, aes(ymin = lower.CL, ymax = upper.CL, y= NULL), fill= clrs[5], alpha = 0.6) +
#                 geom_line(data= predict, aes(y = emmean), col = "black", linewidth=1.5) + 
#                 geom_point(size = 7, fill = clr_high, alpha = 0.6, col = clr_high) + 
#                 labs(x = expression("z-transformed "*Delta*" body mass"), y = expression(Delta*" methylation %"), 
#                         title = labels_mass_dif$lab[[i]]) -> plot
# 
#         ggsave(plot, file= paste0("plots/model_out/physio/raw_no_pre_mass_dif_", i, ".png"), width=8, height=8)
# 
#         list_plot_mass_dif[[i]] <- plot
# }
# 
# ### microf ####
# 
# ## microf_dif 
# microf_dif <- microf_dif %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))
# labels_microf_dif <- data.frame(chr_pos = c("ScEsiA3_17655__HRSCAF_20552_15684322"),
#                             lab = c("CpG R: unannotated"))
# 
# microf_dif <- left_join(microf_dif, labels_microf_dif)
# 
# ggplot(microf_dif, aes(x = parameter_estimate, y = -log10(parameter_qval))) + 
#     geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
#     scale_color_manual(values=c(clrs[5], clr_sig)) +
#     scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
#     labs(x = expression(paste(beta, " estimate")), y = "-log10(q-value)", title = expression(Delta~" Microfilaria spp.")) +
#     geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
#     geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
#     theme(legend.position="none") +
#     xlim(c(-0.3,0.3))+
#     geom_label_repel(aes(label = lab, x = parameter_estimate, y = -log10(parameter_qval)), 
#               nudge_x = 0.05, nudge_y = -5, size = 6) -> volcano_microf_dif
# 
# ggsave(volcano_microf_dif, file="plots/test.png", width=10, height=10)
# 
# #### fitting model raw data ####
# 
# list_plot_microf_dif <- list()
# for (i in 1:nrow(labels_microf_dif)){
#         sub <- subset(delta_meth, chr_pos == labels_microf_dif$chr_pos[[i]] & !is.na(delta_meth) & !is.na(microf_dif_scl))
#         
#         model <- lmerTest::lmer(delta_meth ~ microf_dif_scl + methperc_pre + (1|site), 
#                             data = sub)
#                             
#         gr <- ref_grid(model, cov.keep= c('microf_dif_scl'))
#         predict <- as.data.frame(emmeans(gr, spec="microf_dif_scl", level=0.95))
#         
#         ggplot(sub, aes(x = microf_dif_scl, y = delta_meth)) + 
#                 geom_ribbon(data= predict, aes(ymin = lower.CL, ymax = upper.CL, y= NULL), fill= clrs[5], alpha = 0.6) +
#                 geom_line(data= predict, aes(y = emmean), col = "black", linewidth=1.5) + 
#                 geom_point(size = 7, fill = clr_high, alpha = 0.6, col = clr_high) + 
#                 labs(x = expression("z-transformed "*Delta*" Microfilaria spp."), y = expression(Delta*" methylation %"), 
#                         title = labels_microf_dif$lab[[i]]) -> plot
# 
#         ggsave(plot, file= paste0("plots/pheno/raw_no_pre_microf_dif_", i, ".png"), width=8, height=8)
# 
#         list_plot_microf_dif[[i]] <- plot
# }
# 
# #### trypa ####
# ## trypa_dif 
# 
# trypa_dif <- trypa_dif %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))
# 
# ggplot(trypa_dif, aes(x = parameter_estimate, y = -log10(parameter_qval))) + 
#     geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
#     scale_color_manual(values=c(clrs[5], clr_sig)) +
#     scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
#     labs(x = expression(paste(beta, " estimate")), y = "-log10(q-value)", title = expression(Delta~" Trypanosoma spp.")) +
#     geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
#     geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
#     theme(legend.position="none") +
#     xlim(c(-0.3,0.3))-> volcano_trypa_dif
# 
# ggsave(volcano_trypa_dif, file="plots/test.png", width=10, height=10)
# 
# ### HCT ####
# 
# ## hct_dif 
# hct_dif <- hct_dif %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))
# 
# ggplot(hct_dif, aes(x = parameter_estimate, y = -log10(parameter_qval))) + 
#     geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
#     scale_color_manual(values=c(clrs[5], clr_sig)) +
#     scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
#     labs(x = expression(paste(beta, " estimate")), y = "-log10(q-value)", title = expression(Delta~" HCT")) +
#     geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
#     geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
#     theme(legend.position="none") +
#     xlim(c(-0.3,0.3))-> volcano_hct_dif
# 
# ggsave(volcano_hct_dif, file="plots/test.png", width=10, height=10)
# 
# #### IgG ####
# 
# ## igg_dif 
# igg_dif <- igg_dif %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))
# labels_igg_dif <- data.frame(chr_pos = c("ScEsiA3_21979__HRSCAF_26929_1282292",
#                                          "ScEsiA3_21979__HRSCAF_26929_1282301"),
#                             lab = c("CpG G: upstream of NFIC",
#                                     "CpG M: upstream of NFIC"))
# 
# igg_dif <- left_join(igg_dif, labels_igg_dif)
# 
# ggplot(igg_dif, aes(x = parameter_estimate, y = -log10(parameter_qval))) + 
#     geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
#     scale_color_manual(values=c(clrs[5], clr_sig)) +
#     scale_fill_manual(values=alpha(c(clrs[5], clr_sig), 0.6)) +
#     labs(x = expression(paste(beta, " estimate")), y = "-log10(q-value)", title = expression(Delta~" IgG")) +
#     geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
#     geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
#     theme(legend.position="none") +
#     xlim(c(-0.3,0.3))+
#     geom_label_repel(aes(label = lab, x = parameter_estimate, y = -log10(parameter_qval)), 
#               nudge_x = .05, nudge_y = 0.5, size = 6) -> volcano_igg_dif
# 
# ggsave(volcano_igg_dif, file="plots/test.png", width=10, height=10)
# 
# #### fitting model raw data (same though as with pre-lekking) ####
# 
# list_plot_igg_dif <- list()
# for (i in 1:nrow(labels_igg_dif)){
#         sub <- subset(delta_meth, chr_pos == labels_igg_dif$chr_pos[[i]] & !is.na(delta_meth) & !is.na(ig_dif_scl))
#         
#         model <- lmerTest::lmer(delta_meth ~ ig_dif_scl + methperc_pre + (1|site), 
#                             data = sub)
#                             
#         gr <- ref_grid(model, cov.keep= c('ig_dif_scl'))
#         predict <- as.data.frame(emmeans(gr, spec="ig_dif_scl", level=0.95))
#         
#         ggplot(sub, aes(x = ig_dif_scl, y = delta_meth)) + 
#                 geom_ribbon(data= predict, aes(ymin = lower.CL, ymax = upper.CL, y= NULL), fill= clrs[5], alpha = 0.6) +
#                 geom_line(data= predict, aes(y = emmean), col = "black", linewidth=1.5) + 
#                 geom_point(size = 7, fill = clr_high, alpha = 0.6, col = clr_high) + 
#                 labs(x = expression("z-transformed "*Delta*" IgG"), y = expression(Delta*" methylation %"), 
#                         title = labels_igg_dif$lab[[i]]) -> plot
# 
#         ggsave(plot, file= paste0("plots/pheno/raw_no_pre_igg_dif_", i, ".png"), width=8, height=8)
# 
#         list_plot_igg_dif[[i]] <- plot
# }
# 
# #### assemble figure
# plot_grid(volcano_mass_dif, list_plot_mass_dif[[1]], 
#           volcano_microf_dif, list_plot_microf_dif[[1]],
#            volcano_igg_dif, list_plot_igg_dif[[1]],
#           ncol=2, labels="auto", label_fontface = "plain", label_size = 22) -> fig3
# 
# ggsave(fig3, file="plots/final/main/fig_physio_no_pre.png", width=16, height=24)
