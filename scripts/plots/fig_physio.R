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
                            lab = c("CpG H: upstream of genes BPI and ID1",
                                    "CpG B: in gene body, unknown gene"))

mass_dif <- mass_dif %>% mutate(sig = case_when(parameter_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))
mass_dif <- left_join(mass_dif, labels_mass_dif)

ggplot(mass_dif, aes(x = parameter_estimate, y = -log10(parameter_qval))) + 
    geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
    scale_color_manual(values=c(clrs[5], clrs[17])) +
    scale_fill_manual(values=alpha(c(clrs[5], clrs[17]), 0.6)) +
    labs(x = expression(paste(beta, " estimate")), y = "-log10(q-value)", title = expression(Delta~" body mass")) +
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
        
        ggplot(sub, aes(x = mass_dif_scl, y = delta_meth)) + 
                geom_ribbon(data= predict, aes(ymin = lower.CL, ymax = upper.CL, y= NULL), fill= clrs[5], alpha = 0.6) +
                geom_line(data= predict, aes(y = emmean), col = "black", linewidth=1.5) + 
                geom_point(size = 7, fill = clr_high, alpha = 0.6, col = clr_high) + 
                labs(x = "z-transformed mass_difance", y = expression(Delta*" methylation %"), 
                        title = labels_mass_dif$lab[[i]]) -> plot

        ggsave(plot, file= paste0("plots/pheno/raw_mass_dif_", i, ".png"), width=8, height=8)

        list_plot_mass_dif[[i]] <- plot
}
   