#### what are the effects of pre-lekking methylation? ####

#### Packages #####
pacman::p_load(tidyverse, cowplot, data.table, ggrepel, emmeans, lme4)

#### Theme ####
source("scripts/plotting_theme.R")

#### Data ####

## attendance
load(file="results/modeloutput/effort/out_attend_with_pre.RData")
attend <- data

## fight 
load(file="results/modeloutput/effort/out_fight_with_pre.RData")
fight <- data

## dist 
load(file="results/modeloutput/effort/out_dist_with_pre.RData")
dist <- data

## mass 
load(file="results/modeloutput/physio/out_mass_dif_with_pre.RData")
mass <- data

## microf 
load(file="results/modeloutput/physio/out_microf_dif_with_pre.RData")
microf <- data

## trypa 
load(file="results/modeloutput/physio/out_trypa_dif_with_pre.RData")
trypa <- data

## hct 
load(file="results/modeloutput/physio/out_hct_dif_with_pre.RData")
hct <- data

## ig 
load(file="results/modeloutput/physio/out_ig_dif_with_pre.RData")
ig <- data

rm(data)

#### fdr correction, n_sig per model ####

models <- list(attend, fight, dist, mass, microf, trypa, hct, ig)
sig_pre <- list()
for (i in 1:length(models)){
    df <- as.data.frame(models[[i]])
    df$pre_qval <- p.adjust(df$pre_pval, method = "fdr", n = nrow(df))
    models[[i]] <- df
    sig <- subset(df, pre_qval < 0.05)
    sig_pre[[i]] <- sig
}

for (i in 1:length(sig_pre)){print(nrow(sig_pre[[i]]))}
for (i in 1:length(sig_pre)){print(nrow(models[[i]]))}

# many are sig! 
### volcano plots pre-meth ####
plots <- list()
for (i in 1:length(models)){
    df <- as.data.frame(models[[i]])
    df <- df %>% mutate(sig = case_when(pre_qval < 0.05 ~ "sig", TRUE ~ "nonsig"))
    
    ggplot(df, aes(x = pre_estimate, y = -log10(pre_qval))) + 
        geom_point(size=7, alpha=0.5, aes(col = sig, fill = sig)) +
        scale_color_manual(values=c(clr_grey, clr_sig)) +
        scale_fill_manual(values=alpha(c(clr_grey, clr_sig), 0.6)) +
        labs(x = expression(paste(beta, " estimate pre-lekking methylation %")), y = "-log10(q-value)", title = as.character(df$parameter[1])) +
        geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dotted", linewidth = 1) +
        geom_vline(xintercept = 0, col = "darkred", linetype = "dotted", linewidth = 1) +
        theme(legend.position="none") -> volcano
    plots[[i]] <- volcano        

    ggsave(volcano, file= paste0("plots/model_out/pre_lekking/volcano_", as.character(df$parameter[1]), ".png"), width=10, height=10)}

### combine in plot
plot_grid(plots[[1]] + labs(title = "Attendance"), 
          plots[[2]] + labs(title = "Fighting"), 
          plots[[3]] + labs(title = "Lek centrality"), 
          plots[[4]] + labs(title = expression(Delta~body~mass)),
          plots[[5]] + labs(title = expression(Delta~italic(Microfilaria~spp.))),
          plots[[6]] + labs(title = expression(Delta~italic(Trypanosoma~spp.))),
          plots[[8]] + labs(title = expression(Delta~IgG)),
          plots[[7]] + labs(title = expression(Delta~HCT)),
          ncol=2, align="hv", axis="lb", labels="auto", label_fontface = "plain", label_size = 22 ) -> all_plots
all_plots
ggsave(all_plots, file="plots/model_out/pre_lekking/combined_volcanoes.png", width=16, height=20)            