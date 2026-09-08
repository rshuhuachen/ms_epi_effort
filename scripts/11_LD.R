### here we check if dynamic cpg sites clustering in single genes are correlated to each other ##

load(file = "results/modeloutput/changing/gene_ids_sig_changing_similar.RData") #annotated changing cpg sites
load(file = "results/modeloutput/all_sites_deltameth.RData") #delta meth data
source("scripts/plotting_theme.R")
pacman::p_load(ggcorrplot, tidyverse)

plot_gene_cor <- function(gene_name, n_cpgs){
  
  gene <- subset(annotated_changing, similar == gene_name)
  
  cpgs <- subset(delta_meth, chr_pos %in% gene$chr_pos)
  
  cpgs_w <- cpgs %>%
    pivot_wider(
      id_cols = c(id, year, age),
      names_from = chr_pos,
      values_from = delta_meth
    )
  
  names(cpgs_w) <- c("id", "year", "age",
                     paste0("CpG ", seq_len(n_cpgs)))
  
  # Select only CpG columns
  cpg_data <- cpgs_w %>%
    dplyr::select(-id, -year, -age)
  
  # Spearman correlations + p-values
  res <- rcorr(
    as.matrix(cpg_data),
    type = "spearman"
  )
  
  ggcorrplot(
    res$r,
    p.mat = res$P,
    type = "upper",
    lab = TRUE,
    lab_size = 4,
    insig = "stars",
    sig.level = 0.05,
    outline.color = "white",
    colors = c("#284651", "#efefef", "#d34e38")) + 
    guides(fill = guide_legend(title = "Correlation")) +
             theme(text=element_text(size=18, family = "Arial"),
                                                       legend.text =  element_text(size = 14, family = "Arial"),
                                                       legend.title = element_text(size = 16, family = "Arial"),
                                                       plot.margin = margin(1,1,1,1, "cm"), 
                                                       panel.background = element_rect(fill = "white", colour = NA),
                                                       plot.background = element_rect(fill = "white", colour = NA),
                                                       panel.grid.major = element_blank(),
                                                       panel.grid.minor = element_blank(),
                                                       panel.border = element_blank()
                                                       )
}

plot_gene_cor("MAB21L2", 6) -> cor_matrix_gene1
plot_gene_cor("BEST1", 4) -> cor_matrix_gene2
plot_gene_cor("HES1-B", 6) -> cor_matrix_gene3

plot_grid(cor_matrix_gene1, cor_matrix_gene2, cor_matrix_gene3, ncol = 2, 
          align="hv", axis="lb", labels="auto", label_fontface = "plain", label_size = 22 ) -> ld_plots

ggsave(ld_plots, file = "plots/final/supp/sfig_5_cor_matrix.png", width = 12, height = 10)

