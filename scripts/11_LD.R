### here we check if dynamic cpg sites clustering in single genes are correlated to each other ##

load(file = "results/modeloutput/changing/gene_ids_sig_changing_similar.RData") #annotated changing cpg sites
load(file = "results/modeloutput/all_sites_deltameth.RData") #delta meth data
source("scripts/plotting_theme.R")
pacman::p_load(ggcorrplot, tidyverse)

# extract gene 1 with six dynamic cpg sites

gene_1 <- subset(annotated_changing, similar == "MAB21L2")
cpgs_gene1 <- subset(delta_meth, chr_pos %in% gene_1$chr_pos)

cpgs_gene1_w <- cpgs_gene1 %>%
  pivot_wider(
    id_cols = c(id, year, age),
    names_from = chr_pos,
    values_from = delta_meth
  )

names(cpgs_gene1_w) <- c("id", "year", "age",
                         "CpG 1", "CpG 2", "CpG 3", "CpG 4", "CpG 5", "CpG 5")

cpg_data_gene1 <- cpgs_gene1_w[, 4:ncol(cpgs_gene1_w)]

cor_gene1 <- cor(cpg_data_gene1, method = "spearman", use = "pairwise.complete.obs")
res_gene1 <- rcorr(
  as.matrix(cpg_data_gene1),
  type = "spearman"
)

ggcorrplot(
  res_gene1$r,
  p.mat = res_gene1$P,
  sig.level = 0.05,
  insig = "blank",     
  type = "upper",
  lab = TRUE,
  lab_size = 3,
  outline.color = "white",
  colors = c("#284651", "white", "#d34e38")
)

# extract gene 2 with four dynamic cpg sites

gene_2 <- subset(annotated_changing, similar == "BEST1")
cpgs_gene2 <- subset(delta_meth, chr_pos %in% gene_2$chr_pos)

cpgs_gene2_w <- cpgs_gene2 %>%
  pivot_wider(
    id_cols = c(id, year, age),
    names_from = chr_pos,
    values_from = delta_meth
  )

names(cpgs_gene2_w) <- c("id", "year", "age",
                         "cpg1", "cpg2", "cpg3", "cpg4")

cpg_data_gene2 <- cpgs_gene2_w[, 4:ncol(cpgs_gene2_w)]

cor_gene2 <- cor(
  cpg_data_gene2,
  method = "spearman",
  use = "pairwise.complete.obs"
)

ggcorrplot(
  cor_gene2,
  type = "upper",
  lab = TRUE,
  lab_size = 3,
  outline.color = "white",
  colors = c("#2166AC", "white", "#B2182B")
)
# extract gene 3 with six dynamic cpg sites

gene_3 <- subset(annotated_changing, similar == "HES1-B")
cpgs_gene3 <- subset(delta_meth, chr_pos %in% gene_3$chr_pos)

cpgs_gene3_w <- cpgs_gene3 %>%
  pivot_wider(
    id_cols = c(id, year, age),
    names_from = chr_pos,
    values_from = delta_meth
  )

names(cpgs_gene3_w) <- c("id", "year", "age",
                         "cpg1", "cpg2", "cpg3", "cpg4", "cpg5", "cpg6")

cpg_data_gene3 <- cpgs_gene3_w[, 4:ncol(cpgs_gene3_w)]

cor_gene3 <- cor(
  cpg_data_gene3,
  method = "spearman",
  use = "pairwise.complete.obs"
)

ggcorrplot(
  cor_gene3,
  type = "upper",
  lab = TRUE,
  lab_size = 3,
  outline.color = "white",
  colors = c("#2166AC", "white", "#B2182B")
)


library(dplyr)
library(tidyr)
library(Hmisc)
library(ggcorrplot)
library(ggplot2)

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
    colors = c("#284651", "white", "#d34e38"))
}

plot_gene_cor("MAB21L2", 6) -> cor_matrix_gene1
plot_gene_cor("BEST1", 4) -> cor_matrix_gene2
plot_gene_cor("HES1-B", 6) -> cor_matrix_gene3

plot_grid(cor_matrix_gene1, cor_matrix_gene2, cor_matrix_gene3, ncol = 1, 
          align="hv", axis="lb", labels="auto", label_fontface = "plain", label_size = 22 ) -> ld_plots

ggsave(ld_plots, file = "plots/final/supp/sfig_5_cor_matrix.png", width = 8, height = 14)

