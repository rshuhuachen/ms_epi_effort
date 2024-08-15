### Annotate all significant CpG sites

### Packages ####
pacman::p_load(tidyverse, data.table, genomation, GenomicFeatures, rtracklayer, 
               GenomicRanges)

### Plotting ###
source("scripts/plotting_theme.R")

### Significants sites ####
# changing
load(file="results/modeloutput/prepost_modeloutput_glmer_min0.75.RData")
cpg_change <- subset(out_glmer, prepost_qval < 0.05 & abs(prepost_estimate) > 0.1)

# AMS
load(file="results/modeloutput/AMS_deltameth_modeloutput_filtered.RData")
cpg_ams <- subset(delta_out_ams, ams_delta_meth_qval < 0.05 & abs(ams_delta_meth_estimate) > 0.1)

# effort
load(file="results/modeloutput/effort_deltameth_modeloutput_filtered.RData")
cpg_effort <- subset(delta_out_all, parameter_qval < 0.05 & abs(parameter_estimate) > 0.1)

#physio
load(file="results/modeloutput/physio_deltameth_modeloutput_filtered.RData")
cpg_physio <- subset(delta_out_all, parameter_qval < 0.05 & abs(parameter_estimate) > 0.1)



### combine
cpg_all <- out_glmer %>% dplyr::select(c(chr_pos, prepost_qval))
names(cpg_all)[2] <- "parameter_qval"
cpg_all$parameter <- "all"

cpg_changing_select <- cpg_change %>% dplyr::select(c(chr_pos, prepost_qval))
names(cpg_changing_select)[2] <- "parameter_qval"
cpg_changing_select$parameter <- "time_period"

cpg_ams_select <- cpg_ams %>% dplyr::select(c(chr_pos, ams_delta_meth_qval))
names(cpg_ams_select)[2] <- "parameter_qval"
cpg_ams_select$parameter <- "AMS"

cpg_effort_select <- cpg_effort %>% dplyr::select(c(chr_pos, parameter, parameter_qval))
cpg_physio_select <- cpg_physio %>% dplyr::select(c(chr_pos, parameter, parameter_qval))

all_models_sig <- rbind(cpg_all, cpg_changing_select, cpg_ams_select, cpg_effort_select, cpg_physio_select)

### Rename chr_pos and divide ###
all_models_sig$chr_pos <- gsub("__", ";", all_models_sig$chr_pos)
all_models_sig$chr_pos <- gsub("HRSCAF_", "HRSCAF=", all_models_sig$chr_pos, )

# Extract the numbers following HRSCAF=XXX_number
# Split the chr_pos column into two columns based on the first "_"
split_chr_pos <- strsplit(all_models_sig$chr_pos, "_", fixed = TRUE)

all_models_sig$chr <- paste0(sapply(split_chr_pos, "[", 1), "_",
                             sapply(split_chr_pos, "[", 2))

all_models_sig$pos <- sapply(split_chr_pos, "[", 3)

all_models_sig <- all_models_sig %>% 
  relocate(chr, .after = chr_pos) %>%
  relocate(pos, .after = chr_pos)

#revert scafnames
all_models_sig$chr_pos <- gsub(";","__", all_models_sig$chr_pos)
all_models_sig$chr_pos <- gsub("HRSCAF=", "HRSCAF_", all_models_sig$chr_pos)

all_models_sig$chr <- gsub(";","__", all_models_sig$chr)
all_models_sig$chr <- gsub("HRSCAF=", "HRSCAF_", all_models_sig$chr)

####### Gene regions ###########

### load annotation data
annotation_dir <- "~/PhD_grouse/grouse-annotation/output"

promoter=unique(gffToGRanges(paste0(annotation_dir, "/promoters.gff3")))
genes=unique(gffToGRanges(paste0(annotation_dir, "/genes.gff3")))
TSS=unique(gffToGRanges(paste0(annotation_dir, "/TSS.gff3")))
exons_gene=unique(gffToGRanges(paste0(annotation_dir, "/exons_gene.gff3")))
introns=unique(gffToGRanges(paste0(annotation_dir, "/introns_transcripts.gff3")))
downstream=unique(gffToGRanges(paste0(annotation_dir, "/downstream.gff3")))
upstream=unique(gffToGRanges(paste0(annotation_dir, "/upstream.gff3")))
threeUTR =unique(gffToGRanges(paste0(annotation_dir, "/threeUTRs.gff3")))
fiveUTR=unique(gffToGRanges(paste0(annotation_dir, "/fiveUTRs.gff3")))

#### Annotate ####
all_models_sig$end <- all_models_sig$pos
all_models_sig$start <- all_models_sig$pos
sig_gr <- as(all_models_sig, "GRanges")

sig_promoter <- subsetByOverlaps(sig_gr, promoter) %>% as.data.frame() %>%
  add_column("region" = "promoter", .after="parameter") %>% 
  dplyr::select(-c(seqnames:strand)) 

sig_gene <- as.data.frame(subsetByOverlaps(sig_gr, genes)) %>% as.data.frame() %>%
  add_column("region" = "gene", .after="parameter") %>% 
  dplyr::select(-c(seqnames:strand)) 

sig_tss <- as.data.frame(subsetByOverlaps(sig_gr, TSS)) %>% as.data.frame() %>%
  add_column("region" = "TSS", .after="parameter") %>%
  dplyr::select(-c(seqnames:strand)) 

sig_exon <- as.data.frame(subsetByOverlaps(sig_gr, exons_gene)) %>% as.data.frame() %>%
  add_column("region" = "exon", .after="parameter") %>% 
  dplyr::select(-c(seqnames:strand)) 

sig_intron <- as.data.frame(subsetByOverlaps(sig_gr, introns))  %>% as.data.frame() %>%
  add_column("region" = "intron", .after="parameter") %>% 
  dplyr::select(-c(seqnames:strand)) 

sig_down <- as.data.frame(subsetByOverlaps(sig_gr, downstream)) %>% as.data.frame() %>%
  add_column("region" = "downstream", .after="parameter") %>% 
  dplyr::select(-c(seqnames:strand)) 

sig_up <- as.data.frame(subsetByOverlaps(sig_gr, upstream))  %>% as.data.frame() %>%
  add_column("region" = "upstream", .after="parameter") %>% 
  dplyr::select(-c(seqnames:strand)) 

sig_threeUTR <- as.data.frame(subsetByOverlaps(sig_gr, threeUTR))  %>% as.data.frame() %>%
  add_column("region" = "threeUTR", .after="parameter") %>% 
  dplyr::select(-c(seqnames:strand)) 

sig_fiveUTR <- as.data.frame(subsetByOverlaps(sig_gr, fiveUTR))  %>% as.data.frame() %>%
  add_column("region" = "fiveUTR", .after="parameter") %>% 
  dplyr::select(-c(seqnames:strand)) 

all_models_sig_annotated <- rbind(sig_promoter, sig_gene,
                                  sig_tss, sig_exon, sig_intron, sig_down,
                                  sig_up, sig_threeUTR,  sig_fiveUTR)

summary(as.factor(all_models_sig_annotated$region))
save(all_models_sig_annotated, file="results/annotated/annotated_modeloutput_sig_annotated.RData")

sum_annotated <- as.data.frame(table(as.factor(all_models_sig_annotated$region), all_models_sig_annotated$parameter))
names(sum_annotated) <- c("region", "model", "n")
sum_annotated$model <- gsub("AMS", "Annual mating success", sum_annotated$model)
sum_annotated$model <- gsub("attend", "Attendance", sum_annotated$model)
sum_annotated$model <- gsub("dist", "Centrality", sum_annotated$model)
sum_annotated$model <- gsub("fight", "Fighting rate", sum_annotated$model)
sum_annotated$model <- gsub("time_period", "Time period", sum_annotated$model)

write.csv(sum_annotated, file="results/tables/summary_regions_sig_cpgs_all.csv", row.names=F, quote=F)

source("scripts/plotting_theme.R")
ggplot(sum_annotated, aes(x = region, y = n)) + geom_bar(stat="identity") + 
  facet_wrap(~model, ncol=2, scales="free") + coord_flip() -> num_region_cpg
ggsave(num_region_cpg, file="plots/summary/n_sig_cpg_per_region.png", width=10, height=18)

#### Based on chicken liftover ####
### load annotation data
gff <- makeTxDbFromGFF(paste0(annotation_dir, "/liftoff_gallus_ltet.gff"), 
                       format="gff3", organism="Lyrurus tetrix") 

promoters_chicken <- promoters(gff, upstream=2000, downstream=200, columns=c("tx_name", "gene_id")) # From NIOO
genes_chicken <- genes(gff)
TSS_chicken <- promoters(gff, upstream=300, downstream=50, columns=c("tx_name", "gene_id")) # TSS as in Laine et al., 2016. Nature Communications
downstream_chicken <- flank(genes(gff), 10000, start=FALSE, both=FALSE, use.names=TRUE)
upstream_chicken <- promoters(genes(gff), upstream=10000, downstream=0)
exons_gene_chicken <- unlist(exonsBy(gff, "gene")) # group exons by genes
introns_chicken <- unlist(intronsByTranscript(gff, use.names=TRUE))

exons <- exonsBy(gff, "gene")
gene <- data.frame()
for (i in 1:length(exons)){
  df <- as.data.frame(exons[[i]])
  for (j in 1:nrow(df)){
    if(grepl("XM_015281696.4-1", df$exon_name[j]) == TRUE){
    gene <- rbind(gene, names(exons[i]))  }
  }}

#rename seqnames which can't be done in gff

all_models_sig$chr_pos <- gsub("__", ";",all_models_sig$chr_pos)
all_models_sig$chr_pos <- gsub( "HRSCAF_", "HRSCAF=",all_models_sig$chr_pos)

all_models_sig$chr <- gsub("__", ";",all_models_sig$chr)
all_models_sig$chr <- gsub("HRSCAF_", "HRSCAF=", all_models_sig$chr)

sig_gr <- as(all_models_sig, "GRanges")

## annotate with chicken
sig_promoter_chicken <- mergeByOverlaps(sig_gr, promoters_chicken) 
sig_promoter_chicken_df <- unique(data.frame(chr_pos = sig_promoter_chicken$chr_pos,
                                      chr = sig_promoter_chicken$sig_gr@seqnames,
                                      pos = sig_promoter_chicken$sig_gr$pos,
                                      parameter = sig_promoter_chicken$parameter,
                                      parameter_qval = sig_promoter_chicken$parameter_qval,
                                      gene_id = sig_promoter_chicken$gene_id@unlistData,
                                      region = "promoter"))

sig_gene_chicken <- mergeByOverlaps(sig_gr, genes_chicken) 
sig_gene_chicken_df <- unique(data.frame(chr_pos = sig_gene_chicken$chr_pos,
                                             chr = sig_gene_chicken$sig_gr@seqnames,
                                             pos = sig_gene_chicken$sig_gr$pos,
                                             parameter = sig_gene_chicken$parameter,
                                      parameter_qval = sig_gene_chicken$parameter_qval,
                                           gene_id = sig_gene_chicken$gene_id,
                                             region = "gene"))

sig_TSS_chicken <- mergeByOverlaps(sig_gr, TSS_chicken) 
sig_TSS_chicken_df <- unique(data.frame(chr_pos = sig_TSS_chicken$chr_pos,
                                         chr = sig_TSS_chicken$sig_gr@seqnames,
                                         pos = sig_TSS_chicken$sig_gr$pos,
                                         parameter = sig_TSS_chicken$parameter,
                                      parameter_qval = sig_TSS_chicken$parameter_qval,
                                       gene_id = unlist(sig_TSS_chicken$gene_id),
                                         region = "TSS"))

sig_exon_chicken <- mergeByOverlaps(sig_gr, exons_gene_chicken) 
sig_exon_chicken_df <- unique(data.frame(chr_pos = sig_exon_chicken$chr_pos,
                                        chr = sig_exon_chicken$sig_gr@seqnames,
                                        pos = sig_exon_chicken$sig_gr$pos,
                                        parameter = sig_exon_chicken$parameter,
                                      parameter_qval = sig_exon_chicken$parameter_qval,
                                        gene_id = sig_exon_chicken$exon_name,
                                        region = "exon"))

sig_intron_chicken <- mergeByOverlaps(sig_gr, introns_chicken) 
sig_intron_chicken_df <- unique(data.frame(chr_pos = sig_intron_chicken$chr_pos,
                                         chr = sig_intron_chicken$sig_gr@seqnames,
                                         pos = sig_intron_chicken$sig_gr$pos,
                                         parameter = sig_intron_chicken$parameter,
                                      parameter_qval = sig_intron_chicken$parameter_qval,
                                        gene_id = NA,
                                         region = "intron"))

sig_down_chicken <- mergeByOverlaps(sig_gr, downstream_chicken) 
sig_down_chicken_df <- unique(data.frame(chr_pos = sig_down_chicken$chr_pos,
                                        chr = sig_down_chicken$sig_gr@seqnames,
                                        pos = sig_down_chicken$sig_gr$pos,
                                        parameter = sig_down_chicken$parameter,
                                      parameter_qval = sig_down_chicken$parameter_qval,
                                         gene_id = sig_down_chicken$gene_id,
                                        region = "downstream"))

sig_up_chicken <- mergeByOverlaps(sig_gr, upstream_chicken) 
sig_up_chicken_df <- unique(data.frame(chr_pos = sig_up_chicken$chr_pos,
                                         chr = sig_up_chicken$sig_gr@seqnames,
                                         pos = sig_up_chicken$sig_gr$pos,
                                        parameter = sig_up_chicken$parameter,
                                      parameter_qval = sig_up_chicken$parameter_qval,
                                       gene_id = sig_up_chicken$gene_id,
                                         region = "upstream"))


all_models_sig_annotated_chicken <- rbind(sig_promoter_chicken_df, sig_gene_chicken_df,
                                          sig_TSS_chicken_df, sig_exon_chicken_df, sig_intron_chicken_df, sig_down_chicken_df,
                                  sig_up_chicken_df)

all_models_sig_annotated_chicken <- subset(all_models_sig_annotated_chicken,
                                           !is.na(gene_id) &
                                             !grepl("LOC", gene_id))

subset(all_models_sig_annotated_chicken, parameter == "all" & region != "exon") %>% arrange(parameter_qval) %>%
  dplyr::select(gene_id) %>% unique() %>% write.csv("results/tables/all_gene_ids.csv", quote=F, row.names=F, col.names = F)

subset(all_models_sig_annotated_chicken, parameter == "time_period" & region != "exon") %>% arrange(parameter_qval) %>%
  dplyr::select(gene_id) %>% unique() %>% write.csv("results/tables/sig_gene_ids_time_period.csv", quote=F, row.names=F, col.names = F)

subset(all_models_sig_annotated_chicken, parameter == "AMS" & region != "exon") %>% arrange(parameter_qval) %>%
  dplyr::select(gene_id) %>% unique() %>% write.csv("results/tables/sig_gene_ids_AMS.csv", quote=F, row.names=F, col.names = F)

subset(all_models_sig_annotated_chicken, parameter == "dist" & region != "exon") %>% arrange(parameter_qval) %>%
  dplyr::select(gene_id) %>% unique() %>% write.csv("results/tables/sig_gene_ids_dist.csv", quote=F, row.names=F, col.names = F)

subset(all_models_sig_annotated_chicken, parameter == "attend" & region != "exon") %>% arrange(parameter_qval) %>%
  dplyr::select(gene_id) %>% unique() %>% write.csv("results/tables/sig_gene_ids_attend.csv", quote=F, row.names=F, col.names = F)

subset(all_models_sig_annotated_chicken, parameter == "hct" & region != "exon") %>% arrange(parameter_qval) %>%
  dplyr::select(gene_id) %>% unique() %>% write.csv("results/tables/sig_gene_ids_htc.csv", quote=F, row.names=F, col.names = F)

subset(all_models_sig_annotated_chicken, parameter == "ig" & region != "exon") %>% arrange(parameter_qval) %>%
  dplyr::select(gene_id) %>% unique() %>% write.csv("results/tables/sig_gene_ids_igg.csv", quote=F, row.names=F, col.names = F)

subset(all_models_sig_annotated_chicken, parameter == "trypa" & region != "exon") %>% arrange(parameter_qval) %>%
  dplyr::select(gene_id) %>% unique() %>% write.csv("results/tables/sig_gene_ids_trypa.csv", quote=F, row.names=F, col.names = F)
