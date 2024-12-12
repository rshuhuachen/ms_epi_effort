### Annotate all significant CpG sites

### Packages ####
pacman::p_load(tidyverse, data.table, genomation, GenomicFeatures, rtracklayer, 
               GenomicRanges)

### Plotting ###
source("scripts/plotting_theme.R")

### Significants sites ####

## attendance
load(file="results/modeloutput/effort/out_attend_with_pre.RData")
cpg_attend_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, parameter, parameter_estimate,
                                                                         parameter_se, parameter_df,
                                                                         parameter_tval, parameter_qval))
cpg_attend_pre$pre_control <- "pre"

# load(file="results/modeloutput/effort/out_attend_no_pre.RData")
# cpg_attend_no_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, parameter, parameter_qval))
# cpg_attend_no_pre$pre_control <- "no_pre"

## fight 
load(file="results/modeloutput/effort/out_fight_with_pre.RData")
cpg_fight_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, parameter, parameter_estimate,
                                                                        parameter_se, parameter_df,
                                                                        parameter_tval, parameter_qval))
cpg_fight_pre$pre_control <- "pre"

# load(file="results/modeloutput/effort/out_fight_no_pre.RData")
# cpg_fight_no_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, parameter, parameter_qval))
# cpg_fight_no_pre$pre_control <- "no_pre"

## dist 
load(file="results/modeloutput/effort/out_dist_with_pre.RData")
cpg_dist_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, parameter, parameter_estimate,
                                                                       parameter_se, parameter_df,
                                                                       parameter_tval, parameter_qval))
cpg_dist_pre$pre_control <- "pre"

# load(file="results/modeloutput/effort/out_dist_no_pre.RData")
# cpg_dist_no_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, parameter, parameter_qval))
# cpg_dist_no_pre$pre_control <- "no_pre"

## mass 
load(file="results/modeloutput/physio/out_mass_dif_with_pre.RData")
cpg_mass_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, parameter, parameter_estimate,
                                                                       parameter_se, parameter_df,
                                                                       parameter_tval, parameter_qval))
cpg_mass_pre$pre_control <- "pre"

# load(file="results/modeloutput/physio/out_mass_dif_no_pre.RData")
# cpg_mass_no_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, parameter, parameter_qval))
# cpg_mass_no_pre$pre_control <- "no_pre"

## microf 
load(file="results/modeloutput/physio/out_microf_dif_with_pre.RData")
cpg_microf_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, parameter, parameter_estimate,
                                                                         parameter_se, parameter_df,
                                                                         parameter_tval, parameter_qval))
cpg_microf_pre$pre_control <- "pre"

# load(file="results/modeloutput/physio/out_microf_dif_no_pre.RData")
# cpg_microf_no_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, parameter, parameter_qval))
# cpg_microf_no_pre$pre_control <- "no_pre"

## trypa -> ns
# load(file="results/modeloutput/physio/out_trypa_dif_with_pre.RData")
# cpg_trypa_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, parameter, parameter_qval))
# cpg_trypa_pre$pre_control <- "pre"

# load(file="results/modeloutput/physio/out_trypa_dif_no_pre.RData")
# cpg_trypa_no_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, parameter, parameter_qval))
# cpg_trypa_no_pre$pre_control <- "no_pre"

## hct 
load(file="results/modeloutput/physio/out_hct_dif_with_pre.RData")
cpg_hct_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, parameter, parameter_estimate,
                                                                      parameter_se, parameter_df,
                                                                      parameter_tval, parameter_qval))
cpg_hct_pre$pre_control <- "pre"

# load(file="results/modeloutput/physio/out_hct_dif_no_pre.RData")
# cpg_hct_no_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, parameter, parameter_qval))
# cpg_hct_no_pre$pre_control <- "no_pre"

## ig 
load(file="results/modeloutput/physio/out_ig_dif_with_pre.RData")
cpg_ig_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, parameter, parameter_estimate,
                                                                     parameter_se, parameter_df,
                                                                     parameter_tval, parameter_qval))
cpg_ig_pre$pre_control <- "pre"

# load(file="results/modeloutput/physio/out_ig_dif_no_pre.RData")
# cpg_ig_no_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, parameter, parameter_qval))
# cpg_ig_no_pre$pre_control <- "no_pre"

## injury: pc1
load(file="results/modeloutput/injury/out_pc1_with_pre.RData")
cpg_pca1_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, parameter, parameter_estimate,
                                                                       parameter_se, parameter_df,
                                                                       parameter_tval, parameter_qval))
cpg_pca1_pre$pre_control <- "pre"

## injury: pc1
load(file="results/modeloutput/injury/out_pc2_with_pre.RData")
cpg_pca2_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, parameter, parameter_estimate,
                                                                       parameter_se, parameter_df,
                                                                       parameter_tval, parameter_qval))
cpg_pca2_pre$pre_control <- "pre"

### combine

all_models_sig <- rbind(cpg_attend_pre,
                        cpg_fight_pre,
                        cpg_dist_pre,
                        cpg_mass_pre,
                        cpg_microf_pre,
                        cpg_ig_pre,
                        cpg_hct_pre, 
                        cpg_pca1_pre,
                        cpg_pca2_pre)

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

save(all_models_sig, file = "results/annotated_model_sig.RData")

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

all_models_sig_unique <- unique(all_models_sig[,c("chr_pos")])
sig_gr <- as(all_models_sig, "GRanges")

promoter <- subset(promoter, ID != "ANN20946-RA") #one cpg matches 2 annotations, error
sig_promoter <- mergeByOverlaps(sig_gr, promoter) %>% as.data.frame() %>%
  add_column("region" = "promoter", .after="parameter") %>%
  mutate(distance = as.numeric(promoter.end) - as.numeric(pos) ) %>% 
  dplyr::select(c(chr_pos, pos, parameter, region, parameter_estimate:parameter_qval, pre_control, ID, distance)) 

sig_gene <- as.data.frame(mergeByOverlaps(sig_gr, genes)) %>% 
  add_column("region" = "gene", .after="parameter") %>% 
  mutate(distance = as.numeric(pos) - as.numeric(genes.start)) %>% 
  dplyr::select(c(chr_pos, pos, parameter, region, parameter_estimate:parameter_qval, pre_control, ID, distance)) 

sig_tss <- as.data.frame(mergeByOverlaps(sig_gr, TSS)) %>% 
  add_column("region" = "TSS", .after="parameter") %>%
  mutate(distance = as.numeric(TSS.end) - as.numeric(pos) ) %>% 
  dplyr::select(c(chr_pos, pos, parameter, region, parameter_estimate:parameter_qval, pre_control, ID, distance)) 

sig_exon <- as.data.frame(mergeByOverlaps(sig_gr, exons_gene)) %>% 
  add_column("region" = "exon", .after="parameter") %>% 
  mutate(distance = as.numeric(pos) - as.numeric(exons_gene.start)) %>% 
  dplyr::select(c(chr_pos, pos, parameter, region, parameter_estimate:parameter_qval, pre_control, ID, distance)) 

sig_intron <- as.data.frame(mergeByOverlaps(sig_gr, introns))  %>%
  add_column("region" = "intron", .after="parameter") %>% 
  mutate(distance = as.numeric(pos) - as.numeric(introns.start)) %>% 
  dplyr::select(c(chr_pos, pos, parameter, region, parameter_estimate:parameter_qval, pre_control, ID, distance)) 

sig_down <- mergeByOverlaps(downstream, sig_gr)
sig_down <- as.data.frame(sig_down@listData)
sig_down <- sig_down %>% add_column("region" = "downstream", .after="parameter") %>% 
  mutate(distance = as.numeric(pos) - as.numeric(downstream.start)) %>% 
  dplyr::select(c(chr_pos, pos, parameter, region, parameter_estimate:parameter_qval, pre_control, ID, distance)) 

sig_up <- mergeByOverlaps(upstream, sig_gr)
sig_up <- as.data.frame(sig_up@listData)
sig_up <- sig_up %>% add_column("region" = "upstream", .after="parameter") %>% 
  mutate(distance =  as.numeric(upstream.end) - as.numeric(pos)) %>% 
  dplyr::select(c(chr_pos, pos, parameter, region, parameter_estimate:parameter_qval, pre_control, ID, distance)) 

sig_threeUTR <- as.data.frame(mergeByOverlaps(sig_gr, threeUTR))  %>%
  add_column("region" = "threeUTR", .after="parameter") %>% 
  mutate(distance =  as.numeric(threeUTR.start) - as.numeric(pos)) %>% 
  dplyr::select(c(chr_pos, pos, parameter, region, parameter_estimate:parameter_qval, pre_control, ID, distance)) 

sig_fiveUTR <- as.data.frame(mergeByOverlaps(sig_gr, fiveUTR))  %>% 
  add_column("region" = "fiveUTR", .after="parameter") %>% 
  mutate(distance =  as.numeric(pos) - as.numeric(fiveUTR.start)) %>%
  dplyr::select(c(chr_pos, pos, parameter, region, parameter_estimate:parameter_qval, pre_control, ID, distance)) 

all_models_sig_annotated_raw <- rbind(sig_promoter, sig_gene,
                                  sig_tss, sig_exon, sig_intron, sig_down,
                                  sig_up, sig_threeUTR,  sig_fiveUTR)

summary(as.factor(all_models_sig_annotated_raw$region))

save(all_models_sig_annotated_raw, file="results/annotated/annotated_modeloutput_sig_annotated.RData")

#### Priority workflow: TSS > promoter > exon/intron > down/up ####

prioritize_region <- function(region) {
  # Create a lookup table for region priorities
  priority_table <- c(TSS = 1, promoter = 2, exon = 3, intron = 3, downstream = 4, upstream = 4, gene = 5)
  
  # Get the priority for the given region
  priority <- priority_table[region]
  
  # Handle missing regions (assign a low priority)
  priority[is.na(priority)] <- 7
  
  priority
}


all_models_sig_annotated_correct <- all_models_sig_annotated_raw %>%
  group_by(chr_pos, parameter) %>% #doesn't really matter to group it by gene id, same priority applies
  mutate(region_priority = prioritize_region(region)) %>%
  # Select the row with the highest priority (lowest numeric value)
  slice_min(region_priority, with_ties=T) %>%
  # If multiple rows have the same priority (the case for up/downstream), select the one with the lowest distance
  slice_min(distance, with_ties=T) %>%
  ungroup() %>%
  dplyr::select(-region_priority)

# manually add unannotated cpg for microf
all_models_sig_annotated_correct <- rbind(all_models_sig_annotated_correct,
                                          data.frame(chr_pos = cpg_microf_pre$chr_pos[1],
                                                     pos = 28300178,
                                                     parameter = "microf_dif",
                                                     region = NA,
                                                     parameter_estimate = cpg_microf_pre$parameter_estimate[1],
                                                     parameter_se = cpg_microf_pre$parameter_se[1],
                                                     parameter_df = cpg_microf_pre$parameter_df[1],
                                                     parameter_tval = cpg_microf_pre$parameter_tval[1],
                                                     parameter_qval = cpg_microf_pre$parameter_qval[1],
                                                     pre_control = "pre",
                                                     ID = NA,
                                                     distance = NA))

save(all_models_sig_annotated_correct, file="results/annotated/annotated_modeloutput_sig_annotated_priority.RData")

### summarise per region ####
sum_annotated <- as.data.frame(table(as.factor(all_models_sig_annotated_correct$region), all_models_sig_annotated_correct$parameter, all_models_sig_annotated_correct$pre_control))
names(sum_annotated) <- c("region", "model", "pre_control", "n")

sum_annotated$model <- gsub("attend", "Attendance", sum_annotated$model)
sum_annotated$model <- gsub("dist", "Centrality", sum_annotated$model)
sum_annotated$model <- gsub("fight", "Fighting rate", sum_annotated$model)
sum_annotated$model <- gsub("mass_dif", "Delta body mass", sum_annotated$model)
sum_annotated$model <- gsub("microf_dif", "Delta Microfilaria spp. ", sum_annotated$model)
sum_annotated$model <- gsub("ig_dif", "Delta IgG", sum_annotated$model)
sum_annotated$model <- gsub("hct_dif", "Delta HCT", sum_annotated$model)
sum_annotated$model <- gsub("PC1", "Injury - pecking", sum_annotated$model)
sum_annotated$model <- gsub("PC2", "Injury - kicking", sum_annotated$model)

sum_annotated$pre_control <- gsub("no_pre", "No pre-lekking methylation %", sum_annotated$pre_control )
sum_annotated$pre_control <- gsub("pre", "With pre-lekking methylation %", sum_annotated$pre_control )

sum_annotated$region <- gsub("downstream", "Downstream", sum_annotated$region)
sum_annotated$region <- gsub("upstream", "Upstream", sum_annotated$region)
sum_annotated$region <- gsub("exon", "Exon", sum_annotated$region)
sum_annotated$region <- gsub("fiveUTR", "5' UTR", sum_annotated$region)
sum_annotated$region <- gsub("gene", "Gene body", sum_annotated$region)
sum_annotated$region <- gsub("intron", "Intron", sum_annotated$region)
sum_annotated$region <- gsub("promoter", "Promoter", sum_annotated$region)
sum_annotated$region <- gsub("threeUTR", "3' UTR", sum_annotated$region)

sum_annotated$region <- factor(sum_annotated$region, levels = c("3' UTR", "5' UTR", "Downstream", "Upstream", "Gene body", "Exon", "Intron", "Promoter", "TSS"))

#sum_annotated <- sum_annotated %>% mutate(perc = n / n_total * 100)
sum_annotated <- subset(sum_annotated, n > 0)
write.csv(sum_annotated, file="results/tables/summary_regions_sig_cpgs_all.csv", row.names=F, quote=F)

### plot distribution
source("scripts/plotting_theme.R")
ggplot(sum_annotated, aes(x = region, y = n)) + geom_bar(stat="identity") + 
  facet_wrap(~model, ncol=2, scales="free") + coord_flip() -> num_region_cpg
ggsave(num_region_cpg, file="plots/summary/n_sig_cpg_per_region.png", width=10, height=18)

### merge with 'similar to' column to get gene IDs ####
lookup <- fread("/home/nioo/rebeccash/PhD_grouse/grouse-annotation/data/lookup_ANN_gene_id.txt")
names(lookup) <- c("original", "similar")

all_models_sig_annotated_id <- left_join(all_models_sig_annotated_correct, lookup, by = c("ID" = "original"))

all_models_sig_annotated_id$similar <- toupper(all_models_sig_annotated_id$similar)

unique_cpg <- unique(data.frame(chr_pos = all_models_sig_annotated_id$chr_pos, parameter = all_models_sig_annotated_id$parameter))
unique_cpg <- arrange(unique_cpg, parameter)
unique_cpg <- data.frame(chr_pos = unique(unique_cpg$chr_pos))
unique_cpg$cpg_name <- NA
for (i in 1:nrow(unique_cpg)){
  unique_cpg$cpg_name[i] <- LETTERS[i]
}
unique_cpg$cpg_name[27] <- "AA"
unique_cpg$cpg_name[28] <- "AB"

all_models_sig_annotated_id <- left_join(all_models_sig_annotated_id, unique(unique_cpg[,c("chr_pos", "cpg_name")]), by = "chr_pos")
all_models_sig_annotated_id <- all_models_sig_annotated_id %>% relocate(cpg_name, .before = chr_pos)
all_models_sig_annotated_id <- all_models_sig_annotated_id %>% arrange(cpg_name)

## which ones are also sig for AMS? ###
load(file="results/modeloutput/fitness/out_ams_deltameth_filtered.RData")
cpg_sig_ams <- subset(delta_out_ams, ams_delta_meth_qval < 0.05) #362

all_models_sig_annotated_id <- all_models_sig_annotated_id %>% mutate(sig_ams = case_when(
  chr_pos %in% cpg_sig_ams$chr_pos ~ "sig_AMS",
  TRUE ~ "not_sig_AMS"
))

write.csv(all_models_sig_annotated_id, file="results/tables/summary_trait_sig_cpgs.csv", row.names=F, quote=F)



# #### Based on chicken liftover ####
# ### load annotation data
# gff <- makeTxDbFromGFF(paste0(annotation_dir, "/liftoff_gallus_ltet.gff"), 
#                        format="gff3", organism="Lyrurus tetrix") 
# 
# promoters_chicken <- promoters(gff, upstream=2000, downstream=200, columns=c("tx_name", "gene_id")) # From NIOO
# genes_chicken <- genes(gff)
# TSS_chicken <- promoters(gff, upstream=300, downstream=50, columns=c("tx_name", "gene_id")) # TSS as in Laine et al., 2016. Nature Communications
# downstream_chicken <- flank(genes(gff), 10000, start=FALSE, both=FALSE, use.names=TRUE)
# upstream_chicken <- promoters(genes(gff), upstream=10000, downstream=0)
# exons_gene_chicken <- unlist(exonsBy(gff, "gene")) # group exons by genes
# introns_chicken <- unlist(intronsByTranscript(gff, use.names=TRUE))
# 
# exons <- exonsBy(gff, "gene")
# gene <- data.frame()
# for (i in 1:length(exons)){
#   df <- as.data.frame(exons[[i]])
#   for (j in 1:nrow(df)){
#     if(grepl("XM_015281696.4-1", df$exon_name[j]) == TRUE){
#     gene <- rbind(gene, names(exons[i]))  }
#   }}
# 
# #rename seqnames which can't be done in gff
# 
# all_models_sig$chr_pos <- gsub("__", ";",all_models_sig$chr_pos)
# all_models_sig$chr_pos <- gsub( "HRSCAF_", "HRSCAF=",all_models_sig$chr_pos)
# 
# all_models_sig$chr <- gsub("__", ";",all_models_sig$chr)
# all_models_sig$chr <- gsub("HRSCAF_", "HRSCAF=", all_models_sig$chr)
# 
# sig_gr <- as(all_models_sig, "GRanges")
# 
# ## annotate with chicken
# sig_promoter_chicken <- mergeByOverlaps(sig_gr, promoters_chicken) 
# sig_promoter_chicken_df <- unique(data.frame(chr_pos = sig_promoter_chicken$chr_pos,
#                                       chr = sig_promoter_chicken$sig_gr@seqnames,
#                                       pos = sig_promoter_chicken$sig_gr$pos,
#                                       parameter = sig_promoter_chicken$parameter,
#                                       parameter_qval = sig_promoter_chicken$parameter_qval,
#                                       gene_id = sig_promoter_chicken$gene_id@unlistData,
#                                       region = "promoter"))
# 
# sig_gene_chicken <- mergeByOverlaps(sig_gr, genes_chicken) 
# sig_gene_chicken_df <- unique(data.frame(chr_pos = sig_gene_chicken$chr_pos,
#                                              chr = sig_gene_chicken$sig_gr@seqnames,
#                                              pos = sig_gene_chicken$sig_gr$pos,
#                                              parameter = sig_gene_chicken$parameter,
#                                       parameter_qval = sig_gene_chicken$parameter_qval,
#                                            gene_id = sig_gene_chicken$gene_id,
#                                              region = "gene"))
# 
# sig_TSS_chicken <- mergeByOverlaps(sig_gr, TSS_chicken) 
# sig_TSS_chicken_df <- unique(data.frame(chr_pos = sig_TSS_chicken$chr_pos,
#                                          chr = sig_TSS_chicken$sig_gr@seqnames,
#                                          pos = sig_TSS_chicken$sig_gr$pos,
#                                          parameter = sig_TSS_chicken$parameter,
#                                       parameter_qval = sig_TSS_chicken$parameter_qval,
#                                        gene_id = unlist(sig_TSS_chicken$gene_id),
#                                          region = "TSS"))
# 
# sig_exon_chicken <- mergeByOverlaps(sig_gr, exons_gene_chicken) 
# sig_exon_chicken_df <- unique(data.frame(chr_pos = sig_exon_chicken$chr_pos,
#                                         chr = sig_exon_chicken$sig_gr@seqnames,
#                                         pos = sig_exon_chicken$sig_gr$pos,
#                                         parameter = sig_exon_chicken$parameter,
#                                       parameter_qval = sig_exon_chicken$parameter_qval,
#                                         gene_id = sig_exon_chicken$exon_name,
#                                         region = "exon"))
# 
# sig_intron_chicken <- mergeByOverlaps(sig_gr, introns_chicken) 
# sig_intron_chicken_df <- unique(data.frame(chr_pos = sig_intron_chicken$chr_pos,
#                                          chr = sig_intron_chicken$sig_gr@seqnames,
#                                          pos = sig_intron_chicken$sig_gr$pos,
#                                          parameter = sig_intron_chicken$parameter,
#                                       parameter_qval = sig_intron_chicken$parameter_qval,
#                                         gene_id = NA,
#                                          region = "intron"))
# 
# sig_down_chicken <- mergeByOverlaps(sig_gr, downstream_chicken) 
# sig_down_chicken_df <- unique(data.frame(chr_pos = sig_down_chicken$chr_pos,
#                                         chr = sig_down_chicken$sig_gr@seqnames,
#                                         pos = sig_down_chicken$sig_gr$pos,
#                                         parameter = sig_down_chicken$parameter,
#                                       parameter_qval = sig_down_chicken$parameter_qval,
#                                          gene_id = sig_down_chicken$gene_id,
#                                         region = "downstream"))
# 
# sig_up_chicken <- mergeByOverlaps(sig_gr, upstream_chicken) 
# sig_up_chicken_df <- unique(data.frame(chr_pos = sig_up_chicken$chr_pos,
#                                          chr = sig_up_chicken$sig_gr@seqnames,
#                                          pos = sig_up_chicken$sig_gr$pos,
#                                         parameter = sig_up_chicken$parameter,
#                                       parameter_qval = sig_up_chicken$parameter_qval,
#                                        gene_id = sig_up_chicken$gene_id,
#                                          region = "upstream"))
# 
# 
# all_models_sig_annotated_chicken <- rbind(sig_promoter_chicken_df, sig_gene_chicken_df,
#                                           sig_TSS_chicken_df, sig_exon_chicken_df, sig_down_chicken_df,
#                                   sig_up_chicken_df) # left out intron due to error
# 
# all_models_sig_annotated_chicken <- subset(all_models_sig_annotated_chicken,
#                                            !is.na(gene_id) &
#                                              !grepl("LOC", gene_id))
# 
# subset(all_models_sig_annotated_chicken, parameter == "all" & region != "exon") %>% arrange(parameter_qval) %>%
#   dplyr::select(gene_id) %>% unique() %>% write.csv("results/tables/all_gene_ids.csv", quote=F, row.names=F, col.names = F)
# 
# subset(all_models_sig_annotated_chicken, parameter == "time_period" & region != "exon") %>% arrange(parameter_qval) %>%
#   dplyr::select(gene_id) %>% unique() %>% write.csv("results/tables/sig_gene_ids_time_period.csv", quote=F, row.names=F, col.names = F)
# 
# subset(all_models_sig_annotated_chicken, parameter == "AMS" & region != "exon") %>% arrange(parameter_qval) %>%
#   dplyr::select(gene_id) %>% unique() %>% write.csv("results/tables/sig_gene_ids_AMS.csv", quote=F, row.names=F, col.names = F)
# 
# subset(all_models_sig_annotated_chicken, parameter == "dist" & region != "exon") %>% arrange(parameter_qval) %>%
#   dplyr::select(gene_id) %>% unique() %>% write.csv("results/tables/sig_gene_ids_dist.csv", quote=F, row.names=F, col.names = F)
# 
# subset(all_models_sig_annotated_chicken, parameter == "attend" & region != "exon") %>% arrange(parameter_qval) %>%
#   dplyr::select(gene_id) %>% unique() %>% write.csv("results/tables/sig_gene_ids_attend.csv", quote=F, row.names=F, col.names = F)
# 
# subset(all_models_sig_annotated_chicken, parameter == "fight" & region != "exon") %>% arrange(parameter_qval) %>%
#   dplyr::select(gene_id) %>% unique() %>% write.csv("results/tables/sig_gene_ids_fight.csv", quote=F, row.names=F, col.names = F)
# 
# subset(all_models_sig_annotated_chicken, parameter == "Delta HCT" & region != "exon") %>% arrange(parameter_qval) %>%
#   dplyr::select(gene_id) %>% unique() %>% write.csv("results/tables/sig_gene_ids_htc.csv", quote=F, row.names=F, col.names = F)
# 
# subset(all_models_sig_annotated_chicken, parameter == "Delta IgG" & region != "exon") %>% arrange(parameter_qval) %>%
#   dplyr::select(gene_id) %>% unique() %>% write.csv("results/tables/sig_gene_ids_igg.csv", quote=F, row.names=F, col.names = F)
# 
# subset(all_models_sig_annotated_chicken, parameter == "Delta Trypanosoma spp." & region != "exon") %>% arrange(parameter_qval) %>%
#   dplyr::select(gene_id) %>% unique() %>% write.csv("results/tables/sig_gene_ids_trypa.csv", quote=F, row.names=F, col.names = F)
# 
# subset(all_models_sig_annotated_chicken, parameter == "Delta Microfilaria spp." & region != "exon") %>% arrange(parameter_qval) %>%
#   dplyr::select(gene_id) %>% unique() %>% write.csv("results/tables/sig_gene_ids_microf.csv", quote=F, row.names=F, col.names = F)
# 
# subset(all_models_sig_annotated_chicken, parameter == "Delta body mass" & region != "exon") %>% arrange(parameter_qval) %>%
#   dplyr::select(gene_id) %>% unique() %>% write.csv("results/tables/sig_gene_ids_mass.csv", quote=F, row.names=F, col.names = F)
# 
# 
# ### isolate confusion
# cpg_gr <- as(subset(all_models_sig, parameter == "attend" & pos == 2574288), "GRanges")
# down_sub <- subset(downstream, gene_id == "ANN31167")
# as.data.frame(mergeByOverlaps(cpg_gr, downstream))
