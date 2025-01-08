### load packages
pacman::p_load(tidyverse, data.table, tibble, performance, gaston,
               parallel, performance, lmerTest, tidystats, ggpointdensity)
               
source("scripts/plotting_theme.R")

### load data

load(file="results/modeloutput/changing/changing_sites_glmer.RData")

### load phenotypic data

load(file = "data/phenotypes/fulldata_complete_epi_withdates.RData")

### methylation difference

load(file = "results/modeloutput/all_sites_deltameth.RData")

delta_meth <- subset(delta_meth, chr_pos %in% changing_cpg$chr_pos)

### Plot pre_meth effect on delta_meth ####

lmer_pre_delta <- lmerTest::lmer(delta_meth ~ methperc_pre + (1|id), data = delta_meth)

sum_pre_delta <- as.data.frame(summary(lmer_pre_delta)$coef) # negative relationship

ggplot(delta_meth, aes(methperc_pre, delta_meth)) + 
  geom_pointdensity() + 
  scale_color_viridis_c() + 
  geom_abline(intercept=sum_pre_delta$Estimate[1], slope = sum_pre_delta$Estimate[2], 
                                          color="red", linewidth=1)+
  labs(x = "Methylation % pre-lekking", y = expression(Delta*" methylation %")) -> cor_pre_delta

ggplot(delta_meth, aes(methperc_pre, abs(delta_meth))) + 
  geom_pointdensity() + 
  scale_color_viridis_c() + geom_smooth() + 
  labs(x = "Methylation % pre-lekking", y = expression("Absolute "*Delta*" methylation %")) -> cor_pre_delta_abs

cowplot::plot_grid(cor_pre_delta, cor_pre_delta_abs, labs="auto", align="hv", axis="lb", ncol=2, label_fontface = "plain", label_size = 22) %>%
  ggsave(file = "plots/explore/pre_vs_delta.png", width=14, height=10)

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
                                           
### Run the model per trait ####
source("scripts/function_models.R")

#### attendance ####
## remove repeated samples
delta_meth_attend <- delta_meth %>%
  filter(!is.na(delta_meth)& !is.na(attend))%>%
  group_by(chr_pos, id) %>%
  sample_n(1) %>%
  ungroup()

## only select cpg sites with enough data
delta_meth_n_attend <- delta_meth_attend %>% group_by(chr_pos) %>% tally()
delta_meth_n_attend <- subset(delta_meth_n_attend, n > 20)
 
delta_meth_sub_attend <- subset(delta_meth_attend, chr_pos %in% delta_meth_n_attend$chr_pos)
length(unique(delta_meth_sub_attend$chr_pos)) # 607 sites

delta_meth_attend_ls <- delta_meth_sub_attend %>% group_split(chr_pos)

save(delta_meth_sub_attend, file = "data/processed/delta_meth_ls_attend.RData")

## run the model
m_attend_pre <- parallel::mclapply(delta_meth_attend_ls, function_model_delta_pheno_norepeat, parameter="attend", pre="control", mc.cores=4)
m_attend_pre_out <- function_process_model(m_attend_pre, dir_plots = "plots/model_out/effort", dir_data = "results/modeloutput/effort",
                                            name_file = "attend_with_pre", pretty_name = "Attendance", filter_disp=FALSE) 
nrow(m_attend_pre_out$data) #n=607
nrow(m_attend_pre_out$sig) #n=1

#### fighting ####
## remove repeated samples
delta_meth_fight <- delta_meth %>%
  filter(!is.na(delta_meth)& !is.na(fight))%>%
  group_by(chr_pos, id) %>%
  sample_n(1) %>%
  ungroup()

### only select cpg sites with enough data
delta_meth_n_fight <- delta_meth_fight %>% group_by(chr_pos) %>% tally()
delta_meth_n_fight <- subset(delta_meth_n_fight, n > 20)
 
delta_meth_sub_fight <- subset(delta_meth_fight, chr_pos %in% delta_meth_n_fight$chr_pos)
length(unique(delta_meth_sub_fight$chr_pos)) # 312 sites

delta_meth_fight_ls <- delta_meth_sub_fight %>% group_split(chr_pos)
save(delta_meth_fight_ls, file = "data/processed/delta_meth_ls_fight.RData")

## run the model
m_fight_pre <- parallel::mclapply(delta_meth_fight_ls, function_model_delta_pheno_norepeat, parameter="fight", pre="control", mc.cores=4)
m_fight_pre_out <- function_process_model(m_fight_pre, dir_plots = "plots/model_out/effort", dir_data = "results/modeloutput/effort",
                                            name_file = "fight_with_pre", pretty_name = "Centrality", filter_disp=FALSE) 
nrow(m_fight_pre_out$data) #n=312
nrow(m_fight_pre_out$sig) #n=0

#### centrality ####
## remove repeated samples
delta_meth_dist <- delta_meth %>%
  filter(!is.na(delta_meth)& !is.na(dist))%>%
  group_by(chr_pos, id) %>%
  sample_n(1) %>%
  ungroup()

## only select cpg sites with enough data
delta_meth_n_dist <- delta_meth_dist %>% group_by(chr_pos) %>% tally()
delta_meth_n_dist <- subset(delta_meth_n_dist, n > 20)
 
delta_meth_sub_dist <- subset(delta_meth_dist, chr_pos %in% delta_meth_n_dist$chr_pos)
length(unique(delta_meth_sub_dist$chr_pos)) # 534 sites

delta_meth_dist_ls <- delta_meth_sub_dist %>% group_split(chr_pos)
save(delta_meth_dist_ls, file = "data/processed/delta_meth_ls_dist.RData")

## model
m_dist_pre <- parallel::mclapply(delta_meth_dist_ls, function_model_delta_pheno_norepeat, parameter="dist", pre="control", mc.cores=4)
m_dist_pre_out <- function_process_model(m_dist_pre, dir_plots = "plots/model_out/effort", dir_data = "results/modeloutput/effort",
                                            name_file = "dist_with_pre", pretty_name = "Centrality", filter_disp=FALSE)
nrow(m_dist_pre_out$data) #n=534
nrow(m_dist_pre_out$sig) #n=0

### Annotate sig CpG attend #####
#### Packages ####
pacman::p_load(tidyverse, data.table, genomation, GenomicFeatures, rtracklayer, 
               GenomicRanges)

#### Load significant site ####

## attendance
load(file="results/modeloutput/effort/out_attend_with_pre.RData")
cpg_attend_pre <- subset(data, parameter_qval < 0.05)%>% dplyr::select(c(chr_pos, parameter, parameter_estimate,
                                                                         parameter_se, parameter_df,
                                                                         parameter_tval, parameter_qval))
cpg_attend_pre$pre_control <- "pre"

### Rename chr_pos and divide ###
cpg_attend_pre$chr_pos <- gsub("__", ";", cpg_attend_pre$chr_pos)
cpg_attend_pre$chr_pos <- gsub("HRSCAF_", "HRSCAF=", cpg_attend_pre$chr_pos, )

# Extract the numbers following HRSCAF=XXX_number
# Split the chr_pos column into two columns based on the first "_"
split_chr_pos <- strsplit(cpg_attend_pre$chr_pos, "_", fixed = TRUE)

cpg_attend_pre$chr <- paste0(sapply(split_chr_pos, "[", 1), "_",
                             sapply(split_chr_pos, "[", 2))

cpg_attend_pre$pos <- sapply(split_chr_pos, "[", 3)

cpg_attend_pre <- cpg_attend_pre %>% 
  relocate(chr, .after = chr_pos) %>%
  relocate(pos, .after = chr_pos)

#revert scafnames
cpg_attend_pre$chr_pos <- gsub(";","__", cpg_attend_pre$chr_pos)
cpg_attend_pre$chr_pos <- gsub("HRSCAF=", "HRSCAF_", cpg_attend_pre$chr_pos)

cpg_attend_pre$chr <- gsub(";","__", cpg_attend_pre$chr)
cpg_attend_pre$chr <- gsub("HRSCAF=", "HRSCAF_", cpg_attend_pre$chr)

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
cpg_attend_pre$end <- cpg_attend_pre$pos
cpg_attend_pre$start <- cpg_attend_pre$pos

cpg_attend_pre_unique <- unique(cpg_attend_pre[,c("chr_pos")])
sig_gr <- as(cpg_attend_pre, "GRanges")

#promoter <- subset(promoter, ID != "ANN20946-RA") #one cpg matches 2 annotations, error
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

cpg_attend_pre_annotated_raw <- rbind(sig_promoter, sig_gene,
                                      sig_tss, sig_exon, sig_intron, sig_down,
                                      sig_up, sig_threeUTR,  sig_fiveUTR)

summary(as.factor(cpg_attend_pre_annotated_raw$region))

#save(cpg_attend_pre_annotated_raw, file="results/annotated/annotated_modeloutput_sig_annotated.RData")

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


cpg_attend_pre_annotated_correct <- cpg_attend_pre_annotated_raw %>%
  group_by(chr_pos, parameter) %>% #doesn't really matter to group it by gene id, same priority applies
  mutate(region_priority = prioritize_region(region)) %>%
  # Select the row with the highest priority (lowest numeric value)
  slice_min(region_priority, with_ties=T) %>%
  # If multiple rows have the same priority (the case for up/downstream), select the one with the lowest distance
  slice_min(distance, with_ties=T) %>%
  ungroup() %>%
  dplyr::select(-region_priority)

### merge with 'similar to' column to get gene IDs ####
lookup <- fread("/home/nioo/rebeccash/PhD_grouse/grouse-annotation/data/lookup_ANN_gene_id.txt")
names(lookup) <- c("original", "similar")

cpg_attend_pre_annotated_id <- left_join(cpg_attend_pre_annotated_correct, lookup, by = c("ID" = "original"))

cpg_attend_pre_annotated_id$similar <- toupper(cpg_attend_pre_annotated_id$similar)

unique_cpg <- unique(data.frame(chr_pos = cpg_attend_pre_annotated_id$chr_pos, parameter = cpg_attend_pre_annotated_id$parameter))
unique_cpg <- arrange(unique_cpg, parameter)
unique_cpg <- data.frame(chr_pos = unique(unique_cpg$chr_pos))
unique_cpg$cpg_name <- NA
for (i in 1:nrow(unique_cpg)){
  unique_cpg$cpg_name[i] <- LETTERS[i]
}

cpg_attend_pre_annotated_id <- left_join(cpg_attend_pre_annotated_id, unique(unique_cpg[,c("chr_pos", "cpg_name")]), by = "chr_pos")
cpg_attend_pre_annotated_id <- cpg_attend_pre_annotated_id %>% relocate(cpg_name, .before = chr_pos)
cpg_attend_pre_annotated_id <- cpg_attend_pre_annotated_id %>% arrange(cpg_name)

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

cpg_attend_pre$chr_pos <- gsub("__", ";",cpg_attend_pre$chr_pos)
cpg_attend_pre$chr_pos <- gsub( "HRSCAF_", "HRSCAF=",cpg_attend_pre$chr_pos)

cpg_attend_pre$chr <- gsub("__", ";",cpg_attend_pre$chr)
cpg_attend_pre$chr <- gsub("HRSCAF_", "HRSCAF=", cpg_attend_pre$chr)

sig_gr <- as(cpg_attend_pre, "GRanges")

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


cpg_attend_pre_annotated_chicken <- rbind(sig_promoter_chicken_df, sig_gene_chicken_df,
                                          sig_exon_chicken_df, sig_down_chicken_df,
                                  sig_up_chicken_df) # left out intron due to error



