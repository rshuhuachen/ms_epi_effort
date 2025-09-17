### let's load and compare all 'almost' sig CpG sites ####

# attend
load(file="results/modeloutput/effort/out_attend_with_pre.RData")
m_attend_out <- data

# dist
load(file="results/modeloutput/effort/out_dist_with_pre.RData")
m_dist_out <- data

# ms
load(file="results/modeloutput/effort/out_MS_with_pre.RData")
m_MS_out <- data

# blue
load(file="results/modeloutput/nextyear/out_blue_ny.RData")
m_blue_out <- data

# lyre
load(file="results/modeloutput/nextyear/out_lyre_ny.RData")
m_lyre_out <- data

# eye comb
load(file="results/modeloutput/nextyear/out_eyec_ny.RData")
m_eyec_out <- data

# survival
load(file="results/modeloutput/fitness/out_surv_deltameth_filtered.RData")

almostsig_attend <- subset(m_attend_out, parameter_pval < 0.05)
almostsig_ms <- subset(m_MS_out, parameter_pval < 0.05)
almostsig_dist <- subset(m_dist_out, parameter_pval < 0.05)
almostsig_blue <- subset(m_blue_out, deltameth_pval < 0.05)
almostsig_lyre <- subset(m_lyre_out, deltameth_pval < 0.05)
almostsig_eyec <- subset(m_eyec_out, deltameth_pval < 0.05)
almostsig_surv <- subset(delta_out_surv, surv_delta_meth_pval < 0.05)
almostsig_surv$response = "surv"

almostsig_all <- data.frame(chr_pos = c(almostsig_attend$chr_pos, 
                                        almostsig_ms$chr_pos,
                                        almostsig_dist$chr_pos,
                                        almostsig_blue$chr_pos, 
                                        almostsig_lyre$chr_pos,
                                        almostsig_eyec$chr_pos,
                                        almostsig_surv$chr_pos),
                            trait = c(as.character(almostsig_attend$parameter), 
                                      as.character(almostsig_ms$parameter),
                                      as.character(almostsig_dist$parameter),
                                      as.character(almostsig_blue$response), 
                                      as.character(almostsig_lyre$response),
                                      as.character(almostsig_eyec$response),
                                      as.character(almostsig_surv$response)),
                            estimate = c(almostsig_attend$parameter_estimate, 
                                         almostsig_ms$parameter_estimate,
                                         almostsig_dist$parameter_estimate,
                                         almostsig_blue$deltameth_estimate, 
                                         almostsig_lyre$deltameth_estimate,
                                         almostsig_eyec$deltameth_estimate,
                                         almostsig_surv$surv_delta_meth_estimate),
                            pval = c(almostsig_attend$parameter_pval, 
                                         almostsig_ms$parameter_pval,
                                         almostsig_dist$parameter_pval,
                                         almostsig_blue$deltameth_pval, 
                                         almostsig_lyre$deltameth_pval,
                                         almostsig_eyec$deltameth_pval,
                                         almostsig_surv$surv_delta_meth_pval),
                            qval = c(almostsig_attend$parameter_qval, 
                                     almostsig_ms$parameter_qval,
                                     almostsig_dist$parameter_qval,
                                     almostsig_blue$deltameth_qval, 
                                     almostsig_lyre$deltameth_qval,
                                     almostsig_eyec$deltameth_qval,
                                     almostsig_surv$surv_delta_meth_qval))

dup <- almostsig_all[duplicated(almostsig_all$chr_pos),]
dup <- unique(left_join(data.frame(chr_pos=dup[,c(1)]), almostsig_all[,c(1,2)],by="chr_pos"))

#### Annotate all sites #####
#### Packages ####
pacman::p_load(dplyr, data.table, genomation, GenomicFeatures, rtracklayer, 
               GenomicRanges)


### Rename chr_pos and divide ###
almostsig_all$chr_pos <- gsub("__", ";", almostsig_all$chr_pos)
almostsig_all$chr_pos <- gsub("HRSCAF_", "HRSCAF=", almostsig_all$chr_pos, )

# Extract the numbers following HRSCAF=XXX_number
# Split the chr_pos column into two columns based on the first "_"
split_chr_pos <- strsplit(almostsig_all$chr_pos, "_", fixed = TRUE)

almostsig_all$chr <- paste0(sapply(split_chr_pos, "[", 1), "_",
                      sapply(split_chr_pos, "[", 2))

almostsig_all$pos <- sapply(split_chr_pos, "[", 3)

almostsig_all <- almostsig_all %>% 
  relocate(chr, .after = chr_pos) %>%
  relocate(pos, .after = chr_pos)

#revert scafnames
almostsig_all$chr_pos <- gsub(";","__", almostsig_all$chr_pos)
almostsig_all$chr_pos <- gsub("HRSCAF=", "HRSCAF_", almostsig_all$chr_pos)

almostsig_all$chr <- gsub(";","__", almostsig_all$chr)
almostsig_all$chr <- gsub("HRSCAF=", "HRSCAF_", almostsig_all$chr)

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
almostsig_all$end <- almostsig_all$pos
almostsig_all$start <- almostsig_all$pos

almostsig_all_pre_unique <- unique(almostsig_all[,c("chr_pos")])
sig_gr <- as(almostsig_all, "GRanges")

#promoter <- subset(promoter, ID != "ANN20946-RA") #one cpg matches 2 annotations, error
sig_promoter <- mergeByOverlaps(sig_gr, promoter) %>% as.data.frame() %>%
  add_column("region" = "promoter", .after="trait") %>%
  mutate(distance = as.numeric(promoter.end) - as.numeric(pos) ) %>% 
  dplyr::select(c(chr_pos:qval, ID, distance)) 

sig_gene <- as.data.frame(mergeByOverlaps(sig_gr, genes)) %>% 
  add_column("region" = "gene", .after="trait") %>% 
  mutate(distance = as.numeric(pos) - as.numeric(genes.start)) %>% 
  dplyr::select(c(chr_pos:qval, ID, distance)) 

sig_tss <- as.data.frame(mergeByOverlaps(sig_gr, TSS)) %>% 
  add_column("region" = "TSS", .after="trait") %>%
  mutate(distance = as.numeric(TSS.end) - as.numeric(pos) ) %>% 
  dplyr::select(c(chr_pos:qval, ID, distance)) 

sig_exon <- as.data.frame(mergeByOverlaps(sig_gr, exons_gene)) %>% 
  add_column("region" = "exon", .after="trait") %>% 
  mutate(distance = as.numeric(pos) - as.numeric(exons_gene.start)) %>% 
  dplyr::select(c(chr_pos:qval, ID, distance)) 

sig_intron <- as.data.frame(mergeByOverlaps(sig_gr, introns))  %>%
  add_column("region" = "intron", .after="trait") %>% 
  mutate(distance = as.numeric(pos) - as.numeric(introns.start)) %>% 
  dplyr::select(c(chr_pos:qval, ID, distance)) 

sig_down <- mergeByOverlaps(downstream, sig_gr)
sig_down <- as.data.frame(sig_down@listData)
sig_down <- sig_down %>% add_column("region" = "downstream", .after="trait") %>% 
  mutate(distance = as.numeric(pos) - as.numeric(downstream.start)) %>% 
  dplyr::select(c(chr_pos:qval, ID, distance)) 

sig_up <- mergeByOverlaps(upstream, sig_gr)
sig_up <- as.data.frame(sig_up@listData)
sig_up <- sig_up %>% add_column("region" = "upstream", .after="trait") %>% 
  mutate(distance =  as.numeric(upstream.end) - as.numeric(pos)) %>% 
  dplyr::select(c(chr_pos:qval, ID, distance)) 

sig_threeUTR <- as.data.frame(mergeByOverlaps(sig_gr, threeUTR))  %>%
  add_column("region" = "threeUTR", .after="trait") %>% 
  mutate(distance =  as.numeric(threeUTR.start) - as.numeric(pos)) %>% 
  dplyr::select(c(chr_pos:qval, ID, distance)) 

sig_fiveUTR <- as.data.frame(mergeByOverlaps(sig_gr, fiveUTR))  %>% 
  add_column("region" = "fiveUTR", .after="trait") %>% 
  mutate(distance =  as.numeric(pos) - as.numeric(fiveUTR.start)) %>%
  dplyr::select(c(chr_pos:qval, ID, distance)) 

almostsig_all_annotated_raw <- rbind(sig_promoter, sig_gene,
                               sig_tss, sig_exon, sig_intron, sig_down,
                               sig_up, sig_threeUTR,  sig_fiveUTR)

summary(as.factor(almostsig_all_annotated_raw$region))

save(almostsig_all_annotated_raw, file="results/annotated/annotated_modeloutput_almostsig_annotated.RData")

#### Priority workflow: TSS > promoter > exon/intron > down/up ####

prioritize_region <- function(region) {
  # Create a lookup table for region priorities
  priority_table <- c(TSS = 1, promoter = 2, exon = 3, intron = 3, downstream = 4, upstream = 4, fiveUTR = 5, threeUTR = 5, gene = 6)
  
  # Get the priority for the given region
  priority <- priority_table[region]
  
  # Handle missing regions (assign a low priority)
  priority[is.na(priority)] <- 7
  
  priority
}


almostsig_all_annotated_correct <- almostsig_all_annotated_raw %>%
  group_by(chr_pos, trait) %>% #doesn't really matter to group it by gene id, same priority applies
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

almostsig_all_annotated_id <- left_join(almostsig_all_annotated_correct, lookup, by = c("ID" = "original"))

almostsig_all_annotated_id$similar <- toupper(almostsig_all_annotated_id$similar)

unique_cpg <- unique(data.frame(chr_pos = almostsig_all_annotated_id$chr_pos, trait = almostsig_all_annotated_id$trait))
unique_cpg <- arrange(unique_cpg, trait)
unique_cpg <- data.frame(chr_pos = unique(unique_cpg$chr_pos))
unique_cpg$cpg_name <- NA
for (i in 1:nrow(unique_cpg)){
  unique_cpg$cpg_name[i] <- LETTERS[i]
}

almostsig_all_annotated_id <- left_join(almostsig_all_annotated_id, unique(unique_cpg[,c("chr_pos", "cpg_name")]), by = "chr_pos")
almostsig_all_annotated_id <- almostsig_all_annotated_id %>% relocate(cpg_name, .before = chr_pos)
almostsig_all_annotated_id <- almostsig_all_annotated_id %>% arrange(cpg_name)

save(almostsig_all_annotated_id, file="results/annotated/annotated_modeloutput_almostsig_all_annotated_priority.RData")
almostsig_all_annotated_id$similar %>% write.csv("results/go_almostsig_all.csv", quote=F, row.names=F)

# load scaffold numbers
load("data/scaffold_names_dovetail.RData")

# Split the chr_pos column into two columns based on the first "_"
split_chr_pos <- strsplit(almostsig_all_annotated_id$chr_pos, "_", fixed = TRUE)

# Extract the numbers following HRSCAF=XXX_number
almostsig_all_annotated_id$chr <- paste0(sapply(split_chr_pos, "[", 1), "_",
                             sapply(split_chr_pos, "[", 2), ";", 
                             sapply(split_chr_pos, "[", 4), "=",
                             sapply(split_chr_pos, "[", 5))

almostsig_all_annotated_id$pos <- as.numeric(sapply(split_chr_pos, "[", 6))

# join
almostsig_all_annotated_id <- left_join(almostsig_all_annotated_id, genome[,c("contig", "scaf_nr")], by = c("chr" = "contig"))

#### select only those duplicates ####
almostsig_all_annotated_id %>%
  group_by(chr_pos) %>%
  filter(n() > 1 & trait != "eyec_nextyear") %>%
  ungroup -> almostsig_dups

almostsig_all_annotated_id %>%
  filter(trait == "attend"|trait=="dist"|trait=="MS") -> almostsig_effort

almostsig_all_annotated_id %>%
  filter(trait == "blue_nextyear"|trait=="lyre_nextyear"|trait=="surv") -> almostsig_cost

### for table select all ####
almostsig_all_annotated_id <- almostsig_all_annotated_id %>% ungroup()
almostsig_all_annotated_id %>% dplyr::select(scaf_nr, pos, trait, estimate, pval, qval, region, similar) %>% filter(trait != "eyec_nextyear")-> table_sig
table_sig$estimate <- round(table_sig$estimate, 2)
table_sig$pval <- round(table_sig$pval, 3)
table_sig$qval <- round(table_sig$qval, 2)
table_sig$pos <- prettyNum(table_sig$pos, big.mark = ",", bi.interval = 3L)

table_sig$trait <- gsub("attend", "Attendance", table_sig$trait)
table_sig$trait <- gsub("blue_nextyear", "Blue chroma (next year)", table_sig$trait)
table_sig$trait <- gsub("dist", "Centrality", table_sig$trait)
table_sig$trait <- gsub("lyre_nextyear", "Lyre size (next year)", table_sig$trait)
table_sig$trait <- gsub("MS", "Mating success", table_sig$trait)
table_sig$trait <- gsub("surv", "Survival", table_sig$trait)

table_sig$trait <- factor(table_sig$trait, levels = c("Attendance", "Centrality", "Mating success", "Survival", "Blue chroma (next year)", "Lyre size (next year)"))
table_sig <- table_sig %>% arrange(similar)

write.table(sep = "\t", table_sig, file = "results/annotated/annotated_modeloutput_almostsig_annotated.tsv", quote=F, row.names = F)

### make a manhattan plot ####
### load model output
load(file="results/modeloutput/changing/modeloutput_glmer.RData")

# add col about significance
out_glmer <- out_glmer %>% mutate(sig = as.factor(case_when(abs(mean_delta_meth) >= 0.1 & prepost_qval < 0.05 & chr_pos %in% almostsig_effort$chr_pos & chr_pos %in% almostsig_cost$chr_pos~ "Associated with effort and cost",
                                                            abs(mean_delta_meth) >= 0.1 & prepost_qval < 0.05 & chr_pos %in% almostsig_effort$chr_pos ~ "Associated with effort ",
                                                            abs(mean_delta_meth) >= 0.1 & prepost_qval < 0.05 & chr_pos %in% almostsig_cost$chr_pos~ "Associated with cost",
                                                            abs(mean_delta_meth) >= 0.1 & prepost_qval < 0.05 ~ "Dynamic",
                                                            TRUE ~ "Stable")))

out_changing <- subset(out_glmer, sig != "Stable")

# manhattan plot

# plot 
# lmer
out_changing <- out_changing %>% mutate(shade = case_when(scaf_nr %% 2 == 0 ~ "even",
                                                  TRUE ~ "odd"))

#clrs <- viridisLite::viridis(6)
out_changing %>% subset(scaf_nr <= 10) %>% 
  ggplot(aes(x = pos, y = -log10(as.numeric(prepost_pval)))) + 
  geom_point(size=3, aes(alpha =shade, col = sig, fill = sig)) +
  facet_grid(~scaf_nr,scales = 'free_x', space = 'free_x', switch = 'x') +
  labs(x = "Scaffold number", y = expression(-log[10]*"(p-value)")) +
  # scale_color_manual(values=c(clrs[5], clrs[17])) +
  # scale_fill_manual(values=alpha(c(clrs[5], clrs[17]), 0.5)) +
  scale_alpha_discrete(range=c(0.4,1))+
  theme(axis.text.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        plot.margin = margin(r = 0.5, l = 0.1, b = 0.1, t = 0.1, unit = "cm"),
        axis.line.x = element_blank(),
        axis.title.x = element_text(margin=margin(t=10)),
        axis.title.y = element_text(margin=margin(r=5)),
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank()) 


#### Binomial test #####
# load all changing cpg site numbers
load(file = "results/modeloutput/changing/gene_ids_sig_changing_similar.RData")
sum_annotated <- as.data.frame(table(as.factor(annotated_changing$similar), annotated_changing$parameter)) 
sum_annotated <- subset(sum_annotated,Var2=="time_period")
names(sum_annotated) <- c("gene", "model", "n_changing")

#load all sig cpgs 
load(file="results/annotated/annotated_modeloutput_almostsig_all_annotated_priority.RData")
almostsig_all_annotated_id_unique <- almostsig_all_annotated_id %>% dplyr::select(c(chr_pos, similar)) %>% unique()
sum_sig <- as.data.frame(table(as.factor(almostsig_all_annotated_id_unique$similar))) 
names(sum_sig) <- c("gene", "n_sig")

sum_sig_changing <- left_join(sum_annotated[,c(1,3)], sum_sig, by = "gene")
sum_sig_changing$n_sig[which(is.na(sum_sig_changing$n_sig))] <- 0
sum_sig_changing$prob_null <- sum_sig_changing$n_changing / sum(sum_sig_changing$n_changing)

c("BEST1", "NFIC", "STARD3", "UBTF")

## gene BEST1  
binom.test(x = sum_sig_changing$n_sig[which(sum_sig_changing$gene == "BEST1")], 
           n = sum(sum_sig_changing$n_changing),
           p = sum_sig_changing$prob_null[which(sum_sig_changing$gene == "BEST1")]/100) 

## gene NFIC  
binom.test(x = sum_sig_changing$n_sig[which(sum_sig_changing$gene == "NFIC")], 
           n = sum(sum_sig_changing$n_changing),
           p = sum_sig_changing$prob_null[which(sum_sig_changing$gene == "NFIC")]/100) 

## gene STARD3  
binom.test(x = sum_sig_changing$n_sig[which(sum_sig_changing$gene == "STARD3")], 
           n = sum(sum_sig_changing$n_changing),
           p = sum_sig_changing$prob_null[which(sum_sig_changing$gene == "STARD3")]/100) 

## gene UBTF  
binom.test(x = sum_sig_changing$n_sig[which(sum_sig_changing$gene == "UBTF")], 
           n = sum(sum_sig_changing$n_changing),
           p = sum_sig_changing$prob_null[which(sum_sig_changing$gene == "UBTF")]/100) 

