### load packages
pacman::p_load(dplyr, data.table, tibble, performance, gaston, cowplot,
               parallel, performance, lmerTest, tidystats, glmmTMB, DHARMa)

### load data

load(file="results/modeloutput/changing/changing_sites_glmer.RData")

### load phenotypic data

pheno <- read.csv("data/phenotypes/data_for_nextyear_corrected.csv")

# load("data/phenotypes/pheno_dif_prepost.RData") ## differences in physiology

### load methylation difference

load(file = "results/modeloutput/all_sites_deltameth.RData")

delta_meth_raw <- subset(delta_meth, chr_pos %in% changing_cpg$chr_pos)

### load formulas
source("scripts/function_models.R")

# merge metadata
delta_meth_nextyear <- left_join(delta_meth_raw[,-c(6,8,9)], pheno, by = (c("id", "year")))

delta_meth_nextyear <- delta_meth_nextyear %>%
  group_by(chr_pos, id) %>%
  sample_n(1) %>%
  ungroup()

delta_meth_nextyear_ls <- delta_meth_nextyear %>% group_split(chr_pos)

#### blue chroma ####
ggplot(pheno, aes(x = blue_nextyear)) + geom_histogram() # normal distribution

## remove repeated samples
delta_meth_blue_ny <- delta_meth_nextyear %>%
  filter(!is.na(delta_meth)& !is.na(blue_nextyear))%>%
  group_by(chr_pos, id) %>%
  sample_n(1) %>%
  ungroup()

## only select cpg sites with enough data
delta_meth_n_blue_ny <- delta_meth_blue_ny %>% group_by(chr_pos) %>% tally()
delta_meth_n_blue_ny_min12 <- subset(delta_meth_n_blue_ny, n > 12)

delta_meth_sub_blue_ny <- subset(delta_meth_blue_ny, chr_pos %in% delta_meth_n_blue_ny_min12$chr_pos)
length(unique(delta_meth_sub_blue_ny$chr_pos)) # 125 sites
save(delta_meth_sub_blue_ny, file = "data/processed/delta_meth_sub_blue_ny.RData")

delta_meth_sub_blue_ny_ls <- delta_meth_sub_blue_ny %>% group_split(chr_pos)

# run model
m_blue_ny <- parallel::mclapply(delta_meth_sub_blue_ny_ls, function_lmer_nextyear, parameter="blue_nextyear", mc.cores=4)
m_blue_ny_out <- function_process_model_nextyear(m_blue_ny, dir_plots = "plots/model_out/nextyear", dir_data = "results/modeloutput/nextyear",
                                         name_file = "blue_ny", pretty_name = "Next year blue chroma", filter_disp=FALSE)

nrow(m_blue_ny_out$data) #n=125
nrow(m_blue_ny_out$sig) #n=0 but n=5 without FDR correction

almostsig_blue_ny <- subset(m_blue_ny_out$data, deltameth_pval < 0.05)
summary(almostsig_blue_ny$deltameth_estimate)

#### lyre  ####
ggplot(pheno, aes(x = lyre_nextyear)) + geom_histogram() # normal distribution?

## remove repeated samples
delta_meth_lyre_ny <- delta_meth_nextyear %>%
  filter(!is.na(delta_meth)& !is.na(lyre_nextyear))%>%
  group_by(chr_pos, id) %>%
  sample_n(1) %>%
  ungroup()

## only select cpg sites with enough data
delta_meth_n_lyre_ny <- delta_meth_lyre_ny %>% group_by(chr_pos) %>% tally()
delta_meth_n_lyre_ny_min12 <- subset(delta_meth_n_lyre_ny, n > 12)

delta_meth_sub_lyre_ny <- subset(delta_meth_lyre_ny, chr_pos %in% delta_meth_n_lyre_ny_min12$chr_pos)
length(unique(delta_meth_sub_lyre_ny$chr_pos)) # 125 sites
save(delta_meth_sub_lyre_ny, file = "data/processed/delta_meth_sub_lyre_ny.RData")

delta_meth_sub_lyre_ny_ls <- delta_meth_sub_lyre_ny %>% group_split(chr_pos)

# run model
m_lyre_ny <- parallel::mclapply(delta_meth_sub_lyre_ny_ls, function_lmer_nextyear, parameter="lyre_nextyear", mc.cores=4)
m_lyre_ny_out <- function_process_model_nextyear(m_lyre_ny, dir_plots = "plots/model_out/nextyear", dir_data = "results/modeloutput/nextyear",
                                                 name_file = "lyre_ny", pretty_name = "Next year lyre", filter_disp=FALSE)

nrow(m_lyre_ny_out$data) #n=125
nrow(m_lyre_ny_out$sig) #n=2 but n=16 without FDR correction

almostsig_lyre_ny = subset(m_lyre_ny_out$data, deltameth_pval < 0.05)
summary(almostsig_lyre_ny$deltameth_estimate)

ggplot(subset(delta_meth_sub_lyre_ny, chr_pos == "ScEsiA3_16641__HRSCAF_18713_15489035"), 
       aes(x = delta_meth, y = lyre_nextyear)) + geom_point() + geom_smooth(method='lm')

ggplot(subset(delta_meth_sub_lyre_ny, chr_pos == "ScEsiA3_16858__HRSCAF_19404_134"), 
       aes(x = delta_meth, y = lyre_nextyear)) + geom_point() + geom_smooth(method='lm')

#### red eyec - null test  ####
ggplot(pheno, aes(x = eyec_nextyear)) + geom_histogram() # normal distribution?

## remove repeated samples
delta_meth_eyec_ny <- delta_meth_nextyear %>%
  filter(!is.na(delta_meth)& !is.na(eyec_nextyear))%>%
  group_by(chr_pos, id) %>%
  sample_n(1) %>%
  ungroup()

## only select cpg sites with enough data
delta_meth_n_eyec_ny <- delta_meth_eyec_ny %>% group_by(chr_pos) %>% tally()
delta_meth_n_eyec_ny_min12 <- subset(delta_meth_n_eyec_ny, n > 12)

delta_meth_sub_eyec_ny <- subset(delta_meth_eyec_ny, chr_pos %in% delta_meth_n_eyec_ny_min12$chr_pos)
length(unique(delta_meth_sub_eyec_ny$chr_pos)) # 40 sites

delta_meth_sub_eyec_ny_ls <- delta_meth_sub_eyec_ny %>% group_split(chr_pos)

# run model
m_eyec_ny <- parallel::mclapply(delta_meth_sub_eyec_ny_ls, function_lmer_nextyear, parameter="eyec_nextyear", mc.cores=4)
m_eyec_ny_out <- function_process_model_nextyear(m_eyec_ny, dir_plots = "plots/model_out/nextyear", dir_data = "results/modeloutput/nextyear",
                                                 name_file = "eyec_ny", pretty_name = "Next year eye comb size", filter_disp=FALSE)

nrow(m_eyec_ny_out$data) #n=46
nrow(m_eyec_ny_out$sig) #n=2

almostsig_eyec_ny = subset(m_eyec_ny_out$data, deltameth_pval < 0.05)
summary(almostsig_eyec_ny$deltameth_estimate)

ggplot(subset(delta_meth_sub_eyec_ny, chr_pos == "ScEsiA3_18278__HRSCAF_21663_127143450"), 
       aes(x = delta_meth, y = eyec_nextyear)) + geom_point() + geom_smooth(method='lm')

#### Collect and annotate almost sig sites ####
almostsig_ny <- rbind(almostsig_blue_ny, almostsig_eyec_ny, almostsig_lyre_ny)
dup_ny <- almostsig_ny[duplicated(almostsig_ny$chr_pos),]
dup_ny <- left_join(data.frame(chr_pos=dup_ny[,c(1)]), almostsig_ny[,c(1,2)],by="chr_pos")

save(dup_ny, file = "results/modeloutput/effort/multiple_cpgs_ny_nofdr.RData")
load(file = "results/modeloutput/effort/multiple_cpgs_effort_nofdr.RData")

dup$chr_pos %in% dup_ny$chr_pos # none

#### Collect and annotate all sig sites ####

all_sig_ny <- rbind(subset(m_blue_ny_out$data, m_blue_ny_out$data$deltameth_pval < 0.05),
                    subset(m_lyre_ny_out$data, m_lyre_ny_out$data$deltameth_pval < 0.05),
                    subset(m_eyec_ny_out$data, m_eyec_ny_out$data$deltameth_pval < 0.05))

### Packages ####
pacman::p_load(genomation, GenomicFeatures, rtracklayer, 
               GenomicRanges)

### Combine all sites vs changing sites
all_sig_ny <- all_sig_ny %>% dplyr::select(c(response, chr_pos, sig, deltameth_pval))

### Rename chr_pos and divide ###
all_sig_ny$chr_pos <- gsub("__", ";", all_sig_ny$chr_pos)
all_sig_ny$chr_pos <- gsub("HRSCAF_", "HRSCAF=", all_sig_ny$chr_pos, )

# Extract the numbers following HRSCAF=XXX_number
# Split the chr_pos column into two columns based on the first "_"
split_chr_pos <- strsplit(all_sig_ny$chr_pos, "_", fixed = TRUE)

all_sig_ny$chr <- paste0(sapply(split_chr_pos, "[", 1), "_",
                             sapply(split_chr_pos, "[", 2))

all_sig_ny$pos <- sapply(split_chr_pos, "[", 3)

all_sig_ny <- all_sig_ny %>% 
  relocate(chr, .after = chr_pos) %>%
  relocate(pos, .after = chr_pos)

#revert scafnames
all_sig_ny$chr_pos <- gsub(";","__", all_sig_ny$chr_pos)
all_sig_ny$chr_pos <- gsub("HRSCAF=", "HRSCAF_", all_sig_ny$chr_pos)

all_sig_ny$chr <- gsub(";","__", all_sig_ny$chr)
all_sig_ny$chr <- gsub("HRSCAF=", "HRSCAF_", all_sig_ny$chr)

### Load annotation data
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
all_sig_ny$end <- all_sig_ny$pos
all_sig_ny$start <- all_sig_ny$pos
sig_gr <- as(all_sig_ny, "GRanges")

sig_promoter <- mergeByOverlaps(promoter, sig_gr)
sig_promoter <- as.data.frame(sig_promoter@listData)
sig_promoter <- sig_promoter %>% add_column("region" = "promoter", .after="response") %>% 
  mutate(distance = as.numeric(promoter.end) - as.numeric(pos) ) %>% 
  dplyr::select(c(chr_pos, pos, response, region, deltameth_pval,  ID, distance)) 

sig_gene <- mergeByOverlaps(sig_gr, genes)
sig_gene <- as.data.frame(sig_gene@listData)
sig_gene <- sig_gene %>% add_column("region" = "gene", .after="response") %>% 
  mutate(distance = as.numeric(pos) - as.numeric(genes.start)) %>% 
  dplyr::select(c(chr_pos, pos, response, region, deltameth_pval,  ID, distance)) 

sig_tss <- mergeByOverlaps(sig_gr, TSS)
sig_tss <- as.data.frame(sig_tss@listData)
sig_tss <- sig_tss %>% add_column("region" = "TSS", .after="response") %>% 
  mutate(distance = as.numeric(TSS.end) - as.numeric(pos) ) %>% 
  dplyr::select(c(chr_pos, pos, response, region, deltameth_pval,  ID, distance)) 

sig_exon <- mergeByOverlaps(sig_gr, exons_gene)
sig_exon <- as.data.frame(sig_exon@listData)
sig_exon <- sig_exon %>% add_column("region" = "exon", .after="response") %>% 
  mutate(distance = as.numeric(pos) - as.numeric(exons_gene.start)) %>% 
  dplyr::select(c(chr_pos, pos, response, region, deltameth_pval,  ID, distance)) 

sig_intron <- mergeByOverlaps(sig_gr, introns)
sig_intron <- as.data.frame(sig_intron@listData)
sig_intron <- sig_intron %>% add_column("region" = "intron", .after="response") %>% 
  mutate(distance = as.numeric(pos) - as.numeric(introns.start)) %>% 
  dplyr::select(c(chr_pos, pos, response, region, deltameth_pval,  ID, distance)) 

sig_down <- mergeByOverlaps(downstream, sig_gr)
sig_down <- as.data.frame(sig_down@listData)
sig_down <- sig_down %>% add_column("region" = "downstream", .after="response") %>% 
  mutate(distance = as.numeric(pos) - as.numeric(downstream.start)) %>% 
  dplyr::select(c(chr_pos, pos, response, region, deltameth_pval, ID, distance)) 

sig_up <- mergeByOverlaps(upstream, sig_gr)
sig_up <- as.data.frame(sig_up@listData)
sig_up <- sig_up %>% add_column("region" = "upstream", .after="response") %>% 
  mutate(distance =  as.numeric(upstream.end) - as.numeric(pos)) %>% 
  dplyr::select(c(chr_pos, pos, response, region, deltameth_pval,  ID, distance)) 


all_models_sig_annotated_raw <- rbind(sig_promoter, sig_gene,
                                      sig_tss, sig_exon, sig_intron, sig_down,
                                      sig_up)


summary(as.factor(all_models_sig_annotated_raw$region))

save(all_models_sig_annotated_raw, file="results/modeloutput/nextyear/annotated_sig_cpg_raw.RData")

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


all_models_sig_annotated <- all_models_sig_annotated_raw %>%
  group_by(chr_pos, response) %>% #doesn't really matter to group it by gene id, same priority applies
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

all_models_sig_annotated_id <- left_join(all_models_sig_annotated, lookup, by = c("ID" = "original"))

all_models_sig_annotated_id$similar <- toupper(all_models_sig_annotated_id$similar)

unique_cpg <- unique(data.frame(chr_pos = all_models_sig_annotated_id$chr_pos, response = all_models_sig_annotated_id$response))
unique_cpg <- arrange(unique_cpg, response)
unique_cpg <- data.frame(chr_pos = unique(unique_cpg$chr_pos))
unique_cpg$cpg_name <- NA
for (i in 1:nrow(unique_cpg)){
  unique_cpg$cpg_name[i] <- LETTERS[i]
}

all_models_sig_annotated_id <- left_join(all_models_sig_annotated_id, unique(unique_cpg[,c("chr_pos", "cpg_name")]), by = "chr_pos")
all_models_sig_annotated_id <- all_models_sig_annotated_id %>% relocate(cpg_name, .before = chr_pos)
all_models_sig_annotated_id <- all_models_sig_annotated_id %>% arrange(cpg_name)

save(all_models_sig_annotated_id, file="results/modeloutput/nextyear/annotated_sig_cpg_nextyear.RData")

