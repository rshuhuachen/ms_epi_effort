coverage_raw <- fread("results/all_samples_coverage_CG.txt")
coverage_10X <- fread("results/all_samples_coverage_CG_min10X.txt")
read_quality <- fread("results/all_samples_fastqc_quality.txt")
mapping <- fread("results/all_samples_mapping_efficiency.txt")

load("data/phenotypes/fulldata_complete_epi_withdates.RData")
prepost <- subset(all_pheno_epi, !is.na(prepost))

# raw coverage
coverage_raw <- subset(coverage_raw, sample %in% prepost$epi_nr)


# 10X coverage
coverage_10X$sample <- gsub("_bismark_bt2_pe.CX_CpG_report_merge_10X.txt.gz", "", coverage_10X$sample)
coverage_10X_prepost <- subset(coverage_10X, sample %in% prepost$epi_nr)
summary(coverage_10X_prepost$mean_coverage)

# alignment/mapping
mapping_prepost <- subset(mapping, sample %in% prepost$epi_nr)
summary(mapping_prepost$efficiency)

# read quality
read_quality$sample <- gsub("_trimmed_filt_merged.1", "", read_quality$sample)
read_quality$sample <- gsub("_trimmed_filt_merged.2", "", read_quality$sample)
read_quality_prepost <- subset(read_quality, sample %in% prepost$epi_nr)
summary(read_quality_prepost$mean)
read_quality <- read_quality %>% mutate(Q30 = case_when(mean >= 30 ~ "good", TRUE ~ "bad"))
summary(as.factor(read_quality$Q30))

# raw reads
reads <- fread("results/raw_number_of_reads.txt")
names(reads) <- c("sample", "n_reads")
reads$sample <- gsub("output/lib[0-9]{1,2}/output_demultiplex/fastp_clone_stacks/", "", reads$sample)
reads$sample <- gsub("-Crick*", "", reads$sample)
reads$sample <- gsub("-Watson*", "", reads$sample)
reads$sample <- gsub("*.1.fq.gz", "", reads$sample)
reads$sample <- gsub("*.2.fq.gz", "", reads$sample)

reads_summary <- reads %>%
  group_by(sample) %>%
  summarise(total_reads = sum(n_reads, na.rm = TRUE))

reads_summary_prepost <- subset(reads_summary, sample %in% prepost$epi_nr)
