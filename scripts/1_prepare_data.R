#### Prepare phenotypic data file ####

## load packages
pacman::p_load(dplyr, data.table)

## load data
# pre + post
load("data/phenotypes/phenotypes_long_formodels_withsurv.RData")

# differences in season
load("data/phenotypes/pheno_dif_prepost.RData")

# load in some metadata
load("data/phenotypes/fulldata_complete_epi_withdates.RData")

## process data
# add 'next year' columns and remove some unnecessary ones
prepost_data <- prepost_dif %>% select(c(id:julian, MS:julian_dif)) # drop froh and loads

nextyear <- pheno_long_models_ly %>% select(c(id, year, age, lifespan, eyec, blue, lyre, attend, fight, dist, MS)) 
nextyear$lastyear <- nextyear$year - 1

# merge
prepost_data <- left_join(prepost_data, nextyear, suffix = c("", "_nextyear"), by = c("id", "year" = "lastyear"))
prepost_data <- left_join(prepost_data, unique(all_pheno_epi[,c("id", "Core")]), by = "id")

# rearrange
prepost_data <- prepost_data %>% select(c(id, Core, year, age, lifespan, site:eyec, blue, surv, eyec_nextyear:MS_nextyear))

write.csv(prepost_data, file = "data/phenotypes/data_for_nextyear.csv", quote=F, row.names = F)

### Export subset of data for supplementary table ####
## for sup table: export relevant info ###
prepost <- subset(all_pheno_epi, !is.na(prepost))

table <- prepost %>% dplyr::select(c(id, year, site, fulldate, prepost, attend, dist, MS, surv))
ny <- read.csv("data/phenotypes/data_for_nextyear_corrected.csv")
table <- left_join(table, ny[,c("id", "year", "age_cat", "blue_nextyear", "lyre_nextyear")], by = c("id", "year"))
table <- table %>% arrange(id, year) %>% dplyr::select(c(id, site, age_cat, prepost, fulldate, attend, dist, MS, surv, blue_nextyear, lyre_nextyear))
write.csv(table, file = "data/metadata_samples.csv", quote=F, row.names = F)


### for ncbi upload ####
#for metdata ncbi
data <- subset(all_pheno_epi, !is.na(prepost))
data <- data %>% dplyr::select(epi_nr, id, fulldate)
files_prepost$epi_nr <- as.numeric(files_prepost$epi_nr)

data$multiplex_c1 <- paste0(data$epi_nr, "-Crick.1.fq.gz")
data$multiplex_c2 <- paste0(data$epi_nr, "-Crick.2.fq.gz")
data$multiplex_w1 <- paste0(data$epi_nr, "-Watson.1.fq.gz")
data$multiplex_w2 <- paste0(data$epi_nr, "-Watson.2.fq.gz")
write.csv(data, "upload_ncbi/data_ncbi.csv", quote=F, row.names = F)

# to loop over simlink making for upload
upload <- subset(all_pheno_epi, !is.na(prepost))
upload <- upload %>% dplyr::select(epi_nr, id)
upload <- left_join(upload, files_prepost[,c("ids", "epi_nr")], by = "epi_nr")
upload$lib <- sub("_.*", "", upload$ids)
upload$multiplex_c1 <- paste0("/home/nioo/rebeccash/PhD_grouse/epigbs2-grouse/output/", upload$lib, "/output_demultiplex/fastp_clone_stacks/", upload$epi_nr, "-Crick.1.fq.gz")
upload$multiplex_c2 <- paste0("/home/nioo/rebeccash/PhD_grouse/epigbs2-grouse/output/", upload$lib, "/output_demultiplex/fastp_clone_stacks/", upload$epi_nr, "-Crick.2.fq.gz")
upload$multiplex_w1 <- paste0("/home/nioo/rebeccash/PhD_grouse/epigbs2-grouse/output/", upload$lib, "/output_demultiplex/fastp_clone_stacks/", upload$epi_nr, "-Watson.1.fq.gz")
upload$multiplex_w2 <- paste0("/home/nioo/rebeccash/PhD_grouse/epigbs2-grouse/output/", upload$lib, "/output_demultiplex/fastp_clone_stacks/", upload$epi_nr, "-Watson.2.fq.gz")

upload <- subset(upload, ids != "lib99_1")
upload <- subset(upload, ids != "lib20_119")
upload <- subset(upload, ids != "lib20_191")
upload <- subset(upload, ids != "lib7_250")

write.csv(upload, "upload_ncbi/loop_over_file.csv", quote=F, row.names = F)

# bash code to loop over creating simlinks for ncbi ftp upload