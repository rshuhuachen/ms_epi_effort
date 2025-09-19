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
