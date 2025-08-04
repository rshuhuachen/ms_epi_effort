### load packages
pacman::p_load(tidyverse, data.table, tibble, performance, matrixStats, 
               parallel, performance, lmerTest, tidystats, insight, effects)

### load data
load(file = "data/processed/methylkit_prepost_long_onlyvar_thres0.3_min_0.5_group.RData") #prepost_long only variation
prepost_long <- prepost_long_clean
rm(prepost_long_clean)

## phenotype data ##
load("data/phenotypes/fulldata_complete_epi_withdates.RData")
prepost <- subset(all_pheno_epi, !is.na(prepost))

rm(all_pheno_epi)

### merge with some metadata

prepost_long <- left_join(prepost_long, prepost[,c("id", "year", "Core", "born", "fulldate")], 
                          by = c("id", "year", "fulldate"))

prepost_long <- prepost_long %>% mutate(age_year = as.factor(case_when(Core == "Core" ~ year - born,
                                                        Core == "No core" ~ NA)),
                                        age = as.factor(case_when(Core == "Core" & (year - born > 1) ~ "Adult",
                                                        Core == "Core" & (year - born == 1) ~ "Yearling",
                                                        Core == "No core" ~ "Adult")))

### load phenotypic data

load("data/phenotypes/fulldata_complete_epi_withdates.RData")

### Calculate delta methylation by matching up pre-post ####

delta_meth <- left_join(subset(prepost_long, prepost == "pre"),
                            subset(prepost_long, prepost == "post")[,c("chr_pos", "lib_id", "epi_nr", "lib", "methperc", "cov", "id", "year", "fulldate")],
                            by = c("chr_pos", "id", "year"), suffix = c("_pre", "_post"))

delta_meth <- delta_meth %>% dplyr::select(-c(numC, numT, n_sample, prepost))
delta_meth <- delta_meth %>% relocate(c(id, year, born:age), .before=lib_id_pre)

delta_meth <- delta_meth %>% mutate(delta_meth = methperc_post - methperc_pre, .after =born)
delta_meth <- delta_meth %>% mutate(diff_date = fulldate_post - fulldate_pre)
delta_meth$diff_date <- as.numeric(delta_meth$diff_date)

save(delta_meth, file = "results/modeloutput/all_sites_deltameth.RData")

