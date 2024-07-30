#### Loading methylation data ####

# In this script, we will load in the methylation data, do some exploration, and filter the data accordingly

#### Here we will prepare the prepost files for further processing ####
#### Packages ####

pacman::p_load(tidyverse, data.table, methylKit, tibble)

#### Get phenotype data ####
load("data/phenotypes/fulldata_complete_epi_withdates.RData")
prepost <- subset(all_pheno_epi, !is.na(prepost))

#### Read in filtered data ####
files <- NULL
for (i in c(1:20, 99)){
  files_i <- list.files(path=paste0("/data/processed/AnE/VanOersgroup/2023_Lyrurus_tetrix/epigbs_output/lib", i), pattern = "merge",
                        full.names=T)
  files <- c(files, files_i)
  rm(files_i)
}

## get id's
files <- as.data.frame(files)
files$ids <- gsub("/data/processed/AnE/VanOersgroup/2023_Lyrurus_tetrix/epigbs_output/|_bismark_bt2_pe.CX_CpG_report_merge_10X.txt.gz|_bismark_bt2_pe.CpG_CpG_report_merge_10X.txt.gz", "", files$files)
files$ids <- gsub("/", "_", files$ids)
files$epi_nr <- gsub(".*_", "", files$ids)

## subset based on epinr
files_prepost <- subset(files, epi_nr %in% prepost$epi_nr)

## those with replicate samples across libraries, change to the one with most reads
files_prepost <- subset(files_prepost, ids != "lib99_1")
files_prepost <- subset(files_prepost, ids != "lib20_119")
files_prepost <- subset(files_prepost, ids != "lib20_191")
files_prepost <- subset(files_prepost, ids != "lib7_250")

#as list for methylkit
ids <- as.list(files_prepost$ids)
files <- as.list(files_prepost$files)

ltet_meth <- methRead(files, pipeline = "bismarkCytosineReport",
                      sample.id = ids, assembly = "ltet", 
                      treatment = c(rep(1, each =length(ids))), context = "CpG",
                      sep = " ")

ltet_meth <- filterByCoverage(ltet_meth,lo.count=10,lo.perc=NULL,
                                    hi.count=NULL,hi.perc=99.9)   


ltet_meth_unite <- methylKit::unite(ltet_meth, destrand = TRUE, 
                                    min.per.group = 1L, mc.cores = 8) #1,559,800 CpG sites


save(ltet_meth_unite, file = "data/processed/methylkit_prepost_raw.RData")

#unite and keep those with >75% shared cpg sites

ltet_meth_unite_0.75 <- methylKit::unite(ltet_meth, destrand = TRUE, 
                                    min.per.group = as.integer(0.75*length(ids)), mc.cores = 8)

save(ltet_meth_unite_0.75, file = "data/processed/methylkit_prepost_min0.75.RData") # 274,197

#### Restructure the methylation file #####
load(file = "data/processed/methylkit_prepost_min0.75.RData")
source("scripts/function_convert_methfile.R")

prepost_long <- convert_meth(methfile = ltet_meth_unite_0.75, novar = "remove") #Out of 274,197 CpG sites, kept 274,188 which is 1% removed
save(prepost_long, file = "data/processed/methylkit_prepost_long_onlyvar_min0.75.RData")

