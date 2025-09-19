#### Loading methylation data ####

# In this script, we will load in the methylation data, do some exploration, and filter the data accordingly

#### Here we will prepare the prepost files for further processing ####
#### Packages ####

pacman::p_load(dplyr, data.table, methylKit, tibble)

#### Get phenotype data ####
load("data/phenotypes/fulldata_complete_epi_withdates.RData")
prepost <- subset(all_pheno_epi, !is.na(prepost)) #only select the samples relevant for this paper

#### Read in filtered data ####
#if you want to execute this, make sure to change this absolute path to 
# the relevant folder wher eyou executed the epigbs pipeline

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

## we sequenced a few samples in multiple libraries,
# as these samples contain different cell populations, we
# only use one of the two files and chose the files that have
# the most reads 
# the ones with less reads are deleted here

files_prepost <- subset(files_prepost, ids != "lib99_1")
files_prepost <- subset(files_prepost, ids != "lib20_119")
files_prepost <- subset(files_prepost, ids != "lib20_191")
files_prepost <- subset(files_prepost, ids != "lib7_250")

#as list for methylkit
ids <- as.list(files_prepost$ids)
files <- as.list(files_prepost$files)

# import the data with methylkit
ltet_meth <- methRead(files, pipeline = "bismarkCytosineReport",
                      sample.id = ids, assembly = "ltet", 
                      treatment = c(rep(1, each =length(ids))), context = "CpG",
                      sep = " ")

# filter out those with extremely high coverage
ltet_meth <- filterByCoverage(ltet_meth,lo.count=10,lo.perc=NULL,
                              hi.count=NULL,hi.perc=99.9)   

# unite the two strands
ltet_meth_unite <- methylKit::unite(ltet_meth, destrand = TRUE, 
                                    min.per.group = 1L, mc.cores = 8) #1,559,800 CpG sites

# save this raw, intermediate file. it can be retrieved from XXX
save(ltet_meth_unite, file = "data/processed/methylkit_prepost_raw.RData")

#unite and keep those with >75% shared cpg sites

ltet_meth_unite_0.75 <- methylKit::unite(ltet_meth, destrand = TRUE, 
                                         min.per.group = as.integer(0.75*length(ids)), mc.cores = 8)

save(ltet_meth_unite_0.75, file = "data/processed/methylkit_prepost_min0.75.RData") # 274,197

#### Restructure the methylation file #####
load(file = "data/processed/methylkit_prepost_raw.RData")
source("scripts/function_convert_methfile.R")

# we specify a threshold that indicates how many samples need to be variable in methylation
prepost_long <- convert_meth(methfile = ltet_meth_unite, novar = "remove", threshold = 0) #Out of 1559800 CpG sites, kept 1430526 which is 8.29% removed
prepost_long <- convert_meth(methfile = ltet_meth_unite, novar = "remove", threshold = 0.3) #Out of 1559800 CpG sites, kept 815460 which is 47.72% removed  
save(prepost_long, file = "data/processed/methylkit_prepost_long_onlyvar_thres0.3.RData")

#### Filter for at least 50% of samples in both time points (N>30)

#count number of individuals per CpG per time point
n_per_prepost <- prepost_long %>% 
  group_by(chr_pos, prepost) %>% 
  summarise(count=n())

n_per_prepost_wide <-  spread(n_per_prepost,
                              key=prepost,
                              value=count)

colnames(n_per_prepost_wide)[2] <- "n_post"
colnames(n_per_prepost_wide)[3] <- "n_pre"

#keep only if CpG site is covered in at least 50% of samples at both time points
thres = 0.5
n_per_prepost_wide <- n_per_prepost_wide %>% mutate(keep = as.factor(case_when(n_pre  > thres*(118*0.5) & n_post > thres*(118*0.5) ~ "keep")))

summary(n_per_prepost_wide$keep) # 354,649

prepost_long_clean <- left_join(prepost_long, n_per_prepost_wide, by = c("chr_pos"))

prepost_long_clean <- subset(prepost_long_clean, keep == "keep")
prepost_long_clean$keep <- NULL

# this is the clean datafile that we continue working with in subsequent analyses
save(prepost_long_clean, file = "data/processed/methylkit_prepost_long_onlyvar_thres0.3_min_0.5_group.RData")
