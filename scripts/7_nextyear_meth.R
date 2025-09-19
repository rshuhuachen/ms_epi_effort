### load packages ####
pacman::p_load(dplyr, data.table, tibble, performance, gaston, cowplot,
               parallel, performance, lmerTest, tidystats, glmmTMB, DHARMa)

### load data ####

load(file="results/modeloutput/changing/changing_sites_glmer.RData")

### load phenotypic data

pheno <- read.csv("data/phenotypes/data_for_nextyear_corrected.csv")

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

### Run models ####
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

#### Collect and annotate almost sig sites ####
almostsig_ny <- rbind(almostsig_blue_ny, almostsig_lyre_ny)
dup_ny <- almostsig_ny[duplicated(almostsig_ny$chr_pos),]
dup_ny <- left_join(data.frame(chr_pos=dup_ny[,c(1)]), almostsig_ny[,c(1,2)],by="chr_pos")

save(dup_ny, file = "results/modeloutput/effort/multiple_cpgs_ny_nofdr.RData")
load(file = "results/modeloutput/effort/multiple_cpgs_effort_nofdr.RData")

dup$chr_pos %in% dup_ny$chr_pos # none
