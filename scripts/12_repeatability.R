### Here we calculate repeatability ####
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



pacman::p_load(rptR)

#### All sites ####
## all sites repeatability, randomly subset 1000 CpG sites and do this 100 times

set.seed(1908) 

n_iter <- 100
n_sites <- 1000

# all CpG sites
all_sites <- unique(prepost_long$chr_pos)

# empty df for results
results <- data.frame(
  iteration = 1:n_iter,
  repeatability = NA_real_,
  CI_low = NA_real_,
  CI_high = NA_real_
)

for(i in seq_len(n_iter)){
  
  # randomly sample CpG sites
  sampled_sites <- sample(all_sites, n_sites)
  
  # subset dataset
  dat <- prepost_long %>%
    filter(chr_pos %in% sampled_sites)
  
  # calculate repeatability with rptr
  rpt_fit <- rpt(
    methperc ~ (1|id) + (1|chr_pos),
    grname = "id",
    data = dat,
    datatype = "Gaussian"
  )
  
  # store results
  results$repeatability[i] <- rpt_fit$R["id"]
  results$CI_low[i] <- rpt_fit$CI_emp["id", 1]
  results$CI_high[i]<- rpt_fit$CI_emp["id", 2]
  
  cat("Finished iteration", i, "\n")
}

summary(results$repeatability)

mean(results$repeatability, na.rm = TRUE)
sd(results$repeatability, na.rm = TRUE)

quantile(results$repeatability,
         probs = c(0.025, 0.5, 0.975),
         na.rm = TRUE)

save(results, file = "results/repeatability/results_iterations_all_sites.RData")

#### Only for changing cpg sites ####
load(file="results/modeloutput/changing/changing_sites_glmer.RData")

rpt_changing <- rpt(
  methperc ~ (1|id) + (1|chr_pos),
  grname="id",
  data=changing_cpg,
  datatype="Gaussian"
)

rpt_changing

save(rpt_changing, file = "results/repeatability/results_changing_sites.RData")

