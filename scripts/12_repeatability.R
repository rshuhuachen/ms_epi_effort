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
pacman::p_load(tidyverse, data.table)
pacman::p_load(rptR, future.apply)
prepost_long <- left_join(prepost_long, prepost[,c("id", "year", "Core", "born", "fulldate")], 
                          by = c("id", "year", "fulldate"))

prepost_long <- prepost_long %>% mutate(age_year = as.factor(case_when(Core == "Core" ~ year - born,
                                                                       Core == "No core" ~ NA)),
                                        age = as.factor(case_when(Core == "Core" & (year - born > 1) ~ "Adult",
                                                                  Core == "Core" & (year - born == 1) ~ "Yearling",
                                                                  Core == "No core" ~ "Adult")))





#### All sites ####
## all sites repeatability, randomly subset 1000 CpG sites and do this 100 times

set.seed(1908)
n_iter <- 100
n_sites <- 1000
all_sites <- unique(prepost_long$chr_pos)

# pre-generate site samples sequentially so results are reproducible
site_samples <- lapply(seq_len(n_iter), function(i) sample(all_sites, n_sites))

# pre-subset the data SEQUENTIALLY - keep only needed columns too,
# to shrink what gets sent to each worker
data_subsets <- lapply(site_samples, function(sites) {
  prepost_long %>%
    filter(chr_pos %in% sites) %>%
    dplyr::select(id, chr_pos, numC, cov, methperc)   # drop unused columns if any exist
})

# function each worker will run - only needs a single small subset
run_one <- function(dat) {
  library(rptR)   # must load inside worker for PSOCK clusters
  dat$id_chr_pos <- interaction(dat$id, dat$chr_pos, drop = TRUE)

  rpt_fit <- rpt(
    cbind(numC, cov) ~ (1|id_chr_pos),
    grname = c("id_chr_pos"),
    data = dat,
    datatype = "Proportion", 
    nboot=0,
    npermut = 0
  )
  
  data.frame(
    repeatability = rpt_fit$R["id:chr_pos"]
  )
}

# set up cluster
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)

# export only run_one - NOT prepost_long (avoids the huge transfer)
clusterExport(cl, varlist = "run_one")

# set reproducible parallel RNG stream
clusterSetRNGStream(cl, iseed = 1908)

# run in parallel - data_subsets is passed as the X argument,
# so each worker only receives the chunk it's currently processing
results_list <- parLapply(cl, data_subsets, run_one)

stopCluster(cl)

results <- bind_rows(results_list, .id = "iteration") %>%
  mutate(iteration = as.integer(iteration))

summary(results$chr_pos)
mean(results$chr_pos, na.rm = TRUE)
sd(results$chr_pos, na.rm = TRUE)
quantile(results$chr_pos, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
save(results, file = "results/repeatability/results_iterations_all_sites.RData")

#### Only for changing cpg sites ####
load(file="results/modeloutput/changing/changing_sites_glmer.RData")

rpt_changing <- rpt(
  methperc ~ (1|id) + (1|chr_pos),
  grname="chr_pos",
  data=changing_cpg,
  datatype="Gaussian"
)

rpt_changing

save(rpt_changing, file = "results/repeatability/results_changing_sites.RData")

#### Calculate per CpG site repeatability ####
library(dplyr)
library(purrr)

# ---- 1. Make sure grouping vars are factors ----
prepost_long <- prepost_long %>%
  mutate(id = factor(id), chr_pos = factor(chr_pos))

# ---- 2. Load changing sites list ----
load(file = "results/modeloutput/changing/changing_sites_glmer.RData")
changing_sites <- unique(as.character(changing_cpg$chr_pos))

# ---- 3. Split all_sites into changing vs non-changing ----
all_sites <- unique(as.character(prepost_long$chr_pos))

non_changing_sites <- setdiff(all_sites, changing_sites)
changing_sites_in_data <- intersect(all_sites, changing_sites)

length(all_sites)
length(non_changing_sites)
length(changing_sites_in_data)

# ---- 4. Filter to sites with sufficient repeat structure (do this ONCE, upfront) ----
site_summary <- prepost_long %>%
  group_by(chr_pos, id) %>%
  summarise(n_reps = n(), .groups = "drop") %>%
  group_by(chr_pos) %>%
  summarise(n_id_with_reps = sum(n_reps >= 2), .groups = "drop")

valid_sites <- site_summary %>%
  filter(n_id_with_reps >= 2) %>%
  pull(chr_pos) %>%
  as.character()

valid_non_changing <- intersect(non_changing_sites, valid_sites)
valid_changing     <- intersect(changing_sites_in_data, valid_sites)

length(valid_non_changing)
length(valid_changing)  # sanity check - how many changing sites survive filtering

# ---- 5. Take ALL valid changing sites, and a random 1000 from non-changing ----
set.seed(1908)
selected_changing     <- valid_changing               # all of them, no sampling
selected_non_changing <- sample(valid_non_changing, 1000)

# ---- 6. Function: fit repeatability for ONE site (sequential, no parallel) ----
run_site <- function(dat, site_name) {
  library(rptR)
  
  dat$id <- factor(dat$id)
  
  fit <- tryCatch(
    rpt(
      cbind(numC, cov) ~ (1|id),
      grname = "id",
      data = dat,
      datatype = "Proportion",
      nboot = 0,
      npermut = 0
    ),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  
  if (is.null(fit)) {
    data.frame(chr_pos = site_name, R_org = NA_real_, R_link = NA_real_,
               n_id = length(unique(dat$id)), n_obs = nrow(dat), status = "failed")
  } else {
    data.frame(chr_pos = site_name,
               R_org  = fit$R["R_org",  "id"],
               R_link = fit$R["R_link", "id"],
               n_id = length(unique(dat$id)), n_obs = nrow(dat), status = "ok")
  }
}

# ---- 7. Sequential loop over a set of sites, with progress printing ----
run_sites_sequential <- function(sites, data) {
  n <- length(sites)
  results <- vector("list", n)
  
  for (i in seq_along(sites)) {
    site <- sites[i]
    dat <- data %>%
      filter(chr_pos == site) %>%
      dplyr::select(id, chr_pos, numC, cov, methperc)
    
    results[[i]] <- run_site(dat, site)
    
    if (i %% 50 == 0) cat("Processed", i, "/", n, "sites\n")
  }
  
  bind_rows(results)
}

# ---- 8. Run for non-changing sites ----
cat("Running non-changing sites...\n")
rpt_non_changing <- run_sites_sequential(selected_non_changing, prepost_long)

# ---- 9. Run for changing sites (all of them) ----
cat("Running changing sites (all)...\n")
rpt_changing <- run_sites_sequential(selected_changing, prepost_long)

rownames(rpt_changing) <- NULL
rpt_changing <- unique(rpt_changing)

# ---- 10. Save ----
save(rpt_non_changing, file = "results/repeatability/results_non_changing_1000_per_cpg.RData")
save(rpt_changing,     file = "results/repeatability/results_changing_all_per_cpg.RData")

### plot ####
source("scripts/plotting_theme.R")
library(scales)
rpt_non_changing$changing <- "Non-changing sites"
rpt_changing$changing <- "Changing sites"
allsites <- rbind(rpt_non_changing, rpt_changing)

ggplot(allsites, aes(x = R_link, fill = changing, col = changing)) + 
  geom_histogram(position="dodge") + labs(fill = "Category", y = "Count", x = "R") + 
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  scale_fill_manual(values=alpha(c(clr_sig, clrs[5]), 0.7)) +
  scale_color_manual(values=c(clr_sig, clrs[5])) +
  guides(col="none") -> hist_r

ggsave(hist_r, file = "plots/final/supp/histogram_repeatability.png", width=12, height=8)
