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
prepost_long <- left_join(prepost_long, prepost[,c("id", "year", "Core", "born", "fulldate")], 
                          by = c("id", "year", "fulldate"))

prepost_long <- prepost_long %>% mutate(age_year = as.factor(case_when(Core == "Core" ~ year - born,
                                                                       Core == "No core" ~ NA)),
                                        age = as.factor(case_when(Core == "Core" & (year - born > 1) ~ "Adult",
                                                                  Core == "Core" & (year - born == 1) ~ "Yearling",
                                                                  Core == "No core" ~ "Adult")))



pacman::p_load(rptR, future.apply)

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
library(parallel)

# ---- 1. Make sure grouping vars are factors ----
prepost_long <- prepost_long %>%
  mutate(id = factor(id), chr_pos = factor(chr_pos))

# ---- 2. Pre-filter sites with insufficient repeat structure ----
site_summary <- prepost_long %>%
  group_by(chr_pos, id) %>%
  summarise(n_reps = n(), .groups = "drop") %>%
  group_by(chr_pos) %>%
  summarise(n_id_with_reps = sum(n_reps >= 2), .groups = "drop")

valid_sites <- site_summary %>%
  filter(n_id_with_reps >= 2) %>%
  pull(chr_pos)

length(valid_sites)
length(unique(prepost_long$chr_pos)) - length(valid_sites)  # how many dropped

## dynamic cpg sites only
dynamic <- subset(prepost_long, chr_pos %in% valid_sites & chr_pos %in% changing_cpg$chr_pos)
dynamic_site <- dynamic %>% dplyr::select(id, chr_pos, numC, cov, methperc) %>%
  arrange(chr_pos) %>%
  group_split(chr_pos, .keep = TRUE)

`%!in%` = Negate(`%in%`)
nondynamic <- subset(prepost_long, chr_pos %in% valid_sites & chr_pos %!in% changing_cpg$chr_pos)
nondynamic_random <- sample(nondynamic$chr_pos, size = 1000)
nondynamic_site <- subset(prepost_long, chr_pos %in% nondynamic_random)%>%
  dplyr::select(id, chr_pos, numC, cov, methperc) %>%
  arrange(chr_pos) %>%
  group_split(chr_pos, .keep = TRUE)

# ---- 3. Split into per-site data, then chunk for parallel dispatch ----
n_sites_dynamic <- length(dynamic_site)
chunk_size <- 150
chunk_idx_dynamic <- split(seq_len(n_sites_dynamic), ceiling(seq_len(n_sites_dynamic) / chunk_size))
data_chunks_dynamic <- lapply(chunk_idx_dynamic, function(idx) dynamic_site[idx])

n_sites_nondynamic <- length(dynamic_site)
chunk_idx_nondynamic <- split(seq_len(n_sites_nondynamic), ceiling(seq_len(n_sites_nondynamic) / chunk_size))
data_chunks_nondynamic <- lapply(chunk_idx_nondynamic, function(idx) nondynamic_site[idx])

# ---- 4. Worker function: loop over sites within one chunk ----
run_chunk <- function(chunk) {
  library(rptR)
  library(purrr)
  
  map_dfr(chunk, function(dat) {
    site_name <- as.character(dat$chr_pos[1])
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
      data.frame(chr_pos = site_name, R = NA_real_,
                 n_id = length(unique(dat$id)), n_obs = nrow(dat), status = "failed")
    } else {
      data.frame(chr_pos = site_name, R = unname(fit$R[1]),
                 n_id = length(unique(dat$id)), n_obs = nrow(dat), status = "ok")
    }
  })
}

# ---- 5. Run in parallel ----
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
clusterExport(cl, varlist = "run_chunk")

results_list <- parLapply(cl, data_chunks, run_chunk)
stopCluster(cl)

# ---- 6. Combine and save ----
rpt_all_sites <- bind_rows(results_list)

save(rpt_all_sites, file = "results/repeatability/results_all_sites_per_cpg.RData")

rpt_all_sites <- rpt_all_sites %>% mutate(changing = case_when(
  chr_pos %in% changing_cpg$chr_pos ~ "changing",
  TRUE ~ "non-changing"))


### plot ####

ggplot(rpt_all_sites, aes(x = R, fill = changing)) + geom_histogram() + scale_log10() + labs(title = "All sites")