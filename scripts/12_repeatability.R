
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


# ---------------------------------------------------------------
# Per-CpG-site repeatability across pre/post timepoints
# Accounts for individuals with repeated pre/post measures
# across multiple years (non-independence / pseudoreplication)
# ---------------------------------------------------------------


library(tidyverse)
library(lme4)

# ---- 1. Set your column names here ------------------------------
INDIV_COL <- "id"
YEAR_COL  <- "year"
SITE_COL  <- "chr_pos"
VALUE_COL <- "methperc"

# ---- 2. Function to compute ICC-based repeatability for one site -
# Random effects: individual (captures pre/post consistency) and
# year (captures the fact that the same individual's repeat visits
# in different years are not independent of each other).
# If your year values never overlap between individuals (i.e. each
# individual only has their own years), consider nesting instead:
# (1 | individual_id/year)

calc_repeatability <- function(data, value_col, indiv_col, year_col) {
  
  form <- as.formula(
    paste0(value_col, " ~ 1 + (1|", indiv_col, ") + (1|", year_col, ")")
  )
  
  model <- tryCatch(
    suppressWarnings(suppressMessages(
      lmer(form, data = data, REML = TRUE)
    )),
    error = function(e) NULL
  )
  
  if (is.null(model)) return(tibble(R = NA_real_, var_indiv = NA_real_,
                                    var_year = NA_real_, var_resid = NA_real_,
                                    n_obs = nrow(data), converged = FALSE))
  
  vc <- as.data.frame(VarCorr(model))
  
  var_indiv <- vc$vcov[vc$grp == indiv_col]
  var_year  <- vc$vcov[vc$grp == year_col]
  var_resid <- vc$vcov[vc$grp == "Residual"]
  
  # guard against missing components (e.g. singular fit collapsing a term)
  var_indiv <- ifelse(length(var_indiv) == 0, 0, var_indiv)
  var_year  <- ifelse(length(var_year) == 0, 0, var_year)
  
  total_var <- var_indiv + var_year + var_resid
  R <- var_indiv / total_var
  
  tibble(
    R = R,
    var_indiv = var_indiv,
    var_year = var_year,
    var_resid = var_resid,
    n_obs = nrow(data),
    converged = is.null(model@optinfo$conv$lme4$messages)
  )
}

# ---- 3. Apply across all CpG sites --------------------------------
# nest() groups the data per site; map() fits a model to each nested
# tibble. This avoids bootstrapping (unlike rptR), so it scales to
# many thousands of sites.

repeatability_results <- prepost_long %>%                      # <-- your dataframe
  dplyr::rename(indiv = !!INDIV_COL, yr = !!YEAR_COL,
         site = !!SITE_COL, val = !!VALUE_COL) %>%
  group_by(site) %>%
  nest() %>%
  mutate(res = map(data, ~calc_repeatability(.x, "val", "indiv", "yr"))) %>%
  dplyr::select(site, res) %>%
  unnest(res) %>%
  dplyr::rename(!!SITE_COL := site)

# ---- 4. Inspect results -------------------------------------------
repeatability_results %>%
  arrange(desc(R)) %>%
  print(n = 20)

# Flag sites that failed to converge or had too few observations
repeatability_results %>%
  filter(!converged | is.na(R)) %>%
  nrow()

# ---- 5. (Optional) Get bootstrapped CIs for specific sites of interest
# rptR gives you confidence intervals and significance tests, but is
# much slower — use it only on a shortlist of top hits from step 4.
#
# library(rptR)
# site_data <- ap1 %>% filter(chr_pos == "chr1_12345")
# rpt(meth_pct ~ (1|individual_id) + (1|year), grname = c("individual_id","year"),
#     data = site_data, datatype = "Gaussian", nboot = 1000, npermut = 0)