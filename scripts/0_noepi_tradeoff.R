### load packages ####
pacman::p_load(tidyverse, data.table, lmerTest, glmmTMB, brms, bayesplot)

### load data ####

load("data/phenotypes/fulldata_complete_epi_withdates.RData")
load("data/phenotypes/pheno_dif_prepost.RData")

### test for trade-offs ####

data <- left_join(all_pheno_epi, prepost_dif[,c("id", "year", "mass_dif", "microf_dif", "trypa_dif", "ig_dif", "hct_dif")], by = c("year", "id"))
data <- data %>% mutate(age_year = as.factor(case_when(Core == "Core" ~ year - born,
                                                        Core == "No core" ~ NA)),
                                        age = as.factor(case_when(Core == "Core" & (year - born > 1) ~ "Adult",
                                                        Core == "Core" & (year - born == 1) ~ "Yearling",
                                                        Core == "No core" ~ "Adult")))

data$mass <- as.numeric(data$mass)    

#### formal models in brms ####

#### model 1: full dataset - attendance

## 1) reproductive effort affected by available resources?
m_1_a <- bf(attend ~ scale(mass) + (1|site/id))

## 2) fitness affected by loss in resources and/or reproductive effort?
m_2_su_a <- bf(surv ~ scale(attend) + scale(mass) + (1|site/id), family = "bernoulli")
m_2_ms_a <- bf(MS ~ scale(attend) + scale(mass) + (1|site/id), family = "poisson")

sem_attend <- m_1_a + m_2_su_a + m_2_ms_a

fit_attend <- brm(sem_attend, data = data, cores = 8, control = list(adapt_delta = 0.99, max_treedepth = 15),
           prior = prior(normal(0,10), class = b), iter = 200000, thin = 500, warmup = 50000)

save(fit_attend, file="results/modeloutput/brms_fit_attend.RData")

### model 2: full dataset - fighting
## 1) reproductive effort affected by available resources?
m_1_f <- bf(fight ~ scale(mass) + scale(mass) + (1|site/id))

## 2) fitness affected by loss in resources and/or reproductive effort?
m_2_su_f <- bf(surv ~ scale(fight) + scale(mass) + (1|site/id), family = "bernoulli")
m_2_ms_f <- bf(MS ~ scale(fight) + scale(mass) + (1|site/id), family = "poisson")

sem_fight <- m_1_f + m_2_su_f + m_2_ms_f

fit_fight <- brm(sem_fight, data = data, cores = 8, control = list(adapt_delta = 0.99, max_treedepth = 15),
           prior = prior(normal(0,10), class = b), iter = 200000, thin = 500, warmup = 50000)

save(fit_fight, file="results/modeloutput/brms_fit_fight.RData")

### model 3: full dataset - dist 
## 1) reproductive effort affected by available resources?
m_1_d <- bf(dist ~ scale(mass) + (1|site/id))

## 2) fitness affected by loss in resources and/or reproductive effort?
m_2_su_d <- bf(surv ~ scale(dist) + scale(mass) + (1|site/id), family = "bernoulli")
m_2_ms_d <- bf(MS ~ scale(dist) + scale(mass) + (1|site/id), family = "poisson")

sem_dist <- m_1_d + m_2_su_d + m_2_ms_d

fit_dist <- brm(sem_dist, data = data, cores = 8, control = list(adapt_delta = 0.99, max_treedepth = 15),
           prior = prior(normal(0,10), class = b), iter = 200000, thin = 500, warmup = 50000)

save(fit_dist, file="results/modeloutput/brms_fit_dist.RData")

#### model 2: changes

## 1) loss in resources affected by reproductive effort?
m_1 <- bf(mass_dif ~ scale(attend) + scale(fight) + scale(dist) + (1|site/id))

## 2) loss in resources affects fitness?
m_2_su_m <- bf(surv ~ scale(mass_dif) + (1|site/id), family = "bernoulli")
m_2_ms_m <- bf(MS ~ scale(mass_dif) + (1|site/id), family = "poisson")

sem_mass <- m_1 + m_2_su_m + m_2_ms_m 

fit_mass <- brm(sem_mass, data = data, cores = 8, control = list(adapt_delta = 0.99, max_treedepth = 15),
           prior = prior(normal(0,10), class = b), iter = 200000, thin = 500, warmup = 50000)

save(fit_mass, file="results/modeloutput/brms_fit_mass.RData")