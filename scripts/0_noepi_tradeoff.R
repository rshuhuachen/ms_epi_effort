### load packages
pacman::p_load(tidyverse, data.table, lmerTest, glmmTMB)

### load data

load("data/phenotypes/fulldata_complete_epi_withdates.RData")
load("data/phenotypes/pheno_dif_prepost.RData")

### test for trade-offs

data <- left_join(all_pheno_epi, prepost_dif[,c("id", "year", "mass_dif", "microf_dif", "trypa_dif", "ig_dif", "hct_dif")], by = c("year", "id"))
data <- data %>% mutate(age_year = as.factor(case_when(Core == "Core" ~ year - born,
                                                        Core == "No core" ~ NA)),
                                        age = as.factor(case_when(Core == "Core" & (year - born > 1) ~ "Adult",
                                                        Core == "Core" & (year - born == 1) ~ "Yearling",
                                                        Core == "No core" ~ "Adult")))

data$mass <- as.numeric(data$mass)    

### reproduction
## reproduction vs survival
summary(glmer(surv ~ attend + age + (1|year) + (1|site/id), data = data, family = "binomial")) #higher attendance, lower survival probability
summary(glmer(surv ~ fight + age + (1|year) + (1|site/id), data = data, family = "binomial")) #ns
summary(glmer(surv ~ dist + age + (1|year) + (1|site/id), data = data, family = "binomial")) #ns

## quantifying the 'cost' of reproductive effort
summary(lmerTest::lmer(mass_dif ~ attend + age + (1|year) + (1|site/id), data = data)) #ns
summary(lmerTest::lmer(mass_dif ~ fight + age + (1|year) + (1|site/id), data = data)) #ns
summary(lmerTest::lmer(mass_dif ~ dist + age + (1|year) + (1|site/id), data = data)) #more central, higher loss of body mass

## resources vs reproductive investment
summary(lmerTest::lmer(attend ~ mass + age + (1|year) + (1|site/id), data = data)) # higher starting mass, higher attendance
summary(lmerTest::lmer(fight ~ mass + age + (1|year) + (1|site/id), data = data))
summary(lmerTest::lmer(dist ~ mass + age + (1|year) + (1|site/id), data = data)) # higher starting mass, more central

## is the 'cost' of reproduction predictive of survival?
summary(glmer(surv ~ mass_dif + age + (1|year) + (1|site/id), data = data, family = "binomial")) #ns but lower sample size
summary(glmer(surv ~ scale(mass)*scale(dist) + age + (1|year) + (1|site/id), data = data, family = "binomial")) 
summary(glmer(MS ~ scale(mass)*scale(dist) + age + (1|year) + (1|site/id), data = data, family = "poisson")) 


### maintenance
## survival vs reproduction
summary(glmer(MS ~ scale(ig_dif) + age + (1|year) + (1|site/id), data = data, family = "poisson")) #ns
summary(glmer(MS ~ scale(microf_dif) + age + (1|year) + (1|site/id), data = data, family = "poisson")) #more parasites, higher mating success
summary(glmer(MS ~ scale(trypa_dif) + age + (1|year) + (1|site/id), data = data, family = "poisson")) #more parasites, higher mating success
summary(glmer(MS ~ scale(hct_dif) + age + (1|year) + (1|site/id), data = data, family = "poisson")) #almost sig negative

## loss of maintenance because of reproduction?
summary(lmerTest::lmer(ig_dif ~ scale(dist) + age + (1|year) + (1|site/id), data = data)) #ns
summary(lmerTest::lmer(ig_dif ~ scale(attend) + age + (1|year) + (1|site/id), data = data)) #ns

summary(lmerTest::lmer(microf_dif ~ scale(dist) + age + (1|year) + (1|site/id), data = data)) #
summary(lmerTest::lmer(microf_dif ~ scale(attend) + age + (1|year) + (1|site/id), data = data)) #

summary(lmerTest::lmer(trypa_dif ~ scale(dist) + age + (1|year) + (1|site/id), data = data)) #
summary(lmerTest::lmer(trypa_dif ~ scale(attend) + age + (1|year) + (1|site/id), data = data)) #

summary(lmerTest::lmer(hct_dif ~ scale(dist) + age + (1|year) + (1|site/id), data = data)) # almost sig positive, more central, lower increase in HCT
summary(lmerTest::lmer(hct_dif ~ scale(attend) + age + (1|year) + (1|site/id), data = data)) #

## quantifying the 'cost' of maintenance
summary(lmerTest::lmer(mass_dif ~ scale(ig_dif) + age + (1|year) + (1|site/id), data = data)) #sig negative; higher increase in IgG, more loss in mass
ggplot(data, aes(ig_dif, mass_dif)) + geom_point() + geom_smooth(method="lm")
summary(lmerTest::lmer(mass_dif ~ scale(ig_dif) + age + (1|year) + (1|site/id), data = subset(data, ig_dif > -2142964 & ig_dif < 2549625)) ) #ns minus 2 outliers

summary(lmerTest::lmer(mass_dif ~ microf_dif + age + (1|year) + (1|site/id), data = data)) #sig negative; higher increase in parasites, more loss in mass
ggplot(data, aes(microf_dif, mass_dif)) + geom_point() + geom_smooth(method="lm")

summary(lmerTest::lmer(mass_dif ~ trypa_dif + age + (1|year) + (1|site/id), data = data)) #sig negative;  higher increase in parasites, more loss in mass
ggplot(data, aes(trypa_dif, mass_dif)) + geom_point() + geom_smooth(method="lm")

summary(lmerTest::lmer(mass_dif ~ hct_dif + age + (1|year) + (1|site/id), data = data)) #sig positive
ggplot(data, aes(hct_dif, mass_dif)) + geom_point() + geom_smooth(method="lm")

## quantifying the 'cost' of maintenance vs reproduction
summary(lmerTest::lmer(mass_dif ~ scale(dist) + scale(ig_dif) + age + (1|year) + (1|site/id), data = data)) #both ish but dist more
summary(lmerTest::lmer(mass_dif ~ scale(dist) + scale(microf_dif) + age + (1|year) + (1|site/id), data = data)) #only dist
summary(lmerTest::lmer(mass_dif ~ scale(dist) + scale(trypa_dif) + age + (1|year) + (1|site/id), data = data)) #both almost sig
summary(lmerTest::lmer(mass_dif ~ scale(dist) + scale(hct_dif) + age + (1|year) + (1|site/id), data = data)) #only hct sig (positive)

## affording cost of maintenance
summary(lmerTest::lmer(ig_dif ~ scale(mass) + age + (1|year) + (1|site/id), data = data)) #ns
summary(lmerTest::lmer(microf_dif ~ scale(mass) + age + (1|year) + (1|site/id), data = data)) #ns
summary(lmerTest::lmer(trypa_dif ~ scale(mass) + age + (1|year) + (1|site/id), data = data)) #ns
summary(lmerTest::lmer(hct_dif ~ scale(mass) + age + (1|year) + (1|site/id), data = data)) #ns

## relationship igg and parasites?
summary(lmerTest::lmer(scale(ig_dif) ~ scale(microf_dif) + age + (1|year) + (1|site/id), data = data)) #sig positive
ggplot(data, aes(ig_dif, microf_dif)) + geom_point() + geom_smooth(method="lm")

summary(lmerTest::lmer(ig_dif ~ scale(trypa_dif) + age + (1|year) + (1|site/id), data = data)) #ns

ggplot(data, aes(ig_dif, trypa_dif)) + geom_point() + geom_smooth(method="lm")

## physio and survival
summary(glmer(surv ~ scale(ig_dif) + age + (1|year) + (1|site/id), data = data, family = "binomial")) #ns
summary(glmer(surv ~ scale(microf_dif) + age + (1|year) + (1|site/id), data = data, family = "binomial")) #ns
summary(glmer(surv ~ scale(trypa_dif) + age + (1|year) + (1|site/id), data = data, family = "binomial")) #ns
summary(glmer(surv ~ scale(hct_dif) + age + (1|year) + (1|site/id), data = data, family = "binomial")) #ns


