### packages ###
pacman::p_load(lmerTest, lme4, tidyverse, data.table, glmmTMB)

### data ###
load("data/phenotypes/fulldata_complete_epi_withdates.RData") # all 450 captures
load("data/phenotypes/pheno_dif_prepost.RData") # prepost data only

load("/Users/vistor/Documents/Work/GitHub/PhD/grouse-phenotypes/data/cleandata/fulldata_exactdates/lekdata_fulldata_fromeyecomb.RData")

### trade off models ###
summary(lmerTest::lmer(attend ~ as.numeric(mass) + (1|site/id), data = all_pheno_epi)) #positive
summary(glmer(surv ~ as.numeric(attend) + (1|site/id), data = all_pheno_epi, family = "binomial")) #negative
summary(glmmTMB(MS ~ as.numeric(attend) + (1|site/id), data = all_pheno_epi, family = "poisson", ziformula=~1)) #positive

ggplot(all_pheno_epi, aes(x = as.numeric(mass), y = attend)) + geom_point(aes(col = id)) + 
  geom_smooth(method='lm') + theme(legend.position="none")

# with age 
summary(lmerTest::lmer(attend ~ as.numeric(mass) + age_cat + (1|site/id), data = all_pheno_epi)) #ns
summary(glmer(surv ~ as.numeric(attend) + age_cat + (1|site/id), data = all_pheno_epi, family = "binomial")) #negative
summary(glmmTMB(MS ~ as.numeric(attend) + age_cat + (1|site/id), data = all_pheno_epi, family = "poisson", ziformula=~1)) #positive

# on fitness
summary(glmer(surv  ~ as.numeric(attend) + as.numeric(mass) + (1|site/id), data = all_pheno_epi, family = "binomial")) #attend only negative
summary(glmmTMB(MS ~ as.numeric(attend) + as.numeric(mass)+ (1|site/id), data = all_pheno_epi, family = "poisson", ziformula=~1)) #both positive

### test power with prepost samples only ###
summary(lmerTest::lmer(attend ~ as.numeric(mass_pre) + (1|site/id), data = prepost_dif)) #ns
summary(glmer(surv ~ as.numeric(attend) + (1|site/id), data = prepost_dif, family = "binomial")) #ns
summary(glmmTMB(MS ~ as.numeric(attend) + (1|site/id), data = prepost_dif, family = "poisson", ziformula=~1)) #positive

summary(glmer(surv  ~ as.numeric(attend) + as.numeric(mass_dif) + (1|site/id), data = prepost_dif, family = "binomial")) #ns
summary(glmmTMB(MS ~ as.numeric(attend) + as.numeric(mass_dif)+ (1|site/id), data = prepost_dif, family = "poisson", ziformula=~1)) #na

summary(lmerTest::lmer(dist ~ as.numeric(mass) + (1|site/id), data = all_pheno_epi)) #neg
summary(lmerTest::lmer(fight ~ as.numeric(mass) + (1|site/id), data = all_pheno_epi)) #pos
summary(lmerTest::lmer(attend ~ as.numeric(mass) + (1|site/id), data = all_pheno_epi)) #pos

summary(glmmTMB(MS ~ as.numeric(attend) + as.numeric(mass)+ (1|site/id), data = all_pheno_epi, family = "poisson", ziformula=~1)) #na
summary(glmmTMB(MS ~ as.numeric(dist) + as.numeric(mass)+ (1|site/id), data = all_pheno_epi, family = "poisson", ziformula=~1)) #na
summary(glmmTMB(MS ~ as.numeric(fight) + as.numeric(mass)+ (1|site/id), data = all_pheno_epi, family = "poisson", ziformula=~1)) #na

