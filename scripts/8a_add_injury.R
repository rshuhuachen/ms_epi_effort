## load packages
pacman::p_load(readxl, tidyverse, lme4, data.table, ggbiplot)

## load data 
injury <- read_excel("data/phenotypes/Injury data.xlsx")
load(file = "data/phenotypes/fulldata_complete_epi_withdates.RData")
source("scripts/plotting_theme.R")

### pca for injury data

pca <- prcomp(injury[,c(4:8)], center = T, scale = T)
attributes(pca)

pca$scale
pca$rotation # lower PC
summary(pca)

pc <- as.data.frame(pca$x)

injury <- cbind(injury, pc)

### PC1 = eyec, beak, neck (pecking injuries)
### PC2 = belly (kicking injuries)

ggplot(injury, aes(x = injury_right_ey, y = PC1)) + geom_point(size=3) # lower pc1 value, higher injuries
ggplot(injury, aes(x = injury_left_ec, y = PC1)) + geom_point(size=3) # lower pc1 value, higher injuries
ggplot(injury, aes(x = injury_beak, y = PC1)) + geom_point(size=3) # lower pc1 value, higher injuries
ggplot(injury, aes(x = injury_neck, y = PC1)) + geom_point(size=3) # lower pc1 value, higher injuries
ggplot(injury, aes(x = injury_belly, y = PC2)) + geom_point(size=3) # higher pc2 value, higher injuries

## add injury data to pheno data

data <- left_join(all_pheno_epi, injury[,c("Date", "Ring", "weight", "injury_right_ey",
                                           "injury_left_ec", "injury_beak", "injury_neck", "injury_belly", "PC1", "PC2")], 
                  by = c("id" = "Ring", "fulldate" = "Date"))


prepost <- subset(data, !is.na(prepost))

save(prepost, file = "data/phenotypes/prepost_with_injury.RData")

### are injuries related to effort, mass, etc???

# effort is recorded in pre-lekking, but injuries in post-lekking
pre <- prepost %>% select(id, site, Core, prepost, year, mass:surv) %>% filter(prepost=="pre") %>% select(-c(prepost))
post <- prepost %>% select(id, site, Core, prepost, year, mass,froh1,ig,hct:PC2) %>% filter(prepost=="post") %>% select(-c(prepost))

merge <- left_join(pre, post, by = c("id", "site", "Core", "year"), suffix=c("_pre", "_post"))

## effort
ggplot(merge, aes(x = attend, y = PC1)) + geom_point(size=3) + geom_smooth(method='lm')
ggplot(merge, aes(x = attend, y = PC2)) + geom_point(size=3) + geom_smooth(method='lm')

summary(lmerTest::lmer(PC1 ~ attend + (1|year) + (1|site), data = merge))
summary(lmerTest::lmer(PC2 ~ attend + (1|year) + (1|site), data = merge))

ggplot(merge, aes(x = fight, y = PC1)) + geom_point(size=3) + geom_smooth(method='lm')
ggplot(merge, aes(x = fight, y = PC2)) + geom_point(size=3) + geom_smooth(method='lm') 

summary(lmerTest::lmer(PC1 ~ fight + (1|year) + (1|site), data = merge))
summary(lmerTest::lmer(PC2 ~ fight + (1|year) + (1|site), data = merge)) # sig, more fighting, more belly injuries

ggplot(merge, aes(x = dist, y = PC1)) + geom_point(size=3) + geom_smooth(method='lm')
ggplot(merge, aes(x = dist, y = PC2)) + geom_point(size=3) + geom_smooth(method='lm')

summary(lmerTest::lmer(PC1 ~ dist + (1|year) + (1|site), data = merge))
summary(lmerTest::lmer(PC2 ~ dist + (1|year) + (1|site), data = merge)) 

## mass
merge$mass_pre <- as.numeric(as.character(merge$mass_pre))
merge$mass_post <- as.numeric(as.character(merge$mass_post))

# pre
ggplot(merge, aes(x = mass_pre, y = PC1)) + geom_point(size=3) + geom_smooth(method='lm')
ggplot(merge, aes(x = mass_pre, y = PC2)) + geom_point(size=3) + geom_smooth(method='lm')

summary(lmerTest::lmer(PC1 ~ mass_pre + (1|year) + (1|site), data = merge)) # sig neg, higher mass, more face injuries
summary(lmerTest::lmer(PC2 ~ mass_pre + (1|year) + (1|site), data = merge)) # ns

# post
ggplot(merge, aes(x = mass_post, y = PC1)) + geom_point(size=3) + geom_smooth(method='lm')
ggplot(merge, aes(x = mass_post, y = PC2)) + geom_point(size=3) + geom_smooth(method='lm')

summary(lmerTest::lmer(PC1 ~ mass_post + (1|year) + (1|site), data = merge)) # ns
summary(lmerTest::lmer(PC2 ~ mass_post + (1|year) + (1|site), data = merge)) # ns

## igg
ggplot(merge, aes(x = ig, y = PC1)) + geom_point(size=3) + geom_smooth(method='lm')
ggplot(merge, aes(x = ig, y = PC2)) + geom_point(size=3) + geom_smooth(method='lm')

summary(lmerTest::lmer(PC1 ~ ig + (1|year) + (1|site), data = merge)) # sig neg 0.05, higher igg post, more face injuries
summary(lmerTest::lmer(PC2 ~ ig + (1|year) + (1|site), data = merge)) # ns

# surv
summary(glmer(surv ~ PC1 + (1|year) + (1|site), family = "binomial", data = merge)) # ns
summary(glmer(surv ~ PC2 + (1|year) + (1|site), family = "binomial", data = merge)) # ns

# ms
ggplot(merge, aes(x = PC1, y = MS)) + geom_point(size=3) + geom_smooth(method='lm')
summary(glmer(MS ~ PC1 + (1|year) + (1|site), family = "poisson", data = merge)) # very neg sig, more injuries, more MS
summary(glmmTMB::glmmTMB(MS ~ PC1 + (1|year) + (1|site), family = "poisson", ziformula = ~1,data = merge)) # very neg sig, more injuries, more MS
summary(glmmTMB::glmmTMB(MS ~ PC2, family = "poisson", data = merge)) # ns
summary(glmmTMB::glmmTMB(MS ~ PC1 + (1|year) + (1|site), family = "poisson", ziformula = ~1,data = merge)) # very neg sig, more injuries, more MS
