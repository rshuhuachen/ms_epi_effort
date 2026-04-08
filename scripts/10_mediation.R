# load sig cpg
load(file="results/annotated/annotated_modeloutput_almostsig_all_annotated_priority.RData")

# the three antagonistic pleiotropy ones are:
cpg1 <- "ScEsiA3_18290__HRSCAF_21698_9737771"
cpg2 <- "ScEsiA3_16771__HRSCAF_19097_2512626"
cpg3 <- "ScEsiA3_21978__HRSCAF_26928_3652515"

# load data
### load packages
pacman::p_load(dplyr, data.table, tibble, performance, gaston, cowplot,
               parallel, performance, lmerTest, tidystats, glmmTMB, DHARMa)

### load data

load(file="results/modeloutput/changing/changing_sites_glmer.RData")

### load phenotypic data

load(file = "data/phenotypes/fulldata_complete_epi_withdates.RData")

load("data/phenotypes/pheno_dif_prepost.RData") ## differences in physiology

### methylation difference

load(file = "results/modeloutput/all_sites_deltameth.RData")

delta_meth_raw <- subset(delta_meth, chr_pos %in% changing_cpg$chr_pos)

#combine with site and fitness data
pheno_pre <- subset(all_pheno_epi, prepost=="pre")

delta_meth <- left_join(delta_meth_raw, unique(pheno_pre[,c("id", "year", "MS", "surv", "site", "attend", "fight", "dist")]), by = c("id", "year"))
delta_meth <- left_join(delta_meth, unique(prepost_dif[,c("id", "year", "mass_dif", "trypa_dif", "ig_dif", "hct_dif", "microf_dif")]), by = c("id", "year"))

# per one of the three cpgs, make model

data_1 <- subset(delta_meth, chr_pos == cpg1)
data_2 <- subset(delta_meth, chr_pos == cpg2)
data_3 <- subset(delta_meth, chr_pos == cpg3)


### mediation analysis ####

## cpg 1
m_1_indirect_a <- lmerTest::lmer(delta_meth ~ scale(attend) + age + methperc_pre + (1|site), data = data_1)
m_1_indirect_b <- glmmTMB(surv ~ delta_meth + age + (1|site), data = data_1, family = "binomial", REML=FALSE)
m_1_direct <- glmmTMB(surv ~ scale(attend) + delta_meth + age + (1|site), data = data_1, family = "binomial", REML=FALSE)

summary(m_1_indirect_a) # attend sig
summary(m_1_indirect_b) # delta meth sig
summary(m_1_direct) #almost sig delta meth

coef_m1_in_a <- as.data.frame(summary(m_1_indirect_a)$coef)
coef_m1_in_b <- as.data.frame(summary(m_1_indirect_b)$coefficients$cond)
coef_m1_in_a$Estimate[2]*coef_m1_in_b$Estimate[2] # indirect

as.data.frame(summary(m_1_direct)$coefficients$cond)$Estimate[2] # direct

## cpg 2
load(file = "data/processed/delta_meth_MS_sub.RData")
data_2 <- subset(delta_meth_sub_MS, chr_pos == cpg2)
data_2 <- left_join(data_2, unique(pheno_pre[,c("id", "year", "surv")]), by = c("id", "year"))

m_2_indirect_a <- lmerTest::lmer(delta_meth ~ scale(MS) + age + methperc_pre + (1|site), data = data_2)
m_2_indirect_b <- glmmTMB(surv ~ delta_meth + age + (1|site), data = data_2, family = "binomial", REML=FALSE)
m_2_direct <- glmmTMB(surv ~ scale(MS) + delta_meth + age + (1|site), data = data_2, family = "binomial", REML=FALSE)

summary(m_2_indirect_a) # MS non sig
summary(m_2_indirect_b) # delta meth sig
summary(m_2_direct) # delta meth sig

coef_m1_in_a <- as.data.frame(summary(m_2_indirect_a)$coef)
coef_m1_in_b <- as.data.frame(summary(m_2_indirect_b)$coefficients$cond)
coef_m1_in_a$Estimate[2]*coef_m1_in_b$Estimate[2] # indirect

as.data.frame(summary(m_2_direct)$coefficients$cond)$Estimate[2] # direct
