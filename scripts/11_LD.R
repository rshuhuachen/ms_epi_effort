table_sig <- fread("results/annotated/annotated_modeloutput_almostsig_annotated.tsv")
load(file="results/annotated/annotated_modeloutput_almostsig_all_annotated_priority.RData")

load(file = "results/modeloutput/all_sites_deltameth.RData")

# extract cpg number 1 with antagonistic pleiotropy pattern

ap1_id <- subset(almostsig_all_annotated_id, similar == "BEST1") 
ap1 <- subset(delta_meth, chr_pos == ap1_id$chr_pos)
ap1_w <- ap1 %>%
  pivot_wider(
    id_cols = c(id, year, age),  # columns to keep as-is per row
    names_from = chr_pos,
    values_from = delta_meth
  )
names(ap1_w)
names(ap1_w) <- c("id", "year", "age", "cpg1", "cpg2", "cpg3")

summary(lmerTest::lmer(cpg1 ~ cpg2 + (1|year), data=ap1_w))
summary(lmerTest::lmer(cpg2 ~ cpg3 + (1|year), data=ap1_w))
summary(lmerTest::lmer(cpg1 ~ cpg3 + (1|year), data=ap1_w))


# extract cpg number 2 with antagonistic pleiotropy pattern

ap2_id <- subset(almostsig_all_annotated_id, similar == "NFIC") 
ap2 <- subset(delta_meth, chr_pos == ap2_id$chr_pos)
ap2_w <- ap2 %>%
  pivot_wider(
    id_cols = c(id, year, age),  # columns to keep as-is per row
    names_from = chr_pos,
    values_from = delta_meth
  )
names(ap2_w)
names(ap2_w) <- c("id", "year", "age", "cpg1", "cpg2", "cpg3", "cpg4")

summary(lmerTest::lmer(cpg1 ~ cpg2 + (1|year), data=ap2_w))
summary(lmerTest::lmer(cpg1 ~ cpg3 + (1|year), data=ap2_w))
summary(lmerTest::lmer(cpg1 ~ cpg4 + (1|year), data=ap2_w))
summary(lmerTest::lmer(cpg2 ~ cpg3 + (1|year), data=ap2_w))
summary(lmerTest::lmer(cpg2 ~ cpg4 + (1|year), data=ap2_w))
summary(lmerTest::lmer(cpg3 ~ cpg4 + (1|year), data=ap2_w))
