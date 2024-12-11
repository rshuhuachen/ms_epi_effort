## load packages
pacman::p_load(readxl, tidyverse, lme4, data.table)

## load data 
injury <- read_excel("data/phenotypes/Injury data.xlsx")
load(file = "data/phenotypes/fulldata_complete_epi_withdates.RData")

## format injury data

data <- left_join(all_pheno_epi, injury[,c("Date", "Ring", "weight", "injury_right_ey",
                                           "injury_left_ec", "injury_beak", "injury_neck", "injury_belly")], 
                  by = c("id" = "Ring", "fulldate" = "Date"))


prepost <- subset(data, !is.na(prepost))
