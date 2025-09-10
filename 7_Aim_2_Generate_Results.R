# ============================================================
# 7_Aim_2_Generate_Results.R
# Aim 2 analysis: Survival analyses
# Author: James Groves
# Date: 2025-08-30
# ============================================================

#### SET-UP ####

# Packages

library(tidyverse)
library(readxl)
library(broom)
library(haven)
library(survival)
library(openxlsx)

# Directory

inp <- file.path(getwd(), "6. Aim 2 - Survival")
out <- file.path(getwd(), "6. Aim 2 - Survival")

# Load data 

load("dataset_inc.RData")

organs <- c(
"Conventional", "Brain", 
"Heart", "Lung", 
"Liver", "Kidney", 
"Immune", "Artery")

#### COX MODELS ####

# Organ ages

for (organ in organs) {
  
  agevariable <- paste0(organ, "Age60high")
  
  formula1 <- as.formula(paste(
  "Surv(death_yearssincesample, death) ~ ", 
  agevariable, " + Sex_F + Age60"))
  model <- coxph(formula1, data=dataset)
  
  tidy_model <- tidy(model)
  conf_intervals <- confint(model)
  tidy_model$lci <- exp(conf_intervals[,1])
  tidy_model$uci <- exp(conf_intervals[,2])
  tidy_model$hazardratio <- exp(tidy_model$estimate)
  
  model_name <- paste0("mort_", organ)
  ageoutput <- tidy_model[1,]
  assign(model_name, ageoutput)
  
}

mort_results <- bind_rows(
  "Conventional" = mort_Conventional,
  "Brain" = mort_Brain,
  "Heart" = mort_Heart,
  "Lung" = mort_Lung,
  "Liver" = mort_Liver,
  "Kidney" = mort_Kidney,
  "Immune" = mort_Immune,
  "Artery" = mort_Artery,
  .id="Organ")

mort_results <- mort_results %>% 
  subset(select=-c(term))

# Multi-organ extreme ageing

cox_model <- coxph(
Surv(death_yearssincesample, death) ~ as.factor(
multiorgan60high4max) + Sex_F + Age60, data=dataset)

tidy_model <- tidy(cox_model)
conf_intervals <- confint(cox_model)
tidy_model$lci <- exp(conf_intervals[,1])
tidy_model$uci <- exp(conf_intervals[,2])
tidy_model$hazardratio <- exp(tidy_model$estimate)
tidy_model <- tidy_model %>% slice(1:4)

mort_multiorgan <- data.frame(tidy_model)

mort_multiorgan <-  mort_multiorgan %>%
  mutate(term=ifelse(term=="as.factor(multiorgan60high4max)1",
"Extreme ageing in 1 organ",term))

mort_multiorgan <-  mort_multiorgan %>%
  mutate(term=ifelse(term=="as.factor(multiorgan60high4max)2",
"Extreme ageing in 2 organs",term))

mort_multiorgan <-  mort_multiorgan %>%
  mutate(term=ifelse(term=="as.factor(multiorgan60high4max)3",
"Extreme ageing in 3 organs",term))

mort_multiorgan <-  mort_multiorgan %>%
  mutate(term=ifelse(term=="as.factor(multiorgan60high4max)4",
"Extreme ageing in 4 or more organs",term))

multiorgan_death_counts <-
  dataset %>%
  group_by(multiorgan60high4max, death) %>% 
  summarise(n=n(), .groups="drop") %>%
  pivot_wider(names_from = death, values_from=n, 
  names_prefix = "death")

multiorgan_death_counts <-  multiorgan_death_counts %>%
rename("Number of organs with extreme ageing"="multiorgan60high4max",
"Survived"="death0", "Died"="death1")

#### OUTPUT FILES ####

# Aim 2 Supplementary Tables

res <- createWorkbook()

# ST6: Extreme organ ageing & mortality 
addWorksheet(res, "Extreme ageing & mortality")
writeData(res, "Extreme ageing & mortality", mort_ext_results)

# ST7: Multi-organ extreme ageing & mortality 
addWorksheet(res, "Multi-organ ageing & mortality")
writeData(
res, "Multi-organ ageing & mortality", mort_multiorgan)
writeData(
res, "Multi-organ ageing & mortality", 
multiorgan_death_counts, startRow=10, startCol=1)

# Save
saveWorkbook(
res, 
file.path(out, "Aim_2_ST.xlsx", overwrite=TRUE)
)

# Source code files

sc <- createWorkbook()

# Fig 2c
addWorksheet(sc, "Fig 2c")
writeData(sc, "Fig 2c", mort_ext_results)

# Fig 2e
addWorksheet(sc, "Fig 2e")
writeData(sc, "Fig 2e", mort_multiorgan)
writeData(sc, "Fig 2e", multiorgan_death_counts, 
startRow=10, startCol=1)

# Save
saveWorkbook(
sc, 
file.path(out, "SC_Aim_2.xlsx", overwrite=TRUE)
)


