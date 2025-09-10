# ============================================================
# 6_Aim_1_Generate_Results.R
# Aim 1 analysis: Proteomic organ age distributions
# Author: James Groves
# Date: 2025-08-30
# ============================================================

#### SET-UP ####

# Packages

library(tidyverse)
library(readxl)
library(broom)
library(haven)
library(openxlsx)

# Directory

inp <- file.path(getwd(), "6. Aim 1 - Organ age distributions")
out <- file.path(getwd(), "6. Aim 1 - Organ age distributions")

#### LOAD DATA ####

# Load data 

load(file.path(inp,"dataset_inc.RData"))

# Adjust columns for display

organs <- c(
  "Conventional", "Brain", "Heart", "Lung",
  "Liver", "Kidney", "Immune", "Artery"
  )

dataset <- dataset %>% 
  rename(Conventional=ConventionalAge60,
  Brain=BrainAge60,
  Heart=HeartAge60,
  Lung=LungAge60,
  Liver=LiverAge60,
  Kidney=KidneyAge60,
  Immune=ImmuneAge60,
  Artery=ArteryAge60)

#### DESCRIPTIVES ####

# Organ age descriptive statistics 

organs_dist <- dataset %>%
  select(all_of(organs)) %>%
  pivot_longer(cols=everything(), 
  names_to="variable", values_to="value") %>%
  mutate(variable=factor(variable, levels=organs)) %>%
  group_by(variable) %>%
  summarise(mean=mean(value),
  sd=sd(value),min=min(value),max=max(value),
  .groups="drop")

organs_dist <- organs_dist %>% 
  rename(
    Organ=variable,
    "Age gap, mean" = mean,
    "Age gap, sd" = sd, 
    "Age gap, minimum" = min,
    "Age gap, maximum" = max
    )

#### CORRELATION & COMBINATIONS ####

# Correlation matrix

organ_cor_matrix <- as.data.frame(
  cor(
  dataset[, organs], 
  use="complete.obs", 
  method="pearson"))

# Extreme ageing combinations: 2 organs

comb_2 <- dataset %>% 
  filter(as.numeric(multiorgan60high)>2)

comb_2 <- comb_2 %>% 
  subset(select=c(
  id, BrainAge60high, ImmuneAge60high, ArteryAge60high,
  HeartAge60high, KidneyAge60high, LiverAge60high, 
  LungAge60high))

comb_2 <- comb_2 %>% rename(
  Brain = BrainAge60high, 
  Immune = ImmuneAge60high, 
  Artery = ArteryAge60high,
  Heart = HeartAge60high, 
  Kidney = KidneyAge60high, 
  Liver = LiverAge60high,
  Lung = LungAge60high)

ext_org <- c(
  "Brain", "Heart", "Lung", "Liver", 
  "Kidney", "Immune", "Artery"
  )

c_2 <- combn(ext_org, 2, simplify=FALSE)

for (c in c_2) {
  
  col <- paste(c, collapse="_")
  
  comb_2[[col]] <- ifelse(
  rowSums(comb_2[c]==1)==2,1,0)

  }

sums_2 <- colSums(
  comb_2[ , !names(comb_2) %in% ext_org])

sums2_df <- data.frame(
  Combination = names(sums_2),
  Individuals = sums_2
)

sums2_df <- sums2_df[-1, ]

rownames(sums2_df) <- NULL

sums2_df <- sums2_df %>%
  arrange(desc(Individuals))

# Extreme ageing combinations: 2 organs

comb_3 <- dataset %>% filter(as.numeric(multiorgan60high)>3)

comb_3 <- comb_3 %>% subset(
  select=c(id, BrainAge60high, ImmuneAge60high, ArteryAge60high,
  HeartAge60high, KidneyAge60high, LiverAge60high, LungAge60high))

comb_3 <- comb_3 %>% rename(
  Brain = BrainAge60high, 
  Immune = ImmuneAge60high, 
  Artery = ArteryAge60high,
  Heart = HeartAge60high, 
  Kidney = KidneyAge60high, 
  Liver = LiverAge60high,
  Lung = LungAge60high)

c_3 <- combn(ext_org, 3, simplify=FALSE)

for (c in c_3) {
  
  col <- paste(c, collapse="_")
  
  comb_3[[col]] <- ifelse(
rowSums(comb_3[c]==1)==3,1,0)
  
}

sums_3 <- colSums(
  comb_3[,!names(comb_3) %in% ext_org])

sums3_df <- data.frame(
  Combination = names(sums_3),
  Individuals = sums_3)

sums3_df <- sums3_df[-1, ]

rownames(sums3_df) <- NULL

sums3_df <- sums3_df %>%
  arrange(desc(Individuals))

# Aim 1 Supplementary Tables

res <- createWorkbook()

# ST3: Distributions
addWorksheet(res, "Organ age stats")
writeData(res, "Organ age stats", organs_dist)

# ST4: Correlation
addWorksheet(res, "Organ age correlations")
writeData(res, "Organ age correlations", organ_cor_matrix)

# ST5: Extreme ageing combinations
addWorksheet(res, "Extreme ageing combinations")
writeData(res, "Extreme ageing combinations", sums2_df)
writeData(res, "Extreme ageing combinations", 
sums3_df, startRow=25, startCol=1)

# Save
saveWorkbook(res, file.path(out, "Aim_1_ST.xlsx"))

# Source code files

sc <- createWorkbook()

# Ext Data Fig 2a
addWorksheet(sc, "Ext Data Fig 2a")
writeData(sc, "Ext Data Fig 2a", organ_cor_matrix)

# Ext Data Fig 2b
addWorksheet(sc, "Ext Data Fig 2b")
writeData(sc, "Ext Data Fig 2b", sums2_df)
writeData(sc, "Ext Data Fig 2b", 
sums3_df, startRow=25, startCol=1)

saveWorkbook(sc,"SC_Aim_1.xlsx", overwrite=TRUE)




