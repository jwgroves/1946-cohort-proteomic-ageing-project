# ============================================================
# 4_Load_variables_and_clean_data.R
# Data cleaning for life course, organ age and mortality data
# Author: James Groves
# Date: 2025-08-30
# ============================================================

#### SET-UP ####

# Packages

library(tidyverse)
library(readxl)
library(ggsci)
library(broom)
library(haven)
library(purrr)

# Directory

inp <- file.path(getwd(), "3. Load variables & clean data")

#### ORGAN AGES ####

# Load organ ages

ages_long <- read.csv(
  file.path(inp, "organages.csv"))

# Manage columns 

ages_long <- ages_long %>% 
  rename(
    id=ID
    ) 

ages_long <- ages_long %>% 
  subset(
  select=c(
    id, Sex_F, Age, AgeGap, Organ)) %>%
  rename(
    
  )

ages_long <- ages_long %>% 
  filter(
    Organ=="Conventional" | 
    Organ=="Brain" | 
    Organ=="Heart" | 
    Organ=="Kidney" | 
    Organ=="Liver" | 
    Organ=="Immune" | 
    Organ=="Artery" | 
    Organ=="Lung")

ages_long <- ages_long %>%
  rename(
    AgeGap60 = AgeGap,
    Age60 = Age
  )

# Wide format

organs <- c("Conventional", "Brain", "Immune", "Artery", 
"Heart", "Kidney", "Liver", "Lung")

for (organ in organs) {
  
col <- paste(organ, "Age60", sep="")

ages_long[[col]] <- 
ifelse(ages_long$Organ==organ, 
       ages_long$AgeGap60, 
       NA)

}

ages <- ages_long %>%
group_by(id) %>%
summarize_all(~first(na.omit(.))) %>%
  subset(select=-c(
    AgeGap60, Organ
  ))

# Classify extreme ageing

for (organ in organs) {

  col <- paste(organ, "Age60", sep="")
    
  threshold <- quantile(ages[[col]], 0.9)
  
  extreme <- ifelse(
    ages[[col]]>=threshold,
    1,
    0)
  
  ages[[paste0(col, "high")]] <- extreme
  
}

# Extreme ageing counts

ext_cols <- c("BrainAge60high", "ArteryAge60high", 
              "HeartAge60high", "ImmuneAge60high", 
              "KidneyAge60high", "LiverAge60high", 
              "LungAge60high")

ages$multiorgan60high <- rowSums(ages[ext_cols])

ages$multiorgan60high3max <- ages$multiorgan60high

ages$multiorgan60high3max[ages$multiorgan60high3max==6] <- 3
ages$multiorgan60high3max[ages$multiorgan60high3max==5] <- 3
ages$multiorgan60high3max[ages$multiorgan60high3max==4] <- 3

ages$multiorgan60high4max <- ages$multiorgan60high

ages$multiorgan60high4max[ages$multiorgan60high4max==6] <- 4
ages$multiorgan60high4max[ages$multiorgan60high4max==5] <- 4

ages$multiorganageing <- ifelse(ages$multiorgan60high>1,1,0)

# Set classes

ages$multiorgan60high <- as.factor(ages$multiorgan60high)

ages$multiorgan60high3max <- as.factor(ages$multiorgan60high3max)

ages$multiorgan60high4max <- as.factor(ages$multiorgan60high4max)

# Any missing data?

anyNA(ages)

# Record included participants

inc <- ages$id

#### LIFE COURSE EXPOSURES ####


# Load data

lc1 <- read.csv("[Path1]")
lc2 <- read.csv("[Path2]")
lc3 <- read.csv("[Path3]")
lc4 <- read.csv("[Path4]")

# Remove duplicate columns

lc1 <- lc1 %>% subset(select=-c(inf))
lc2 <- lc2 %>% subset(select=-c(sex, inf, cmd6clas))
lc3 <- lc3 %>% subset(select=-c(sex, inf))
lc4 <- lc4 %>% subset(select=-c(sex, inf, cmd53or4))

# Merge files

lc <- Reduce(function(x,y) merge(
x,y, by = "nshdid_NTAG1", all = TRUE), 
list(lc1, lc2, lc3, lc4))

# Subset data

lc <- lc %>% 
  rename(
    id=nshdid_NTAG1)

vars <- c(
  "id", "lowbwt", "chsc", "lhqr", "bmi61u",
  "packyr20", "packyr2025", "packyr2531", "packyr3136", 
  "packyr3643", "packyr4353", "packyr5363", "avalc82u", 
  "avalc89u", "avalc99u", "avalc09u","exer82", 
  "exer89x", "exer99x", "exer09x", "sc1553")

exp <- lc %>% subset(select=c(vars))

# Remove NAs

exp[exp<0] <- NA

exp$chsc[exp$chsc == 9] <- NA

exp$sc1553[exp$sc1553 == 999] <- NA

exp$lhqr[exp$lhqr == 10] <- NA
exp$lhqr[exp$lhqr == 9] <- NA

exp$bmi61u[exp$bmi61u == 9999] <- NA
exp$bmi61u[exp$bmi61u == 7777] <- NA

exp$exer82[exp$exer82 == 9] <- NA

exp$exer89x[exp$exer89x == 9] <- NA

exp$exer99x[exp$exer99x == 88] <- NA
exp$exer99x[exp$exer99x == 7] <- NA

exp$exer09x[exp$exer09x == 999] <- NA
exp$exer09x[exp$exer09x == 998] <- NA
exp$exer09x[exp$exer09x == 777] <- NA

# Check plausible max/mins

maxs <- sapply(exp, max, na.rm=TRUE)

mins <- sapply(exp, min, na.rm=TRUE)

maxmins <- data.frame(
  Var=names(exp), 
  min=mins, 
  max=maxs)

# Add adversity measures

adv <- read.csv("[Path5]")

adv <- adv[ , -1]

# Merge 

exp_final <- merge(exp, adv, by="id")

# Set classes

adv_var <- c("CROW15y", "amenities", "childSES", 
  "matsep", "divorce16", "peers", "phy16",
  "finance36", "employed36", "amenities26", 
  "crow26", "social36", "divorce36", "work36", 
  "work4353", "finance4353", "housecondition43", 
  "support4353", "lostcontact4353","chdconflict4353", 
  "divorce4353", "social4353","finance63", 
  "support63", "lostcontact63", "chdconflict63", 
  "divorce63", "social63", "work63")

cont <- c(
  "bmi61u", "packyr20", "packyr2025", 
  "packyr2531", "packyr3136", "packyr3643", 
  "packyr4353", "packyr5363", "avalc82u", 
  "avalc89u", "avalc99u", "avalc09u") 

cat <- c(
  "lowbwt", "chsc", "lhqr", 
  "exer82", "exer89x", "exer99x", 
  "exer09x", "sc1553", 
  adv_var)

exp_final[cat] <- lapply(
  exp_final[cat], 
  as.factor)

# Remove excluded

exp_final <- exp_final %>%
  filter(id %in% inc)

#### MORTALITY ####

# Load

mort <- read.csv("[Path6]")

# Manage columns

mort <- mort %>% 
  rename(
  death=DTH2CEN2, 
  id=ï..id, 
  months=DTH2DTM2)

mort <- mort %>% 
  subset(
    select=c(id, death, months))

# Set classes

mort$months <- as.numeric(mort$months)

# Reset max months

mort$months[mort$months==995] <- 940

# Remove excluded

mort_final <- mort %>%
  filter(id %in% inc)

#### OTHER VARIABLES ####

oth_vars <- c(
  "id", "disa15x", "chron19tot15x", 
  "summhealth15x", "bmi66u", "bmi72u", 
  "bmi82u", "bmi89u", "bmi99u", 
  "bmi09")

oth_exp <- lc %>% 
  subset(select=oth_vars)

# Remove NAs

oth_exp[oth_exp<0] <- NA

oth_exp$bmi66u[oth_exp$bmi66u == 9999] <- NA
oth_exp$bmi66u[oth_exp$bmi66u == 7777] <- NA

oth_exp$bmi72u[oth_exp$bmi72u == 9999] <- NA
oth_exp$bmi72u[oth_exp$bmi72u == 7777] <- NA

oth_exp$bmi82u[oth_exp$bmi82u == 9999] <- NA
oth_exp$bmi82u[oth_exp$bmi82u == 7777] <- NA

oth_exp$bmi89u[oth_exp$bmi89u == 9999] <- NA
oth_exp$bmi89u[oth_exp$bmi89u == 7777] <- NA

oth_exp$bmi99u[oth_exp$bmi99u == 9999] <- NA
oth_exp$bmi99u[oth_exp$bmi99u == 7777] <- NA

oth_exp$bmi09[oth_exp$bmi09 == 7777] <- NA
oth_exp$bmi09[oth_exp$bmi09 == 8888] <- NA
oth_exp$bmi09[oth_exp$bmi09 == 9999] <- NA

# Verify maxs and mins

mins <- sapply(oth_exp, min, na.rm=TRUE)

maxs <- sapply(oth_exp, max, na.rm=TRUE)

maxmins <- data.frame(
Var=names(oth_exp),
min=mins, 
max=maxs)

# Set classes

cat <- c("disa15x", "chron19tot15x", "summhealth15x")

oth_exp[cat] <- lapply(
oth_exp[cat], as.factor)

# Remove excluded

oth_exp_final <- oth_exp %>%
  filter(id %in% inc)

#### CREATE DATASET ####

# Merge

dataset <- Reduce(function(x,y) merge(x,y, by = "id", all=T), 
list(ages, exp_final, mort_final, oth_exp_final))

# Add follow-up length

dataset$proteomicsampleage_months <- dataset$Age60*12

dataset$death_monthssincesample <- 
  dataset$months - dataset$proteomicsampleage_months

dataset$death_yearssincesample <- 
  dataset$death_monthssincesample/12

#### OUTPUT FILES ####

# Save dataset
rm(list=setdiff(ls(), "dataset"))

out1 <- file.path(getwd(), "4. Sample description & missing data")
out2 <- file.path(getwd(), "5. Multiple Imputation")
out3 <- file.path(getwd(), "6. Aim 1 - Organ age distributions")
out4 <- file.path(getwd(), "7. Aim 2 - Mortality prediction")
out5 <- file.path(getwd(), "9. Aim 4 - Specific proteins connected to mortality & modifiable risk")

# All participants

save.image(file = file.path(out1,"dataset_all.RData"))

# Included participants

save.image(file = file.path(out1,"dataset_inc.RData"))

save.image(file = file.path(out2,"dataset_inc.RData"))

save.image(file =  file.path(out3,"dataset_inc.RData"))

save.image(file =  file.path(out4,"dataset_inc.RData"))

save.image(file =  file.path(out5,"dataset_inc.RData"))

