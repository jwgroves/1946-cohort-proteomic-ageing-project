# ============================================================
# 1_Loading_and_cleaning_proteomic_data.R
# Data cleaning for SomaScan v5.0  proteomic data
# Author: James Groves
# Date: 2025-08-30
# ============================================================

#### SET-UP ####

# Packages

library(cli)
library(crayon)
library(tidyverse)
library(lifecycle)
library(usethis)
library(readxl)
library(SomaDataIO)

# Directory

inp <- file.path(getwd(), "1. Loading & cleaning proteomic data")
out1 <- file.path(getwd(), "2. Calculating organ ages")
out2 <- file.path(getwd(), "8. Aim 3 - Modifiable risk factor associations")
out3 <- file.path(getwd(), "9. Aim 4 - Specific proteins connected to mortality & modifiable risk")

#### LOAD & CLEAN PROTEOMIC DATA ####

# Load proteomic data 

prot <- read_adat(file.path(inp, "SS-2343652_v5.0_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.20231219.adat"))

# Subset age 63 samples
prot_63 <- subset(prot, grepl("EB", SampleId))

# Convert v5.0 data to v4.1 format

getSomaScanVersion(prot_63)
prot_lift <- lift_adat(prot_63, bridge = c("11k_to_7k")) 
is_lifted(prot_lift)

# Write & read CSV for cleaning

write.csv(prot_lift, file=file.path(inp, "raw_prot.csv"), row.names=F)

prot_lift <- read.csv(file.path(inp, "raw_prot.csv"))

# Remove flagged samples

prot_valid <- prot_lift %>% filter(RowCheck=="PASS")

# Remove QC columns

prot_fin <- subset(
  prot_valid, select = -c(
    PlateId, PlateRunDate, ScannerID, 
    PlatePosition, SlideId, Subarray, 
    SampleType, PercentDilution, SampleMatrix, 
    Barcode, Barcode2d, SampleName, 
    SampleNotes, AliquotingNotes, SampleDescription, 
    AssayNotes, TimePoint, ExtIdentifier, 
    SsfExtId, SampleGroup, SiteId, 
    SubjectID, CLI, RMA, 
    HybControlNormScale, RowCheck, NormScale_20, 
    NormScale_0_005, NormScale_0_5, ANMLFractionUsed_20, 
    ANMLFractionUsed_0_005, ANMLFractionUsed_0_5))

# Adjust seqids names

colnames(prot_fin) <- gsub(
  "^seq\\.", 
  "", 
  colnames(prot_fin))

#### LOAD & CLEAN METADATA ####

# Load NSHD ids

ids <- read_xlsx(file.path(inp, "1946_NSHD_EDTA_plasma_Somalogic_key_updated_Apr25.xlsx"))

# Subset columns

ids_only <- ids %>% 
  subset(select=c(ntag1, SampleId)) %>% 
  rename(ID=ntag1)

# Adjust incorrect ids

ids_only[which(ids_only$SampleId=="[anonymsed]"), "ID"] <- "[anonymsed]"
ids_only[which(ids_only$SampleId=="[anonymsed]"), "ID"] <- "[anonymsed]"
ids_only[which(ids_only$SampleId=="[anonymsed]"), "ID"] <- "[anonymsed]"
ids_only[which(ids_only$SampleId=="[anonymsed]"), "ID"] <- "[anonymsed]"
ids_only[which(ids_only$SampleId=="[anonymsed]"), "ID"] <- "[anonymsed]"
ids_only[which(ids_only$SampleId=="[anonymsed]"), "ID"] <- "[anonymsed]"
ids_only <- rbind(ids_only, data.frame(SampleId="[anonymsed]", ID="[anonymsed]"))

# Add NSHD id to proteomics

prot_merged <- merge(ids_only, prot_fin, by="SampleId")

# Remove SampleId 

prot_merged <- prot_merged %>% 
  subset(select=-c(SampleId))

# Any missing data?

anyNA(prot_merged)

# Load  chronological age & sex

age_sex <- read.csv(
  file.path(inp, 
  "agen09.csv"))

# Change column names

age_sex <- age_sex %>% 
  rename(
  ID=ntag1,
  Sex_F=sex,
  Age=AGEN09)

# Chronological age as years

age_sex$Age <- as.numeric(age_sex$Age)
age_sex$Age <- age_sex$Age / 12

# Reset sex variable

age_sex$Sex_F[age_sex$Sex_F == 1] <- 0
age_sex$Sex_F[age_sex$Sex_F == 2] <- 1

# Only individuals with proteomics

prot_ids <- prot_merged$ID

ids_age_sex <- merge(ids_only, age_sex, by="ID")

ids_age_sex <- ids_age_sex %>% 
  filter(ID %in% prot_ids) %>%
distinct(ID, .keep_all = TRUE)

# Drop Sample ID

ids_age_sex <- subset(ids_age_sex, select = -SampleId)

# Any missing data?

anyNA(ids_age_sex)

# IDs match?

setequal(prot_merged$ID, ids_age_sex$ID)

# Remove excluded participants

mort <- read.csv(file.path(inp, "mortality_derived_130824-DSH_sharing.csv"))

mort <- mort %>% 
  rename(
  death=DTH2CEN2, 
  id=ï..id)

mort_inc <- mort %>% 
  filter(
    !is.na(death)
  )

has_m_data <- mort_inc$id

prot_merged <- prot_merged %>%
  filter(ID %in% has_m_data)

ids_age_sex <- ids_age_sex %>%
  filter(ID %in% has_m_data)

#### OUTPUT FILES ####

# CSVs for organ age estimation

write.csv(prot_merged,
          file = file.path(out1, "prot_final.csv"),
          row.names = F)

write.csv(prot_merged,
          file = file.path(out2, "prot_final.csv"),
          row.names = F)

write.csv(prot_merged,
          file = file.path(out3, "prot_final.csv"),
          row.names = F)

write.csv(ids_age_sex, 
          file = file.path(out1, "meta_final.csv"),
          row.names = F)
