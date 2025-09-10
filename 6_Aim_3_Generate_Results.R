
# ============================================================
# 8_Aim_3_Generate_Results.R
# Aim 3 analysis: Life course exposures & organ ageing
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
library(glmnet)
library(mice)
library(mediation)
library(survival)
library(car)
library(mitml)
library(miceadds)

# Directory

inp <- file.path(getwd(), "8. Aim 3 - Modifiable risk factor associations")
out <- file.path(getwd(), "8. Aim 3 - Modifiable risk factor associations")

#### LOAD & PREPARE DATA ####

load(file.path(inp, "imp_data.RData"))

organs <- c("Conventional", "Brain", 
            "Heart", "Lung", 
            "Liver", "Kidney", 
            "Immune", "Artery")

#### EXPOSURES & ORGAN AGEING ####

# Life course exposures & organ ages

for (organ in organs) {
  
agevariable <- paste0(organ, "Age60")
  
formula1 <- as.formula(paste(agevariable, " ~  lowbwt + chsc + Sex_F + Age60"))

models <- with(imputed_data, {
  f <- formula1
  environment(f) <- environment()
  lm(f)
})

  results <- pool(models)
  tidy_model <- summary(results, conf.int=TRUE)
  
  model_name <- paste0("bwt_", organ)
  ageoutput <- tidy_model[2,]
  assign(model_name, ageoutput)
  
  # Examine Associations with Childhood Social Class
  formula1 <- as.formula(paste(agevariable, " ~  high_ses + Sex_F + Age60"))
  
  models <- with(imputed_data, {
    f <- formula1
    environment(f) <- environment()
    lm(f)
  })
  
  results <- pool(models)
  tidy_model <- summary(results, conf.int=TRUE)
  
  model_name <- paste0("socio_", organ)
  ageoutput <- tidy_model[2,]
  assign(model_name, ageoutput)
  
  
  # Examine Associations with Education
  formula1 <- as.formula(paste(agevariable, "~  obtained_gcse + chsc + Sex_F + Age60"))
  
  models <- with(imputed_data, {
    f <- formula1
    environment(f) <- environment()
    lm(f)
  })
  
  results <- pool(models)
  tidy_model <- summary(results, conf.int=TRUE)
  
  model_name <- paste0("gcse_", organ)
  ageoutput <- tidy_model[2,]
  assign(model_name, ageoutput)
  
  # Examine Associations with Adolescent BMI
  formula1 <- as.formula(paste(agevariable," ~  adolescent_overweight + chsc + Sex_F + Age60"))

  models <- with(imputed_data, {
    f <- formula1
    environment(f) <- environment()
    lm(f)
  })
  
  results <- pool(models)
  tidy_model <- summary(results, conf.int=TRUE)
  
  model_name <- paste0("cbmi_", organ)
  ageoutput <- tidy_model[2,]
  assign(model_name, ageoutput)
  
  # Examine Associations with Smoking
  formula1 <- as.formula(paste(agevariable," ~ packyrs + chsc + sc1553 + Sex_F + Age60"))
 
  models <- with(imputed_data, {
    f <- formula1
    environment(f) <- environment()
    lm(f)
  })
  
  results <- pool(models)
  tidy_model <- summary(results, conf.int=TRUE)
  
  model_name <- paste0("smok_", organ)
  ageoutput <- tidy_model[2,]
  assign(model_name, ageoutput)
  
  # Examine Associations with Alcohol
  formula1 <- as.formula(paste(agevariable, "~ alcohol_intake + chsc + sc1553 + Sex_F + Age60"))
  
  models <- with(imputed_data, {
    f <- formula1
    environment(f) <- environment()
    lm(f)
  })
  
  results <- pool(models)
  tidy_model <- summary(results, conf.int=TRUE)
  
  model_name <- paste0("alc_", organ)
  ageoutput <- tidy_model[2,]
  assign(model_name, ageoutput)
  
  # Examine Associations with Activity
  formula1 <- as.formula(paste(agevariable, "~ active + chsc + sc1553 + Sex_F + Age60"))
 
  models <- with(imputed_data, {
    f <- formula1
    environment(f) <- environment()
    lm(f)
  })
  
  results <- pool(models)
  tidy_model <- summary(results, conf.int=TRUE)
  
  model_name <- paste0("act_", organ)
  ageoutput <- tidy_model[2,]
  assign(model_name, ageoutput)
  
  # Examine Associations with Adversity
  formula1 <- as.formula(paste(agevariable, "~ adall + Sex_F + Age60"))
  
  models <- with(imputed_data, {
    f <- formula1
    environment(f) <- environment()
    lm(f)
  })
  
  results <- pool(models)
  tidy_model <- summary(results, conf.int=TRUE)
  
  model_name <- paste0("adv_", organ)
  ageoutput <- tidy_model[2,]
  assign(model_name, ageoutput)
  
}

bwt_res <- bind_rows(
  "Conventional" = bwt_Conventional,
  "Brain" = bwt_Brain,
  "Heart" = bwt_Heart,
  "Lung" = bwt_Lung,
  "Liver" = bwt_Liver,
  "Kidney" = bwt_Kidney,
  "Immune" = bwt_Immune,
  "Artery" = bwt_Artery,
  .id="Organ"
)

bwt_res$padj <- p.adjust(bwt_res$p.value, method="BH")

bwt_res <- bwt_res %>% subset(select=-c(term, df))
bwt_res <- bwt_res %>% rename(lci=`2.5 %`, uci=`97.5 %`)

socio_res <- bind_rows(
  "Conventional" = socio_Conventional,
  "Brain" = socio_Brain,
  "Heart" = socio_Heart,
  "Lung" = socio_Lung,
  "Liver" = socio_Liver,
  "Kidney" = socio_Kidney,
  "Immune" = socio_Immune,
  "Artery" = socio_Artery,
  .id="Organ"
)

socio_res$padj <- p.adjust(socio_res$p.value, method="BH")
socio_res <- socio_res %>% subset(select=-c(term, df))
socio_res <- socio_res %>% rename(lci=`2.5 %`, uci=`97.5 %`)

gcse_res <- bind_rows(
  "Conventional" = gcse_Conventional,
  "Brain" = gcse_Brain,
  "Heart" = gcse_Heart,
  "Lung" = gcse_Lung,
  "Liver" = gcse_Liver,
  "Kidney" = gcse_Kidney,
  "Immune" = gcse_Immune,
  "Artery" = gcse_Artery,
  .id="Organ"
)

gcse_res$padj <- p.adjust(gcse_res$p.value, method="BH")

gcse_res <- gcse_res %>% subset(select=-c(term, df))
gcse_res <- gcse_res %>% rename(lci=`2.5 %`, uci=`97.5 %`)

cbmi_res <- bind_rows(
  "Conventional" = cbmi_Conventional,
  "Brain" = cbmi_Brain,
  "Heart" = cbmi_Heart,
  "Lung" = cbmi_Lung,
  "Liver" = cbmi_Liver,
  "Kidney" = cbmi_Kidney,
  "Immune" = cbmi_Immune,
  "Artery" = cbmi_Artery,
  .id="Organ"
)

cbmi_res$padj <- p.adjust(cbmi_res$p.value, method="BH")

cbmi_res <- cbmi_res %>% subset(select=-c(term, df))
cbmi_res <- cbmi_res %>% rename(lci=`2.5 %`, uci=`97.5 %`)

smoking_res <- bind_rows(
  "Conventional" = smok_Conventional,
  "Brain" = smok_Brain,
  "Heart" = smok_Heart,
  "Lung" = smok_Lung,
  "Liver" = smok_Liver,
  "Kidney" = smok_Kidney,
  "Immune" = smok_Immune,
  "Artery" = smok_Artery,
  .id="Organ"
)

smoking_res$padj <- p.adjust(smoking_res$p.value, method="BH")

smoking_res <- smoking_res %>% subset(select=-c(term, df))
smoking_res <- smoking_res %>% rename(lci=`2.5 %`, uci=`97.5 %`)

alc_res <- bind_rows(
  "Conventional" = alc_Conventional,
  "Brain" = alc_Brain,
  "Heart" = alc_Heart,
  "Lung" = alc_Lung,
  "Liver" = alc_Liver,
  "Kidney" = alc_Kidney,
  "Immune" = alc_Immune,
  "Artery" = alc_Artery,
  .id="Organ"
)

alc_res$padj <- p.adjust(alc_res$p.value, method="BH")

alc_res <- alc_res %>% subset(select=-c(term, df))
alc_res <- alc_res %>% rename(lci=`2.5 %`, uci=`97.5 %`)

act_res <- bind_rows(
  "Conventional" = act_Conventional,
  "Brain" = act_Brain,
  "Heart" = act_Heart,
  "Lung" = act_Lung,
  "Liver" = act_Liver,
  "Kidney" = act_Kidney,
  "Immune" = act_Immune,
  "Artery" = act_Artery,
  .id="Organ"
)

act_res$padj <- p.adjust(act_res$p.value, method="BH")
act_res <- act_res %>% subset(select=-c(term, df))
act_res <- act_res %>% rename(lci=`2.5 %`, uci=`97.5 %`)

adv_res <- bind_rows(
  "Conventional" = adv_Conventional,
  "Brain" = adv_Brain,
  "Heart" = adv_Heart,
  "Lung" = adv_Lung,
  "Liver" = adv_Liver,
  "Kidney" = adv_Kidney,
  "Immune" = adv_Immune,
  "Artery" = adv_Artery,
  .id="Organ"
)

adv_res$padj <- p.adjust(adv_res$p.value, method="BH")

adv_res <- adv_res %>% subset(select=-c(term, df))
adv_res <- adv_res %>% rename(lci=`2.5 %`, uci=`97.5 %`)

# Adolescent overweight - Conventional age with adult BMI controlled

fit <- with(imputed_data,lm(
ConventionalAge60 ~ adolescent_overweight + chsc + Sex_F + sc1553 + ifelse(bmi09>25,1,0) + Age60))

adolescent_overweight_63 <- pool(fit)

adolescent_overweight_63 <- tidy(adolescent_overweight_63, conf.int=T)

adolescent_overweight_63 <- adolescent_overweight_63 %>% 
subset(select=-c(
  statistic, df, dfcom, fmi, lambda, riv, ubar, m, b)) 

fit <- with(imputed_data,lm(
ConventionalAge60 ~ adolescent_overweight + chsc + Sex_F  
+ sc1553 + + ifelse(bmi09>25,1,0) + ifelse(bmi99u>25,1,0)))

adolescent_overweight_53 <- pool(fit)

adolescent_overweight_53 <- tidy(adolescent_overweight_53, conf.int=T)

adolescent_overweight_53 <- adolescent_overweight_53 %>% 
subset(select=-c(
  statistic, df, dfcom, fmi, lambda, riv, ubar, m, b)) 

fit <- with(imputed_data,lm(
ConventionalAge60 ~ adolescent_overweight + chsc + Sex_F + 
sc1553 + ifelse(bmi99u>25,1,0) + ifelse(bmi89u>25,1,0) + ifelse(bmi09>25,1,0) + Age60))

adolescent_overweight_43 <- pool(fit)

adolescent_overweight_43 <- tidy(adolescent_overweight_43, conf.int=T)

adolescent_overweight_43 <- adolescent_overweight_43 %>% 
subset(select=-c(
  statistic, df, dfcom, fmi, lambda, riv, ubar, m, b)) 

fit <- with(imputed_data,lm(
  ConventionalAge60 ~ adolescent_overweight + chsc + Sex_F + 
    sc1553 + ifelse(bmi99u>25,1,0) + ifelse(bmi89u>25,1,0) + 
    ifelse(bmi82u>25,1,0) + ifelse(bmi09>25,1,0) + Age60))

adolescent_overweight_36 <- pool(fit)

adolescent_overweight_36 <- tidy(adolescent_overweight_36, conf.int=T)

adolescent_overweight_36 <- adolescent_overweight_36 %>% 
subset(select=-c(
  statistic, df, dfcom, fmi, lambda, riv, ubar, m, b)) 

# Alcohol - Brain age association controlled

fit <- with(imputed_data,lm(BrainAge60 ~ avalc82u  + chsc + 
obtained_gcse + bmi09 + adall + active + 
social63 + lostcontact63+ work63 + 
finance63 + peers + divorce63 + Sex_F + Age60))

brain_alcohol_controlled_82 <- pool(fit)

brain_alcohol_controlled_82 <- tidy(brain_alcohol_controlled_82)

brain_alcohol_controlled_82 <- brain_alcohol_controlled_82[2,]

fit <- with(imputed_data,lm(BrainAge60 ~ avalc89u  + chsc + 
obtained_gcse + bmi09 + adall + active + 
social63 + lostcontact63+ work63 + finance63 + 
peers + divorce63 + Sex_F + Age60))

brain_alcohol_controlled_89 <- pool(fit)

brain_alcohol_controlled_89 <- tidy(brain_alcohol_controlled_89)

brain_alcohol_controlled_89 <- brain_alcohol_controlled_89[2,]

fit <- with(imputed_data,lm(BrainAge60 ~ avalc99u  + chsc + 
obtained_gcse + bmi09 + adall + active + 
social63 + lostcontact63+ work63 + 
finance63 + peers + divorce63 + Sex_F + Age60))

brain_alcohol_controlled_99 <- pool(fit)

brain_alcohol_controlled_99 <- tidy(brain_alcohol_controlled_99)

brain_alcohol_controlled_99 <- brain_alcohol_controlled_99[2,]

fit <- with(imputed_data,lm(BrainAge60 ~ avalc09u  + chsc + 
obtained_gcse + bmi09 + adall + active + 
social63 + lostcontact63+ work63 + 
finance63 + peers + divorce63 + Sex_F + Age60))

brain_alcohol_controlled_09 <- pool(fit)

brain_alcohol_controlled_09 <- tidy(brain_alcohol_controlled_09)

brain_alcohol_controlled_09 <- brain_alcohol_controlled_09[2,]

brain_alcohol_controlled <- bind_rows(brain_alcohol_controlled_82,
brain_alcohol_controlled_89, brain_alcohol_controlled_99,
brain_alcohol_controlled_09)

brain_alcohol_controlled <- brain_alcohol_controlled %>%
subset(select=-c(statistic, df, dfcom, fmi, lambda, riv, ubar, m, b)) 

# Test dose-dependent association between Alcohol - Brain age

fit <- with(imputed_data,lm(
BrainAge60 ~ as.factor(alcoholcategory)  + chsc +
obtained_gcse + bmi09 + adall + active + 
social63 + lostcontact63+ work63 + 
finance63 + peers + divorce63 + Sex_F + Age60))

brain_alcohol_dose <- pool(fit)

brain_alcohol_dose <- tidy(brain_alcohol_dose)

brain_alcohol_dose <- brain_alcohol_dose[2:4,]

brain_alcohol_dose <- brain_alcohol_dose %>%
subset(select=-c(statistic, df, dfcom, fmi, lambda, riv, ubar, m, b)) 

# Test dose-dependent association between Activity - Chronological age

fit <- with(imputed_data, lm(ConventionalAge60 ~ as.factor(activecategory)  +
chsc + sc1553 + Sex_F + Age60))

activity_dose <- pool(fit)

activity_dose <- tidy(activity_dose)

activity_dose <- activity_dose[2:5,]

activity_dose <- activity_dose %>% subset(
select=-c(statistic, df, dfcom, fmi, 
lambda, riv, ubar, m, b)) 

# Adversity components & accelerated ageing

fit <- with(imputed_data, lm(ConventionalAge60 ~ CROW15y + 
amenities + childSES + phy16 +
finance36 + employed36 + amenities26 + crow26 + 
work4353 + work36 + finance4353 + housecondition43 +
work63 + finance63 + matsep + divorce16 + 
peers + social36 + divorce36 + support4353 + 
lostcontact4353 + chdconflict4353 + divorce4353 + social4353 + 
support63 + lostcontact63 + chdconflict63 + 
divorce63 + social63 + Sex_F + Age60))

adversity_multivariable <- pool(fit)

adversity_multivariable <- tidy(adversity_multivariable)

adversity_multivariable <- adversity_multivariable %>% subset(select=-c(statistic, df, dfcom, fmi, lambda, riv, ubar, m, b)) 

# Stopping smoking & ageing

fit <- with(imputed_data, lm(ConventionalAge60 ~ as.factor(stopstatus) 
+ chsc + sc1553 + Sex_F + Age60))

pooled_smok <- pool(fit)

pooled_smok <- tidy(pooled_smok)

pooled_smok <- pooled_smok[2:5,]

pooled_smok <- pooled_smok %>% subset(
select=-c(statistic, df, dfcom, fmi, 
lambda, riv, ubar, m, b)) 

# Education group & ageing

fit <- with(imputed_data, lm(ConventionalAge60 ~ as.factor(educat) 
+ chsc + Sex_F + Age60))

edu_pooled <- pool(fit)

edu_pooled <- tidy(edu_pooled)

edu_pooled <- edu_pooled[2:5,]

edu_pooled <- edu_pooled %>% subset(
  select=-c(statistic, df, dfcom, fmi, 
            lambda, riv, ubar, m, b)) 

# Protective factor distributions

imp_list <- complete(imputed_data, "all")

pfs_count_list <- lapply(imp_list, function(df){
  df %>% count(protectivefactors_grouped) %>%
    rename(value = protectivefactors_grouped, count=n)
})

pfs_count_df <- bind_rows(pfs_count_list, .id="m")

pf_averagedvalues <- pfs_count_df %>%
group_by(value) %>%
summarise(meancount = mean(count))

pf_averagedvalues <- pf_averagedvalues %>%
rename(pfs_score = value, 
mean_n_across_imputations=meancount)

# Conventional ageing risk

fit <- with(imputed_data, glm(
ConventionalAge60high ~ as.factor(protectivefactors_grouped) + 
Sex_F + Age60, family="binomial"))

pfs_Conv <- pool(fit)

pfs_Conv <- tidy(pfs_Conv, conf.int=TRUE, conf.level=0.95)

pfs_Conv <-  pfs_Conv %>%
mutate(term=ifelse(
term=="as.factor(protectivefactors_grouped)4","4 protective factors",term))
pfs_Conv <-  pfs_Conv %>%
mutate(term=ifelse(
term=="as.factor(protectivefactors_grouped)5","5 protective factors",term))
pfs_Conv <-  pfs_Conv %>%
mutate(term=ifelse(
term=="as.factor(protectivefactors_grouped)6","6 protective factors",term))
pfs_Conv <-  pfs_Conv %>%
mutate(term=ifelse(
term=="as.factor(protectivefactors_grouped)7",
"7 or more protective factors",term))

pfs_Conv$hr <- exp(pfs_Conv$estimate)

pfs_Conv$lower_ci <- exp(pfs_Conv$conf.low)

pfs_Conv$upper_ci <- exp(pfs_Conv$conf.high)

pfs_Conv <- pfs_Conv %>% subset(select=-c(
b, df, dfcom, fmi, lambda, m, riv, ubar, conf.low, conf.high))

pfs_Conv <- pfs_Conv[-1,]

# Mutli-organ extreme ageing risk 

fit <- with(imputed_data, glm(
multiorganageing ~ 
as.factor(protectivefactors_grouped) + Sex_F + Age60, family="binomial"))

pfs_multi <- pool(fit)

pfs_multi <- tidy(pfs_multi, conf.int=TRUE, conf.level=0.95)

pfs_multi <-  pfs_multi %>%
mutate(term=ifelse(
term=="as.factor(protectivefactors_grouped)4",
"4 protective factors",term))
pfs_multi <-  pfs_multi %>%
mutate(term=ifelse(
term=="as.factor(protectivefactors_grouped)5",
"5 protective factors",term))
pfs_multi <-  pfs_multi %>%
mutate(term=ifelse(
term=="as.factor(protectivefactors_grouped)6",
"6 protective factors",term))
pfs_multi <-  pfs_multi %>%
mutate(term=ifelse(term=="as.factor(protectivefactors_grouped)7",
"7 or more protective factors",term))

pfs_multi$hr <- exp(pfs_multi$estimate)

pfs_multi$lower_ci <- exp(pfs_multi$conf.low)

pfs_multi$upper_ci <- exp(pfs_multi$conf.high)

pfs_multi <- pfs_multi %>% subset(
select=-c(b, df, dfcom, fmi, lambda, m, riv, ubar, conf.low, conf.high))

pfs_multi <- pfs_multi[-1,]

# Protective factor multivariable models

model_conv <- with(imputed_data, glm(
ConventionalAge60high ~ 
normal_birthweight + high_ses + 
obtained_gcse + normal_adolescent_bmi + non_smoker +
mod_alcohol + active + no_adversity + Sex_F + Age60, 
family="binomial"))

pooled_conv <- pool(model_conv)

pooled_conv <- tidy(pooled_conv)

pooled_conv <- pooled_conv %>% subset(
select=-c(statistic, df, dfcom, fmi, lambda, riv, ubar, m, b)) 

model_multiorgan <- with(imputed_data, glm(
multiorganageing ~ 
normal_birthweight + high_ses + 
obtained_gcse + normal_adolescent_bmi + non_smoker +
mod_alcohol + active + no_adversity + Sex_F + Age60, 
family="binomial"))

pooled_multiorgan <- pool(model_multiorgan)

pooled_multiorgan <- tidy(pooled_multiorgan)

pooled_multiorgan <- pooled_multiorgan %>% subset(
select=-c(statistic, df, dfcom, fmi, lambda, riv, ubar, m, b)) 

# Protective factors & mortality

pfs_mort <- with(imputed_data, glm(
death ~ protectivefactors + Sex_F + Age60, 
family="binomial"))

pfs_mort <- pool(pfs_mort)
pfs_mort_tidy <- tidy(pfs_mort)

pfs_mort_tidy$oddsratio <- exp(pfs_mort_tidy$estimate)
pfs_mort_tidy <- as.data.frame(pfs_mort_tidy)

# Mediation analysis (~ 25 minutes to run)

# vector of mediators to loop over
mediators <- c("BrainAge60","HeartAge60",
"LungAge60","LiverAge60","KidneyAge60",
"ImmuneAge60","ArteryAge60")

pooled_results <- lapply(mediators, function(med) {
  
  message("Running:", med)
  
  med_list <- lapply(imp_list, function(data) {
    
fit_m <- lm(as.formula(paste(med, "~ protectivefactors + Sex_F + Age60")), 
    data = data)
  
fit_y <- glm(as.formula(paste("death ~", med, 
"+ protectivefactors + Sex_F + Age60")),
data = data, family = "binomial")

mediate(fit_m, fit_y,
            treat    = "protectivefactors",
            mediator = med,
            sims     = 1000)
  })

  acme_est <- sapply(med_list, function(x) x$d.avg)
  ade_est  <- sapply(med_list, function(x) x$z.avg)
  te_est   <- sapply(med_list, function(x) x$tau.coef)
pm_est <- sapply(med_list, function(x) {
    r <- as.numeric(unlist(x$d.avg.sims)) / as.numeric(unlist(x$tau.sims))
    mean(r)})

  m <- length(med_list)
  
  ests <- lapply(seq_len(m), function(i)
    c(ACME=acme_est[i], ADE=ade_est[i], TE=te_est[i], PM=pm_est[i]))
  
  vars <- lapply(seq_len(m), function(i){
    x <- med_list[[i]]
      ACME=as.numeric(unlist(x$d.avg.sims))
      ADE=as.numeric(unlist(x$z.avg.sims)) 
      TE=as.numeric(unlist(x$tau.sims))
    PM=ACME/TE
    cov(cbind(ACME=ACME, ADE=ADE, TE=TE, PM=PM), use="complete.obs")
  })
    
  out <- testEstimates(qhat=ests, uhat=vars, conf.level=0.95)
  
  all <- out$estimates %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    mutate(term = sub("\\..*", "", term),
    lci = Estimate - qt(.975, df) * `Std.Error`,
    uci = Estimate + qt(.975, df) * `Std.Error`,
    p   = `P(>|t|)`) %>%
    dplyr::select(term, estimate = Estimate, se = `Std.Error`, p, lci, uci)
  
   wide <- all %>%
     pivot_wider(names_from = term, 
     values_from = c(estimate, se, p, lci, uci),
     names_glue = "{term}_{.value}") %>%
    mutate(organ=med, .before=1)
   
   wide

})

pooled_results_df <- do.call(rbind, pooled_results)


#### PROTEIN LEVEL ANALYSIS ####

# Add protein level data

prot <- read.csv(file.path(inp,"proteomics_final.csv"))

prot <- prot %>% rename(id=ID)

prot_log = prot %>% mutate(across(-1, ~log10(.)))

prot_norm = prot_log %>% mutate(across(-1, ~scale(.)[,1]))

imp_list <- complete(imputed_data, action="all")

imp_prot <- lapply(imp_list, function(df) {
  merge(df, prot_norm, by="id")
})

imp_p <- as.mitml.list(imp_prot)

# Add protein metadata

prot_metadata <- read_xlsx(file.path(inp, "3.6 SL00000906_SomaScan_11K_v5.0_Plasma_Serum_Annotated_Menu.xlsx"))

colnames(prot_metadata) <- as.character(prot_metadata[4,])

prot_metadata <- prot_metadata[-(1:4),]

prot_metadata$SeqId <- paste0("X", gsub("-",".",prot_metadata$SeqId))

prot_metadata <- prot_metadata %>% 
  subset(select=c("SeqId", "SomaId", "Target Full Name", 
                  "UniProt ID", "Entrez Gene Name"
                  ))

coef <- read_xlsx(file.path(inp, "organ_clock_coef.xlsx"))

soma_seqid <- read_xlsx(file.path(inp, "Somamer_seqids_match.xlsx"))

clock_coef <- merge(soma_seqid, coef, by="Somamer")

clock_coef <- clock_coef %>% subset(select=-c(Somamer))

clock_coef$SeqId <- paste0("X", gsub("-",".",clock_coef$SeqId))

prot_meta_full <- merge(prot_metadata, clock_coef, by="SeqId", all=TRUE)

brain <- prot_meta_full %>% filter(organ=="Brain")
brain_proteins <- brain$SeqId

liver <- prot_meta_full %>% filter(organ=="Liver")
liver_proteins <- liver$SeqId

immune <- prot_meta_full %>% filter(organ=="Immune")
immune_proteins <- immune$SeqId

kidney <- prot_meta_full %>% filter(organ=="Kidney")
kidney_proteins <- kidney$SeqId

# Brain-specific proteins and alcohol

results <- list()

for (protein in brain_proteins) {

formula <- as.formula(
paste0(
protein, " ~ avalc82u + Sex_F + 
chsc + obtained_gcse + bmi09 + 
adall + active + social63 + 
lostcontact63+ work63 + 
finance63 + peers + divorce63 + Age60"))

fit <- lapply(imp_p, function(d) lm(formula, data=d))

pooled <- pool_mi(qhat=lapply(fit, coef), u=lapply(fit,vcov))

results[[protein]] <- summary(pooled)

}

pooled_df <- bind_rows(
lapply(results, function(df) rownames_to_column(df, "term")), .id="SeqId")

pooled_df <- pooled_df %>% filter(term=="avalc82u")

pooled_df_82 <- merge(pooled_df, prot_metadata, by="SeqId")

pooled_df_82$padj <- p.adjust(pooled_df_82$p, method="BH")

pooled_df_82$sig <- ifelse(pooled_df_82$padj<0.05,1,0)

pooled_df_82 <- 
  pooled_df_82 %>%
  rename(
    estimate=results,
    ci_lower=`(lower`,
    ci_upper=`upper)`
  )

pooled_df_82 <- pooled_df_82 %>% subset(select=c(SeqId, estimate,
se, p, ci_lower, ci_upper, SomaId, `Target Full Name`,
`UniProt ID`, `Entrez Gene Name`))

pooled_df_82 <- pooled_df_82[, c("Entrez Gene Name", "Target Full Name", 
"SeqId", "SomaId", 
"UniProt ID", "estimate", "se", 
"ci_lower", "ci_upper", "p")]

results <- list()

for (protein in brain_proteins) {
  
  formula <- as.formula(
  paste0(protein, " ~ avalc89u + Sex_F + chsc + 
  obtained_gcse + bmi09 + adall + active + social63 + 
  lostcontact63+ work63 + finance63 + peers + divorce63 + Age60"))
  
  fit <- lapply(imp_p, function(d) lm(formula, data=d))
  
  pooled <- pool_mi(qhat=lapply(fit, coef), u=lapply(fit,vcov))
  
  results[[protein]] <- summary(pooled)
}
  
pooled_df <- bind_rows(
lapply(results, function(df) rownames_to_column(df, "term")), .id="SeqId")

pooled_df <- pooled_df %>% filter(term=="avalc89u")

pooled_df_89 <- merge(pooled_df, prot_metadata, by="SeqId")

pooled_df_89$padj <- p.adjust(pooled_df_89$p, method="BH")

pooled_df_89$sig <- ifelse(pooled_df_89$padj<0.05,1,0)

pooled_df_89 <- 
  pooled_df_89 %>%
  rename(
    estimate=results,
    ci_lower=`(lower`,
    ci_upper=`upper)`
  )

pooled_df_89 <- pooled_df_89 %>% 
subset(select=c(SeqId, estimate, se, p, ci_lower, ci_upper))

pooled_df_89 <- pooled_df_89[, c("estimate", "se", "ci_lower", "ci_upper", "p")]

# Age 53

results <- list()

for (protein in brain_proteins) {

  formula <- as.formula(
  paste0(protein, " ~ avalc99u + Sex_F + chsc +
  obtained_gcse + bmi09 + adall + active + 
  social63 + lostcontact63+ work63 + finance63 + peers + 
  divorce63 + Age60"))
  
  fit <- lapply(imp_p, function(d) lm(formula, data=d))
  
  pooled <- pool_mi(qhat=lapply(fit, coef), u=lapply(fit,vcov))
  
  results[[protein]] <- summary(pooled)
}
  
pooled_df <- bind_rows(
lapply(results, function(df) rownames_to_column(df, "term")), .id="SeqId")

pooled_df <- pooled_df %>% filter(term=="avalc99u")

pooled_df_99 <- merge(pooled_df, prot_metadata, by="SeqId")

pooled_df_99$padj <- p.adjust(pooled_df_99$p, method="BH")

pooled_df_99$sig <- ifelse(pooled_df_99$padj<0.05,1,0)

pooled_df_99 <- 
  pooled_df_99 %>%
  rename(
    estimate=results,
    ci_lower=`(lower`,
    ci_upper=`upper)`
  )

pooled_df_99 <- pooled_df_99 %>% 
subset(select=c(SeqId, estimate, se, p, ci_lower, ci_upper))

pooled_df_99 <- pooled_df_99[, c("estimate", "se", "ci_lower", "ci_upper", "p")]

# Age 63

results <- list()

for (protein in brain_proteins) {
  
  formula <- as.formula(
    paste0(protein, " ~ avalc09u + Sex_F + chsc + 
    obtained_gcse + bmi09 + adall + active + social63 + 
    lostcontact63+ work63 + finance63 + peers +
    divorce63 + Age60"))
  
  fit <- lapply(imp_p, function(d) lm(formula, data=d))
  
  pooled <- pool_mi(qhat=lapply(fit, coef), u=lapply(fit,vcov))
  
  results[[protein]] <- summary(pooled)
}
  
pooled_df <- bind_rows(
lapply(results, function(df) rownames_to_column(df, "term")), .id="SeqId")

pooled_df <- pooled_df %>% filter(term=="avalc09u")

pooled_df_09 <- merge(pooled_df, prot_metadata, by="SeqId")

pooled_df_09$padj <- p.adjust(pooled_df_09$p, method="BH")

pooled_df_09$sig <- ifelse(pooled_df_09$padj<0.05,1,0)

pooled_df_09 <- 
  pooled_df_09 %>%
  rename(
    estimate=results,
    ci_lower=`(lower`,
    ci_upper=`upper)`
  )

pooled_df_09 <- pooled_df_09 %>% subset(select=c(SeqId, estimate, se, p, ci_lower, ci_upper))

pooled_df_09 <- pooled_df_09[, c("estimate", "se", "ci_lower", "ci_upper", "p")]

# Overall

results <- list()

for (protein in brain_proteins) {
 
   formula <- as.formula(
    paste0(protein, " ~ alcohol_intake + Sex_F + chsc + 
    obtained_gcse + bmi09 + adall + active + 
    social63 + lostcontact63+ work63 + finance63 + peers + divorce63"))
  
   fit <- lapply(imp_p, function(d) lm(formula, data=d))
   
   pooled <- pool_mi(qhat=lapply(fit, coef), u=lapply(fit,vcov))
   
   results[[protein]] <- summary(pooled)
}
  
pooled_df <- bind_rows(
lapply(results, function(df) rownames_to_column(df, "term")), .id="SeqId")

pooled_df_ov <- merge(pooled_df, prot_metadata, by="SeqId")

pooled_df_ov$padj <- p.adjust(pooled_df_ov$p, method="BH")

pooled_df_ov$sig <- ifelse(pooled_df_ov$padj<0.05,1,0)

pooled_df_ov <- 
  pooled_df_ov %>%
  rename(
    estimate=results,
    ci_lower=`(lower`,
    ci_upper=`upper)`
  )

pooled_df_ov <- pooled_df_ov %>% subset(select=c(SeqId, estimate, se, p, ci_lower, ci_upper))

pooled_df_ov <- pooled_df_ov[, c("estimate", "se", "ci_lower", "ci_upper", "p")]

# Liver & Adolescent Overweight 

results <- list()

for (protein in liver_proteins) {
  
  formula <- as.formula(paste0(protein, " ~ adolescent_overweight + Sex_F + chsc"))
  
  fit <- lapply(imp_p, function(d) lm(formula, data=d))
  
  pooled <- pool_mi(qhat=lapply(fit, coef), u=lapply(fit,vcov))
  
  results[[protein]] <- summary(pooled)
}
  
pooled_liver_ow <- bind_rows(
lapply(results, function(df) rownames_to_column(df, "term")), .id="SeqId")

pooled_liver_ow <- pooled_liver_ow %>% filter(term=="adolescent_overweight")

pooled_liver_ow <- merge(pooled_liver_ow, prot_metadata, by="SeqId")

pooled_liver_ow$padj <- p.adjust(pooled_liver_ow$p, method="BH")

pooled_liver_ow$sig <- ifelse(pooled_liver_ow$padj<0.05,1,0)

pooled_liver_ow <- 
  pooled_liver_ow %>%
  rename(
    estimate=results,
    ci_lower=`(lower`,
    ci_upper=`upper)`
  )

pooled_liver_ow <- pooled_liver_ow %>% subset(
select=c(SeqId, estimate, se, p, ci_lower, 
ci_upper, SomaId, `Target Full Name`, 
`UniProt ID`, `Entrez Gene Name`))

pooled_liver_ow <- pooled_liver_ow[, c(
"Entrez Gene Name", "Target Full Name", "SeqId", 
"SomaId", "UniProt ID", "estimate",
"se", "ci_lower", "ci_upper", "p")]

pooled_liver_ow <- pooled_liver_ow %>% arrange(p)

# Kidney & Smoking

results <- list()

for (protein in kidney_proteins) {
  
  formula <- as.formula(
  paste0(protein, " ~ packyrs + Sex_F + chsc + sc1553"))
  
  fit <- lapply(imp_p, function(d) lm(formula, data=d))
  
  pooled <- pool_mi(qhat=lapply(fit, coef), u=lapply(fit,vcov))
  
  results[[protein]] <- summary(pooled)
}
  
pooled_kidney_sm <- bind_rows(
lapply(results, function(df) rownames_to_column(df, "term")), .id="SeqId")

pooled_kidney_sm <- pooled_kidney_sm %>% filter(term=="packyrs")

pooled_kidney_sm <- merge(pooled_kidney_sm, prot_metadata, by="SeqId")

pooled_kidney_sm$padj <- p.adjust(pooled_kidney_sm$p, method="BH")

pooled_kidney_sm$sig <- ifelse(pooled_kidney_sm$padj<0.05,1,0)

pooled_kidney_sm <- 
  pooled_kidney_sm %>%
  rename(
    estimate=results,
    ci_lower=`(lower`,
    ci_upper=`upper)`
  )

pooled_kidney_sm <- pooled_kidney_sm %>% subset(select=c(SeqId, estimate, 
se, p, ci_lower, ci_upper, SomaId, `Target Full Name`, `UniProt ID`, `Entrez Gene Name`))

pooled_kidney_sm <- pooled_kidney_sm[, c("Entrez Gene Name", "Target Full Name",
"SeqId", "SomaId", "UniProt ID", "estimate", "se", "ci_lower", "ci_upper", "p")]

pooled_kidney_sm <- pooled_kidney_sm %>% arrange(p)

# Immune & Activity

results <- list()

for (protein in immune_proteins) {
  
  formula <- as.formula(
  paste0(protein, " ~ active + Sex_F + chsc + sc1553"))
  
  fit <- lapply(imp_p, function(d) lm(formula, data=d))
  
  pooled <- pool_mi(qhat=lapply(fit, coef), u=lapply(fit,vcov))
  
  results[[protein]] <- summary(pooled)
}
  
pooled_immune_act <- bind_rows(
  lapply(results, function(df) rownames_to_column(df, "term")), .id="SeqId")

pooled_immune_act <- pooled_immune_act %>% filter(term=="active")

pooled_immune_act <- merge(pooled_immune_act, prot_metadata, by="SeqId")

pooled_immune_act$padj <- p.adjust(pooled_immune_act$p, method="BH")

pooled_immune_act$sig <- ifelse(pooled_immune_act$padj<0.05,1,0)

pooled_immune_act <- 
  pooled_immune_act %>%
  rename(
    estimate=results,
    ci_lower=`(lower`,
    ci_upper=`upper)`
  )

pooled_immune_act <- pooled_immune_act %>% subset(select=c(SeqId, estimate,
se, p, ci_lower, ci_upper, SomaId, `Target Full Name`, 
`UniProt ID`, `Entrez Gene Name`))

pooled_immune_act <- pooled_immune_act[, c("Entrez Gene Name", "Target Full Name",
"SeqId", "SomaId",
"UniProt ID", "estimate", "se", "ci_lower", "ci_upper", "p")]

pooled_immune_act <- pooled_immune_act %>% arrange(p)

# FIBA analysis for Brain - Alcohol (~ 20 minutes to run)

model <- with(
imputed_data, lm(BrainAge60 ~ alcohol_intake + Sex_F + chsc + sc1553 + Age60))

model <- pool(model)

tidy <- tidy(model)
tidy <- tidy[2,]

prev_coef <- tidy$estimate

brain_proteins <- sub("X","",brain_proteins)
brain_proteins <- gsub("\\.","-",brain_proteins)

results <- list()

for (brain_protein in brain_proteins){
  
  filepath <- "S:/LHA_JG0923/Proteomic Ageing Project/Final Scripts & Data/FIBA_run_1/"
  file <- paste0(filepath, brain_protein, "_run1", '.csv')
  tp1 <- read.csv(file)

  names(tp1)[names(tp1) == "ID"] <- "id"
  names(tp1)[names(tp1) == "AgeGap"] <- "FIBA_BrainAge60"
  
  tp1 <- tp1 %>% subset(select=c(id, FIBA_BrainAge60))

  fiba_data <- lapply(imp_list, function(df) {
  fiba_datasets <- Reduce(function(x,y) merge(x,y, by = "id", all = TRUE), list(tp1, df))
     }) 
  
    formula <- as.formula(paste0("FIBA_BrainAge60 ~ alcohol_intake + Sex_F + chsc + sc1553 + Age60"))
    
    fit <- lapply(fiba_data, function(d) lm(formula, data=d))
    
    pooled <- pool_mi(qhat=lapply(fit, coef), u=lapply(fit,vcov))
    
    res <- summary(pooled)
    
    res <-  res[2,]
  
    fiba_estimate <- res$results

  df <- data.frame(SeqId = brain_protein, fiba_estimate = fiba_estimate)
  
  assign(paste0(brain_protein), df)
  
}

run1 <- bind_rows(lapply(brain_proteins, function(x) get(x)))

run1$FIBA_diff_1 <- run1$fiba_estimate - prev_coef

for (brain_protein in brain_proteins){

  filepath <- "S:/LHA_JG0923/Proteomic Ageing Project/Final Scripts & Data/FIBA_run_2/"
  file <- paste0(filepath, brain_protein, "_run2", '.csv')
  tp1 <- read.csv(file)
  
  names(tp1)[names(tp1) == "ID"] <- "id"
  names(tp1)[names(tp1) == "AgeGap"] <- "FIBA_BrainAge60"
  
  tp1 <- tp1 %>% subset(select=c(id, FIBA_BrainAge60))
  
  fiba_data <- lapply(imp_list, function(df) {
    fiba_datasets <- Reduce(function(x,y) merge(x,y, by = "id", all = TRUE), list(tp1, df))
  }) 
  
  formula <- as.formula(paste0("FIBA_BrainAge60 ~ alcohol_intake + Sex_F + chsc + sc1553 + Age60"))
  
  fit <- lapply(fiba_data, function(d) lm(formula, data=d))
  
  pooled <- pool_mi(qhat=lapply(fit, coef), u=lapply(fit,vcov))
  
  res <- summary(pooled)
  
  res <-  res[2,]
  
  fiba_estimate <- res$results
  
  df <- data.frame(SeqId = brain_protein, fiba_estimate = fiba_estimate)
  
  assign(paste0(brain_protein), df)
  
}

run2 <- bind_rows(lapply(brain_proteins, function(x) get(x)))

run2$FIBA_diff_2 <- run2$fiba_estimate - prev_coef

for (brain_protein in brain_proteins){

  filepath <- "S:/LHA_JG0923/Proteomic Ageing Project/Final Scripts & Data/FIBA_run_3/"
  file <- paste0(filepath, brain_protein, "_run3", '.csv')
  tp1 <- read.csv(file)
  
  names(tp1)[names(tp1) == "ID"] <- "id"
  names(tp1)[names(tp1) == "AgeGap"] <- "FIBA_BrainAge60"
  
  tp1 <- tp1 %>% subset(select=c(id, FIBA_BrainAge60))
  
  fiba_data <- lapply(imp_list, function(df) {
    fiba_datasets <- Reduce(function(x,y) merge(x,y, by = "id", all = TRUE), list(tp1, df))
  }) 
  
  formula <- as.formula(paste0("FIBA_BrainAge60 ~ alcohol_intake + Sex_F + chsc + sc1553 + Age60"))
  
  fit <- lapply(fiba_data, function(d) lm(formula, data=d))
  
  pooled <- pool_mi(qhat=lapply(fit, coef), u=lapply(fit,vcov))
  
  res <- summary(pooled)
  
  res <-  res[2,]
  
  fiba_estimate <- res$results
  
  df <- data.frame(SeqId = brain_protein, fiba_estimate = fiba_estimate)
  
  assign(paste0(brain_protein), df)
}

run3 <- bind_rows(lapply(brain_proteins, function(x) get(x)))

run3$FIBA_diff_3 <- run3$fiba_estimate - prev_coef

for (brain_protein in brain_proteins){

  filepath <- "S:/LHA_JG0923/Proteomic Ageing Project/Final Scripts & Data/FIBA_run_4/"
  file <- paste0(filepath, brain_protein, "_run4", '.csv')
  tp1 <- read.csv(file)
  
  names(tp1)[names(tp1) == "ID"] <- "id"
  names(tp1)[names(tp1) == "AgeGap"] <- "FIBA_BrainAge60"
  
  tp1 <- tp1 %>% subset(select=c(id, FIBA_BrainAge60))
  
  fiba_data <- lapply(imp_list, function(df) {
    fiba_datasets <- Reduce(function(x,y) merge(x,y, by = "id", all = TRUE), list(tp1, df))
  }) 
  
  formula <- as.formula(paste0("FIBA_BrainAge60 ~ alcohol_intake + Sex_F + chsc + sc1553 + Age60"))
  
  fit <- lapply(fiba_data, function(d) lm(formula, data=d))
  
  pooled <- pool_mi(qhat=lapply(fit, coef), u=lapply(fit,vcov))
  
  res <- summary(pooled)
  
  res <-  res[2,]
  
  fiba_estimate <- res$results
  
  df <- data.frame(SeqId = brain_protein, fiba_estimate = fiba_estimate)
  
  assign(paste0(brain_protein), df)
}

run4 <- bind_rows(lapply(brain_proteins, function(x) get(x)))
run4$FIBA_diff_4 <- run4$fiba_estimate - prev_coef

for (brain_protein in brain_proteins){

  filepath <- "S:/LHA_JG0923/Proteomic Ageing Project/Final Scripts & Data/FIBA_run_5/"
  file <- paste0(filepath, brain_protein, "_run5", '.csv')
  tp1 <- read.csv(file)
  
  names(tp1)[names(tp1) == "ID"] <- "id"
  names(tp1)[names(tp1) == "AgeGap"] <- "FIBA_BrainAge60"
  
  tp1 <- tp1 %>% subset(select=c(id, FIBA_BrainAge60))
  
  fiba_data <- lapply(imp_list, function(df) {
    fiba_datasets <- Reduce(function(x,y) merge(x,y, by = "id", all = TRUE), list(tp1, df))
  }) 
  
  formula <- as.formula(paste0("FIBA_BrainAge60 ~ alcohol_intake + Sex_F + chsc + sc1553 + Age60"))
  
  fit <- lapply(fiba_data, function(d) lm(formula, data=d))
  
  pooled <- pool_mi(qhat=lapply(fit, coef), u=lapply(fit,vcov))
  
  res <- summary(pooled)
  
  res <-  res[2,]
  
  fiba_estimate <- res$results
  
  df <- data.frame(SeqId = brain_protein, fiba_estimate = fiba_estimate)
  
  assign(paste0(brain_protein), df)
}

run5 <- bind_rows(lapply(brain_proteins, function(x) get(x)))

run5$FIBA_diff_5 <- run5$fiba_estimate - prev_coef

run1 <- run1 %>% subset(select=-c(fiba_estimate))
run2 <- run2 %>% subset(select=-c(fiba_estimate))
run3 <- run3 %>% subset(select=-c(fiba_estimate))
run4 <- run4 %>% subset(select=-c(fiba_estimate))
run5 <- run5 %>% subset(select=-c(fiba_estimate))

fiba_results <- Reduce(function(x,y) merge(x,y, by = "SeqId", all = TRUE), 
list(run1, run2, run3, run4, run5))

fiba_results <- fiba_results %>%
mutate(av_FIBA_diff=rowMeans(
  across(c(FIBA_diff_1, FIBA_diff_2, 
  FIBA_diff_3, FIBA_diff_4, FIBA_diff_5))))

fiba_results$SeqId <- paste0("X", fiba_results$SeqId)
fiba_results$SeqId <- gsub("\\-",".",fiba_results$SeqId)

fiba_results_prot <- merge(fiba_results, prot_meta_full, by="SeqId")

fiba_results_prot <- fiba_results_prot %>% subset(select=c(SeqId, FIBA_diff_1,
FIBA_diff_2, FIBA_diff_3, FIBA_diff_4, FIBA_diff_5, SomaId, `Target Full Name`,
`UniProt ID`, `Entrez Gene Name`, abs_mean_coef, av_FIBA_diff))

fiba_results_prot <- fiba_results_prot %>% rename(
  run1_FIBA_coef = FIBA_diff_1,
  run2_FIBA_coef = FIBA_diff_2,
  run3_FIBA_coef = FIBA_diff_3,
  run4_FIBA_coef = FIBA_diff_4,
  run5_FIBA_coef = FIBA_diff_5,
  avg_FIBA_coef = av_FIBA_diff,
  clock_coefficient = abs_mean_coef
)

fiba_results_prot <- fiba_results_prot[, c(
"Entrez Gene Name", "Target Full Name", "SeqId", 
"SomaId", "UniProt ID", "clock_coefficient", 
"run1_FIBA_coef", "run2_FIBA_coef", "run3_FIBA_coef", 
"run4_FIBA_coef", "run5_FIBA_coef", "avg_FIBA_coef")]

fiba_results_prot <- fiba_results_prot %>% arrange(desc(clock_coefficient))

#### OUTPUT FILES ####

# Aim 3 Supplementary Tables

res <- createWorkbook()

# ST8: Organ ageing & life-course exposures
addWorksheet(res, "Organs & life-course exposures")
writeData(res, "Organs & life-course exposures", bwt_res, startRow=2)
writeData(res, "Organs & life-course exposures", socio_res, startRow=13)
writeData(res, "Organs & life-course exposures", gcse_res, startRow=24)
writeData(res, "Organs & life-course exposures", cbmi_res, startRow=35)
writeData(res, "Organs & life-course exposures", smoking_res, startRow=47)
writeData(res, "Organs & life-course exposures", alc_res, startRow=59)
writeData(res, "Organs & life-course exposures", act_res, startRow=71)
writeData(res, "Organs & life-course exposures", adv_res, startRow=83)

# ST9: Adolescent overweight & liver proteins
addWorksheet(res, "Liver proteins - ad overweight")
writeData(res, "Liver proteins - ad overweight", pooled_liver_ow)

# ST10: Smoking & kidney proteins
addWorksheet(res, "Kidney proteins - smoking")
writeData(res, "Kidney proteins - smoking", pooled_kidney_sm)

# ST11: Activity & immune proteins
addWorksheet(res, "Immune proteins - activity")
writeData(res, "Immune proteins - activity", pooled_immune_act)

# ST12: Alcohol & brain proteins
addWorksheet(res, "Brain proteins - alcohol")
writeData(res, "Brain proteins - alcohol", pooled_df_82, startRow=2, startCol=1)
writeData(res, "Brain proteins - alcohol", pooled_df_89, startRow=2, startCol=12)
writeData(res, "Brain proteins - alcohol", pooled_df_99, startRow=2, startCol=17)
writeData(res, "Brain proteins - alcohol", pooled_df_09, startRow=2, startCol=22)
writeData(res, "Brain proteins - alcohol", pooled_df_ov, startRow=2, startCol=27)

# ST13: FIBA: Alcohol & brain proteins
addWorksheet(res, "FIBA Brain proteins - alcohol")
writeData(res, "FIBA Brain proteins - alcohol", fiba_results_prot)

# ST14: PFS distributions
addWorksheet(res, "PFS distributions")
writeData(res, "PFS distributions", pf_averagedvalues)

# ST15: Extreme ageing & PFS
addWorksheet(res, "PFS risk for extreme ageing")
writeData(res, "PFS risk for extreme ageing", pfs_Conv)
writeData(res, "PFS risk for extreme ageing", pfs_multi, startRow=10, startCol=1)

# ST16: Mediation analyses
addWorksheet(res, "Organ age mediation")
writeData(res, "Organ age mediation", pooled_results_df)

# Save
saveWorkbook(
res,
file.path(out,"Aim_3_ST.xlsx"),
overwrite=TRUE)

# Source code files

sc <- createWorkbook()

# Fig 3a, 3b
addWorksheet(sc, "Fig 3a, 3b")
writeData(sc, "Fig 3a, 3b", bwt_res, startRow=2)
writeData(sc, "Fig 3a, 3b", socio_res, startRow=13)
writeData(sc, "Fig 3a, 3b", gcse_res, startRow=24)
writeData(sc, "Fig 3a, 3b", cbmi_res, startRow=35)
writeData(sc, "Fig 3a, 3b", smoking_res, startRow=47)
writeData(sc, "Fig 3a, 3b", alc_res, startRow=59)
writeData(sc, "Fig 3a, 3b", act_res, startRow=71)
writeData(sc, "Fig 3a, 3b", adv_res, startRow=83)

# Fig 4a
addWorksheet(sc, "Fig 4a")
writeData(sc, "Fig 4a", adolescent_overweight_63, startRow=2)
writeData(sc, "Fig 4a", adolescent_overweight_53, startRow=7)
writeData(sc, "Fig 4a", adolescent_overweight_43, startRow=12)
writeData(sc, "Fig 4a", adolescent_overweight_36, startRow=17)

# Fig 4b
addWorksheet(sc, "Fig 4b")
writeData(sc, "Fig 4b", activity_dose)

# Fig 4c
addWorksheet(sc, "Fig 4c")
writeData(sc, "Fig 4c", adversity_multivariable)

# Fig 4d
addWorksheet(sc, "Fig 4d")
writeData(sc, "Fig 4d", pooled_smok)

# Fig 4e
addWorksheet(sc, "Fig 4e")
writeData(sc, "Fig 4e", edu_pooled)

# Fig 5b
addWorksheet(sc, "Fig 5b")
writeData(sc, "Fig 5b", pf_averagedvalues)

# Fig 5c, 5d
addWorksheet(sc, "Fig 5c, 5d")
writeData(sc, "Fig 5c, 5d", pfs_Conv)
writeData(sc, "Fig 5c, 5d", pfs_multi, startRow=10, startCol=1)

# Fig 5e, 5f
addWorksheet(sc, "Fig 5e, 5f")
writeData(sc, "Fig 5e, 5f", pooled_conv, startRow=2)
writeData(sc, "Fig 5e, 5f", pooled_multiorgan, startRow=7)

# Ext data Fig 3a
addWorksheet(sc, "Ext data Fig 3a")
writeData(sc, "Ext data Fig 3a", pooled_liver_ow)

# Ext data Fig 3b
addWorksheet(sc, "Ext data Fig 3b")
writeData(sc, "Ext data Fig 3b", pooled_kidney_sm)

# Ext data Fig 3c
addWorksheet(sc, "Ext data Fig 3c")
writeData(sc, "Ext data Fig 3c", pooled_immune_act)

# Ext data Fig 4a
addWorksheet(sc, "Ext data Fig 4a")
writeData(sc, "Ext data Fig 4a", brain_alcohol_controlled)

# Ext data Fig 4b
addWorksheet(sc, "Ext data Fig 4b")
writeData(sc, "Ext data Fig 4b", brain_alcohol_dose)

# Ext data Fig 4c
addWorksheet(sc, "Ext data Fig 4c")
writeData(sc, "Ext data Fig 4c", pooled_df_82, startRow=2, startCol=1)
writeData(sc, "Ext data Fig 4c", pooled_df_89, startRow=2, startCol=12)
writeData(sc, "Ext data Fig 4c", pooled_df_99, startRow=2, startCol=17)
writeData(sc, "Ext data Fig 4c", pooled_df_09, startRow=2, startCol=22)
writeData(sc, "Ext data Fig 4c", pooled_df_ov, startRow=2, startCol=27)

# Ext data Fig 4d
addWorksheet(sc, "Ext data Fig 4d")
writeData(sc, "Ext data Fig 4d", fiba_results_prot)

# Save
saveWorkbook(
sc, 
file.path(out,"SC_Aim_3.xlsx"), 
overwrite=TRUE)
