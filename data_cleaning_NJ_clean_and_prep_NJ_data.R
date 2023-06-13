require(haven); require(tidyverse)

# Read in and save raw data -----------------------------------------------

#' An unabridged NJ dataset; contains >1100 variables, most of which we do not need for now.
#nj <- read_dta('//smb-main.ad.bcm.edu/genepi3/TiffanyChambers/GOBACK Data/New Jersey Data/Complete NJ GOBACK Data.dta')

nj.clean <- read_dta('//smb-main.ad.bcm.edu/genepi3/TiffanyChambers/GOBACK Data/New Jersey Data/Complete NJ GOBACK Data_harmonized birth certificate data.dta')

saveRDS(nj.clean, '//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/raw_data/new_jersey_raw_data_20230104.rds')

# Divide in BD, cancer, vital records files -------------------------------

setwd('//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/')

nj <- readRDS('raw_registry_data/new_jersey_raw_data_20230104.rds')

names(nj) <- names(nj) %>% tolower()

nj <- nj %>% 
  select(-c(382,383,385,387)) %>% # These columns have duplicate names.
  filter(birth_yr > 1989) %>%  # For 1985-1989, the prevalence of major birth defects is very low. These data will be excluded due to concern about non-differential misclassification leading to bias towards the null.
  mutate(studyid = paste0('nj',studyid))

#' Separate NJ data into birth cert, cancer, and BD datasets.
nj.bd <- nj %>% 
  select(studyid, dob, starts_with('icd')) %>% 
  mutate(across(starts_with('icd9'), ~ as.character(.x)),
         across(starts_with('icd'), ~ ifelse(is.na(.x), "", .x))) %>% # Convert ICD9 codes to character class and replace NAs with an empty string.
  filter(if_any(icd9_1:icd10_12, ~ .x != "")) # Keep rows where at least one BD code column is not missing.

nj.can <- nj %>% 
  select(studyid, starts_with(c('histology', 'primary_site', 'behavior', 'laterality', 'seer_summ_stage_2000'))) %>% 
  filter(if_any(starts_with('hist'), ~ !is.na(.x)))

#' Contains the names of the vital records columns I need.
sc <- readRDS("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded_datasets/sc_GOBACK_data_v20220816.rds") %>% 
  select(studyid:person.yrs) %>% 
  names()

nj.vit <- nj %>% 
  rename_with(~str_replace_all(.x, '_', '.')) %>% 
  rename(m.edu = m.edu2,
         person.yrs = timeatrisk.nocancer) %>% 
  mutate(cancer = ifelse(is.na(cancer), 0, cancer),
         person.yrs = ifelse(cancer == 1, age.at.diagnosis.n2301, 
                      ifelse((person.yrs/12) > 19, 19, person.yrs/12))) %>%
  select(all_of( intersect(names(.), sc) ), plu.totalbirths)

saveRDS(nj.bd, '//smb-main.ad.bcm.edu/genepi3/JeremySchraw/GOBACK/Datasets/nj.birth.defect.raw.data.rds')
saveRDS(nj.can, '//smb-main.ad.bcm.edu/genepi3/JeremySchraw/GOBACK/Datasets/nj.cancer.raw.data.rds')
saveRDS(nj.vit, '//smb-main.ad.bcm.edu/genepi3/JeremySchraw/GOBACK/Datasets/nj.vital.records.raw.data.rds')

rm(list = ls()); gc()

# Prepare NJ vit records data ---------------------------------------------

nj <- readRDS('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/GOBACK/Datasets/nj.vital.records.raw.data.rds') %>% 
  as_tibble() %>% 
  mutate(sex = factor(ifelse(is.na(sex), 9, sex), labels = c('Male', 'Female', 'Unknown')),
         m.race = factor(ifelse(is.na(m.race), 7, m.race), labels = c('Hispanic', 'NHW', 'NHB', 'Asian', 'Other', 'AIAN', 'unknown')),
         f.race = factor(ifelse(is.na(f.race), 7, f.race), labels = c('Hispanic', 'NHW', 'NHB', 'Asian', 'Other', 'AIAN', 'unknown')),
         m.edu = factor(ifelse(is.na(m.edu), 7, m.edu), labels = c('< High school', 'High School', '> High school', 'Unknown')),
         birth.wt = ifelse(birth.wt < 300, NA, 
                    ifelse(birth.wt >=10000, NA, birth.wt)),
         birth.wt.cat = ifelse(birth.wt >= 2500 & birth.wt < 4000, 0, 
                        ifelse(birth.wt >= 4000, 2, 
                        ifelse(birth.wt < 2500, 1, 3))),
         plu = plu.totalbirths,
         plu.tmp = ifelse(is.na(plu.totalbirths), 99, plu.totalbirths),
         singleton = factor(ifelse(plu.tmp == 1, 0, 
                            ifelse(plu.tmp %in% 2:98, 1, 2)),
                            labels = c('Singleton', 'Multiple', 'Unknown'))) %>% 
  select(-plu.tmp, -plu.totalbirths)

#' Add a level for missing values of birth weight category.
nj <- nj %>% 
  mutate(birth.wt.cat = factor(ifelse(is.na(birth.wt.cat), 3, birth.wt.cat),
                               labels = c('NBW','LBW','HBW', 'Unknown')))

#' Load in and harmonize GOBACK data.
goback <- readRDS('//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/goback_20230216.rds')

goback <- goback %>% 
  filter(state != 'NJ') %>% 
  mutate(birth.wt.cat = as.character(birth.wt.cat)) %>% 
  mutate(birth.wt.cat = ifelse(is.na(birth.wt.cat), 'Unknown', 
                        ifelse(birth.wt.cat == '', 'Unknown', birth.wt.cat))) %>% 
  mutate(birth.wt.cat = fct_relevel(factor(birth.wt.cat, labels = c('HBW', 'LBW', 'NBW', 'Unknown')), 'NBW', 'LBW', 'HBW', 'Unknown'),
         plu.tmp = ifelse(is.na(plu), 99, plu),
         singleton = factor(ifelse(plu.tmp == 1, 0, 
                            ifelse(plu.tmp %in% 2:98, 1, 2)), labels = c('Singleton', 'Multiple', 'Unknown'))) %>% 
  select(all_of(names(nj) ))

# Append NJ cancer & BD data ----------------------------------------------

setwd('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded_datasets/')

#' Load NJ cancer data.
nj.can <- readRDS('computed.cancer.variables.v20230104.rds') %>% 
  filter(substr(studyid,1,2) == 'nj')

#' Load NJ birth defects data.
bd <- readRDS("birth_defects_variables_20230119.rds") %>% 
  filter(state == 'NJ')

#' Join cancer and vital records data.
nj <- nj %>% 
  left_join(nj.can, by = 'studyid') %>% 
  mutate(cancer = ifelse(is.na(cancer), 0, cancer),
         across(all:num.diagnoses, ~ ifelse(cancer == 0, 0, .x)))

#' Join birth defects and vital records data.
nj <- nj %>% 
  left_join(select(bd, -c(state, coding.system)), by = 'studyid') %>% 
  mutate(birth.defect = ifelse(is.na(birth.defect), 0, birth.defect)) %>% 
  mutate(across(cns.anomaly:major.musculoskeletal.anomaly, ~ ifelse(is.na(.x) & birth.defect == 0, 0, .x)))

saveRDS(nj, 'NJ_GOBACK_data_20230405.rds')

rm(list = ls()); gc()

#' TODO: Compute number.defects variable in NJ.