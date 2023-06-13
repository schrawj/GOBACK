#' In general, the purpose of this script is to document operations made to the harmonized GOBACK analytic file.
#' Most often, this will be the process of combining birth certificate, birth defect, and cancer variables.
#' Birth defects and cancer variables are calculated using other scripts. Birth certificate variables will generally be cleaned here.

#' TODO: Rename AR variables to fit new structure.
require(tidyverse); require(magrittr); require(haven)

setwd('/mount/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/')

goback <- readRDS('Expanded_datasets/goback.tx.ma.mi.nc.ok.v20220620.rds') %>% 
  as_tibble() %>% 
  mutate(sex = factor(sex, labels = c('Male', 'Female', 'Unknown')))

vr.variables <- goback %>% select(studyid:person.yrs) %>% names()

#' Read in and clean SC vital records data.
#' Paternal data were not collected.
#' Time variables are measured on different scales.
sc <- readRDS('Old_datasets/South_Carolina/south.carolina.raw.data.rds') %>% 
  rename(m.edu = m.edu2) %>% 
  mutate(f.race = factor(NA),
         f.age = as.numeric(NA),
         timeatrisk.nocancer = timeatrisk.nocancer/12,
         age.days = age.days/365.2,
         studyid = paste0('sc', studyid),
         sex = factor(sex, labels = c('Male','Female')),
         m.edu = factor(m.edu, labels = levels(goback$m.edu)),
         m.race = factor(m.race, labels = c('Hispanic', 'NHW', 'NHB', 'Other', 'unknown')))

sc$cancer <- ifelse(is.na(sc$cancer), 0, sc$cancer)

sc %<>% mutate(person.yrs = ifelse(cancer == 1, age.days, timeatrisk.nocancer)) 

#' These are children without cancer who have a calculated person-time of zero. 
#' Some were not even noted to have died during infancy.
possible.stillbirths <- filter(sc, person.yrs == 0, cancer == 0) %>% pull(studyid)

sc %<>% 
  filter(!studyid %in% possible.stillbirths) %>% 
  select(all_of(vr.variables)) %>% 
  mutate(birth.wt = ifelse(birth.wt < 300, NA, birth.wt))

#' Read in cancer and birth defects data.
sc.can <- readRDS('Expanded_datasets/computed.cancer.variables.v20220811.rds') %>% 
  filter(substr(studyid,1,2) == 'sc', !duplicated(studyid))

sc.bd <- readRDS('Expanded_datasets/birth_defects_variables_20230119.rds') %>% 
  filter(state == 'SC') %>% 
  as_tibble()

#' Join cancer and vital records data.
sc %<>% 
  left_join(sc.can, by = 'studyid') %>% 
  mutate(cancer = ifelse(is.na(cancer), 0, cancer))

sc %<>% mutate(across(all:num.diagnoses, ~ ifelse(cancer == 0, 0, .x)))

#' Join birth defects and vital records data.
sc %<>% 
  left_join(select(sc.bd, -c(state, coding.system)), by = 'studyid') %>% 
  mutate(birth.defect = ifelse(is.na(birth.defect), 0, birth.defect),
         runif = runif(nrow(.))) %>% 
  mutate(across(cns.anomaly:major.musculoskeletal.anomaly, ~ ifelse(is.na(.x) & birth.defect == 0, 0, .x)))

saveRDS(sc, 'Expanded_datasets/SC_GOBACK_data_20230105.rds')

rm(sc.bd, sc.can, vr.variables, possible.stillbirths); gc()

goback %<>% bind_rows(sc)

goback <- goback %>% mutate(across(number.major.defects:major.musculoskeletal.anomaly, ~ ifelse(is.na(.x) & birth.defect == 0, 0, .x)))

saveRDS(goback, 'goback_20221117.rds')

rm(list = ls()); gc()

#' Proceed to data_cleaning_
#' Proceed to data_cleaning_modify_combined_GOBACK_file.R.