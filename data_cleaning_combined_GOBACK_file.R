#' In general, the purpose of this script is to document operations made to the harmonized GOBACK analytic file.
#' Most often, this will be the process of combining birth certificate, birth defect, and cancer variables.
#' Birth defects and cancer variables are calculated using other scripts. Birth certificate variables will generally be cleaned here.

#' TODO: Rename AR variables to fit new structure.
#' TODO: Review NC person-time. Max of 9 years in kids w/o cancer and 11 years in kids with. Max of 11 years in kids w/o birth defects and 10 years in kids with.

require(tidyverse); require(magrittr); require(haven)

setwd('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/')

# Append 2022 Texas data --------------------------------------------------

goback <- readRDS('goback.v20220215.rds')

#' AR data are not compatible and need to updated separately.
#' TX and MI data are updated below.
bc.data <- goback %>% 
  filter(!state %in% c('TX','MI', 'AR')) %>% 
  select(studyid:runif, person.yrs)

rm(goback); gc()

tx.bcs <- readRDS('Old_datasets/Texas/birth.cert.data.2022.update.raw.rds')

names(tx.bcs) %<>% tolower() %<>% str_replace_all('_', '.')

#' Update missing value codes.
tx.bcs %<>% mutate(m.edu2 = ifelse(is.na(m.edu2), 7, m.edu2),
                   f.race = ifelse(is.na(f.race), 7, f.race),
                   plu = ifelse(plu == 9, NA, plu))

#' Drop implausible birth weights.
tx.bcs %<>% mutate(birth.wt = ifelse(birth.wt < 300, NA, birth.wt))

#' Harmonize variable structures.
tx.bcs <- tx.bcs %>% 
  mutate(runif = runif(nrow(tx.bcs)),
         studyid = paste0('tx',studyid),
         m.race = fct_explicit_na(factor(m.race, labels = levels(bc.data$m.race)), na_level = 'unknown'),
         m.edu2 = factor(m.edu2, labels = levels(bc.data$m.edu)),
         f.race = factor(f.race, labels = levels(bc.data$f.race))) %>% 
  rename(m.edu = m.edu2,
         person.yrs = timeatrisk.nocaner)

tx.bcs %<>% select(all_of(names(bc.data)))

#' Harmonized birth certificate data, sans Arkansas and Michigan.
bc.data <- bind_rows(bc.data, tx.bcs)

rm(tx.bcs); gc()

defects <- readRDS('Expanded_datasets/birth_defects_variables_20230119.rds') %>% 
  filter(substr(studyid,1,2) %in% c('ma','nc','ok','tx')) %>% 
  select(studyid, birth.defect, number.major.defects, cns.anomaly:major.musculoskeletal.anomaly)

cancers <- readRDS('Expanded_datasets/computed.cancer.variables.v20220603.rds') %>% 
  filter(substr(studyid,1,2) %in% c('ma','nc','ok','tx')) %>% 
  select(studyid, cancer, first.primary, all:num.diagnoses)

bc.data %<>% 
  left_join(cancers, by = 'studyid') %>% 
  mutate(cancer = ifelse(is.na(cancer), 0, cancer))

bc.data %<>% mutate(across(all:num.diagnoses, ~ ifelse(cancer == 0, 0, .x)))

bc.data %<>% 
  left_join(defects, by = 'studyid') %>% 
  mutate(birth.defect = ifelse(is.na(birth.defect), 0, birth.defect))

#' Harmonized data file for TX, MA, NC, and OK.
bc.data %<>% mutate(across(number.major.defects:major.musculoskeletal.anomaly, ~ ifelse(birth.defect == 0, 0, .x)))

saveRDS(bc.data, 'Expanded_datasets/goback.tx.ma.nc.ok.v20220608.rds')

rm(list = ls()); gc()

# Append 2022 MI data -----------------------------------------------------

goback <- readRDS('Expanded_datasets/goback.tx.ma.nc.ok.v20220608.rds')

#' Read in new MI data. Select BC columns.
mi <- readRDS('Old_Datasets/Michigan/mi.birth.cert.and.cancer.data.v20220517.rds') %>% 
  mutate(m.race = fct_recode(m.race, unknown = 'Unknown'),
         m.edu = fct_recode(m.edu, "< High school" = '< High School', "> High school" = '> High School'))

bc.cols <- select(mi, studyid:person.yrs) %>% names()

mi.bc <- select(mi, all_of(bc.cols))

rm(mi); gc()

#' Read in MI defects and cancer data.
defects <- readRDS('Expanded_datasets/birth_defects_variables_20230119.rds') %>% 
  filter(substr(studyid,1,2) == 'mi') %>% 
  select(studyid, birth.defect, number.major.defects, cns.anomaly:major.musculoskeletal.anomaly)

cancers <- readRDS('Expanded_datasets/computed.cancer.variables.v20220603.rds') %>% 
  filter(substr(studyid,1,2) == 'mi') %>% 
  select(studyid, cancer, first.primary, all:num.diagnoses)

mi.bc %<>% 
  left_join(cancers, by = 'studyid') %>% 
  mutate(cancer = ifelse(is.na(cancer), 0, cancer))

mi.bc %<>% mutate(across(all:num.diagnoses, ~ ifelse(cancer == 0, 0, .x)))

mi.bc %<>% 
  left_join(defects, by = 'studyid') %>% 
  mutate(birth.defect = ifelse(is.na(birth.defect), 0, birth.defect),
         runif = runif(nrow(.)))

#' Harmonized data file for MI.
mi.bc %<>% mutate(across(number.major.defects:major.musculoskeletal.anomaly, ~ ifelse(birth.defect == 0, 0, .x)))

rm(cancers, defects); gc()

#bc.data <- goback %>% 
#  select(all_of(bc.cols)) #%>% 
  #mutate(m.race = fct_recode(m.race, unknown = 'Unknown'))

#mi.bc %<>% select(all_of(bc.cols))

#bc.data %<>% bind_rows(mi.bc)

#rm(goback, mi.bc)

#' Read in defects and cancer data.
#defects <- readRDS('Expanded_datasets/birth_defects_variables_20230119.rds') %>% 
#  filter(substr(studyid,1,2) %in% tolower(unique(bc.data$state))) %>% 
#  select(studyid, birth.defect, number.major.defects, cns.anomaly:major.musculoskeletal.anomaly)

#cancers <- readRDS('Expanded_datasets/computed.cancer.variables.v20220603.rds') %>% 
#  filter(substr(studyid,1,2) %in% tolower(unique(bc.data$state))) %>% 
#  select(studyid, cancer, first.primary, all:num.diagnoses)

#bc.data %<>% left_join(defects, by = 'studyid')
#bc.data %<>% left_join(cancers, by = 'studyid') 

#' Set birth.defect and cancer variables to 0 in children w/o diagnoses.
#bc.data %<>% mutate(across(c(birth.defect, cancer), ~ ifelse(is.na(.x), 0, .x)))

#bc.data %<>% 
#  mutate(across(number.major.defects:major.musculoskeletal.anomaly, ~ ifelse(birth.defect == 0, 0, .x)),
#         across(all:num.diagnoses, ~ ifelse(cancer == 0, 0, .x)))

#rm(cancers, defects); gc()

mi.bc <- mi.bc %>% bind_rows(goback)

saveRDS(#bc.data, 
        mi.bc,
        'Expanded_datasets/goback.tx.ma.mi.nc.ok.v20220620.rds')

rm(list = ls()); gc()

#' Proceed to data_cleaning_SC_clean_and_append_SC_data.R.