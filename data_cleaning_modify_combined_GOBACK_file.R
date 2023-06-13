require(tidyverse); require(magrittr)

setwd('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/')

# Some updates to demographic variables -----------------------------------

goback <- readRDS('goback.v20211026.rds')

goback %<>% 
  rename(m.edu = m.edu2) %>% 
  mutate(m.race = fct_collapse(m.race, unknown = c('Unknown', 'UNK')),
         m.edu = factor(ifelse(is.na(goback$m.edu2), 7, goback$m.edu2), 
                        labels = c('< High school', 'High School', '> High school', 'Unknown')),
         birth.wt = ifelse(birth.wt <300, NA, 
                    ifelse(birth.wt >9999, NA, birth.wt)),
         birth.wt.cat = factor(ifelse(birth.wt < 2500, 1, 
                               ifelse(birth.wt >= 2500 & birth.wt <4000, 0,
                               ifelse(birth.wt >= 4000, 2, 3))),
                        labels = c('NBW', 'LBW', 'HBW')),
         f.race = factor(ifelse(is.na(goback$f.race), 7, goback$f.race),
                  labels = c('Hispanic', 'NHW', 'NHB', 'Asian', 'Other', 'AIAN', 'Unknown'))
         )

saveRDS(goback, file = 'goback.v20220204.rds')

# Update cancer variables based on new algorithm --------------------------

goback <- readRDS('goback.v20220204.rds')

goback %<>% 
  select(-laterality1, -cancer1) %>% 
  rename(other.tumor = other.any)

cancer <- readRDS('Expanded datasets/computed.cancer.variables.v20220211.rds') %>% 
  mutate(cancer = 1)

#' SPLIT BY CANCER STATUS. USE STUDYIDS IN CANCER DATA FRAME TO DETERMINE.
cancer.ids <- pull(cancer, studyid)

goback.cancer <- filter(goback, studyid %in% cancer.ids)

goback %<>% filter(!studyid %in% cancer.ids)

#' NO CANCER - CREATE FIRST PRIMARY VARIABLE WITH VALUE 'NONE'. Set the other cancer variables to zero.
goback %<>% mutate(first.primary = 'None')

set.to.zero <- function(x){ x = 0}

goback %<>% mutate(across(cancer:gct.any, set.to.zero))

#' Cancer - Set num.diagnoses to missing in AR, since they only provided data on first primary tumor.

goback.cancer %<>% mutate(num.diagnoses = ifelse(state == 'AR', NA, num.diagnoses))

#' CANCER - DROP OLD CANCER COLUMNS.
goback.cancer %<>% select(!(cancer:gct.any))

#' CANCER - BIND IN NEW CANCER COLUMNS.
goback.cancer %<>% left_join(cancer, by = 'studyid')

#' BIND_ROWS TO MERGE THE TWO AND SAVE.
goback %<>% bind_rows(goback.cancer)

saveRDS(goback, file = 'goback.v20220211.rds')

rm(list = ls())

# Update BD variables based on new CDC-BPA algorithm ----------------------

goback <- readRDS('goback.v20220214.rds') %>% 
  select(studyid, state, any.birth.defect:deletion.22q)

birth.defects.variables <- readRDS('Expanded datasets/computed.birth.defects.variabls.cdc.bpa.states.v20220214.rds')

cases <- goback %>% filter(studyid %in% birth.defects.variables$studyid)

goback  %<>% filter(!studyid %in% birth.defects.variables$studyid)

#' Set all birth defects variables for all MA, NC, OK, and TX children in GOBACK to 0.
#' This may already be how the data are structured, but let's make sure.
for (i in 3:ncol(goback)){
  
  goback[, i] <- ifelse(goback$state %in% c('MA', 'NC', 'OK', 'TX'), 0, goback[,i])
  
}

#' PUll out analytic BD variables from birth defects variable file using cases as a reference.
birth.defects.variables %<>% 
  mutate(state = toupper(substr(studyid, 1, 2))) %>% 
  select(all_of(names(goback)))

#' Drop analytic variables from cases file and bind in updated ones from birth defects variables DF.
cases %<>% 
  select(studyid) %>% 
  left_join(birth.defects.variables, by = 'studyid')
  
#' Bind GOBACK and cases data frames back together.
goback.bd.vars <- bind_rows(goback, cases)

rm(birth.defects.variables, cases, goback, i); gc()

#' Join the new birth defects variables to the demographic and cancer variables.
goback <- readRDS('goback.v20220214.rds') %>% 
  select(!(any.birth.defect:deletion.22q))

#goback <- readRDS('goback.v20220211.rds') %>% 
#  select(!(any.birth.defect:conotruncal.defects))

goback %<>% left_join(goback.bd.vars, by = c('studyid', 'state'))

rm(goback.bd.vars); gc()

#' Rearrange columns: demographics, then birth defects, then cancer.
goback %<>% select(studyid:f.age, runif, 
                   any.birth.defect:deletion.22q,
                   cancer, first.primary, person.yrs, num.diagnoses, all:gct.any) %>% 
  rename(tapvr = total.anomalous.pulmonary.venous.return)

saveRDS(goback, 
        file = 'goback.v20220215.rds')

write_csv(goback, ('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/R outputs/goback.v20220215.csv'))

# Generate syndromic and non-syndromic files ------------------------------

setwd('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/')

goback <- readRDS('goback.v20220215.rds')

goback %<>% filter(any.genetic.anomaly == 0)

saveRDS(goback, 'goback.non.sydromic.v20220520.rds')

# Add variable for singleton or not ---------------------------------------

goback <- readRDS("goback_20221117.rds")

#' Corrections to missing values in the categorical birthweight variable.
#' Compute a varible denoting singleton or multiple birth.
goback <- goback %>% 
  mutate(birth.wt.cat = as.character(birth.wt.cat)) %>% 
  mutate(birth.wt.cat = ifelse(is.na(birth.wt.cat), 'Unknown', 
                        ifelse(birth.wt.cat == '', 'Unknown', birth.wt.cat))) %>% 
  mutate(birth.wt.cat = fct_relevel(factor(birth.wt.cat, labels = c('HBW', 'LBW', 'NBW', 'Unknown')), 'NBW', 'LBW', 'HBW', 'Unknown'),
         plu.tmp = ifelse(is.na(plu), 99, plu),
         singleton = factor(ifelse(plu.tmp == 1, 0, 
                            ifelse(plu.tmp %in% 2:98, 1, 2)), labels = c('Singleton', 'Multiple', 'Unknown'))) %>% 
  select(-plu.tmp)

saveRDS(goback, 'goback_20230104.rds')

rm(list = ls()); gc()

#' Proceed to data_cleaning_NJ_append_NJ_data.R.
#' 

# Corrected OK BD variables -----------------------------------------------

#' Done on 4/17/2023.

#' Updated computed birth defects variables.
bds <- readRDS('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded_datasets/birth_defects_variables_20230417.rds')

#' Existing GOBACK data.
goback <- readRDS('//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/goback_20230411.rds')

#' Split; update BD variables in affected kids; recombine.
goback.bd <- filter(goback, studyid %in% bds$studyid)

goback <- goback %>% filter(!studyid %in% bds$studyid)

cols <- setdiff(names(goback), names(bds))

#' select all the columns that AREN'T in the new BD data frame, and add in the new variables.
#' Some OK kids have NA values; they have codes for conditions we are not considering as birth defects.
goback.bd <- goback.bd %>% 
  select(studyid, all_of(cols)) %>% 
  left_join(bds, by = 'studyid')

goback <- goback %>% 
  bind_rows(goback.bd) %>% 
  mutate(birth.defect = ifelse(is.na(birth.defect), 0, birth.defect),
         number.major.defects = ifelse(birth.defect == 0, 0, number.major.defects)) %>% 
  select(-coding.system)

saveRDS(goback, 
        '//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/goback_20230411.rds')



# Remove old, questionable NJ data ----------------------------------------

#' Reported birth defects prevalence was extremely low in NJ prior to 1990.
goback <- readRDS('//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/goback_20230411.rds')

goback <- goback %>% 
  filter(state != 'NJ' | (state == 'NJ' & birth.yr >= 1990) )

saveRDS(goback,
        '//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/goback_20230418.rds')

# Updated organ system major BD vars --------------------------------------

#' Performed 5/23/2023.
#' New BD variables.
bd <- readRDS('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded_datasets/birth_defects_variables_20230523.rds') %>% 
  filter(!is.na(birth.defect))

#' Current GOBACK file.
goback <- readRDS('//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/goback_20230418.rds')

#' Remove old BD variables.
goback <- goback %>% select(-c(birth.defect:number.major.defects))

#' Join in new bd variables.
goback <- goback %>% left_join(select(bd, -c(state, coding.system)), by = 'studyid')

goback <- goback %>% 
  mutate(birth.defect = ifelse(is.na(birth.defect), 0, birth.defect),
         number.major.defects = ifelse(birth.defect == 0, 0, number.major.defects),
         across(anencephalus:major.musculoskeletal.anomaly, ~ifelse(birth.defect == 0, 0, .x)))

saveRDS(goback,
        '//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/goback_20230523.rds')

# Subset the newest dataset  ----------------------------------------------

#' For testing code and to initiate people to the GOBACK file, I create a version that 
#' includes a random 10% of the data.
goback <- readRDS("//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/goback_20230523.rds")

goback <- goback %>% 
  filter(runif <= 0.1)

saveRDS(goback,
        '//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/goback_random_subset_20230523.rds')
