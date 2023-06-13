#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2021.09.07.
#' 
#' Bind columns for the MA demographics, cancer, and birth defects variables.
#' 
#' Reformat so that they can be appended to the GOBACK analytic file.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(tidyverse); require(haven); require(magrittr)

setwd('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/')

ma.cancer <- readRDS('Expanded datasets/ma.cancer.variables.v20210907.rds') %>% 
  mutate(KIDUID = paste0('ma', KIDUID)) %>% 
  select(-cancer, -laterality1)

load('Expanded datasets/computed.bd.vars.mancoktx.v20210809.rdata')

ma.bd <- filter(birth.defects.results, substr(studyid, 1, 2) == 'ma')

cols <- grep('canonical', names(ma.bd), invert = T)

ma.bd <- ma.bd[ , cols]

ma.demo <- read_dta('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/BD-CC Projects/GOBACK/Massachusetts Data/Data (June 2021)/Full dataset_harmonized birth certificate variables_092421.dta') 

#' This demographics file seems to have 2 variables for sex, one of which is basically just blank.
ma.demo <- ma.demo[ , -131]

names(ma.demo) %<>% tolower() %>% str_replace_all('_','.')

ma.demo %<>%  
  mutate(studyid = paste0('ma', studyid),
         cancer = ifelse(is.na(cancer), 0, 1)) %>% 
  filter(!duplicated(studyid)) 

#tmp <- filter(ma.demo, timeatrisk.nocancer == 0) %>% 
#  mutate(death.year = as.numeric(substr(dod2,1,4)))

ma <- left_join(ma.demo, ma.cancer, by = c('studyid' = 'KIDUID')) %>% 
  filter(!duplicated(studyid)) %>% 
  left_join(ma.bd, by = 'studyid') %>% 
  mutate(m.race = ifelse(is.na(m.race), 7, m.race)) %>% 
  mutate(f.race = 7,
         f.age = as.numeric(NA),
         runif = runif(nrow(.), 0, 1),
         m.race = factor(m.race, 
                         labels = c('Hispanic','NHW','NHB','Asian','Other', 'UNK')),
         birth.wt.cat = factor(ifelse(birth.wt.cat == 'NBW',0, 
                               ifelse(birth.wt.cat == 'HBW', 1, 2)),
                           labels = c('NBW','HBW','LBW')),
         person.yrs = ifelse(is.na(person.yrs) & cancer == 0, timeatrisk.nocancer/12, person.yrs),
         num.diagnoses = ifelse(cancer == 0, 0, num.diagnoses),
         any.birth.defect = ifelse(is.na(any.birth.defect), 0, any.birth.defect))

rm(ma.cancer, ma.bd, ma.demo, birth.defects.results); gc()

#' The only values for specific birth defects are 1 (the expected value when the index child has the index birth defect) and
#' NA (the expected value when the index child has one or more other birth defects). Children should have the value 0 if they have no birth 
#' defects. Set birth defect variables to zero if they have none,
update <- function(col){ ifelse(ma$studyid %in% no.bd.ids, 0, col) }

no.bd.ids <- ma %>% filter(any.birth.defect == 0) %>% select(studyid) %>% unlist()

ma  %<>% mutate(across(conganomalies.cns:conotruncal.defects, update))

#' TODO: Harmonize factor variable structures to GOBACK.
#' TODO: Convert f.race to a factor variable in GOBACK.

goback <- readRDS("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20210616.rds") %>% 
  select(-f.edu2, -birth.char.flag)

goback  %<>% 
  bind_rows(select(ma, all_of(names(goback)))) %>% 
  select(-state.num)

saveRDS(goback, file = 'goback.v20211026.rds')

rm(list = ls()); gc()
