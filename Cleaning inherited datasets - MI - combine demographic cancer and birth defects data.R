require(tidyverse); require(magrittr)

setwd('//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/')

goback <- readRDS('goback_20230405.rds')

goback %<>% filter(state != 'MI')

#' Retrieve MI birth cert data.
setwd('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/')

mi <- readRDS('Old_Datasets/Michigan/mi.v20220516.rds')

names(mi) %<>% tolower() %>% str_replace_all('_', '.')

mi %<>% 
  select(studyid, state, sex, m.race, m.edu2, m.age, birth.wt, birth.wt.cat, gest.age, birth.yr, plu, f.race, f.age, person.yrs) %>% 
  rename(m.edu = m.edu2) %>% 
  mutate(birth.wt.cat = ifelse(birth.wt.cat == "", NA, birth.wt.cat)) %>% 
  mutate(runif = runif(nrow(.)),
         studyid = paste0('mi',studyid),
         m.race = factor(m.race, labels = c('Hispanic','NHW','NHB','Asian','Other','AIAN', 'Unknown')),
         f.race = factor(f.race, labels = c('Hispanic','NHW','NHB','Asian','Other','AIAN', 'Unknown')),
         m.edu = factor(m.edu, labels = c('< High School', 'High School', '> High School', 'Unknown')),
         birth.wt.cat = relevel(factor(birth.wt.cat, labels = c('HBW','LBW','NBW')), ref = 'NBW'))

#' Join with MI cancer data.
mi.cancer <- readRDS('Expanded_datasets/computed.cancer.variables.v20220516.rds')
mi.cancer %<>% filter(substr(studyid, 1, 2) == 'mi')

mi %<>% 
  left_join(mi.cancer, by = 'studyid') %>% 
  mutate(num.diagnoses = ifelse(is.na(num.diagnoses), 0, num.diagnoses),
         cancer = ifelse(num.diagnoses >= 1, 1, 0))

#' Set birthweights below 300g to missing.
mi %<>% mutate(birth.wt = ifelse(birth.wt < 300, NA, birth.wt),
               birth.wt.cat = ifelse(is.na(birth.wt), NA, birth.wt.cat))

#' Back to a factor with the same levels as GOBACK.
mi %<>% mutate(birth.wt.cat = factor(ifelse(birth.wt.cat == 1, 1, 
                                     ifelse(birth.wt.cat == 2, 3, 2)),
                               labels = c('NBW','LBW','HBW')))

saveRDS(mi, 
        'Old_Datasets/Michigan/mi.birth.cert.and.cancer.data.v20230406.rds')
