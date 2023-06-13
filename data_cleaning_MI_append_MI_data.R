require(tidyverse)

#' TODO: Figure out why no MI gastroschisis cases and fix.

#' Load current GOBACK data and remove old MI data.
goback <- readRDS("//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/old_GOBACK_datasets/goback_20230216.rds")

#' Fix an issue with the paternal race variable's format.
goback <- goback %>% 
  filter(state != 'MI') %>% 
  mutate(f.race = fct_collapse(fct_na_value_to_level(f.race, 'unknown'),
                               unknown = c('Unknown', 'unknown')))

vital.cols <- goback %>% slice(1) %>% select(studyid:person.yrs) %>% names()

#' COmbine MI datasets.
mi.bd <- readRDS('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded_datasets/birth_defects_variables_20230411.rds') %>% 
  filter(state == 'MI') %>% 
  select(-state)

mi.vital <- readRDS('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old_datasets/Michigan/mi.birth.cert.and.cancer.data.v20230406.rds') 

tmp <- mi.vital$sex %>% as.character()
  
mi.vital <- mi.vital %>% 
  mutate(singleton = factor(case_when(plu == 1 ~ 0,
                                      plu > 1 ~ 1,
                                      is.na(plu) ~ 2,
                                      .default = as.numeric(plu)),
                            labels = c('Singleton', 'Multiple', 'Unknown')),
         sex = factor(tmp, labels = c('Male', 'Female', 'Unknown')),
         m.race = factor(m.race, labels = levels(goback$m.race))) %>% 
  select(all_of(vital.cols))

mi <- mi.vital %>% 
  left_join(mi.bd, by = 'studyid')

mi.cancer <- readRDS('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded_datasets/computed_cancer_variables_20230216.rds') %>% 
  filter(substr(studyid,1,2) == 'mi')

mi <- left_join(mi, mi.cancer, by = 'studyid')

#' Update incorrect missing values for kids without birth defects or cancer.
mi <- mi %>% 
  mutate(across(c(birth.defect, cancer), ~ ifelse(is.na(.x), 0, .x)),
         across(anencephalus:major.musculoskeletal.anomaly, ~ ifelse(birth.defect == 0, 0, .x)),
         across(all:other.tumor, ~ ifelse(cancer == 0, 0, .x)),
         number.major.defects = ifelse(birth.defect == 0, 0, number.major.defects),
         num.diagnoses = ifelse(cancer == 0, 0, num.diagnoses))

#' Append to GOBACK.
goback <- goback %>% 
  bind_rows(mi)

#' Save updated GOBACK files.
saveRDS(goback,
        '//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/goback_20230411.rds')

sub <- goback %>% filter(runif < 0.1)

saveRDS(sub,
        '//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/goback_random_subset_20230411.rds')
