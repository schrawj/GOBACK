require(tidyverse)

#setwd('/mount/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/')
setwd('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/')

#' Current GOBACK and NJ datasets.
goback <- readRDS("goback_20230105.rds")

nj <- readRDS('Expanded_datasets/NJ_GOBACK_data_20230405.rds') %>% 
  select(studyid:plu, singleton, person.yrs, all:number.major.defects) %>%  # Re-order columns so plu and singleton are adjacent.
  mutate(runif = runif(nrow(.)))

goback <- nj %>% bind_rows(filter(goback, state != 'NJ'))

goback <- goback %>% mutate(across(number.major.defects:major.musculoskeletal.anomaly, ~ ifelse(is.na(.x) & birth.defect == 0, 0, .x)))

saveRDS(goback, 'goback_20230405.rds')
