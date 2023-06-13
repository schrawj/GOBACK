require(tidyverse); require(magrittr)

setwd('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded_datasets/')

lookup.table <- readRDS('birth.defects.lookup.table.v20230119.rds')

# Compute variables for specific BDs and organ systems --------------------

#' Load birth defects codes.
bpa.subjects <- readRDS('birth_defects_codes_bpa_states_20230417.rds')
icd9.subjects <- readRDS('birth_defects_codes_icd9_states_20230410.rds')
icd10.subjects <- readRDS('birth_defects_codes_icd10_states_20230410.rds')
mixed.system.subjects <- readRDS('birth_defects_codes_mixed_coding_states_20230523.rds')

#' Compute BD variables from BPA codes.
#' Variable assignment makes use of interpolation syntax from the 'glue' package: https://dplyr.tidyverse.org/reference/dplyr_data_masking.html
#' If any code falls within the range for the index defect, set to 1. All other children will be set to NA because they have at least one other birth defect.
#' Later, these variables will be set to zero for children who do not link to a record in the birth defects codes data frames.

#' Compute BD variables from BPA codes.
for (i in seq_along(lookup.table$diagnosis)){
  
  print(paste0('Computing ', lookup.table$diagnosis[i], '.'))
  
  new.var <- lookup.table$diagnosis[i]
  
  #bpa.subjects %<>% mutate( "{ new.var }" := if_any(starts_with('bpa'), ~ str_detect(.x, lookup.table$bpa.codes[i]), 1, NA ))
  bpa.subjects <- bpa.subjects %>% 
    mutate( "{ new.var }" := ifelse(if_any(starts_with('bpa'), ~ str_detect(.x, lookup.table$bpa.codes[i])) == T, 1, NA)) 
  #tmp <- bpa.subjects %>% 
   # mutate( "{ new.var }" := if_any(starts_with('bpa'), ~ length(str_subset(.x, lookup.table$bpa.codes[i])) >= 1      ))
  
}

#' Compute BD variables from ICD9 codes.
for (i in seq_along(lookup.table$diagnosis)){
  
  print(paste0('Computing ', lookup.table$diagnosis[i], '.'))
  
  new.var <- lookup.table$diagnosis[i]
  
  #icd9.subjects %<>% mutate( "{ new.var }" := if_any(starts_with('icd9'), ~ str_detect(.x, lookup.table$icd9.codes[i]), 1, NA ))
  #icd9.subjects %<>% mutate( "{ new.var }" := if_any(starts_with('icd9'), ~ str_detect(.x, lookup.table$icd9.codes[i])))
  icd9.subjects <- icd9.subjects %>% 
    mutate( "{ new.var }" := ifelse(if_any(starts_with('icd'), ~ str_detect(.x, lookup.table$icd9.codes[i])) == T, 1, NA)) 
  
}

#' Compute BD variables from ICD10 codes.
for (i in seq_along(lookup.table$diagnosis)){
  
  print(paste0('Computing ', lookup.table$diagnosis[i], '.'))
  
  new.var <- lookup.table$diagnosis[i]
  
  #icd10.subjects %<>% mutate( "{ new.var }" := if_any(starts_with('icd10'), ~ str_detect(.x, lookup.table$icd10.codes[i]), 1, NA ))
  #icd10.subjects %<>% mutate( "{ new.var }" := if_any(starts_with('icd10'), ~ str_detect(.x, lookup.table$icd10.codes[i])))
  icd10.subjects <- icd10.subjects %>% 
    mutate( "{ new.var }" := ifelse(if_any(starts_with('icd'), ~ str_detect(.x, lookup.table$icd10.codes[i])) == T, 1, NA)) 
  
}

#' Compute BD variables from a mixture of ICD9 and ICD10 codes.
for (i in seq_along(lookup.table$diagnosis)){
  
  print(paste0('Computing ', lookup.table$diagnosis[i], '.'))
  
  new.var <- lookup.table$diagnosis[i]
  
  mixed.system.subjects %<>% 
    mutate( "{ new.var }" := if_any(starts_with('icd'), 
                                   ~str_detect(.x, lookup.table$icd9.codes[i]) | str_detect(.x, lookup.table$icd10.codes[i]), 1, NA)) %>% 
    mutate(across(starts_with('bpa'), ~ ifelse(.x == T, 1, NA)))
    
}

tmp <- bind_rows(select(bpa.subjects, studyid:coding.system, birth.defect:deletion.22q),
                 select(icd9.subjects, studyid:coding.system, birth.defect:deletion.22q),
                 select(icd10.subjects, studyid:coding.system, birth.defect:deletion.22q),
                 select(mixed.system.subjects, studyid:coding.system, birth.defect:deletion.22q)) %>% 
  select(-c(cns.anomaly, ear.face.neck.anomaly, heart.circulatory.anom, respiratory.anomaly,
            digestive.anomaly, `genitourinary,anomaly`, musculoskeletal.anomaly, eye.anomaly))

# Compute major defects variables -----------------------------------------

new.vars <- c('major.cns.anomaly', 'major.eye.anomaly', 'major.ear.face.neck.anomaly', 'major.heart.circulatory.anomaly',
              'major.respiratory.anomaly', 'major.digestive.anomaly', 'major.genitourinary.anomaly', 'major.musculoskeletal.anomaly')

organ.system.vars <- c('cns', 'eye', 'ear.face.neck', 'heart.circ.sys', 
                       'resp.sys', 'digestive.sys', 'genitourinary', 'musculoskeletal.sys')

for (i in seq_along(new.vars)){
  
  new.var <- new.vars[i]
  
  contributing.vars <- lookup.table %>% filter(organ.system == organ.system.vars[i]) %>% pull(diagnosis)
  
  tmp <- tmp %>% 
    mutate('{ new.var }' := if_any(all_of(contributing.vars), ~ . == 1))
  
}

#' The computed major anomaly variables are logical. Set to 1 if true, NA otherwise.
tmp %<>% mutate(across(starts_with('major'), ~ ifelse(.x == F, NA, 
                                               ifelse(.x == T, 1, .x))))

#' Compute a variable for number of major defects.
vars <- lookup.table %>% filter(!is.na(organ.system), organ.system != 'syndromes') %>% pull(diagnosis)

#' Sum across individual birth defects variables to calculate the number of major defects.
#' Then, set to missing in children with syndromes, so that they are not accidentally included in these analyses moving forward.
syndrome.ids <- tmp %>% filter(genetic.anomaly == 1) %>% pull(studyid)

tmp <- tmp %>% 
  mutate(number.major.defects = rowSums(tmp[,vars], na.rm = T)) %>% 
  mutate(number.major.defects = ifelse(studyid %in% syndrome.ids, NA, number.major.defects))

saveRDS(tmp, 
        'birth_defects_variables_20230523.rds')

rm(list = ls()); gc()
