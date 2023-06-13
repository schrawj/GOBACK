
# Prepare environment -----------------------------------------------------

require(tidyverse); require(magrittr)

compute.cancer.vars <- function(input.data, output.data){
  
  ids <-input.data %>% filter(iccc.extended.classification %in% codes[[i]]) %>% pull(studyid)
  
  cancer.name <- cancers[[i]]
  
  output.data %<>% mutate('{cancer.name}' := ifelse(studyid %in% ids, 1, 0))
  
}

compute.rms.vars <- function(input.data, output.data){
  
  erms.ids <- filter(input.data, iccc.extended.classification == 55, histology == 8910) %>% pull(studyid) #%>% unlist() %>% as.numeric()
  
  arms.ids <- filter(input.data, iccc.extended.classification == 55, histology == 8920) %>% pull(studyid) #%>% unlist() %>% as.numeric()
  
  rms.other.ids <- setdiff(unlist(select(filter(input.data, iccc.extended.classification == 55), studyid)) , 
                           union(erms.ids, arms.ids))
  
  output.data %<>% mutate(erms = ifelse(studyid %in% erms.ids, 1, 0),
                          arms = ifelse(studyid %in% arms.ids, 1, 0),
                          rms.other = ifelse(studyid %in% rms.other.ids, 1, 0)
  )
  
  return(output.data)
  
}

#' Extended recodes for each tumor type, per https://seer.cancer.gov/iccc/iccc-iarc-2017.html.  
codes <- list(1, 5, 1:8, c(2,3,4,6:8),
              9, 10:13, 9:16, 14:16,
              19, 20, 17, 21, 17:32, c(18,22:32),
              33, 33:34, 34,
              35,
              36, 36:40, 37:40,
              41, 41:45, 42:45,
              46, c(48,60), 46:54, c(47,49:54),
              55, 
              54:70, c(56:59, 61:71),
              72:77, 78:83, 84:92, 72:92,
              93:108,
              109:115)

#' The names of those tumor types, as they will appear in the modified GOBACK file.
cancers <- list('all', 'aml', 'leu.any', 'leu.other',
                'hl', 'nhl', 'lym.any', 'lym.other',
                'astro','medullo','ependymoma','pnet','cns.any','cns.other',
                'neuro','pns.any','pns.other',
                'retino',
                'nephro', 'renal.any', 'renal.other',
                'hepato', 'hepatic.any', 'hepatic.other',
                'osteo', 'ewing', 'bone.any', 'bone.other',
                'rms.any',
                'soft.any', 'soft.other',
                'gct.intra', 'gct.extra', 'gct.gonad', 'gct.any',
                'epithe',
                'other.any')

setwd('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/')

# Load datasets and encode variables --------------------------------------

cancer.codes <- readRDS('Datasets/Expanded_datasets/cancer_codes_20230216.rds')

load('Datasets/Expanded_datasets/iccc.lookup.table.v20210826.rdata')

#' Convert cancer codes from wide to long.
histologies <- cancer.codes %>%  
  select(studyid, starts_with('hist')) %>% 
  pivot_longer(cols = starts_with('hist'),
               names_to = 'cancer',
               names_prefix = 'histology.',
               values_to = 'histology') 

sites <- cancer.codes %>% 
  select(studyid, starts_with('site')) %>% 
  pivot_longer(cols = starts_with('site'),
               names_to = 'cancer',
               names_prefix = 'site.',
               values_to = 'site') %>% 
  mutate(site = ifelse(!is.na(site), paste0('C',site), NA))

behavior <- cancer.codes %>% 
  select(studyid, starts_with('behav')) %>% 
  pivot_longer(cols = starts_with('behav'),
               names_to = 'cancer',
               names_prefix = 'behavior.',
               values_to = 'behavior')

#' Join the long data frames together.
cancer.long <- 
  left_join(histologies, sites, by = c('studyid','cancer')) %>% 
  left_join(behavior, by = c('studyid','cancer')) %>% 
  mutate(site = ifelse(str_length(site) == 3, 
                       paste0(substr(site,1,1), '0', substr(site,2,3)),
                       site))

#' Join the ICCC encodings.
cancer.long %<>% left_join(select(iccc.lookup, iccc.extended.classification, iccc.name, icdo3.hist, icdo3.site.code), 
                           by = c('histology' = 'icdo3.hist', 'site' = 'icdo3.site.code')) %>% 
  mutate(iccc.extended.classification = as.numeric(iccc.extended.classification))

#' Compute calculated cancer variables.
cancer.variables <- select(cancer.long, studyid)

for (i in 1:length(codes)){
  
  cancer.variables <- compute.cancer.vars(cancer.long, cancer.variables)
  
}

cancer.variables <- compute.rms.vars(cancer.long, cancer.variables)

#' Each row contains identical information; there are five rows per child regardless of the number of cancers.
#' Delete the 2nd:5th.
cancer.variables %<>% filter(!duplicated(studyid))

cancer.variables %<>% rename(other.tumor = other.any) 

#' Select columns for specific cancers and sum across them.
tmp <- cancer.variables %>% select(!ends_with('.any')) 

tmp$num.diagnoses <- rowSums(tmp[,-1], na.rm = T)

#' A small number of children (~0.1% as of this writing) have histology/site codes that don't align with an ICCC extended class.
#' Set number of diagnoses == 1 and 'other.tumor == 1'.
tmp %<>% mutate(other.tumor = ifelse(num.diagnoses == 0, 1, other.tumor),
                num.diagnoses = ifelse(num.diagnoses == 0, 1, num.diagnoses))

#' Add these columns to the cancer.variables data frame.
cancer.variables %<>% 
  select(-other.tumor) %>% 
  left_join(select(tmp, studyid, other.tumor, num.diagnoses), by = 'studyid')

#' Add a variable for first primary, by pulling it from the 'wide' codes data frame.
iccc.lookup %<>% mutate(icdo3.site.code = as.numeric(str_remove_all(icdo3.site.code, 'C')))

cancer.codes  %<>% 
  left_join(select(iccc.lookup, icdo3.hist, icdo3.site.code, iccc.name), 
            by = c('histology.1' = 'icdo3.hist', 'site.1' = 'icdo3.site.code')) %>% 
  rename(first.primary = iccc.name)
  
cancer.variables %<>% left_join(select(cancer.codes, studyid, first.primary),
                                by = 'studyid')

#' Add an indicator variable for whether any cancer. Presumed to be yes in all cases.
cancer.variables$cancer <- 1

saveRDS(cancer.variables, 'Datasets/Expanded_datasets/computed_cancer_variables_20230216.rds')

rm(list = ls()); gc()
