require(tidyverse); require(lubridate); require(haven)

#' TODO: Clean birth defects variables.
#' 5. Compute BD variables.
#' 6. Join computed BD variables with other states.

# Read and save raw data --------------------------------------------------

fl <- read_dta('//smb-main.ad.bcm.edu/genepi3/TiffanyChambers/GOBACK Data/Florida/all years-harmonized(clean).dta')

saveRDS(fl, '//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/raw_data/florida_raw_data_20230215.rds')

rm(list = ls())

# Clean birth certificate data --------------------------------------------

setwd('//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/')

#' For referencing existing column names and structures.
goback <- readRDS('goback_20230105.rds')
goback <- goback %>% 
  slice(1:1000)

vr <- goback %>% 
  select(studyid:person.yrs) %>% 
  names()

goback.vr <- goback %>% 
  select(all_of(vr))

fl <- readRDS('raw_data/florida_raw_data_20230215.rds') %>% 
  as_tibble() %>% 
  rename_with(~ tolower(str_replace_all(.x,'_','.')))

setdiff(vr, names(fl))

#' Add missing vital records columns.
fl <- fl %>% 
  mutate(studyid = paste0('fl',state.file.number),
         dob = as_date(dob),
         sex = factor(sex, labels = c('Male', 'Female', 'Unknown')),
         birth.wt = ifelse(birth.wt < 300 | birth.wt > 10000, NA, birth.wt),
         m.race = factor(m.race, labels = c('Hispanic','NHW', 'NHB', 'Asian', 'Other', 'AIAN', 'unknown')),
         f.race = factor(f.race, labels = c('Hispanic','NHW', 'NHB', 'Asian', 'Other', 'AIAN', 'unknown')),
         f.age= ifelse(f.age < 10 | f.age > 100, NA, f.age),
         m.edu = factor(m.edu2, labels = c('< High school', 'High School', '> High school', 'Unknown')),
         birth.wt.cat = factor(ifelse(is.na(birth.wt.cat), 4, birth.wt.cat), labels = c('NBW', 'LBW', 'HBW', 'Unknown')),
         infantdied = ifelse(is.na(infantdied), 0, infantdied),
         cancer = ifelse(is.na(cancer), 0, 
                  ifelse(timetodiag.n3901 >= 5840, 0, cancer)),
         time.to.death.missing = ifelse(is.na(timetodeath.mon), 1, 0),
         state = 'FL',
         plu = ifelse(plu %in% c(9,99), NA, plu),
         singleton = case_when(plu == 1 ~ 'Singelton',
                               plu %in% 2:6 ~ 'Multiple',
                               is.na(plu) ~ 'Unknown'),
         person.yrs = case_when(cancer == 0 & time.to.death.missing == 0 & (timetodeath.mon/12) <= 16 ~ timetodeath.mon/12,
                                cancer == 0 & time.to.death.missing == 0 & (timetodeath.mon/12) > 16 ~ 16,
                                cancer == 0 & time.to.death.missing == 1 ~ as.numeric(as_date('2013-12-31') - dob)/365,
                                cancer == 1 & timetodiag.n3901 <= 0 ~ 0,
                                cancer == 1 & timetodiag.n3901 > 0 & timetodiag.n3901 < 5840 ~ timetodiag.n3901/365))

fl <- fl %>% 
  select(all_of(vr))

tmp <- bind_rows(goback, fl)

saveRDS(fl, 
        '//smb-main.ad.bcm.edu/genepi3/JeremySchraw/GOBACK/Datasets/fl.vital.records.raw.data.rds')

rm(list = ls()); gc()

# Clean birth defects data ------------------------------------------------

#' The existing final data structure, for reference.
bd.ref <- readRDS('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded_datasets/birth_defects_codes_icd9_states.v20220609.rds')

setwd('//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/')

fl <- readRDS('raw_data/florida_raw_data_20230215.rds') %>% 
  as_tibble() %>% 
  rename_with(~ tolower(str_replace_all(.x,'_','.')))

#' Pull and rename, or create, FL birth defect columns.
bd <- fl %>% 
  select(state.file.number, dob, starts_with('bd.code')) %>% 
  mutate(studyid = paste0('fl',state.file.number),
         state = 'FL',
         coding.system = 'ICD-9',
         #across(starts_with('bd.code'), ~ as.character(.x)),
         across(starts_with('bd.code'), ~ ifelse(nchar(as.character(.x)) == 3, paste0(.x,'.0'), as.character(.x)))
         ) %>% 
  rename_with(~str_replace(.x, 'bd.code','icd9code'), starts_with('bd.code'))

for (i in 20:30){
  
  new.col <- paste0('icd9code',i)
  
  bd <- bd %>% 
    mutate('{new.col}' := as.character(NA))
  
}

bd <- bd %>% 
  select(all_of(names(bd.ref))) %>% 
  filter(if_any(starts_with('icd9'), ~ !is.na(.x)))

saveRDS(bd, 
        '//smb-main.ad.bcm.edu/genepi3/JeremySchraw/GOBACK/Datasets/fl.birth.defect.raw.data.rds')

rm(list = ls()); gc()

#' Proceed to data_cleaning_generate_a_data_frame_of_birth_defects_codes.Rs

# Clean cancer data -------------------------------------------------------

#' Final data structure, for reference.
codes <- readRDS('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded_datasets/cancer.codes.v20230104.rds')

#' Load FL data.
setwd('//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/')

fl <- readRDS('raw_data/florida_raw_data_20230215.rds') %>% 
  as_tibble() %>% 
  rename_with(~ tolower(str_replace_all(.x,'_','.'))) %>% 
  filter(!is.na(histologic.type.icd.o.3.n5221))

#' Select and rename FL cancer columns.
fl <- fl %>% 
  select(state.file.number,
         starts_with('histologic'),
         starts_with('seer.summary.stage.2000'),
         starts_with('laterality'),
         starts_with('primary.site'),
         starts_with('behavior')) %>% 
  rename_with(~ str_replace_all(.x, 'histologic.type.icd.o.3.n522', 'histology.'), starts_with('histologic.type.icd.o.3.n522')) %>% 
  rename_with(~ str_replace_all(.x, 'seer.summary.stage.2000.n759', 'seer.stage.'), starts_with('seer.summary.stage.2000.n759')) %>% 
  rename_with(~ str_remove_all(.x, 'n410'), starts_with('lateral')) %>% 
  rename_with(~ str_replace_all(.x, 'primary.site.n400', 'site.'), starts_with('primary.site.n400')) %>% 
  rename_with(~ str_remove(.x, 'code.icd.o.3.n523'), starts_with('behavior.code.icd.o.3.n523')) %>% 
  mutate(studyid = paste0('fl', state.file.number),
         across(starts_with('site'), ~ as.numeric(str_remove(.x, '^C')))) %>% 
  select(-state.file.number)

#' Add missing columns for fifth primary.
col.types <- fl %>% 
  select(-studyid) %>% 
  names() %>% 
  str_remove_all('[:digit:]') %>% 
  unique()

for (i in col.types){
  
  new.col <- paste0(i,5)
  
  fl <- fl %>% 
    mutate('{new.col}' := as.double(NA))
  
}

#' Select and sort columns.
fl <- fl %>% 
  select(all_of(names(codes)))

saveRDS(fl, 
        '//smb-main.ad.bcm.edu/genepi3/JeremySchraw/GOBACK/Datasets/fl.cancer.raw.data.rds')

rm(list = ls()); gc()
# Combine FL datasets -----------------------------------------------------

setwd('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded_datasets/')

fl.vit <- readRDS('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/GOBACK/Datasets/fl.vital.records.raw.data.rds')

#' Load FL cancer data. Variables computed elsewhere from raw data created above.
fl.can <- readRDS('computed_cancer_variables_20230216.rds') %>% 
  filter(substr(studyid,1,2) == 'fl')

#' Load FL birth defects data. Variables computed elsewhere from raw data created above.
fl.bd <- readRDS("birth_defects_variables_20230216.rds") %>% 
  filter(state == 'FL')

#' Join cancer and vital records data.
fl <- fl.vit %>% 
  left_join(fl.can, by = 'studyid') %>% 
  mutate(cancer = ifelse(is.na(cancer), 0, cancer),
         across(all:num.diagnoses, ~ ifelse(cancer == 0, 0, .x)))

#' Join birth defects and vital records data.
fl <- fl %>% 
  left_join(select(fl.bd, -c(state, coding.system)), by = 'studyid') %>% 
  mutate(birth.defect = ifelse(is.na(birth.defect), 0, birth.defect),
         runif = runif(nrow(.),0,1)) %>% 
  mutate(across(cns.anomaly:major.musculoskeletal.anomaly, ~ ifelse(is.na(.x) & birth.defect == 0, 0, .x)))

saveRDS(fl, 'FL_GOBACK_data_20230216.rds')

rm(list = ls()); gc()

# Append FL data to GOBACK ------------------------------------------------

setwd('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/')

goback <- readRDS('goback_20230105.rds')

fl <- readRDS('Expanded_datasets/FL_GOBACK_data_20230216.rds')

goback <- bind_rows(goback, fl)

saveRDS(goback, 
        '//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/goback_20230216.rds')
