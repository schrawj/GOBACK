
# README and TODO comments ------------------------------------------------

#' Questions for Tiffany:
#' 1. 12/31/2018 was the censoring date for individuals without cancer. Is ascertainment of vital status by the cancer registry complete up through 2018?
#' 2. How does your computed race/ethnicity variable work for kids born before 2007.

#' TODO:
#' 1. Compute runif.
#' 2. Check how race-ethnicity was handled for MI kids in the existing GOBACK file; update m_race and f_race variables in new data as necessary.
#' 3. Change state_num to 3.
#' 4. Compute cancer variables.
#' 5. Harmonize the structure of the ID variable.

require(tidyverse); require(haven); require(magrittr); require(icd)

mi <- read_dta('//smb-main.ad.bcm.edu/genepi3/TiffanyChambers/GOBACK Data/Michigan Data 2021/Updated Data (June 2021)/Stata/all years (harmonized).dta')

saveRDS(mi, '//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/raw_registry_data/michigan.2021.raw.data.rds')

rm(list = ls()); gc()

# Read in and prepare raw data --------------------------------------------

#' Load the raw MI data, compute race/ethnicity, and standardize the formatting of the column names.
mi <- readRDS('//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/raw_registry_data/michigan.2021.raw.data.rds') 

names(mi) %<>% str_replace_all('_','.') %>% tolower()

#' Compute month of birth.
mi %<>% mutate(birth.mo = as.numeric(substr(dob, 6, 7))) %>% 
  filter(birth.yr > 1991) # Prior to this year, birth defects prevalence is suspiciously low, indicating possible incomplete ascertainment. Remove these years due to concern for bias towards the null caused by non-differential misclassification.
  select(studyid, m.race ,birth.yr, birth.mo, bdefect.icd9cod1:bdefect.icd10cod30)

# Clean up ICD10 codes in newer data --------------------------------------

#' Birth defects were recorded using ICD 10 codes beginning in October 2015.
mi.recent <- filter(mi, birth.yr > 2015 | 
                       (birth.yr == 2015 & birth.mo >= 10) ) %>% 
  select(studyid:birth.mo, contains('icd10'))

#' A function that removes ICD10 codes outside the range for evaluated birth defects (Q00-Q99).
is.icd10.defect <- function(x) { ifelse(str_detect(x, '^Q'), x, '') }

#' A vector of positions for the ICD 10 columns.
icd10.cols <- 5:34

#' Remove irrelevant codes.
mi.recent %<>% mutate(across(c(icd10.cols), is.icd10.defect))

#' Add a decimal to each code at the appropriate place.
mi.recent[ , icd10.cols] %<>% apply(2, function(x){ ifelse(nchar(x) > 3, paste0(substr(x, 1, 3), '.', substr(x, 4, nchar(x))), x) })

#' Remove ICD9 codes from new data.
mi.recent <- mi.recent %>% 
  filter(if_any(starts_with('bdefect'), ~ .x != ""))

saveRDS(mi.recent, 
        file = '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old_Datasets/Michigan/2021 update/michigan.icd10.birth.defect.codes.v20230405.rds')

# Clean up ICD9 codes in older data ---------------------------------------

#' Prior to that, they were recorded using ICD9 codes.
mi.old <- filter(mi, birth.yr < 2015 | 
                   (birth.yr == 2015 & birth.mo < 10) ) %>% 
  select(studyid:birth.mo, contains('icd9'))

#' A function that removes ICD9 codes outside the range for evaluated birth defects (216, 228.0, 237.7, 740-759).
is.icd9.defect <- function(x) { ifelse(str_detect(x, valid.ranges), x, '') }

#' Vector of positions for ICD9 columns.
icd9.cols <- 5:34

#' Valid birth defect ICD9 code prefixes.
valid.ranges <- str_c(c('^216', '^2280', '^2377', '^7[45]'), collapse = '|')

table(mi.old$bdefect.icd9cod1, useNA = 'ifany')

mi.old[,icd9.cols] %<>% apply(2, function(x) { str_remove_all(x, '\\Q.\\E') })

mi.old %<>% mutate(across(c(icd9.cols), is.icd9.defect)) %>% 
  mutate(across(all_of(icd9.cols), function (x) { ifelse(nchar(x) > 3, paste0(substr(x, 1, 3), '.', substr(x, 4, nchar(x))), x) }))

#' Remove ICD10 codes from old data and remove kids without birth defects.
mi.old <- mi.old %>% 
  filter(if_any(starts_with('bdefect'), ~ .x != ""))

saveRDS(mi.old, 
        file = '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old_Datasets/Michigan/2021 update/michigan.icd9.birth.defect.codes.v20230405.rds')