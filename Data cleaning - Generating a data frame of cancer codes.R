require(tidyverse); require(haven); require(magrittr)

# Load in state-level data ------------------------------------------------

#' An intermediate version the AR data: has the computed studyid along with cancer variables.
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/Arkansas/arkansas.v20170828.2.rdata")

load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/North Carolina/nc.cancer.data.rdata")
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/Michigan/mi.raw.data.v20170818.rdata")

tx.can <- as.data.frame(read_dta(file = '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Birth Defects-Childhood Cancer Projects/GOBACK Project (Tiffany)/For Jeremy/Datasets with SEER Summary Stages/tcr_17_007_july20 (STATA).dta'))

# Texas -------------------------------------------------------------------

tx.can$birthID <- paste0('tx',as.character(tx.can$birthid))

tx.can <- arrange(select(tx.can, birthID, histtypeicdo3, psite, behavioricdo3, lateral, seersumstgbest), 
                  birthID)

#' Strip leading "C" from site codes, then convert all non-ID variables to numeric.
tx.can$psite <- substr(tx.can$psite, 2, 4)

for (i in 2:ncol(tx.can)){
  
  tx.can[,i] <- as.numeric(tx.can[,i])
  
}

#' Need to convert from long to wide.
#' Create a vector equal to the number of rows in the cancer data frame. Each entry represents the ith time that particular ID has occurred.
out <- c()

for (i in unique(tx.can$birthID)){
  
  stop <- sum(tx.can$birthID == i)
  
  out <- c(out, c(1:stop))
  
}

#' Pass the new vector to the reshape function in the form of the count column. This will tell R how many columns to create for each measured
#' variable when converting from long to wide.
tx.can$count <- out
tx.can <- reshape(tx.can, idvar = 'birthID', timevar = 'count', direction = 'wide')

tx.can[,(ncol(tx.can)+1):(ncol(tx.can)+5)] <- as.numeric(NA)

names(tx.can)[22:26] <- c('histtypeicdo3.5','psite.5','behavioricdo3.5','lateral.5','seersumstgbest.5')

names(tx.can) <- str_replace(names(tx.can), 'seersumstgbest','seer.stage2000')

tx.can <- rename(tx.can, studyid = birthID)

# North Carolina ----------------------------------------------------------

nc.cancer <- filter(nc.cancer, cancer == 1)

nc.trim <- select(nc.cancer, ncid, 
                  histologic.type.icdo3.1:histologic.type.icdo3.5, 
                  primary.site.1:primary.site.5, 
                  behavior.code.icdo3.1:behavior.code.icdo3.5,
                  seer.summary.stage.2000.1:seer.summary.stage.2000.5)

nc.trim$ncid <- paste0('nc',nc.trim$ncid)
nc.trim <- rename(nc.trim, studyid = ncid)

#' Removes the leading "C" from site codes.
for (i in 7:11){
  nc.trim[,i] <- substr(nc.trim[,i], 2, 4)
}

#' Remove empty strings and convert character variables to numeric.
for (i in 2:21){
  nc.trim[,i] <- ifelse(nc.trim[,i] == "", NA, nc.trim[,i])
  nc.trim[,i] <- as.numeric(nc.trim[,i])
}

#' Add 5 empty variables for laterality, which I don't think we got from NC.
first <- ncol(nc.trim) + 1
last <- ncol(nc.trim) + 5

for (i in first:last){
  nc.trim[,i] <- as.numeric(NA)
}

names(nc.trim)[first:last] <- paste0(rep('lateral.',5),1:5)

#' Vectors of current and desired column names, along with code to do the replacements and merge the data.
strings <- c('histologic.type.icdo3','primary.site','behavior.code.icdo3','seer.summary.stage.')

replacements <- c('histtypeicdo3', 'psite','behavioricdo3', 'seer.stage')

for (i in seq_along(strings)){
    
    names(nc.trim) <- stringr::str_replace(names(nc.trim), strings[i], replacements[i])
    
}

tmp <- bind_rows(tx.can, nc.trim)

# Michigan ----------------------------------------------------------------

mi.trim <- select(filter(mi, any_cancer == 1), 
                  studyid, morph31:morph34,site_code1:site_code4, behave31:behave34, 
                  laterality1:laterality4, sumstage21:sumstage24)

mi.trim$studyid <- paste0('mi',mi.trim$studyid)

#' Add 5 empty variables for 5th set of cancer variables. MI only provided info on up to four diagnoses.
first <- ncol(mi.trim) + 1
last <- ncol(mi.trim) + 5

mi.trim[,first:last] <- as.numeric(NA)

names(mi.trim)[first:last] <- c('morph35','site_code5','behave35','laterality5','sumstage25')

strings <- c('morph3','site_code','behave3','sumstage2','laterality')

replacements <- c('histtypeicdo3.', 'psite.','behavioricdo3.', 'seer.stage2000.', 'lateral.')

#' Vectors of current and desired column names, along with code to do the replacements and merge the data.
for (i in seq_along(strings)){
  
  names(mi.trim) <- stringr::str_replace(names(mi.trim), strings[i], replacements[i])
  
}

tmp <- bind_rows(tmp, mi.trim)

# Arkansas ----------------------------------------------------------------

ar <- filter(ar, cancer == 1)

#' Actually, these variables are missing at almost 100% in the AR data.
ar.trim <- select(rename(ar, psite.1 = primarysit, histtypeicdo3.1 = histology3, behavioricdo3.1 = behavior3, 
                             lateral.1 = laterality, seer.stage2000.1 = grade),
                  studyid, histtypeicdo3.1, psite.1, behavioricdo3.1, seer.stage2000.1, lateral.1)

ar.trim$psite.1 <- substr(ar.trim$psite.1, 2, 4)

for(i in 2:6){
  ar.trim[,i] <- as.numeric(ar.trim[,i])
}

#' Add 20 empty variables for the missing cancer variables. AR only provided this info for the first cancer diagnosis, and only in a few patients.
first <- ncol(ar.trim) + 1
last <- ncol(ar.trim) + 20

ar.trim[, first:last] <- as.numeric(NA)

names(ar.trim)[first:last] <- paste0(rep(c('histtypeicdo3.', 'psite.','behavioricdo3.', 'seer.stage2000.', 'lateral.'), 4), rep(2:5, each = 5))

cancer.codes <- bind_rows(tmp, ar.trim)

save(cancer.codes, file = 'W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/cancer.codes.v20190625.1.rdata')

# Oklahoma ----------------------------------------------------------------

#' OK cancer registry data.
load('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/Oklahoma/oklahoma.raw.data.rdata')

#' Merged AR, MI, NC, TX cancer registry data has already been generated.
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/cancer.codes.v20190625.1.rdata")

ok <- filter(ok, cancer == 1)

#' Remove leading C's from site variables
site.vars <- c(paste0('site', 1:7))

for (i in site.vars){
  
  ok[, i] <- str_remove(ok[, i], 'C')
  
}

#' OK histology codes have 5 digits. The fifth is the same as behavior.
hist.vars <- c(paste0('hist',1:7))

for (i in hist.vars){
  
  ok[ , i] <- as.numeric(substr(as.character(ok[, i]), 1, 4))
  
}

#' By NAACR conventions, malignant primaries are recorded in cols 1-3, non-malignant neoplasms in cols 4-7.
#' In practice no child has more than 5 diagnoses. Therefore we do not necessarily need to carry forward all 7 columns.
histology.codes <- as.data.frame(matrix(nrow = 0, ncol = 6))
names(histology.codes) <- c('studyid',paste0('histtypeicdo3.',1:5))

primary.site.codes <- as.data.frame(matrix(nrow = 0, ncol = 6))
names(primary.site.codes) <- c('studyid',paste0('psite.',1:5))

behavior.codes <- as.data.frame(matrix(nrow = 0, ncol = 6))
names(behavior.codes) <- c('studyid',paste0('behavioricdo3.',1:5))

laterality.codes <- as.data.frame(matrix(nrow = 0, ncol = 6))
names(laterality.codes) <- c('studyid',paste0('lateral.',1:5))

#' Other states didn't provide data in this same format. 
#' Concatenate codes from all cancers in each child; we will use behavior codes to separate the unclassifiable ones later.
for (i in 1:nrow(ok)){

  site.codes <- as.numeric(ok[i, c(77, 2:8)])
  site.codes <- subset(site.codes, !is.na(site.codes))
  site.codes <- c(site.codes, rep(NA, times = 6 - length(site.codes)))
  
  primary.site.codes[i,] <- site.codes
  
  hist.codes <- as.numeric(ok[i, c(77,9:15)])
  hist.codes <- subset(hist.codes, !is.na(hist.codes))
  hist.codes <- c(hist.codes, rep(NA, times = 6 - length(hist.codes)))
  
  histology.codes[i,] <- hist.codes
  
  lat.codes <- as.numeric(ok[i, c(77, 16:22)])
  lat.codes <- subset(lat.codes, !is.na(lat.codes))
  lat.codes <- c(lat.codes, rep(NA, times = 6 - length(lat.codes)))
  
  laterality.codes[i,] <- lat.codes
  
  beh.codes <- as.numeric(ok[i, c(77, 30:36)])
  beh.codes <- subset(beh.codes, !is.na(beh.codes))
  beh.codes <- c(beh.codes, rep(NA, times = 6 - length(beh.codes)))
  
  behavior.codes[i,] <- beh.codes
  
}

#' SEER summary stage was not provided.
seer.stage <- as.data.frame(matrix(nrow = 971, ncol = 6))
names(seer.stage) <- c('studyid',paste0('seer.stage2000.',1:5))
seer.stage$studyid <- c(unique(ok$randID))

for (i in 2:ncol(seer.stage)){
  
  seer.stage[,i] <- as.numeric(seer.stage[,i])
  
}

ok.cancer <- left_join(histology.codes, primary.site.codes, by = 'studyid')
ok.cancer <- left_join(ok.cancer, laterality.codes, by = 'studyid')
ok.cancer <- left_join(ok.cancer, behavior.codes, by = 'studyid')
ok.cancer <- left_join(ok.cancer, seer.stage, by = 'studyid')

ok.cancer$studyid <- paste0('ok', ok.cancer$studyid)

cancer.codes <- bind_rows(cancer.codes, ok.cancer)

save(cancer.codes,
     file = '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/cancer.codes.v20191008.1.rdata')

rm(list = ls()); gc()

# Michigan 2022 updates ---------------------------------------------------

require(tidyverse); require(magrittr)

load('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded_datasets/cancer.codes.v20191008.1.rdata')

#' Remove old MI data.
cancer.codes %<>% filter(substr(studyid, 1, 2) != 'mi')

#' Load updated MI data.
mi.cancer <- readRDS('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded_datasets/mi.cancer.codes.v20220516.rds') %>% 
  mutate(STUDYID = paste0('mi', STUDYID)) %>% 
  select(-Cancer_Histology6, -Cancer_BehaviorICDO36, -Cancer_PrimarySite6, -Cancer_Laterality6, -Cancer_SEERSumStg20006)

#' Rename MI variables to match combined dataset.
cancer.names <- c('studyid',
                  paste0('histtypeicdo3.', 1:5),
                  paste0('behavioricdo3.', 1:5),
                  paste0('psite.', 1:5),
                  paste0('lateral.', 1:5),
                  paste0('seer.stage2000.', 1:5))

names(mi.cancer) <- cancer.names

#' Change MI site variables to numeric.
mi.cancer %<>%
  mutate(across(psite.1:psite.5, ~ str_remove(.x, 'C')  )) %>% 
  mutate(across(psite.1:psite.5, as.numeric)) #' Warning messages attributable to coercion of "" and "." to NA values.

cancer.codes %<>% bind_rows(mi.cancer)

saveRDS(cancer.codes,
        '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded_datasets/cancer.codes.v20220516.rds')

rm(list = ls()); gc()

# Texas 2022 update -------------------------------------------------------

cancer.codes <- readRDS('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded_datasets/cancer.codes.v20220516.rds')

names(cancer.codes) %<>% str_remove_all('\\Q.\\E')

#' Remove old TX data; change names of SEER stage variables to reflect that we are no longer just using 2000 staging.
cancer.codes %<>% 
  filter(!str_detect(studyid, '^tx')) %>% 
  rename_with(~paste0('seerstagebest',substr(.x, 14, 14)), .cols = starts_with('seer'))

#' Load new Texas data.
tx <- read_dta('//smb-main.ad.bcm.edu/genepi3/TiffanyChambers/GOBACK Data/Texas/Texas 2021 Data/Stata/TCR Data.dta')

names(tx) %<>% tolower() %>% str_replace_all('_','.')

#' Add placeholder columns for fifth primary tumor, which was recorded in a few kids in the merged file but none in the new Texas file.
new.cols <- paste0(c('histtype','psite','behavior','lateral','seersumstgbest'),5)

for (i in new.cols){ tx %<>% mutate('{ i }' := as.numeric(NA)) }

#' Select relevant columns; add state prefix to study ID; remove leading C from site codes and convert to numeric; rename variables.
tx %<>% 
  select(studyid, starts_with(c('histtype','psite','behavior','lateral','seersum'))) %>% 
  mutate(studyid = paste0('tx',studyid),
         across(psite1:psite4, ~ as.numeric(str_remove(.x, '^C')))) %>% 
  rename_with(~paste0('seerstagebest',substr(.x,15,15)), .cols = starts_with('seer'))

#' Append new TX data.
cancer.codes %<>% bind_rows(tx)

saveRDS(cancer.codes,
        '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded_datasets/cancer.codes.v20220601.rds')

rm(list = ls()); gc()

# Massachusetts -----------------------------------------------------------

merged <- readRDS('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded_datasets/cancer.codes.v20220601.rds') %>% 
  select(-histtype5, -behavior5) # empty columns that hitched a ride in this data frame somehow.

ma <- read_dta('//smb-main.ad.bcm.edu/genepi3/TiffanyChambers/GOBACK Data/Massachusetts Data/Data (June 2021)/All cancer_092421.dta')

names(ma) %<>% tolower() %>% str_replace_all('_', '.')

ma <- ma %>% 
  select(kiduid, starts_with(c('hist', 'primary.site', 'behav', 'lateral', 'seer'))) %>% 
  mutate(studyid = paste0('ma',kiduid),
         across(starts_with('seer'), ~as.numeric(.x))) %>% 
  select(-kiduid) %>% 
  rename_with(~paste0('histtypeicdo3', substr(.x,11,11)), .cols = starts_with('hist.')) %>% 
  rename_with(~paste0('psite', substr(.x,13,13)), .cols = starts_with('primary')) %>% 
  rename_with(~paste0('behavioricdo3',1:3), .cols = starts_with('behav')) %>% 
  rename_with(~paste0('lateral',1:3), .cols = starts_with('lateral')) %>% 
  rename_with(~paste0('seerstagebest',1:3), .cols = starts_with('seer')) 

#' Add placeholder columns for cancers 4 and 5.
new.vars <- ma  %>% 
  select(-studyid) %>% 
  names() %>% 
  str_remove_all('\\d{1}$') # remove exactly 1 number from the end of the string. Ensures the 3 in icdo3 isn't dropped. 

new.vars <- c(paste0(unique(new.vars),4), 
              paste0(unique(new.vars), 5))

for (i in new.vars){ ma %<>% mutate('{ i }' := as.numeric(NA)) }

ma %<>% select(names(merged))

merged %<>% bind_rows(ma)

saveRDS(merged,
        '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded_datasets/cancer.codes.v20220602.rds')

rm(list = ls()); gc()

# Simplify naming conventions ---------------------------------------------

setwd('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/')

#' Existing cancer codes data frame.
merged <- readRDS('Expanded_datasets/cancer.codes.v20220602.rds')

new.names <- c('studyid', rep(c('histology.', 'site.', 'behavior.', 'laterality.', 'seer.stage.'), 5))

new.names <- c(new.names[1], paste0(new.names[2:26], rep(1:5, each = 5)) )

names(merged) <- new.names

saveRDS(merged, 'Expanded_datasets/cancer.codes.v20220810.rds')

rm(list = ls()); gc()

# South Carolina ----------------------------------------------------------

#' Existing cancer codes data frame.
merged <- readRDS('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded_datasets/cancer.codes.v20220810.rds')

#' Raw SC data.
sc <- readRDS('//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/raw_registry_data/south_carolina_raw_data_20220810.rds')

names(sc) %<>% tolower() %>% str_replace_all('_', '.')

sc.can <- sc %>% 
  filter(naaccr523.behavior != '') %>% 
  mutate(studyid = paste0('sc',studyid)) %>% 
  select(studyid, naaccr522.histology, naaccr400.primsite, naaccr523.behavior, naaccr410.laterality, ss.stage) %>% 
  rename(histology.1 = naaccr522.histology,
         site.1 = naaccr400.primsite,
         behavior.1 = naaccr523.behavior,
         laterality.1 = naaccr410.laterality,
         seer.stage.1 = ss.stage) %>% 
  mutate(site.1 = as.numeric(str_remove(.$site.1, '^[:alpha:]')))

min.col <- ncol(sc.can) + 1
max.col <- ncol(merged)

new.cols <- length(seq(ncol(sc.can) + 1, ncol(merged)))
new.col.names <- names(merged)[min.col:max.col]

for (i in 1:new.cols){
  
  new.col <- new.col.names[i]
  
  sc.can %<>% mutate('{ new.col }' := as.numeric(NA))
  
}

merged %<>% bind_rows(sc.can) %>% as_tibble()

saveRDS(merged, '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded_datasets/cancer.codes.v20220811.rds')

rm(list = ls()); gc()

# New Jersey --------------------------------------------------------------

setwd('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/')

cancer.codes <- readRDS('Expanded_datasets/cancer.codes.v20220811.rds')

nj <- readRDS('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/GOBACK/Datasets/nj.cancer.raw.data.rds')

nj.clean <- nj %>% 
  rename_with(~paste0('histology.',1:5), starts_with('hist')) %>% 
  rename_with(~paste0('site.',1:5), starts_with('primary_site')) %>% 
  rename_with(~paste0('behavior.',1:5), starts_with('behavior')) %>% 
  rename_with(~paste0('laterality.',1:5), starts_with('laterality')) %>% 
  rename_with(~paste0('seer.stage.',1:5), starts_with('seer_summ')) %>% 
  select(all_of(names(cancer.codes))) %>% 
  mutate(across(starts_with('site'), ~ as.numeric(str_remove(.x, '^C'))))

cancer.codes <- cancer.codes %>% bind_rows(nj.clean)

saveRDS(cancer.codes, 'Expanded_datasets/cancer.codes.v20230104.rds')

rm(list = ls()); gc()

# Florida -----------------------------------------------------------------

setwd('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded_datasets/')

#' Current merged cancer codes data frame.
codes <- readRDS('cancer.codes.v20230104.rds')

#' These data have already been cleaned and formatted. See 'data_cleaning_clean_and_prep_FL_data.R for details.
fl <- readRDS('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/GOBACK/Datasets/fl.cancer.raw.data.rds')

#' Combine and save.
codes <- codes %>% 
  bind_rows(fl)

saveRDS(codes, 'cancer_codes_20230216.rds')