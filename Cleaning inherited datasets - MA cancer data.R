#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2021.08.25.
#' 
#' Assign cancer diagnoses to children in Massachusetts.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(haven); require(tidyverse); require(magrittr)

compute.cancer.vars <- function(input.data, output.data, iccc.extended, new.var){
  
  input.data$iccc.extended.classification <- as.numeric(input.data$iccc.extended.classification)
  
  ids <- filter(input.data, iccc.extended.classification %in% iccc.extended) %>% select(KIDUID) %>% unlist()
  
  output.data[ , new.var] <- 0
    
  output.data[ , new.var] <- ifelse(output.data$KIDUID %in% ids, 1, 0)
  
  return(output.data)
  
}

compute.rms.vars <- function(input.data, output.data){
  
  input.data$iccc.extended.classification <- as.numeric(input.data$iccc.extended.classification)
  
  erms.ids <- filter(input.data, iccc.extended.classification == 55 & Hist_ICDO3 == 8910) %>% select(KIDUID) %>% unlist() %>% as.numeric()
  
  arms.ids <- filter(input.data, iccc.extended.classification == 55 & Hist_ICDO3 == 8920) %>% select(KIDUID) %>% unlist() %>% as.numeric()
  
  rms.other.ids <- setdiff(unlist(select(filter(input.data, iccc.extended.classification == 55), KIDUID)) , 
                           union(erms.ids, arms.ids))
  
  output.data %<>% mutate(erms = ifelse(KIDUID %in% erms.ids, 1, 0),
                          arms = ifelse(KIDUID %in% arms.ids, 1, 0),
                          rms.other = ifelse(KIDUID %in% rms.other.ids, 1, 0)
                          )
  return(output.data)
  
}

#' An ICCC-3 extended classification lookup table.
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/iccc.lookup.table.v20210826.rdata")

#' Massachusetts cancer registry data.
ma <- read_dta('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/BD-CC Projects/GOBACK/Massachusetts Data/Data (June 2021)/All cancer.dta') %>% 
  mutate(Primary_Site = ifelse(nchar(Primary_Site) == 2, paste0('0',Primary_Site), as.character(Primary_Site))) %>% 
  mutate(Primary_Site = paste0('C',Primary_Site))

#' Append ICCC-3 classifications; drop non-malignant tumors unless CNS or GCT.
ma  %<>% left_join(iccc.lookup, by = c('Hist_ICDO3' = 'icdo3.hist', 'Primary_Site' = 'icdo3.site.code')) %>% 
  filter(Behav_ICDO3 == 3 | (iccc.extended.classification %in% 
                                      c('017','018','019','020','021','022','023','024','025','026','027','028',
                                       '029','030','031','032','071','072','073','074','075','076') & Behav_ICDO3 %in% c(0,1,3)) )

#' Sort by year of diagnosis so tumors will be in chronological order.
ma  %<>% arrange(KIDUID, YearDx)

#' Add columns to hold histology, site, and behavior code for 2nd and 3rd malignant neoplasms.
start <- ncol(ma)+1
stop <- ncol(ma)+20

ma[,start:stop] <- as.numeric(NA)

names(ma)[start:stop] <- c(paste0(c('behavioricdo3.', 'psite.', 'histtypeicdo3.', 'lateral.', 'seer.stage2000.'),2), 
                           paste0(c('behavioricdo3.', 'psite.', 'histtypeicdo3.', 'lateral.', 'seer.stage2000.'),3),
                           paste0(c('behavioricdo3.', 'psite.', 'histtypeicdo3.', 'lateral.', 'seer.stage2000.'),4),
                           paste0(c('behavioricdo3.', 'psite.', 'histtypeicdo3.', 'lateral.', 'seer.stage2000.'),5))

#' Iterate over rows: 
#' If this is the second time an ID appears, update 2nd neoplasm codes in preceding rows with values from this row.
#' If this is the third time an ID appears, update 3rd neoplasm codes 2 rows above with corresponding values from this row.
for (i in 1:nrow(ma)){
  
  index.id <- ma$KIDUID[i]
  
  ids <- ma$KIDUID[1:i]
  ids %<>% subset(ids %in% index.id)
  
  if (length(ids) == 1){
    
    next()
    
  }
  
  else if (length(ids) == 2){
    
    ma$histtypeicdo3.2[i-1] <- ma$Hist_ICDO3[i]
    ma$psite.2[i-1] <- ma$Primary_Site[i]
    ma$behavioricdo3.2[i-1] <- ma$Behav_ICDO3[i]
    ma$lateral.2[i-1] <- ma$Laterality[i]
    
  }
  
  else if (length(ids) == 3){
    
    ma$histtypeicdo3.3[i-2] <- ma$Hist_ICDO3[i]
    ma$psite.3[i-2] <- ma$Primary_Site[i]
    ma$behavioricdo3.3[i-2] <- ma$Behav_ICDO3[i]
    ma$lateral.3[i-2] <- ma$Laterality[i] 
    
  }
  
  else {
    
    break()
    
  }
  
}

#' Create data frame housing computed cancer variables.
ma.cancer <- select(ma, KIDUID)

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

for (i in 1:length(codes)){
  
  ma.cancer <- compute.cancer.vars(ma, ma.cancer, codes[[i]], cancers[[i]])
  
}

ma.cancer <- compute.rms.vars(ma, ma.cancer) 

first.primaries <- ma %>% 
  filter(!duplicated(KIDUID)) %>% 
  rename(cancer1 = iccc.name,
         laterality1 = Laterality)

ma.cancer %<>% left_join(select(first.primaries, KIDUID, cancer1, laterality1), by = 'KIDUID') %>% 
  mutate(cancer = ifelse(rowSums(.[2:41]) >= 1, 1, 0),
         num.diagnoses = rowSums(.[c(2,3,5:7,9:13, 15:16, 18:20, 22:23, 25:27, 29,32:35,37:41)])) %>% 
  filter(num.diagnoses > 0)

#' Compute a person-years variable for kids with cancer.
ma  %<>% mutate(person.yrs = as.numeric(Age_Dx))

ma.cancer  %<>% left_join(select(ma, KIDUID, person.yrs), by = 'KIDUID')

saveRDS(ma.cancer, '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/ma.cancer.variables.v20210907.rds')

# Append MA cancer codes to the existing states ---------------------------

#' First, load in the existing cancer codes data frame.
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/cancer.codes.v20191008.1.rdata")

#' Remove duplicated rows.
ma  %<>% filter(!duplicated(KIDUID))

#' Then, rename and select relevant columns from 'ma' data frame. 
#' Convert site variables from character to numeric after removing leading C.
ma.cancer.codes <- ma %>% 
  mutate(studyid = paste0('ma', KIDUID)) %>% 
  rename(histtypeicdo3.1 = Hist_ICDO3, psite.1 = Primary_Site, behavioricdo3.1 = Behav_ICDO3, lateral.1 = Laterality, seer.stage2000.1 = SEER_Sum_Stg_2000) %>% 
  select(studyid, starts_with(c('histtype', 'behavior', 'psite','seer.stage2000','lateral'))) %>% 
  mutate(across(psite.1:psite.5, ~str_remove(.x, 'C'))) %>% 
  mutate(across(psite.1:psite.5, ~as.numeric(.x)))

#' Bind together.
cancer.codes %<>% bind_rows(ma.cancer.codes)

saveRDS(cancer.codes, 
        '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/cancer.codes.v20220210.rds')

# Save first primaries ----------------------------------------------------

first.primaries %<>% mutate(studyid = paste0('ma', KIDUID)) %>% 
  rename(histtypeicdo3.1 = Hist_ICDO3, psite.1 = Primary_Site, behavioricdo3.1 = Behav_ICDO3, lateral.1 = laterality1) %>% 
  select(all_of(selection.vars))

saveRDS(first.primaries,
        '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/ma.cancer.codes.v20210907.rds')

# Scratch paper -----------------------------------------------------------

load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/iccc.lookup.table.v20210826.rdata")

#' Massachusetts cancer registry data.
ma <- read_dta('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/BD-CC Projects/GOBACK/Massachusetts Data/Data (June 2021)/All cancer.dta') %>% 
  mutate(Primary_Site = ifelse(nchar(Primary_Site) == 2, paste0('0',Primary_Site), as.character(Primary_Site))) %>% 
  mutate(Primary_Site = paste0('C',Primary_Site))

#' Append ICCC-3 classifications; drop non-malignant tumors unless CNS or GCT.
ma  %<>% left_join(iccc.lookup, by = c('Hist_ICDO3' = 'icdo3.hist', 'Primary_Site' = 'icdo3.site.code')) %>% 
  filter(Behav_ICDO3 == 3 | (iccc.extended.classification %in% 
                               c('017','018','019','020','021','022','023','024','025','026','027','028',
                                 '029','030','031','032','071','072','073','074','075','076') & Behav_ICDO3 %in% c(0,1,3)) )

start <- ncol(ma)+1
stop <- ncol(ma)+12

ma[,start:stop] <- as.numeric(NA)

names(ma)[start:stop] <- c(paste0(c('behavioricdo3.', 'psite.', 'histtypeicdo3.'),2), 
                              paste0(c('behavioricdo3.', 'psite.', 'histtypeicdo3.'),3),
                              paste0(c('behavioricdo3.', 'psite.', 'histtypeicdo3.'),4),
                              paste0(c('behavioricdo3.', 'psite.', 'histtypeicdo3.'),5))

for (i in 1:nrow(ma)){
  
  index.id <- ma$KIDUID[i]

  ids <- ma$KIDUID[1:i]
  ids %<>% subset(ids %in% index.id)
  
  if (length(ids) == 1){
    
    next()
    
  }
  
  else if (length(ids) == 2){
    
    ma$histtypeicdo3.2[i-1] <- ma$Hist_ICDO3[i]
    ma$psite.2[i-1] <- ma$Primary_Site[i]
    ma$behavioricdo3.2[i-1] <- ma$Behav_ICDO3[i]
    
  }
  
  else if (length(ids) == 3){
  
    ma$histtypeicdo3.3[i-2] <- ma$Hist_ICDO3[i]
    ma$psite.3[i-2] <- ma$Primary_Site[i]
    ma$behavioricdo3.3[i-2] <- ma$Behav_ICDO3[i]
    
  }
  
  else {
    
    break()
    
  }
    
}

ma  %<>% filter(!duplicated(KIDUID))



ma.cancer.codes <- ma %>% 
  mutate(studyid = paste0('ma', KIDUID)) %>% 
  rename(histtypeicdo3.1 = Hist_ICDO3, psite.1 = Primary_Site, behavioricdo3.1 = Behav_ICDO3, lateral.1 = Laterality, seer.stage2000.1 = SEER_Sum_Stg_2000) %>% 
  select(studyid, starts_with(c('histtype', 'behavior', 'psite','seer.stage2000','lateral'))) %>% 
  mutate(across(psite.1:psite.5, ~str_remove(.x, 'C'))) %>% 
  mutate(across(psite.1:psite.5, ~as.numeric(.x)))

duplicated.ids <- ma[which(duplicated(ma$KIDUID)), ] %>% pull(KIDUID)

ma.long.cancer.codes <- ma %>% 
  filter(KIDUID %in% duplicated.ids) %>% 
  select(KIDUID, Primary_Site, Hist_ICDO3, Behav_ICDO3) %>% 
  mutate(studyid = paste0('ma', KIDUID)) 

start <- ncol(ma.long)+1
stop <- ncol(ma)+12

ma[,start:stop] <- as.numeric(NA)

names(ma)[start:stop] <- c(paste0(c('behavioricdo3.', 'psite.', 'histtypeicdo3.'),2), 
                           paste0(c('behavioricdo3.', 'psite.', 'histtypeicdo3.'),3),
                           paste0(c('behavioricdo3.', 'psite.', 'histtypeicdo3.'),4),
                           paste0(c('behavioricdo3.', 'psite.', 'histtypeicdo3.'),5))

for (i in 1:nrow(ma.long.cancer.codes)){
  
  index.id <- ma.long.cancer.codes$KIDUID[i]
  
  ids <- ma.long.cancer.codes$KIDUID[1:i]
  ids %<>% subset(ids %in% index.id)
  
  if (length(ids) == 1){
    
    next()
    
  }
  
  else if (length(ids) == 2){
    
    ma.long.cancer.codes$histtypeicdo3.2[i-1] <- ma.long.cancer.codes$Hist_ICDO3[i]
    ma.long.cancer.codes$psite.2[i-1] <- ma.long.cancer.codes$Primary_Site[i]
    ma.long.cancer.codes$behavioricdo3.2[i-1] <- ma.long.cancer.codes$Behav_ICDO3[i]
    
  }
  
  else if (length(ids) == 3){
    
    ma.long.cancer.codes$histtypeicdo3.3[i-2] <- ma.long.cancer.codes$Hist_ICDO3[i]
    ma.long.cancer.codes$psite.3[i-2] <- ma.long.cancer.codes$Primary_Site[i]
    ma.long.cancer.codes$behavioricdo3.3[i-2] <- ma.long.cancer.codes$Behav_ICDO3[i]
    
  }
  
  else {
    
    break()
    
  }
  
}


ma.wide.cancer.codes <- ma.long.cancer.codes %>% 
  pivot_wider()



%>% 
  rename(histtypeicdo3.1 = Hist_ICDO3, psite.1 = Primary_Site, behavioricdo3.1 = Behav_ICDO3, lateral.1 = Laterality, seer.stage2000.1 = SEER_Sum_Stg_2000) %>% 
  select(studyid, starts_with(c('histtype', 'behavior', 'psite','seer.stage2000','lateral'))) %>% 
  mutate(across(psite.1:psite.5, ~str_remove(.x, 'C'))) %>% 
  mutate(across(psite.1:psite.5, ~as.numeric(.x)))


