# TODO comments -----------------------------------------------------------

#' TODO: Clean covariates for Cox model: mother's age and child's sex.
#' TODO: Compute verified variables for any birth defect and any chromosomal defect.
#' TODO: COmpute time to event variable (last follow-up, death within first year, or cancer DX)
#' TODO: Run a Cox model comparing cancer risk in those with and without birth defects.
#' TODO: Run Cox models comparing cancer risk in those with and without birth defects,
#'       stratifying on whether they're chromosomal or not.

# Prep environment --------------------------------------------------------

require(tidyverse); require(magrittr); require(haven); require(tictoc)

setwd('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets')

# Load and save raw data --------------------------------------------------

wa <- read_dta('Old Datasets/Washington State/baylor2020_data_malf_20220111.dta')

saveRDS(wa, 'Old Datasets/Washington State/washington.raw.data.rds')

# Data cleaning -----------------------------------------------------------

#' Add periods to ICD9 codes.
recode <- function(x) (ifelse(nchar(x) >= 4, 
                              paste0(substr(x, 1, 3),'.',substr(x, 4, nchar(x))),
                              x))

#codes <- c('^7[45]', '^237.7')

wa <- readRDS('Old Datasets/Washington State/washington.state.raw.data.rds')

wa %<>% mutate(across(cbdiag1:cbdiag25, recode))

#' Remove codes that are neither structural birth defects (740-759) or neurofibromatosis (237.7).
wa %<>% mutate(across(cbdiag1:cbdiag25, function(x) (ifelse(str_detect(x, '^7[45]|^237.7'), x, "None" ))))

#' Create a list with as many elements as children. Each element's name corresponds to a study ID. 
#' Each element will be all the BD codes for that child.
#' Takes ~3 minutes.
master.list <- wa %>% 
  split(.$idno) %>% 
  map(~str_c(.[, 187:211]))

#' A single TRUE/FALSE vlaue indicating whether the child has any ICD9 code indicative of a BD.
#' FALSE indicates one or more BDs.
tic()
cases <- lapply(master.list, function(x) { all(str_detect(x, "None")) })
toc()

cases <- data.frame(ids = names(master.list),
                    cases = unlist(cases))

cases <- filter(cases, cases == F) %>% pull(ids) %>% as.character()

wa %<>% mutate(any.birth.defect = ifelse(idno %in% cases, 1, 0))

#' A single TRUE/FALSE value for every child indicating whether they have DS, and update the WA data. 
ds <- '758.0'

cases <- lapply(master.list, function(x){ any(str_detect(x, ds)) })

cases <- data.frame(ids = names(master.list),
                    ds = unlist(cases))

ids <- cases %>% filter(ds == T) %>% pull(ids) %>% as.character()

wa <- wa %>% mutate(ds = ifelse(idno %in% ids, 1, 0))

# DS-leukemia counts for Philip -------------------------------------------

#' Count number of children with ALL or AML, excluding children with DS (1,333).
#' Count number of children with DS-ALL or DS-AML (29 and 14, respectively).
#' Count number of children with DS and no cancer diagnosis (34).

table(wa$ds, useNA = 'ifany')
table(wa$iccc_ext, wa$ds, useNA = 'ifany')

# Scratch paper -----------------------------------------------------------
