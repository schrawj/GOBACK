
require(survival); require(tidyverse); require(magrittr); require(tictoc)

setwd('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/')

source('R scripts/GOBACK-R-Scripts/Analysis - user-defined functions.R')

goback <- readRDS('Datasets/goback.v20210616.rds')

#' Buyer beware: this line might take a day or more to execute.
goback.results <- generate.models(goback = goback, bd.index = 67:143, cancer.index = 20:59)

saveRDS(goback.results, 'Datasets/Expanded datasets/goback.coxph.results.v20210616.rds')

rm(goback, goback.results); gc()

goback.non.syn <- readRDS('Datasets/goback.nochrom.v20210616.rds')

goback.non.syn.results <- generate.models(goback.non.syn, c(67:128,141:143), cancer.index = 20:59)

saveRDS(goback.non.syn.results, 'Datasets/Expanded datasets/goback.non.syndromic.coxph.results.v20210616.rds')

rm(goback.non.syn, goback.non.syn.results); gc()

# Bind results ------------------------------------------------------------

overall <- readRDS('Datasets/Expanded datasets/goback.coxph.results.v20210616.rds') %>% 
  mutate(set = 'all.defects')

non.syndomic <- readRDS('Datasets/Expanded datasets/goback.non.syndromic.coxph.results.v20210616.rds') %>% 
  mutate(set = 'non.syndromic.defects')

combined <- bind_rows(overall, non.syndomic) %>% 
  mutate(q.value.coef = p.adjust(p.value.coef, method = 'fdr'),
         fdr.significant = q.value.coef < 0.05) %>% 
  select(defect, cancer, set, num.comorbid, HR:p.value.coef, q.value.coef, fdr.significant)

write_csv(combined, 'R outputs/goback.cox.models.v20210621.csv')
