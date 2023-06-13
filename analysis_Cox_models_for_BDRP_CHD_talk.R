require(tidyverse); require(survival); require(survminer); require(ggsci)

#' Load current GOBACK file.
goback <- readRDS("//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/goback_20230523.rds")

#' Organ system-level major defects variables, apart from the heart.
#' Will be used to identify children with CHD complicated by the co-occurrence of other major anomalies.
other.systems <- goback %>% 
  select(-c(major.heart.circulatory.anomaly, number.major.defects)) %>% 
  names() %>% 
  str_subset('major')

#' Remove children with genetic or chromosomal syndromes or unknown sex, right censor at 18, and identify more complex cases.
goback <- goback %>% 
  filter(genetic.anomaly == 0 | is.na(genetic.anomaly), sex != 'Unknown') %>% 
  select(birth.defect, cancer, all:gct.any, person.yrs, sex, birth.wt, gest.age, state, m.age, major.heart.circulatory.anomaly, all_of(other.systems), number.major.defects) %>% 
  mutate(cancer = ifelse(cancer == 1 & person.yrs > 18, 0, cancer),
         person.yrs = ifelse(person.yrs > 18, 18, person.yrs),
         non.chd.defects = rowSums(.[, other.systems], na.rm = T),
         chd.plus = factor(ifelse(major.heart.circulatory.anomaly == 1 & non.chd.defects >= 1, 2,
                           ifelse(major.heart.circulatory.anomaly == 1 & non.chd.defects == 0, 1, 
                           ifelse(birth.defect == 0, 0, NA))),
                           labels = c('No birth defect', 'CHD only', 'CHD plus')))

#' Run Cox model for all cancers combined using CHD plus variable and covariates.
cox <- coxph(Surv(person.yrs, cancer) ~ chd.plus, data = goback) %>% 
  summary()
cox

cox.adj <- coxph(Surv(person.yrs, cancer) ~ chd.plus + sex + m.age + birth.wt + gest.age + state, data = goback) %>% 
  summary() 
cox.adj

table(goback$chd.plus, goback$cancer, useNA = 'ifany')

#' Cox model for any major CHD.
cox.adj <- coxph(Surv(person.yrs, cancer) ~ major.heart.circulatory.anomaly + sex + m.age + birth.wt + gest.age + state, data = goback) %>% 
  summary() 
cox.adj

table(goback$major.heart.circulatory.anomaly, goback$cancer, useNA = 'ifany')

# Identify specific cancers with at least 10 CHD+ cases -------------------

cancers <- goback %>% select(all:gct.any) %>% names()

counts <- data.frame()

for (i in cancers){
  
  var.name <- i
  
  data <- goback %>% 
    select(all_of(i), chd.plus)
  
  count <- data %>% 
    count(.data[[var.name]], chd.plus) %>% 
    filter(.data[[var.name]] == 1, !is.na(chd.plus)) %>% 
    mutate('{var.name}' := i)
  
  names(count) <- c('cancer','chd.plus','n')
  
  counts <- rbind(counts, count)
  
}

tmp <- counts %>% filter(chd.plus != 'No birth defect')
tmp <- filter(tmp, n >= 10, duplicated(cancer))

cancers <- tmp %>% pull(cancer)

# Run Cox models for those cancers ----------------------------------------

estimates <- data.frame()

#' Estimates for my CHD alone versus CHD with extracardiac defect versus no birth defect variable.
for (i in 1:length(cancers)){
  
  print(paste0('Computing adjusted model for risk of ', cancers[i], '.'))
  
  goback.surv <- goback %>% 
    select(person.yrs, 
           all_of(cancers[i]), 
           chd.plus, sex, m.age, birth.wt, gest.age, state) 
  
  names(goback.surv) <- c('time',
                          'cancer',
                          'chd.plus', 'sex', 'm.age', 'birth.wt', 'gest.age', 'state')
  
  cox <- coxph(Surv(time, cancer) ~ chd.plus + sex + m.age + birth.wt + state, data = goback.surv) %>% 
    summary() %>% 
    coefficients()
  
  tmp <- cox[1:2,]
  
  new.estimate <- data.frame(cancer = cancers[i], 
                             
                             chd.only.hr = cox[1,2], 
                             chd.only.ci.lower = exp( cox[1,1]-(1.96*cox[1,3]) ), 
                             chd.only.ci.upper = exp( cox[1,1]+(1.96*cox[1,3]) ),
                             chd.only.p.value.coef = cox[1,5],
                             
                             chd.plus.hr = cox[2,2], 
                             chd.plus.ci.lower = exp( cox[2,1]-(1.96*cox[1,3]) ), 
                             chd.plus.ci.upper = exp( cox[2,1]+(1.96*cox[1,3]) ),
                             chd.plus.p.value.coef = cox[2,5])
  
  estimates <- rbind(estimates, new.estimate)
  
  saveRDS(estimates,
          '//smb-main.ad.bcm.edu/genepi/GOBACK/GOBACK_association_analyses/3_results/goback_2_cox_results_isolated_versus_nonisolated_CHD_20230530.rds')

}

estimates %>% 
  mutate(across(chd.only.hr:chd.plus.ci.upper, ~ round(.x, 2))) %>% 
  write_csv('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/R outputs/goback_2_cox_results_isolated_versus_nonisolated_CHD_20230530.csv')

#' Estimates for any major CHD versus no birth defect.
estimates <- data.frame()

for (i in 1:length(cancers)){
  
  print(paste0('Computing adjusted model for risk of ', cancers[i], '.'))
  
  goback.surv <- goback %>% 
    select(person.yrs, 
           all_of(cancers[i]), 
           major.heart.circulatory.anomaly, sex, m.age, birth.wt, gest.age, state) 
  
  names(goback.surv) <- c('time',
                          'cancer',
                          'major.heart.circulatory.anomaly', 'sex', 'm.age', 'birth.wt', 'gest.age', 'state')
  
  cox <- coxph(Surv(time, cancer) ~ major.heart.circulatory.anomaly + sex + m.age + birth.wt + state, data = goback.surv) %>% 
    summary() %>% 
    coefficients()
  
  tab <- table(goback.surv$major.heart.circulatory.anomaly, goback.surv$cancer)
  
  new.estimate <- data.frame(cancer = cancers[i], 
                             
                             n = tab[2,2], # number of co-occurring cases
                             
                             hr = cox[1,2], 
                             ci.lower = exp( cox[1,1]-(1.96*cox[1,3]) ), 
                             ci.upper = exp( cox[1,1]+(1.96*cox[1,3]) ),
                             p.value.coef = cox[1,5])
  
  estimates <- rbind(estimates, new.estimate)
  
  saveRDS(estimates,
          '//smb-main.ad.bcm.edu/genepi/GOBACK/GOBACK_association_analyses/3_results/goback_2_cox_results_CHD_kids_over_five_20230531.rds')
  
}

estimates %>% mutate(across(hr:ci.upper, ~ round(.x, 2)))

# Cox models in children 5 or older ---------------------------------------

#' Remove children with cancers diagnosed at less than 5 years.
tmp <- goback %>%
  select(cancer, all:gct.any, person.yrs, 
         major.heart.circulatory.anomaly, chd.plus, number.major.defects,
         sex, m.age, birth.wt, gest.age, state) %>%
  filter(cancer == 0 | (cancer == 1 & person.yrs >= 5))

#' Any cancer.
cox.adj <- coxph(Surv(person.yrs, cancer) ~ major.heart.circulatory.anomaly + sex + m.age + birth.wt + gest.age + state, data = tmp) %>% 
  summary() 
cox.adj

#' Generate an adjusted HR for each specific cancer where there were at least 10 co-occurring cases in the entire dataset.
estimates <- data.frame()

for (i in 1:length(cancers)){
  
  print(paste0('Computing adjusted model for risk of ', cancers[i], '.'))
  
  goback.surv <- tmp %>% 
    select(person.yrs, 
           all_of(cancers[i]), 
           major.heart.circulatory.anomaly, sex, m.age, birth.wt, gest.age, state) 
  
  names(goback.surv) <- c('time',
                          'cancer',
                          'major.heart.circulatory.anomaly', 'sex', 'm.age', 'birth.wt', 'gest.age', 'state')
  
  cox <- coxph(Surv(time, cancer) ~ major.heart.circulatory.anomaly + sex + m.age + birth.wt + state, data = goback.surv) %>% 
    summary() %>% 
    coefficients()
  
  tab <- table(goback.surv$major.heart.circulatory.anomaly, goback.surv$cancer)
  
  new.estimate <- data.frame(cancer = cancers[i], 
                             
                             n = tab[2,2], # number of co-occurring cases
                             
                             hr = cox[1,2], 
                             ci.lower = exp( cox[1,1]-(1.96*cox[1,3]) ), 
                             ci.upper = exp( cox[1,1]+(1.96*cox[1,3]) ),
                             p.value.coef = cox[1,5])
  
  estimates <- rbind(estimates, new.estimate)
  
  saveRDS(estimates,
          '//smb-main.ad.bcm.edu/genepi/GOBACK/GOBACK_association_analyses/3_results/goback_2_cox_results_CHD_kids_over_five_20230531.rds')
  
}

estimates %>% mutate(across(hr:ci.upper, ~ round(.x, 2)))
