require(tidyverse); require(survival)

setwd('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/GOBACK/Datasets/')

#' This is a data frame that includes birth-defect cancer combinations that have at least five co-occurring cases throughout the dataset.
counts <- readRDS('goback_co_occurring_case_counts_20230224.rds')

all.defects <- counts$defect %>% unique()

syndromes <- all.defects[64:75]

not.syndromes <- setdiff(all.defects, syndromes)

model.defects <- counts %>% 
  filter(n >= 5,
        (defect %in% syndromes & phenotype == 'all.bd.cases') | (defect %in% not.syndromes & phenotype == 'non.syndromic.cases')
         ) %>% 
  arrange(defect)

#' Keep only birth defects for which there were five or more cases with cancer.
not.syndromes <- subset(not.syndromes, not.syndromes %in% model.defects$defect)

#' Initialize a data frame to hold the results.
estimates <- data.frame()

#' Generates Cox models FOR NON-SYNDROMIC CASES and stores the results in a data frame.
for (j in not.syndromes){
  
  goback <- readRDS('//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/goback_20230216.rds') %>% 
    filter(genetic.anomaly != 1 | is.na(genetic.anomaly), !is.na(.data[[j]])) %>% 
    select(sex, m.age, state, birth.wt, person.yrs, all:other.tumor,cancer, all_of(j)) %>% 
    mutate(state = factor(state))
  
  gc()
  
  cancers <- model.defects %>% 
    filter(defect == j) %>% 
    pull(cancer)

  for (i in 1:length(cancers)){
    
    print(paste0('Computing adjusted model for risk of ', cancers[i], ' in kids with ', j,'.'))
    
    goback.surv <- goback %>% 
      select(person.yrs, 
             all_of(cancers[i]), 
             all_of(j), 
             sex, m.age, birth.wt, state) 
    
    names(goback.surv) <- c('time',
                            'cancer',
                            'defect', 'sex', 'm.age', 'birth.wt', 'state')
    
    cox <- coxph(Surv(time, cancer) ~ defect + sex + m.age + birth.wt + state, data = goback.surv) %>% 
      summary() %>% 
      coefficients()
    
    new.estimate <- data.frame(defect = j, 
                               cancer = cancers[i], 
                               hr = cox[1,2], 
                               ci.lower = exp( cox[1,1]-(1.96*cox[1,3]) ), 
                               ci.upper = exp( cox[1,1]+(1.96*cox[1,3]) ),
                               p.value.coef = cox[1,5])
    
    estimates <- rbind(estimates, new.estimate)
    
    saveRDS(estimates,
            '//smb-main.ad.bcm.edu/genepi/GOBACK/GOBACK_association_analyses/3_results/goback_2_cox_results_non_syndromic_defects_20230227.rds')
    
    rm(goback.surv, cox, new.estimate)
    
    gc()
    
  }
  
}

estimates <- estimates %>% 
  mutate(across(hr:ci.upper, ~ round(.x, 2)))

estimates
