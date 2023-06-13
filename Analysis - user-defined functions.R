#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' GOBACK Cox proportional hazards modeling.
#' 
#' Liftover of code used to generate logistic regression models in GOBACK.
#' 
#' Will generate Cox PH models for all cancer-birth defect associations 
#' with at least 5 cormorbid cases.
#' 
#' Two sets of tables: one for kids with chromosomal or genetic anomalies, 
#' one for kids with non-chromosomal birth defects.
#' 
#' Information to record in each: 
#'    HR and 95% CI for the birth defect 
#'    P-value for the birth defect
#'    P-value for the test of the PH assumption (not performed in more 
#'    recent versions, in the interest of time)
#'    Number of comorbid cases
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

#' Only models for BD-CC associations with 5 or more co-occurring cases are computed.
#' Models for ALL, Wilms, and hepatoblastoma are adjusted for birthweight.
#' I am NOT performing tests of the proportional hazards assumption this time.  
#' Testing this assumption accounts for more than 75% of the elapsed time in every iteration of this loop.  

generate.models <- function(goback, bd.index, cancer.index){
  
  estimates <- data.frame()
  
  for (i in bd.index){
    
    for (j in cancer.index){
      
      data <- subset(goback, goback[, i] == 1 & goback[ , j] == 1)
      
      if (nrow(data) > 4 & names(goback)[j] %in% c('all','hepato','nephro')){
        
        print(paste('Computing BW-adjusted model for risk of', names(goback)[j], 'in kids with', names(goback[i])))
        
        goback.surv <- goback %>% 
          select(person.yrs, all_of(names(goback)[j]), all_of(names(goback)[i]), sex, m.age, state.num, birth.wt) %>% 
          mutate(sex = factor(sex, levels = 1:2, labels = c('Male', 'Female')))
        names(goback.surv) <- c('time',
                                'cancer',
                                'defect', 'sex', 'm.age', 'state', 'birth.wt')
        
        cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state + birth.wt, data = goback.surv)
        
        cox.coef <- summary(cox)$coefficients
        #'    test.ph <- cox.zph(cox)
        #'    test.ph <- test.ph$table['defect','p']
        
        new.estimate <- data.frame(defect = names(goback)[i], 
                                   cancer = names(goback)[j], 
                                   HR = exp(cox.coef[1,1]), 
                                   ci.lower = exp( cox.coef[1,1]-(1.96*cox.coef[1,3]) ), 
                                   ci.upper = exp( cox.coef[1,1]+(1.96*cox.coef[1,3]) ),
                                   p.value.coef = cox.coef[1,5],
                                   #' p.value.zph = test.ph,
                                   num.comorbid = as.numeric(table(goback[,i], goback[,j])[2,2]))
        
        estimates <- rbind(estimates, new.estimate)
        
      }
      
      else if (nrow(data) > 4 & !(names(goback)[j] %in% c('all','hepato','nephro'))){
        
        print(paste('Computing adjusted model for risk of', names(goback)[j], 'in kids with', names(goback[i])))
        
        goback.surv <- goback %>% 
          select(person.yrs, all_of(names(goback)[j]), all_of(names(goback)[i]), sex, m.age, state.num) %>% 
          mutate(sex = factor(sex, levels = 1:2, labels = c('Male', 'Female')))
        names(goback.surv) <- c('time',
                                'cancer',
                                'defect', 'sex', 'm.age', 'state')
        
        cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
        
        cox.coef <- summary(cox)$coefficients
        #'    test.ph <- cox.zph(cox)
        #'    test.ph <- test.ph$table['defect','p']
        
        new.estimate <- data.frame(defect = names(goback)[i], 
                                   cancer = names(goback)[j], 
                                   HR = exp(cox.coef[1,1]), 
                                   ci.lower = exp( cox.coef[1,1]-(1.96*cox.coef[1,3]) ), 
                                   ci.upper = exp( cox.coef[1,1]+(1.96*cox.coef[1,3]) ),
                                   p.value.coef = cox.coef[1,5],
                                   #' p.value.zph = test.ph,
                                   num.comorbid = as.numeric(table(goback[,i], goback[,j])[2,2]))
        
        estimates <- rbind(estimates, new.estimate)
        
      }
      
      else {
        
        next
        
      }
      
    }
    
  }
  
  return(estimates)
  
}
