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

# Cox PH models in kids without chromosomal defects -----------------------

require(survival); require(dplyr)

load('Z:/Jeremy/GOBACK/Datasets/goback.nochrom.v20180611.rdata')

#' Legacy code demonstrating how state.num was computed is commented out.
#' codes <- data.frame(state = c('NC','TX','MI','AR'),
#'                     state.num = 1:4)
#' goback.nochrom <- filter(goback.nochrom, sex != 9)
#' goback.nochrom <- left_join(goback.nochrom, codes, by = 'state')

#' Models for specific cancers in children with non-chromosomal defects.
#' Only models for BD-CC associations with 5 or more co-occurring cases are computed.
#' Models for ALL, Wilms, and hepatoblastoma are adjusted for birthweight.
#' I am NOT performing tests of the proportional hazards assumption this time.  
#' Testing this assumption accounts for more than 75% of the elapsed time in every iteration of this loop.  
#' This code is commented out.  
#' Results are available in prior versions of this output and should have changed negligibly.
for (i in 22:94){
  
  tmp <- table(goback.nochrom[,i], goback.nochrom$cancer1)
  tmp <- which(tmp[2, ] >= 5)
  tmp <- tmp + 111
  
  if (length(tmp) > 0){
    
    for (j in tmp){
      
      if (j %in% c(112,126,131)){
        
        print(paste('Computing BW-adjusted model for risk of',names(goback.nochrom[j]), 'in kids with', names(goback.nochrom[i])))
        
        goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                                  cancer = goback.nochrom[,j],
                                  defect = goback.nochrom[,i],
                                  sex = factor(goback.nochrom$sex,
                                               levels = c(1,2),
                                               labels = c('Male','Female')),
                                  m.age = goback.nochrom$m.age,
                                  state = goback.nochrom$state.num,
                                  birth.wt = goback.nochrom$birth.wt)
        
        cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state + birth.wt, data = goback.surv)
        cox.coef <- summary(cox)$coefficients
        #'    test.ph <- cox.zph(cox)
        #'    test.ph <- test.ph$table['defect','p']
        
        estimates <- data.frame(defect = names(goback.nochrom[i]), 
                                cancer = names(goback.nochrom[j]), 
                                HR = exp(cox.coef[1,1]), 
                                ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                                ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                                p.value.coef = cox.coef[1,5],
                                #' p.value.zph = test.ph,
                                num.comorbid = as.numeric(table(goback.nochrom[,i], goback.nochrom[,j])[2,2]))
        
        write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.coxph.models.v20180612.csv', 
                    sep=',', append = TRUE, row.names = FALSE, col.names = FALSE)
        
      }
      
      else{
        
        print(paste('Computing model for risk of',names(goback.nochrom[j]), 'in kids with', names(goback.nochrom[i])))
        
        goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                                  cancer = goback.nochrom[,j],
                                  defect = goback.nochrom[,i],
                                  sex = factor(goback.nochrom$sex,
                                               levels = c(1,2),
                                               labels = c('Male','Female')),
                                  m.age = goback.nochrom$m.age,
                                  state = goback.nochrom$state.num)
        
        cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
        cox.coef <- summary(cox)$coefficients
        #'    test.ph <- cox.zph(cox)
        #'    test.ph <- test.ph$table['defect','p']
        
        estimates <- data.frame(defect = names(goback.nochrom[i]), 
                                cancer = names(goback.nochrom[j]), 
                                HR = exp(cox.coef[1,1]), 
                                ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                                ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                                p.value.coef = cox.coef[1,5],
                                #' p.value.zph = test.ph,
                                num.comorbid = as.numeric(table(goback.nochrom[,i], goback.nochrom[,j])[2,2]))
        
        write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.coxph.models.v20180612.csv', 
                    sep=',', append = TRUE, row.names = FALSE, col.names = FALSE)
      }
      
    }
    
  }
  
  else{
    
    next
    
  }
}

rm(cox.coef, estimates, goback.surv, cox, i, j, tmp); gc()

#' models for '[cancer].any' variables.
for (i in 22:94){
  
  for (j in 142:151){

    print(paste('Computing model for risk of',names(goback.nochrom[j]), 'in kids with', names(goback.nochrom[i])))
    
    comorbid.cases <- table(goback.nochrom[,i], goback.nochrom[,j])[2,2]
    
    if (comorbid.cases >= 5){

      goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                                cancer = goback.nochrom[,j],
                                defect = goback.nochrom[,i],
                                sex = factor(goback.nochrom$sex,
                                             levels = c(1,2),
                                             labels = c('Male','Female')),
                                m.age = goback.nochrom$m.age,
                                state = goback.nochrom$state.num)      

      cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
      cox.coef <- summary(cox)$coefficients
#'    test.ph <- cox.zph(cox)
#'    test.ph <- test.ph$table['defect','p']

      estimates <- data.frame(defect = names(goback.nochrom[i]), 
                              cancer = names(goback.nochrom[j]), 
                              HR = exp(cox.coef[1,1]), 
                              ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                              ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                              p.value.coef = cox.coef[1,5],
#'                              p.value.zph = test.ph,
                              num.comorbid = comorbid.cases)
      
      write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.coxph.models.v20180612.csv', 
                  sep=',', append = TRUE, row.names = FALSE, col.names = FALSE)
    }
    
    else{
      
      next
      
    }
  }
}

rm(list = ls()); gc()



# Cox PH models in kids with chromosomal defects --------------------------

require(survival)

load('Z:/Jeremy/GOBACK/Datasets/goback.chrom.v20180611.rdata')

for (i in c(97:104)){

  tmp <- table(goback.chrom[,i], goback.chrom$cancer1)
  tmp <- which(tmp[2, ] >= 5)
  tmp <- tmp + 111

  if (length(tmp) > 0){
    
    for (j in tmp){
      
      if (j %in% c(112,126,131)){
        
        print(paste('Computing BW-adjusted model for risk of',names(goback.chrom[j]), 'in kids with', names(goback.chrom[i])))
        
        goback.surv <- data.frame(time = goback.chrom$person.yrs,
                                  cancer = goback.chrom[,j],
                                  defect = goback.chrom[,i],
                                  sex = factor(goback.chrom$sex,
                                               levels = c(1,2),
                                               labels = c('Male','Female')),
                                  m.age = goback.chrom$m.age,
                                  state = goback.chrom$state.num,
                                  birth.wt = goback.chrom$birth.wt)
        
        cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state + birth.wt, data = goback.surv)
        cox.coef <- summary(cox)$coefficients
        #'    test.ph <- cox.zph(cox)
        #'    test.ph <- test.ph$table['defect','p']
        
        estimates <- data.frame(defect = names(goback.chrom[i]), 
                                cancer = names(goback.chrom[j]), 
                                HR = exp(cox.coef[1,1]), 
                                ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                                ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                                p.value.coef = cox.coef[1,5],
                                #' p.value.zph = test.ph,
                                num.comorbid = as.numeric(table(goback.chrom[,i], goback.chrom[,j])[2,2]))
        
        write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.coxph.models.v20180612.csv', 
                    sep=',', append = TRUE, row.names = FALSE, col.names = FALSE)
        
      }
      
      else{
        
        print(paste('Computing model for risk of',names(goback.chrom[j]), 'in kids with', names(goback.chrom[i])))
        
        goback.surv <- data.frame(time = goback.chrom$person.yrs,
                                  cancer = goback.chrom[,j],
                                  defect = goback.chrom[,i],
                                  sex = factor(goback.chrom$sex,
                                               levels = c(1,2),
                                               labels = c('Male','Female')),
                                  m.age = goback.chrom$m.age,
                                  state = goback.chrom$state.num)
        
        cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
        cox.coef <- summary(cox)$coefficients
        #'    test.ph <- cox.zph(cox)
        #'    test.ph <- test.ph$table['defect','p']
        
        estimates <- data.frame(defect = names(goback.chrom[i]), 
                                cancer = names(goback.chrom[j]), 
                                HR = exp(cox.coef[1,1]), 
                                ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                                ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                                p.value.coef = cox.coef[1,5],
                                #' p.value.zph = test.ph,
                                num.comorbid = as.numeric(table(goback.chrom[,i], goback.chrom[,j])[2,2]))
        
        write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.coxph.models.v20180612.csv', 
                    sep=',', append = TRUE, row.names = FALSE, col.names = FALSE)
      }
      
    }
    
  }
  
  else{
    
    next
    
  }
}

rm(cox.coef, estimates, goback.surv, cox, i, j, tmp); gc()

#' models for '[cancer].any' variables.
for (i in c(97:104)){
  
  for (j in 142:151){
    
    print(paste('Computing model for risk of',names(goback.chrom[j]), 'in kids with', names(goback.chrom[i])))

    comorbid.cases <- table(goback.chrom[,i], goback.chrom[,j])[2,2]
    
    if (comorbid.cases >= 5){
      
      goback.surv <- data.frame(time = goback.chrom$person.yrs,
                                cancer = goback.chrom[,j],
                                defect = goback.chrom[,i],
                                sex = factor(goback.chrom$sex,
                                             levels = c(1,2),
                                             labels = c('Male','Female')),
                                m.age = goback.chrom$m.age,
                                state = goback.chrom$state.num)      
      
      cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
      cox.coef <- summary(cox)$coefficients
#'    test.ph <- cox.zph(cox)
#'    test.ph <- test.ph$table['defect','p']
      
      
      estimates <- data.frame(defect = names(goback.chrom[i]), 
                              cancer = names(goback.chrom[j]), 
                              HR = exp(cox.coef[1,1]), 
                              ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                              ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                              p.value.coef = cox.coef[1,5],
                              #' p.value.zph = test.ph,
                              num.comorbid = as.numeric(table(goback.chrom[,i], goback.chrom[,j])[2,2]))
      
      write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.coxph.models.v20180612.csv', 
                  sep=',', append = TRUE, row.names = FALSE, col.names = FALSE)
    }
    
    else{
      
      next
      
    }
  }
}

rm(list = ls()); gc()
