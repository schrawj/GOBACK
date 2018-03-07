#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.03.07.
#' 
#' Analyses for BTEC abstract.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------



# Prep environment --------------------------------------------------------

require(dplyr); require(survival)

setwd('Z:/Jeremy/GOBACK/Datasets/')

cbts <- c('astro','medullo','ependymoma','pnet')

defects <- c('any.birthdefect', 'conganomalies.cns')



# Generate Cox model: non-chromosomal -------------------------------------

#' There are no striking associations of chromosomal defects with CBT, 
#' so we will only perform this analysis in the non-chromosomal set.
load('./goback.no.chrom.v20180122.1.rdata')

#' Compute a variable 'any.cbt', i.e. any childhood brain tumor.
#' Set to 1 if child DX'd with astro, medullo, ependymoma or PNET.
goback.nochrom$any.cbt <- 0

for (i in cbts){
  goback.nochrom[,'any.cbt'] <- ifelse(goback.nochrom[,i] == 1, 1, goback.nochrom[,'any.cbt'])
}

for (i in defects){
  
  tab <- as.numeric(table(goback.nochrom[,i], goback.nochrom$any.cbt)[2,2])
  
  goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                            cancer = goback.nochrom$any.cbt,
                            defect = goback.nochrom[,i],
                            sex = factor(goback.nochrom$sex,
                                         levels = c(1,2),
                                         labels = c('Male','Female')),
                            m.age = goback.nochrom$m.age,
                            state = goback.nochrom$state.num)
  
  cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
  cox.coef <- summary(cox)$coefficients
  test.ph <- cox.zph(cox)
  test.ph <- test.ph$table['defect','p']
  
  estimates <- data.frame(defect = as.character(i), 
                          cancer = 'any.cbt', 
                          HR = exp(cox.coef[1,1]), 
                          ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                          ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                          p.value.coef = cox.coef[1,5],
                          p.value.zph = test.ph,
                          num.comorbid = tab, 
                          set = 'goback.nochrom')
  
  write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/btec.models.csv', sep=',', append = TRUE, 
              row.names = FALSE, col.names = FALSE)
  
}

rm(cox, cox.coef, estimates, goback.nochrom, goback.surv, tab, test.ph, i, cbts, defects)




