#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.07.02.
#' 
#' Risk of any cancer, astrocytoma, medulloblastoma, nephroblastoma, 
#' hepatoblastoma and neuroblastoma by number of birth defects.
#' 
#' Parameterize as categorical variable: 0 vs. 1, 2, 3, 4 or more defects.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(survival); require(survminer)

setwd('Z:/Jeremy/GOBACK/')
load('./Datasets/goback.nochrom.v20180611.rdata')

goback.nochrom$majordefect.cat <- factor(ifelse(goback.nochrom$majordefect.total == 0, 0,
                                                ifelse(goback.nochrom$majordefect.total == 1, 1,
                                                       ifelse(goback.nochrom$majordefect.total == 2, 2, 
                                                              ifelse(goback.nochrom$majordefect.total == 3, 3,
                                                                     ifelse(goback.nochrom$majordefect.total >= 4, 4, NA))))),
                                         levels = c(0:4),
                                         labels = c('0', '1', '2', '3', '4 or more'))

#' A vector of the cancer diagnoses we are interested in.
outcomes <- c('cancer','astro','hepato','medullo','nephro','neuro')

#' Cox PH models for outcomes of interest by number of birth defects.
for (i in outcomes){
  
  goback.surv <- data.frame(studyid = goback.nochrom$studyid,
                            time = goback.nochrom$person.yrs, 
                            cancer = goback.nochrom[,i], 
                            defect = goback.nochrom$majordefect.cat,
                            sex = factor(goback.nochrom$sex, 
                                         levels = c(1,2),
                                         labels = c('Male','Female')),
                            m.age = goback.nochrom$m.age,
                            state = goback.nochrom$state)
  
  cox <- cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
  cox.coef <- as.data.frame(summary(cox)$coefficients)
  
  estimates <- data.frame(var = rownames(cox.coef), 
                          hr = cox.coef$`exp(coef)`, 
                          ci.lower = exp(cox.coef$coef-(1.96*cox.coef$`se(coef)`)), 
                          ci.upper = exp(cox.coef$coef+(1.96*cox.coef$`se(coef)`)))
  
  write.csv(estimates, file = paste0('Z:/Jeremy/GOBACK/R outputs/Cancer risk by number of birth defects/',i,'.risk.by.num.defects.csv'), row.names = FALSE)

}

rm(cox.coef, estimates, goback.surv, cox, i); gc()

#' Generate plots and store as elements of a list.
the.plot.thickens <- list()

for (i in 1:length(outcomes)){
  
  goback.surv <- data.frame(studyid = goback.nochrom$studyid,
                            time = goback.nochrom$person.yrs, 
                            cancer = goback.nochrom[,outcomes[i]], 
                            defect = goback.nochrom$majordefect.cat,
                            sex = factor(goback.nochrom$sex, 
                                         levels = c(1,2),
                                         labels = c('Male','Female')),
                            m.age = goback.nochrom$m.age,
                            state = goback.nochrom$state)
  
  fit <- survfit(Surv(time, cancer) ~ defect, data = goback.surv)
  
  #' The plot for any cancer needs different y-axis limits.
  if (i == 1){
  
  new.plot <- ggsurvplot(fit, 
                         conf.int = FALSE, 
                         ylim = c(0.98, 1), 
                         ylab = 'Survival Probability', 
                         xlim = c(0,18), 
                         xlab = 'Time in Years', 
                         linetype = 'strata', 
                         legend.labs = c('No birth defect', '1 defect', '2 defects', '3 defects', '4 or more defects'))
  
  }
  
  else {
    
    new.plot <- ggsurvplot(fit, 
                           conf.int = FALSE, 
                           ylim = c(0.995, 1), 
                           ylab = 'Survival Probability', 
                           xlim = c(0,18), 
                           xlab = 'Time in Years', 
                           linetype = 'strata', 
                           legend.labs = c('No birth defect', '1 defect', '2 defects', '3 defects', '4 or more defects'))
  }
  
  the.plot.thickens[i] <- new.plot

}

names(the.plot.thickens) <- outcomes




