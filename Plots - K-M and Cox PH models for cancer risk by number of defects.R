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

# Generate Cox PH estimates -----------------------------------------------

require(survival)

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
outcomes <- c('cancer','cns.any','neuro')

#' Cox PH models for 3 cancers of interest by number of birth defects.
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

#' BW-adjusted model for hepatoblastoma.
goback.surv <- data.frame(studyid = goback.nochrom$studyid,
                          time = goback.nochrom$person.yrs, 
                          cancer = goback.nochrom$hepato, 
                          defect = goback.nochrom$majordefect.cat,
                          sex = factor(goback.nochrom$sex, 
                                       levels = c(1,2),
                                       labels = c('Male','Female')),
                          m.age = goback.nochrom$m.age,
                          state = goback.nochrom$state,
                          birth.wt = goback.nochrom$birth.wt)

cox <- cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state + birth.wt, data = goback.surv)
cox.coef <- as.data.frame(summary(cox)$coefficients)

estimates <- data.frame(var = rownames(cox.coef), 
                        hr = cox.coef$`exp(coef)`, 
                        ci.lower = exp(cox.coef$coef-(1.96*cox.coef$`se(coef)`)), 
                        ci.upper = exp(cox.coef$coef+(1.96*cox.coef$`se(coef)`)))

write.csv(estimates, file = paste0('Z:/Jeremy/GOBACK/R outputs/Cancer risk by number of birth defects/hepato.risk.by.num.defects.csv'), row.names = FALSE)

rm(list = ls()); gc()



# Generate K-M curves -----------------------------------------------------

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
outcomes <- c('cancer','cns.any','hepato','neuro')

the.plots.thicken <- list()

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
  #' All plots have axis and legend labels silenced.
  if (i == 1){
  
    new.plot <- ggsurvplot(fit, 
                           conf.int = FALSE, 
                           ylab = NULL,
                           ylim = c(0.98, 1), 
                           xlab = NULL,
                           xlim = c(0,18), 
                           font.tickslab = c(15, 'bold', 'black'),
                           linetype = 'strata', 
                           legend = "none")
  
  }
  
  else {
    
    new.plot <- ggsurvplot(fit, 
                           conf.int = FALSE, 
                           ylab = NULL,
                           ylim = c(0.995, 1), 
                           xlab = NULL,
                           xlim = c(0,18), 
                           font.tickslab = c(15, 'bold', 'black'),
                           linetype = 'strata', 
                           legend = "none")
  }
  
  the.plots.thicken[i] <- new.plot

}

names(the.plots.thicken) <- outcomes

#' We need one plot with the strata in the legend so we can steal it for the figure.
print(ggsurvplot(fit, 
                       conf.int = FALSE, 
                       ylab = NULL,
                       ylim = c(0.995, 1), 
                       font.tickslab = c(15, 'bold', 'black'),
                       xlab = NULL,
                       xlim = c(0,18), 
                       linetype = 'strata', 
                       legend.labs = c('No defect', '1 defect', '2 defects', '3 defects', '4 or more defects')))



# Generate Cox models for new variables -----------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.07.17.
#' 
#' Sharon suggested we use any cancer, any hematologic cancer, any cns
#' tumor and any non-cns solid tumor as the categories.
#' 
#' Parameterize as categorical variable: 0 vs. 1, 2, 3, 4 or more defects.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(survival)

load('Z:/Jeremy/GOBACK/Datasets/goback.nochrom.v20180711.rdata')

goback.nochrom$majordefect.cat <- factor(ifelse(goback.nochrom$majordefect.total == 0, 0,
                                                ifelse(goback.nochrom$majordefect.total == 1, 1,
                                                       ifelse(goback.nochrom$majordefect.total == 2, 2, 
                                                              ifelse(goback.nochrom$majordefect.total == 3, 3,
                                                                     ifelse(goback.nochrom$majordefect.total >= 4, 4, NA))))),
                                         levels = c(0:4),
                                         labels = c('0', '1', '2', '3', '4 or more'))


#' Need new variables for any heme cancer (any leukemia or any lymphoma) + 
#' any non-CNS solid tumor.
goback.nochrom$any.heme.cancer <- ifelse(goback.nochrom$leu.any == 1 | goback.nochrom$lym.any == 1, 1, 0)

solid.tumors <- unique(goback.nochrom$cancer1)
solid.tumors <- subset(solid.tumors, !(solid.tumors %in% c(NA, 'all','leu.other','aml','hl','nhl','cns.other','medullo','pnet','lym.other',
                                                           'gct.intra','astro','ependymoma')))

goback.nochrom$any.non.cns.solid.tumor <- ifelse(goback.nochrom$cancer1 %in% solid.tumors, 1, 0)

outcomes <- c('any.heme.cancer','any.non.cns.solid.tumor')

#' Cox PH models for cancers of interest by number of birth defects.
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



# Generate K-M curves for new variables -----------------------------------

require(survival); require(survminer)

load('Z:/Jeremy/GOBACK/Datasets/goback.nochrom.v20180711.rdata')

goback.nochrom$majordefect.cat <- factor(ifelse(goback.nochrom$majordefect.total == 0, 0,
                                                ifelse(goback.nochrom$majordefect.total == 1, 1,
                                                       ifelse(goback.nochrom$majordefect.total == 2, 2, 
                                                              ifelse(goback.nochrom$majordefect.total == 3, 3,
                                                                     ifelse(goback.nochrom$majordefect.total >= 4, 4, NA))))),
                                         levels = c(0:4),
                                         labels = c('0', '1', '2', '3', '4 or more'))


#' Need new variables for any heme cancer (any leukemia or any lymphoma) + 
#' any non-CNS solid tumor.
goback.nochrom$any.heme.cancer <- ifelse(goback.nochrom$leu.any == 1 | goback.nochrom$lym.any == 1, 1, 0)

solid.tumors <- unique(goback.nochrom$cancer1)
solid.tumors <- subset(solid.tumors, !(solid.tumors %in% c(NA, 'all','leu.other','aml','hl','nhl','cns.other','medullo','pnet','lym.other',
                                                           'gct.intra','astro','ependymoma')))

goback.nochrom$any.non.cns.solid.tumor <- ifelse(goback.nochrom$cancer1 %in% solid.tumors, 1, 0)

outcomes <- c('any.heme.cancer','cns.any','any.non.cns.solid.tumor')

the.plots.thicken <- list()

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
  
  new.plot <- ggsurvplot(fit, 
                         conf.int = FALSE, 
                         ylab = NULL,
                         ylim = c(0.99, 1), 
                         xlab = NULL,
                         xlim = c(0,18), 
                         font.tickslab = c(15, 'bold', 'black'),
                         linetype = 'strata', 
                         legend = "none")
  
  the.plots.thicken[i] <- new.plot
  
}

names(the.plots.thicken) <- outcomes



