#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.05.07.
#' 
#' Generate a model for the risk of any cancer according to number of 
#' defects.
#' 
#' Parameterize as categorical variable: 0 vs. 1, 2, 3, 4 or more defects.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(survival); require(survminer)

setwd('Z:/Jeremy/GOBACK/')
load('./Datasets/goback.nochrom.v20180507.rdata')

goback.nochrom$majordefect.cat <- factor(ifelse(goback.nochrom$majordefect.total == 0, 0,
                                                  ifelse(goback.nochrom$majordefect.total == 1, 1,
                                                         ifelse(goback.nochrom$majordefect.total == 2, 2, 
                                                                ifelse(goback.nochrom$majordefect.total == 3, 3,
                                                                       ifelse(goback.nochrom$majordefect.total >= 4, 4, NA))))),
                                           levels = c(0:4),
                                           labels = c('0', '1', '2', '3', '4 or more'))

goback.surv <- data.frame(studyid = goback.nochrom$studyid,
                          time = goback.nochrom$person.yrs, 
                          cancer = goback.nochrom$cancer, 
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

write.csv(estimates, file = './R outputs/Number of defects - model estimates.csv', row.names = FALSE)

# Diagnostics -------------------------------------------------------------

#' Schoenfeld residuals.
test.ph <- cox.zph(cox)
ggcoxzph(test.ph)  

#' Kaplan-Meier curves.
fit <- survfit(Surv(time, cancer) ~ defect, data = goback.surv)

#' With confidence bands.
ggsurvplot(fit, 
           conf.int = TRUE, 
           ylim = c(0.985, 1), ylab = 'Survival Probability',
           xlim = c(0,18), xlab = 'Time in Years', 
           linetype = 'strata', 
           risk.table = TRUE,
           legend.labs = c('No birth defect', '1 defect', '2 defects', '3 defects', '4 or more defects'))

#' Without confidence bands.
ggsurvplot(fit, 
           conf.int = FALSE, 
           ylim = c(0.985, 1), ylab = 'Survival Probability', 
           xlim = c(0,18), xlab = 'Time in Years', 
           linetype = 'strata', 
           risk.table = TRUE,
           legend.labs = c('No birth defect', '1 defect', '2 defects', '3 defects', '4 or more defects'))
