#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2017.12.19.
#' 
#' Convert GOBACK data into a form that's suited for survival analysis.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------



# Prep environment --------------------------------------------------------

load("Z:/Jeremy/GOBACK/Datasets/goback.no.chrom.v20171211.1.rdata")

require(dplyr)
require(survival)
require(survminer)
require(SurvRegCensCov)

setwd('Z:/Jeremy/GOBACK/')




# User-defined functions --------------------------------------------------

#' A function to print upper and lower bounds for 95% CIs on backtransformed
#' log logistic model parameters.
print.ci <- function(estimate, se, scale){
  print(paste('lower bound:', exp(((-1*estimate)-(1.96*se)) * (1/scale))
  ))
  print(paste('upper bound:', exp(((-1*estimate)+(1.96*se)) * (1/scale))
  ))
}



# Reshape data ------------------------------------------------------------

goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                          cancer = goback.nochrom$cancer,
                          defect.status = goback.nochrom$any.birthdefect,
                          sex = factor(goback.nochrom$sex,
                                       levels = c(1,2,9),
                                       labels = c('Male','Female','Unknown')),
                          m.age = goback.nochrom$m.age,
                          state = goback.nochrom$state)

tmp <- Surv(goback.surv$time, goback.surv$cancer)



# Kaplan Meier and Cox models and plots -----------------------------------

#' Fit Kaplan Meier estimators.
fit <- survfit(tmp ~ defect.status, data = goback.surv)
fit <- survfit(Surv(time, cancer) ~ defect.status, data = goback.surv)
summary(fit)
ggsurvplot(fit, conf.int = TRUE, ylim = c(0.995, 1), linetype = 'strata', risk.table = TRUE)
plot(fit, main = 'Kaplan-Meier estimate with 95% confidence bounds', ylab = 'Survival function', xlab = 'Time', xlim = c(0,15), ylim = c(0.99,1))
ggsave('kaplan meier curves any cancer by any birth defect.pdf', scale = 2)

#' Fit Cox PH estimators.
goback.cox <- coxph(tmp ~ defect.status + sex + m.age + state, data = goback.surv)
summary(goback.cox)

ggsurvplot(survfit(goback.cox), data = goback.surv,  ylim = c(0.995, 1))

defect.df <- data.frame(defect.status = c(0,1),
                        m.age = rep(mean(goback.nochrom$m.age, na.rm = TRUE),2),
                        sex = 2,
                        state = 'TX')
defect.df

fit <- survfit(goback.cox, newdata = defect.df)
ggsurvplot(fit, data = defect.df, ylim = c(0.995,1))




# Test assumptions --------------------------------------------------------

#' Test that proportional hazards holds.
test.ph <- cox.zph(goback.cox)
test.ph
ggcoxzph(test.ph)  
ggsave('schoenfeld residuals.pdf', scale = 2)

#' Influential observations.
ggcoxdiagnostics(goback.cox, type = 'dfbeta', linear.predictions = FALSE, ggtheme = theme_bw())
  
  

# Parametric survival models ----------------------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' We discussed a few parametric survival models we might use when we
#' met with Sue Hilsenbeck.  She mentioned Weibull, log normal and log
#' logistic.
#' 
#' Exponential and Weibull are identical when the scale parameter in the 
#' model statement is set to 1.
#' 
#' The survival package supports each of these, plus exponential.
#' 
#' Build models using each distribution and compare their efficacy.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                          cancer = goback.nochrom$cancer,
                          defect.status = goback.nochrom$any.birthdefect,
                          sex = factor(goback.nochrom$sex,
                                       levels = c(1,2,9),
                                       labels = c('Male','Female','Unknown')),
                          m.age = goback.nochrom$m.age,
                          state = goback.nochrom$state)

goback.surv$time <- goback.surv$time + 0.00001

goback.exp <- survreg(Surv(goback.surv$time, goback.surv$cancer) ~ goback.surv$defect.status + goback.surv$sex + goback.surv$m.age + goback.surv$state, 
                      dist = 'exponential')
summary(goback.exp)

#' It is more convenient to perform Weibull AFT regression via the WeibullReg 
#' function in the SurvRegCensCov package, which requires survival.
#' In this way we automatically back-transform estimates to HRs and associated CIs.
goback.wei <- WeibullReg(Surv(time, cancer) ~ defect.status + m.age + sex + state, data = goback.surv)
goback.wei

goback.loglog <- survreg(Surv(goback.surv$time, goback.surv$cancer) ~ goback.surv$defect.status + goback.surv$sex + goback.surv$m.age + goback.surv$state, 
                      dist = 'loglogistic')
summary(goback.loglog)

print.ci(0.3854,0.05230,1.59) # TX
print.ci(0.1398,0.0666,1.59) # NC
print.ci(0.3913,0.0555,1.59) # MI
print.ci(-0.0163,0.00215,1.59) # maternal age
print.ci(0.0951,0.05630,1.59) # Female sex
print.ci(-1.5105,0.04230,1.59) # Defect status

goback.lognorm <- survreg(Surv(goback.surv$time, goback.surv$cancer) ~ goback.surv$defect.status + goback.surv$sex + goback.surv$m.age + goback.surv$state, 
                      dist = 'lognormal')
summary(goback.lognorm)

defect.df <- data.frame(defect.status = c(0,1),
                        m.age = rep(mean(goback.nochrom$m.age, na.rm = TRUE),2),
                        sex = 2,
                        state = 'TX')
defect.df

fit <- survfit(goback.wei, newdata = defect.df)
ggsurvplot(fit, data = defect.df, ylim = c(0.995,1))




# Compare different parametric survival models via AIC --------------------

#' Obtain AIC values for each model.
l <- list(goback.cox, goback.wei, goback.loglog, goback.lognorm, goback.exp)

lapply(l, function(x) {extractAIC(x)})

rm(l)



# Weibull diagnostic plots ------------------------------------------------

#' If the Weibull model is a good fit, plots should be linear and parallel.
#' There is evidence that the Weibull model is a poor fit, particular as 
#' survival time increases.
WeibullDiag(Surv(time, cancer) ~ defect.status, data = goback.surv)



# Censoring at 18 ---------------------------------------------------------

#' Both the Cox and Weibull models appear to fail at or after age 18.
#' What would happen if we right-censored all individuals at 18?
goback.nochrom$dxby18 <- ifelse(goback.nochrom$cancer == 1 & goback.nochrom$person.yrs <= 18, 1, 0)

#' Set max follow up to 18.
#' Only consider children to have cancer if they are flagged as ever cancer and person years is <= 18.
goback.surv.pedi <- data.frame(time = ifelse(goback.nochrom$person.yrs > 18, 18, goback.nochrom$person.yrs),
                              cancer = ifelse(goback.nochrom$dxby18 == 1, 1, 0),
                              defect.status = goback.nochrom$any.birthdefect,
                              sex = goback.nochrom$sex,
                              m.age = goback.nochrom$m.age,
                              state = goback.nochrom$state)

#' Plot Kaplan-Meier survival curves.
goback.km <- survfit(Surv(time, cancer) ~ defect.status, data = goback.surv.pedi)
summary(goback.km)
ggsurvplot(goback.km, conf.int = TRUE, ylim = c(0.99, 1), linetype = 'strata', risk.table = TRUE)

#' Cox PH model and diagnostics.
goback.cox <- coxph(Surv(time, cancer) ~ defect.status + sex + m.age + state, data = goback.surv.pedi)
summary(goback.cox)

test.ph <- cox.zph(goback.cox)
test.ph
ggcoxzph(test.ph)  

goback.cox.interact <- coxph(Surv(time, cancer) ~ defect.status + sex + m.age + state + defect.status*time, data = goback.surv.pedi)
summary(goback.cox.interact)

test.ph <- cox.zph(goback.cox.interact)
test.ph
ggcoxzph(test.ph)  

#' Re-assess Weibull model fit.
goback.wei <- WeibullReg(Surv(time, cancer) ~ defect.status + m.age + sex + state, data = goback.surv.pedi)
goback.wei
WeibullDiag(Surv(time, cancer) ~ defect.status, data = goback.surv.pedi)

gobacksurv <- Surv(goback.surv.pedi$time, goback.surv.pedi$cancer)
plot(survfit(gobacksurv ~ goback.surv.pedi$defect.status), col = c('black','red'), fun = 'cloglog')
plot(survfit(gobacksurv ~ goback.surv.pedi$sex), col = c('black','red'), fun = 'cloglog')

# Repeat, excluding MI ----------------------------------------------------

goback.nochrom <- filter(goback.nochrom, state != 'MI')

#' Both the Cox and Weibull models appear to fail at or after age 18.
#' What would happen if we right-censored all individuals at 18?
goback.nochrom$dxby18 <- ifelse(goback.nochrom$cancer == 1 & goback.nochrom$person.yrs <= 18, 1, 0)

#' Set max follow up to 18.
#' Only consider children to have cancer if they are flagged as ever cancer and person years is <= 18.
goback.surv.pedi <- data.frame(time = ifelse(goback.nochrom$person.yrs > 18, 18, goback.nochrom$person.yrs),
                               cancer = ifelse(goback.nochrom$dxby18 == 1, 1, 0),
                               defect.status = goback.nochrom$any.birthdefect,
                               sex = goback.nochrom$sex,
                               m.age = goback.nochrom$m.age,
                               state = goback.nochrom$state)

#' Plot Kaplan-Meier survival curves.
goback.km <- survfit(Surv(time, cancer) ~ defect.status, data = goback.surv.pedi)
ggsurvplot(goback.km, conf.int = TRUE, ylim = c(0.99, 1), linetype = 'strata', risk.table = TRUE)

#' Cox PH model and diagnostics.
goback.cox <- coxph(Surv(time, cancer) ~ defect.status + sex + m.age + state, data = goback.surv.pedi)
summary(goback.cox)

test.ph <- cox.zph(goback.cox)
test.ph
ggcoxzph(test.ph)  

goback.cox.interact <- coxph(Surv(time, cancer) ~ defect.status + sex + m.age + state + defect.status*time, data = goback.surv.pedi)
summary(goback.cox.interact)

test.ph <- cox.zph(goback.cox.interact)
test.ph
ggcoxzph(test.ph)  

#' Re-assess Weibull model fit.
goback.surv.pedi$time <- goback.surv.pedi$time + 0.000001
goback.wei <- WeibullReg(Surv(time, cancer) ~ defect.status + m.age + sex + state, data = goback.surv.pedi)
goback.wei
WeibullDiag(Surv(time, cancer) ~ defect.status, data = goback.surv.pedi)

gobacksurv <- Surv(goback.surv.pedi$time, goback.surv.pedi$cancer)
plot(survfit(gobacksurv ~ goback.surv.pedi$defect.status), col = c('black','red'), fun = 'cloglog')
plot(survfit(gobacksurv ~ goback.surv.pedi$sex), col = c('black','red'), fun = 'cloglog')
# Log logistic model diagnostic plots -------------------------------------

tmp <- filter(goback.surv, defect.status == 1)
fit <- survfit(Surv(time, cancer) ~ 1, data = tmp)
ggsurvplot(fit, conf.int = TRUE, ylim = c(0.995, 1), linetype = 'strata', risk.table = TRUE)

tmp3 <- data.frame(time = fit$time, surv = fit$surv)

tmp2 <- filter(goback.surv, defect.status == 0)
fit <- survfit(Surv(time, cancer) ~ 1, data = tmp2)
ggsurvplot(fit, conf.int = TRUE, ylim = c(0.995, 1), linetype = 'strata', risk.table = TRUE)

tmp4 <- data.frame(time = fit$time, surv = fit$surv)

print(ggplot(data = tmp2, aes(x = time, y = surv)) + geom_point())
print(ggplot(data = tmp3, aes(x = log(time), y = log(-log(surv)))) + geom_point() + lims(y = c(-10, 1)))
print(ggplot() + 
        geom_smooth(data = tmp3, aes(x = log(time), y = log(-log(surv)))) + 
        geom_smooth(data = tmp4, aes(x = log(time), y = log(-log(surv)))) +  
        lims(y = c(-10, 1)))

#' Log logistic regression diagnostic plot.
print(ggplot() + 
        geom_smooth(data = tmp3, aes(x = log(time), y = log(surv)), color = 'red') + 
        geom_smooth(data = tmp4, aes(x = log(time), y = log(surv)), color = 'blue') + 
        ggtitle('Log Logistic Diagnostic Plot', subtitle = 'Children with non-chromosomal defects (red) and without (blue)'))



# Censoring at 18 ---------------------------------------------------------

goback.nochrom$dxby18 <- ifelse(goback.nochrom$cancer == 1 & goback.nochrom$person.yrs <= 18, 1, 0)

#' Set max follow up to 18.
#' Only consider children to have cancer if they are flagged as ever cancer and person years is <= 18.
goback.surv.pedi <- data.frame(time = ifelse(goback.nochrom$person.yrs > 18, 18, goback.nochrom$person.yrs),
                               cancer = ifelse(goback.nochrom$dxby18 == 1, 1, 0),
                               defect.status = goback.nochrom$any.birthdefect,
                               sex = goback.nochrom$sex,
                               m.age = goback.nochrom$m.age,
                               state = goback.nochrom$state)

tmp <- filter(goback.surv.pedi, defect.status == 1)
fit <- survfit(Surv(time, cancer) ~ 1, data = tmp)
ggsurvplot(fit, conf.int = TRUE, ylim = c(0.995, 1), linetype = 'strata', risk.table = TRUE)

tmp3 <- data.frame(time = fit$time, surv = fit$surv)

tmp2 <- filter(goback.surv.pedi, defect.status == 0)
fit <- survfit(Surv(time, cancer) ~ 1, data = tmp2)
ggsurvplot(fit, conf.int = TRUE, ylim = c(0.995, 1), linetype = 'strata', risk.table = TRUE)

tmp4 <- data.frame(time = fit$time, surv = fit$surv)

print(ggplot() + 
      geom_point(data = tmp3, aes(x = time, y = surv)) + 
      geom_point(data = tmp4, aes(x = time, y = surv)))

#' Log logistic regression diagnostic plot.
print(ggplot() + 
        geom_smooth(data = tmp3, aes(x = log(time), y = log(surv)), color = 'red') + 
        geom_smooth(data = tmp4, aes(x = log(time), y = log(surv)), color = 'blue') + 
        ggtitle('Log Logistic Diagnostic Plot', subtitle = 'Children with non-chromosomal defects (red) and without (blue).  Data right-censored at 18 years.') +
        lims(y=c(-0.1,0.1)))





print(ggplot() + 
        geom_smooth(data = tmp3, aes(x = time, y = -log(surv)), color = 'red') + 
        geom_smooth(data = tmp4, aes(x = time, y = -log(surv)), color = 'blue') + 
        ggtitle('Log Logistic Diagnostic Plot', subtitle = 'Children with non-chromosomal defects (red) and without (blue)'))

print(ggplot() + 
        geom_smooth(data = tmp3, aes(x = log(time), y = log(-log(surv))), color = 'red') + 
        geom_smooth(data = tmp4, aes(x = log(time), y = log(-log(surv))), color = 'blue') + 
        ggtitle('Log Logistic Diagnostic Plot', subtitle = 'Children with non-chromosomal defects (red) and without (blue)'))
