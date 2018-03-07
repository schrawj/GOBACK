#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.01.03.
#' 
#' Philip wants to compare different survival models in the context of a 
#' specific BD-CC association.  Somewhat arbitrarily, we chose atrial 
#' septal defect and hepatoblastoma.
#' 
#' Compare model AIC for Cox, Weibull and Log logistic models.
#' 
#' Compre parameter estimates from the three models.
#' 
#' Generate test statistics for the PH assumption and Schoenfeld residual 
#' plots for the Cox model.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------



# Prep environment --------------------------------------------------------

load("Z:/Jeremy/GOBACK/Datasets/goback.no.chrom.v20171211.1.rdata")

require(dplyr); require(survival); require(survminer); require(SurvRegCensCov)

setwd('Z:/Jeremy/GOBACK/')



# User-defined functions --------------------------------------------------

#' A function to print upper and lower bounds for 95% CIs on backtransformed
#' log logistic model parameters.
print.ci <- function(estimate, se, scale){
  point.est <- exp((-1*estimate)*(1/scale))
  lower.bound <- exp(((-1*estimate)-(1.96*se)) * (1/scale))
  upper.bound <- exp(((-1*estimate)+(1.96*se)) * (1/scale))
  print(paste(point.est, '(', lower.bound, '-', upper.bound, ')'))
}



# Reshape data ------------------------------------------------------------

#' If no cancer and hepato is NA, set hepato to 0.
#' goback.nochrom$hepato2 <- ifelse(is.na(goback.nochrom$hepato) & goback.nochrom$cancer == 0, 0, goback.nochrom$hepato)

#' Make sure that kids with cancer and NA values of hepato do not have hepato.
#' tmp <- filter(goback.nochrom, cancer == 1 & is.na(hepato2))
#' tmp <- c(tmp$studyid)

#' load("Z:/Jeremy/GOBACK/Datasets/cancer.codes.rdata")

#' tmp <- cancer.codes[cancer.codes$studyid %in% tmp, ]
#' tmp2 <- filter(tmp, morph31 == 8970)
#' tmp2 <- filter(tmp, morph32 == 8970)
#' tmp2 <- filter(tmp, morph33 == 8970)
#' tmp2 <- filter(tmp, morph34 == 8970)
#' tmp2 <- filter(tmp, morph35 == 8970)

#' None do.
#' rm(tmp, tmp2, cancer.codes)

goback.nochrom$hepato <- ifelse(is.na(goback.nochrom$hepato) & goback.nochrom$cancer == 0, 0, goback.nochrom$hepato)
goback.nochrom$hepato <- ifelse(is.na(goback.nochrom$hepato) & goback.nochrom$cancer == 1, 0, goback.nochrom$hepato)

table(goback.nochrom$atrialseptaldefect, goback.nochrom$hepato, useNA = 'always')

#' We are interested in children with ASD and hepatoblastoma.
goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                          hepatoblastoma = goback.nochrom$hepato,
                          asd = goback.nochrom$atrialseptaldefect,
                          sex = factor(goback.nochrom$sex,
                                       levels = c(1,2,9),
                                       labels = c('Male','Female','Unknown')),
                          m.age = goback.nochrom$m.age,
                          state = goback.nochrom$state)



# Fit Kaplan-Meier survival curve -----------------------------------------

fit <- survfit(Surv(time, hepatoblastoma) ~ asd, data = goback.surv)
summary(fit)
ggsurvplot(fit, conf.int = TRUE, ylim = c(0.995, 1), linetype = 'strata', risk.table = TRUE)
plot(fit, main = 'Kaplan-Meier estimate with 95% confidence bounds', ylab = 'Survival function', xlab = 'Time', xlim = c(0,15), ylim = c(0.99,1))




# Fit Cox model: no ASD*time interaction ----------------------------------

cox <- coxph(Surv(time, hepatoblastoma) ~ asd + sex + m.age + state, data = goback.surv)
summary(cox)
cox$loglik

test.ph <- cox.zph(cox)
test.ph
ggcoxzph(test.ph)  
ggsave('schoenfeld residuals for hepatoblastoma by ASD.pdf', scale = 2)

ggcoxdiagnostics(cox, type = 'dfbeta', linear.predictions = FALSE, ggtheme = theme_bw())

defect.df <- data.frame(asd = c(0,1),
                        m.age = rep(mean(goback.nochrom$m.age, na.rm = TRUE),2),
                        sex = 2,
                        state = 'TX')
defect.df

fit <- survfit(cox, newdata = defect.df)
ggsurvplot(fit, data = defect.df, ylim = c(0.999,1))




# Fit Cox model: ASD*time interaction -------------------------------------

#' Not working...
cox.int <- coxph(Surv(time, hepatoblastoma) ~ asd + sex + m.age + state + asd*time, data = goback.surv)
summary(cox.int)

defect.df <- data.frame(asd = c(0,1),
                        m.age = rep(mean(goback.nochrom$m.age, na.rm = TRUE),2),
                        sex = 2,
                        state = 'TX')
defect.df

test.ph.int <- cox.zph(cox.int)
test.ph.int
ggcoxzph(test.ph.int)  
ggsave('schoenfeld residuals for hepatoblastoma by ASD.pdf', scale = 2)

# Parametric survival models ----------------------------------------------

goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                          hepatoblastoma = goback.nochrom$hepato,
                          asd = goback.nochrom$atrialseptaldefect,
                          sex = factor(goback.nochrom$sex,
                                       levels = c(1,2,9),
                                       labels = c('Male','Female','Unknown')),
                          m.age = goback.nochrom$m.age,
                          state = goback.nochrom$state)

goback.surv$time <- goback.surv$time + 0.00001

#' It is more convenient to perform Weibull AFT regression via the WeibullReg 
#' function in the SurvRegCensCov package, which requires survival.
#' In this way we automatically back-transform estimates to HRs and associated CIs.
goback.wei <- WeibullReg(Surv(time, hepatoblastoma) ~ asd + m.age + sex + state, data = goback.surv)
goback.wei

WeibullDiag(Surv(time, hepatoblastoma) ~ asd, data = goback.surv)

goback.loglog <- survreg(Surv(goback.surv$time, goback.surv$hepatoblastoma) ~ goback.surv$asd + goback.surv$sex + goback.surv$m.age + goback.surv$state, 
                         dist = 'loglogistic')
summary(goback.loglog)

print.ci(-7.2568,0.798,2.7) # ASD
print.ci(0.9572,0.407,2.7) # Female sex
print.ci(-0.0372,0.032,2.7) # Maternal age
print.ci(0.4948,1.153,2.7) # MI
print.ci(-2.2686,1.12,2.7) # NC
print.ci(-1.3475,1.057,2.7) # TX
