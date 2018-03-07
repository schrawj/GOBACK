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



# Generate Cox model: non-chromosomal -------------------------------------

load('./goback.no.chrom.v20180122.1.rdata')

#' Compute a variable 'any.cbt', i.e. any childhood brain tumor.
#' Set to 1 if child DX'd with astro, medullo, ependymoma or PNET.
goback.nochrom$any.cbt <- 0

for (i in cbts){
  goback.nochrom[,'any.cbt'] <- ifelse(goback.nochrom[,i] == 1, 1, goback.nochrom[,'any.cbt'])
}

goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                          cancer = goback.nochrom$any.cbt,
                          defect = goback.nochrom$any.birthdefect,
                          sex = factor(goback.nochrom$sex,
                                       levels = c(1,2),
                                       labels = c('Male','Female')),
                          m.age = goback.nochrom$m.age,
                          state = goback.nochrom$state.num)


tab <- as.numeric(table(goback.nochrom$any.birthdefect, goback.nochrom$any.cbt)[2,2])

cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
cox.coef <- summary(cox)$coefficients
test.ph <- cox.zph(cox)
test.ph <- test.ph$table['defect','p']

estimates <- data.frame(defect = 'any.nonchromosomal.defect', 
                        cancer = 'any.cbt', 
                        HR = exp(cox.coef[1,1]), 
                        ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                        ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                        p.value.coef = cox.coef[1,5],
                        p.value.zph = test.ph,
                        num.comorbid = tab)

write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/btec.models.csv', sep=',', append = TRUE, 
            row.names = FALSE, col.names = FALSE)

rm(cox, cox.coef, estimates, goback.nochrom, goback.surv, tab, test.ph, i)




# Generate Cox model: chromosomal -----------------------------------------

load('./goback.chrom.v20180122.1.rdata')

#' Compute a variable 'any.cbt', i.e. any childhood brain tumor.
#' Set to 1 if child DX'd with astro, medullo, ependymoma or PNET.
goback.chrom$any.cbt <- 0

for (i in cbts){
  goback.chrom[,'any.cbt'] <- ifelse(goback.chrom[,i] == 1, 1, goback.chrom[,'any.cbt'])
}

goback.surv <- data.frame(time = goback.chrom$person.yrs,
                          cancer = goback.chrom$any.cbt,
                          defect = goback.chrom$any.birthdefect,
                          sex = factor(goback.chrom$sex,
                                       levels = c(1,2),
                                       labels = c('Male','Female')),
                          m.age = goback.chrom$m.age,
                          state = goback.chrom$state.num)

tab <- as.numeric(table(goback.chrom$any.birthdefect, goback.chrom$any.cbt)[2,2])

cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
cox.coef <- summary(cox)$coefficients
test.ph <- cox.zph(cox)
test.ph <- test.ph$table['defect','p']

estimates <- data.frame(defect = 'any.chromosomal.defect', 
                        cancer = 'any.cbt', 
                        HR = exp(cox.coef[1,1]), 
                        ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                        ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                        p.value.coef = cox.coef[1,5],
                        p.value.zph = test.ph,
                        num.comorbid = tab)

write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/btec.models.csv', sep=',', append = TRUE, 
            row.names = FALSE, col.names = FALSE)

rm(cox, cox.coef, estimates, goback.chrom, goback.surv, tab, test.ph, i)
