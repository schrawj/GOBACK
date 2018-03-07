#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' GOBACK Cox proportional hazards modeling
#' 
#' Liftover of code used to generate logistic regression models in GOBACK.
#' 
#' Will generate Cox PH models for all cancer-birth defect associations 
#' with at least 5 cormorbid cases.
#' 
#' Two sets of tables: one for kids with chromosomal birth defects, one 
#' for kids with non-chromosomal birth defects.
#' 
#' Information to record in each: 
#'    HR and 95% CI for the birth defect 
#'    P-value for the birth defect
#'    P-value for the test of the PH assumption
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------



# prep environment --------------------------------------------------------

require(dplyr); require(survival); require(ggplot2)

setwd('Z:/Jeremy/GOBACK/Datasets/')



# Cox PH models in kids without chromosomal defects -----------------------

load('goback.no.chrom.v20180117.1.rdata')

codes <- data.frame(state = c('NC','TX','MI','AR'),
                    state.num = 1:4)

goback.nochrom <- filter(goback.nochrom, sex != 9)
goback.nochrom <- left_join(goback.nochrom, codes, by = 'state')

for (i in 22:94){
  #' Generate a vector holding the column indices of the cancers for which there
  #' are at least 5 comorbid cases in connection with the index defect.  
  tmp <- table(goback.nochrom[,i], goback.nochrom$cancer1)
  tmp <- which(tmp[2, ] >= 5)
  tmp <- tmp + 107
  #' If there is at least one cancer with 5 cases...  
  if (length(tmp) > 0){
    #'...then for all such cancers..    
    for (j in tmp){
      #' get the name of the defect and the cancer.      
      z <- names(goback.nochrom[i])
      y <- names(goback.nochrom[j])
      #' generate a new data frame suitable for survival analysis, 
      #' where the outcome is the jth cancer and the exposure is the ith birth defect.
      goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                                cancer = goback.nochrom[,j],
                                defect = goback.nochrom[,i],
                                sex = factor(goback.nochrom$sex,
                                             levels = c(1,2),
                                             labels = c('Male','Female')),
                                m.age = goback.nochrom$m.age,
                                state = goback.nochrom$state.num)
      #' Run a Cox PH model on that combination, and extract the HR and CI, 
      #' p value for coefficient and test of PH assumption.       
      cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
      cox.coef <- summary(cox)$coefficients
      tab <- as.numeric(table(goback.nochrom[,i], goback.nochrom[,j])[2,2])
      test.ph <- cox.zph(cox)
      test.ph <- test.ph$table['defect','p']
      #' Write this info to a data frame with one row.      
      estimates <- data.frame(defect = z, 
                              cancer = y, 
                              HR = exp(cox.coef[1,1]), 
                              ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                              ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                              p.value.coef = cox.coef[1,5],
                              p.value.zph = test.ph,
                              num.comorbid = tab)
      #' Append contents of that dataframe to a csv file.      
      write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.coxph.models.csv', sep=',', append = TRUE, 
                  row.names = FALSE, col.names = FALSE)
    }
    
  }
  
  else{
    
    sink(file = 'Z:/Jeremy/GOBACK/R Outputs/list of unmodeled defects.txt', append = TRUE)
    
    print(paste('There were no cancers with 5 or more comorbid cases for',names(goback.nochrom[i])))
    
    sink()
  }
}

rm(codes, cox.coef, estimates, goback.surv, cox, i, j, tab, test.ph, tmp, y, z); gc()

#' models for '[cancer].any' variables.
for (i in 22:94){
  
  for (j in 138:147){
    
    z <- names(goback.nochrom[i])
    y <- names(goback.nochrom[j])
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
      tab <- as.numeric(table(goback.nochrom[,i], goback.nochrom[,j])[2,2])
      test.ph <- cox.zph(cox)
      test.ph <- test.ph$table['defect','p']
      
      estimates <- data.frame(defect = z, 
                              cancer = y, 
                              HR = exp(cox.coef[1,1]), 
                              ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                              ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                              p.value.coef = cox.coef[1,5],
                              p.value.zph = test.ph,
                              num.comorbid = tab)
      
      write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.coxph.models.csv', sep=',', append = TRUE, 
                  row.names = FALSE, col.names = FALSE)
    }
    
    else{
      
      sink(file = 'Z:/Jeremy/GOBACK/R Outputs/list of unmodeled defects.txt', append = TRUE)
      
      print(paste('There were less than five comorbid instances of',z,'and',y))
      
      sink()
      
    }
  }
}

rm(goback.nochrom, estimates, cox, cox.coef, i, j, y, z, goback.surv, tab, test.ph); gc()

rm(codes, cox.coef, estimates, goback.surv, cox, i, j, tab, test.ph, tmp, y, z); gc()

# Cox PH models in kids with chromosomal defects --------------------------

load('goback.chrom.v20180117.1.rdata')

codes <- data.frame(state = c('NC','MI','TX','AR'),
                    state.num = 1:4)

goback.chrom <- filter(goback.chrom, sex != 9)
goback.chrom <- left_join(goback.chrom, codes, by = 'state')

for (i in 95:101){

  tmp <- table(goback.chrom[,i], goback.chrom$cancer1)
  tmp <- which(tmp[2, ] >= 5)
  tmp <- tmp + 107

  if (length(tmp) > 0){

    for (j in tmp){
 
      z <- names(goback.chrom[i])
      y <- names(goback.chrom[j])

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
      tab <- as.numeric(table(goback.chrom[,i], goback.chrom[,j])[2,2])
      test.ph <- cox.zph(cox)
      test.ph <- test.ph$table['defect','p']

      estimates <- data.frame(defect = z, 
                              cancer = y, 
                              HR = exp(cox.coef[1,1]), 
                              ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                              ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                              p.value.coef = cox.coef[1,5],
                              p.value.zph = test.ph,
                              num.comorbid = tab)
    
      write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.coxph.models.csv', sep=',', append = TRUE, 
                  row.names = FALSE, col.names = FALSE)
    }
    
  }
  
  else{
    
    sink(file = 'Z:/Jeremy/GOBACK/R Outputs/list of unmodeled defects.txt', append = TRUE)
    
    print(paste('There were no cancers with 5 or more comorbid cases for',names(goback.chrom[i])))
    
    sink()
  }
}

rm(codes, cox.coef, estimates, goback.surv, cox, i, j, tab, test.ph, tmp, y, z, comorbid.cases); gc()

#' models for '[cancer].any' variables.
for (i in 95:101){
  
  for (j in 138:147){
    
    z <- names(goback.chrom[i])
    y <- names(goback.chrom[j])
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
      tab <- as.numeric(table(goback.chrom[,i], goback.chrom[,j])[2,2])
      test.ph <- cox.zph(cox)
      test.ph <- test.ph$table['defect','p']
      
      estimates <- data.frame(defect = z, 
                              cancer = y, 
                              HR = exp(cox.coef[1,1]), 
                              ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                              ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                              p.value.coef = cox.coef[1,5],
                              p.value.zph = test.ph,
                              num.comorbid = tab)
      
      write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.coxph.models.csv', sep=',', append = TRUE, 
                  row.names = FALSE, col.names = FALSE)
    }
    
    else{
      
      sink(file = 'Z:/Jeremy/GOBACK/R Outputs/list of unmodeled defects.txt', append = TRUE)
      
      print(paste('There were less than five comorbid instances of',z,'and',y))
      
      sink()
      
    }
  }
}

rm(goback.chrom, goback.surv, cox.coef, estimates, i, j, y, z, comorbid.cases, cox, tab, test.ph); gc()



# Load in model outputs, compute FDR --------------------------------------

setwd('Z:/Jeremy/GOBACK/R outputs/')

goback.coxmodels <- read.csv(file = 'goback.coxph.models.csv', header = FALSE)
names(goback.coxmodels) <- c('defect','cancer','HR','CI.lower','CI.upper','p.val.coef','p.val.zph','n.comorbid')

#' An older version exists, dated 1-16, from before I fixed the NC cancer variables.
save(goback.coxmodels, file = 'goback.coxph.results.v20180119.1.rdata')

#' P-values < 1.11 * 10^-16 are listed as zero.  
#' Impute largest number smaller than this instead.
goback.coxmodels$p.val.coef <- ifelse(goback.coxmodels$p.val.coef < 1.11e-16, 1.10e-16, goback.coxmodels$p.val.coef)

#' Compute Benjamini-Hochberg FDR.
goback.coxmodels <- arrange(goback.coxmodels, p.val.coef) #' Arrange smallest to largest based on p-value.
goback.coxmodels$j <- 1:606 #' Ranks for p-values
m <- 606 #' Total number of tests
delta <- 0.05 #' Desired family-wise error rate
goback.coxmodels$bh.fdr <- (goback.coxmodels$j/m)*delta #' Compute value for FDR statistic at delta = 0.05
#' Find greatest rank for which p-value is less than (j/m)*delta.
#' This test and all tests ranking below are considered significant.
goback.coxmodels$fdr.flag <- goback.coxmodels$bh.fdr - goback.coxmodels$p.val.coef 

write.csv(goback.coxmodels, file = 'goback.coxph.models.with.fdr.csv', row.names = FALSE)

#' Add columns comparing these results to the crude ORs.
setwd('Z:/Jeremy/GOBACK/R outputs/Preliminary Analyses - AR MI TX/')

old.assoc <- read.csv('BD-CC associations.csv')
old.assoc <- rename(old.assoc, defect = birth.defect, crude.or = or, or.ci.lower = ci.lower, or.ci.upper = ci.upper)

goback.coxmodels <- left_join(goback.coxmodels, select(old.assoc, defect, cancer, crude.or, or.ci.lower, or.ci.upper), by = c('defect','cancer'))
goback.coxmodels$delta.hr.minus.or <- goback.coxmodels$HR - goback.coxmodels$crude.or
goback.coxmodels$relative.delta <- (goback.coxmodels$delta.hr.minus.or/goback.coxmodels$HR)*100

#' On average, estimates are not much different under the new models.
p <- ggplot(data = goback.coxmodels) + geom_histogram(aes(x=delta.hr.minus.or), color = 'red', fill = 'white') + ggtitle('Distribution of Differences in Adjusted Hazard Ratio and Crude Odds Ratio')
print(p)

p2 <- ggplot(data = goback.coxmodels) + geom_histogram(aes(x=relative.delta), color = 'red', fill = 'white') + 
  ggtitle('Distribution of Differences in Adjusted Hazard Ratio and Crude Odds Ratio',
          subtitle = 'As a Percentage of the Hazard Ratio')
print(p2)

p3 <- ggplot(data = goback.coxmodels) + geom_histogram(aes(x=p.val.zph), color = 'red', fill = 'white') +
  ggtitle('Distribution of P-values for Tests of Proportional Hazards for Defects Variables in Adjusted Cox Models')
print(p3)

setwd('Z:/Jeremy/GOBACK/R outputs/'); write.csv(goback.coxmodels, file = 'goback.cox.models.v20180119.csv', row.names = FALSE)

save(goback.coxmodels, file = 'Z:/Jeremy/GOBACK/Datasets/goback.cox.ph.results.v20180119.2.rdata')

# Compute person.years stats by cancer DX and state -----------------------

setwd('Z:/Jeremy/GOBACK/Datasets/'); load ('goback.no.chrom.v20180117.1.rdata')

setwd('Z:/Jeremy/GOBACK/R outputs/')

sink(file = 'person.years.by.cancer.and.state.for.nonchromosomal.dataset.txt', append = TRUE)
for (i in 108:147){
  print(names(goback.nochrom[i]))
  print(aggregate(person.yrs ~ goback.nochrom[,i] + goback.nochrom$state, data = goback.nochrom, summary))
}
sink()

rm(goback.nochrom); gc()

setwd('Z:/Jeremy/GOBACK/Datasets/'); load('goback.chrom.v20180117.1.rdata')

setwd('Z:/Jeremy/GOBACK/R outputs/')

sink(file = 'person.years.by.cancer.and.state.for.chromosomal.dataset.txt', append = TRUE)
for (i in 108:147){
  print(names(goback.chrom[i]))
  print(aggregate(person.yrs ~ goback.chrom[,i] + goback.chrom$state, data = goback.chrom, summary))
}
sink()

rm(goback.chrom, i); gc()
# Sensitivity analysis: Hepato ~ ASD + BW ---------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Bob Meyer from the group in NC raised a concern that birthweight may be 
#' a confounder.  He mentioned this in the context of ASD and 
#' hepatoblastoma, thought it might apply to other BD-CC associations as
#' well.  In this instance, there is uncertainty whether the registries
#' can differentiate ASD2 from PFO, which is common in premature infants.
#' 
#' He suggested adjusting for gestational age.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

load('goback.no.chrom.v20180117.1.rdata')

codes <- data.frame(state = c('NC','TX','MI','AR'),
                    state.num = 1:4)

goback.nochrom <- filter(goback.nochrom, sex != 9)
goback.nochrom <- left_join(goback.nochrom, codes, by = 'state')

#' I realize I never updated the birthweight category variable after I 
#' recovered the original continuous AR birthweight data.
goback.nochrom$birth.wt.cat <- factor(
                                      ifelse(goback.nochrom$birth.wt >= 2500 & goback.nochrom$birth.wt < 4000, 0,
                                             ifelse(goback.nochrom$birth.wt >= 4000, 1, 2)),
                                      levels = c(0:2),
                                      labels = c('NBW','HBW','LBW'))

#' Fair number of birthweights missing in TX, but mostly they're in non-BD and non-cancer
#' kids.  Probably not consequential for modeling if we exclude them.
table(goback.nochrom$birth.wt.cat, useNA = 'ifany')
table(is.na(goback.nochrom$birth.wt), goback.nochrom$state)
with(subset(goback.nochrom, state == 'TX'), table(is.na(birth.wt), any.birthdefect))
with(subset(goback.nochrom, state == 'TX'), table(is.na(birth.wt), cancer))

#' Generate a variable for pre-term birth.
goback.nochrom$preterm.birth <- factor(ifelse(goback.nochrom$gest.age < 37, 1, 0),
                                      levels = c(0,1),
                                      labels = c("Not Preterm", "Preterm"))

goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                          cancer = goback.nochrom$hepato,
                          defect = goback.nochrom$atrialseptaldefect,
                          sex = factor(goback.nochrom$sex,
                                       levels = c(1,2),
                                       labels = c('Male','Female')),
                          m.age = goback.nochrom$m.age,
                          state = goback.nochrom$state.num,
                          birth.wt = goback.nochrom$birth.wt,
                          birth.wt.cat = goback.nochrom$birth.wt.cat,
                          gest.age = goback.nochrom$gest.age,
                          preterm.birth = goback.nochrom$preterm.birth)

#' Models stratified on NBW vs. LBW.
cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = subset(goback.surv, birth.wt.cat == 'LBW'))
cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = subset(goback.surv, birth.wt.cat == 'NBW'))

#' Different adjusted models for comparison.
cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state + gest.age, data = goback.surv)
cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state + preterm.birth, data = goback.surv)
cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state + birth.wt, data = goback.surv)
cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state + birth.wt.cat, data = goback.surv)


# BW-adjusted models for Wilms/ALL/hepatoblastoma: non-chromosomal --------

load('goback.no.chrom.v20180122.1.rdata')

tmp <- filter(goback.nochrom, cancer == 0)
tmp2 <- filter(goback.nochrom, cancer == 1 & (cancer1 == 'all' | cancer1 == 'hepato' | cancer1 == 'nephro'))

rm(goback.nochrom); gc()
goback.select <- rbind(tmp, tmp2); rm(tmp, tmp2); gc()

#' Remove unnecessary columns. 
#' This task is causing repeated issues with RAM usage.
goback.select <- goback.select[,c(1,152,3,6,7,15,22:94,108,122,127)]

for (i in 7:79){
  
  for (j in 80:82){
    
    tab <- table(goback.select[,i], goback.select[,j])[2,2]
    
    if (tab >= 5){
      
      goback.surv <- data.frame(time = goback.select$person.yrs,
                                cancer = goback.select[,j],
                                defect = goback.select[,i],
                                sex = factor(goback.select$sex,
                                             levels = c(1,2),
                                             labels = c('Male','Female')),
                                m.age = goback.select$m.age,
                                state = goback.select$state.num,
                                birth.wt = goback.select$birth.wt)
      
      cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state + birth.wt, data = goback.surv)
      cox.coef <- summary(cox)$coefficients
      test.ph <- cox.zph(cox)
      test.ph <- test.ph$table['defect','p']
      rm(cox, goback.surv); gc()
      
      estimates <- data.frame(defect = names(goback.select[i]), 
                              cancer = names(goback.select[j]), 
                              HR = exp(cox.coef[1,1]), 
                              ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                              ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                              p.value.coef = cox.coef[1,5],
                              p.value.zph = test.ph,
                              num.comorbid = tab)
   
      write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.coxph.models.adjust.birthwt.csv', sep=',', append = TRUE, 
                  row.names = FALSE, col.names = FALSE)
      
    }
    
    else{
      
      next
      
    }
  }
}

rm(goback.select, cox.coef, estimates, i, j, tab, test.ph); gc()



# BW-adjusted models for Wilms/ALL/hepatoblastoma: chromosomal ------------

load('goback.chrom.v20180122.1.rdata')

tmp <- filter(goback.chrom, cancer == 0)
tmp2 <- filter(goback.chrom, cancer == 1 & (cancer1 == 'all' | cancer1 == 'hepato' | cancer1 == 'nephro'))

rm(goback.chrom); gc()
goback.select <- rbind(tmp, tmp2); rm(tmp, tmp2); gc()

goback.select <- goback.select[,c(1,152,3,6,7,15,95:101,108,122,127)]

for (i in 7:13){
  
  for (j in 14:16){
    
    tab <- table(goback.select[,i], goback.select[,j])[2,2]
    
    if (tab >= 5){
      
      goback.surv <- data.frame(time = goback.select$person.yrs,
                                cancer = goback.select[,j],
                                defect = goback.select[,i],
                                sex = factor(goback.select$sex,
                                             levels = c(1,2),
                                             labels = c('Male','Female')),
                                m.age = goback.select$m.age,
                                state = goback.select$state.num,
                                birth.wt = goback.select$birth.wt)
      
      cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state + birth.wt, data = goback.surv)
      cox.coef <- summary(cox)$coefficients
      test.ph <- cox.zph(cox)
      test.ph <- test.ph$table['defect','p']
      rm(cox, goback.surv); gc()
      
      estimates <- data.frame(defect = names(goback.select[i]), 
                              cancer = names(goback.select[j]), 
                              HR = exp(cox.coef[1,1]), 
                              ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                              ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                              p.value.coef = cox.coef[1,5],
                              p.value.zph = test.ph,
                              num.comorbid = tab)
      
      write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.coxph.models.adjust.birthwt.csv', sep=',', append = TRUE, 
                  row.names = FALSE, col.names = FALSE)
      
    }
    
    else{
      
      next
      
    }
  }
}

rm(goback.select, cox.coef, estimates, i, j, tab, test.ph); gc()

# Reconcile old and new Cox model outputs ---------------------------------

setwd('Z:/Jeremy/GOBACK/Datasets/'); load("goback.cox.ph.results.v20180119.2.rdata")

cox.new <- read.csv(file = 'Z:/Jeremy/GOBACK/R outputs/goback.coxph.models.adjust.birthwt.csv', 
                    stringsAsFactors = FALSE, header = FALSE,
                    col.names = c('defect','cancer','HR','CI.lower','CI.upper','p.val.coef','p.val.zph','n.comorbid'))

l <- c('all','hepato','nephro')

goback.coxmodels <- goback.coxmodels[,1:8]
goback.coxmodels <- filter(goback.coxmodels, !(cancer %in% l))
goback.coxmodels <- rbind(goback.coxmodels, cox.new); rm(cox.new, l)

#' Compute FDR at 5% via Benjamini-Hochberg method.
m <- 606; delta <- 0.05 

goback.coxmodels$p.val.coef <- ifelse(goback.coxmodels$p.val.coef < 1.11e-16, 1.10e-16, goback.coxmodels$p.val.coef)
goback.coxmodels <- arrange(goback.coxmodels, p.val.coef) 
goback.coxmodels$j <- 1:606 

goback.coxmodels$bh.fdr <- (goback.coxmodels$j/m)*delta 
goback.coxmodels$fdr.delta <- goback.coxmodels$bh.fdr - goback.coxmodels$p.val.coef 

save(goback.coxmodels, file = 'goback.cox.ph.results.v20180124.1.rdata'); rm(m, delta)



# Filter down to specific BD-CC pairs -------------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' THIS SECTION IS DEPRECATED.  USE THE VERSION REVISED 02/23/2018.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

#' setwd('Z:/Jeremy/GOBACK/Datasets/') 

load('goback.cox.ph.results.v20180124.1.rdata')

#' Compute Bonferroni-Holm adjusted p-values.
#' goback.coxmodels <- arrange(goback.coxmodels, p.val.coef)
#' goback.coxmodels$j.rev <- 606:1
#' goback.coxmodels$p.val.coef.bon <- 0.05/goback.coxmodels$j.rev
#' goback.coxmodels$p.val.coef.bon.delta <- goback.coxmodels$p.val.coef - goback.coxmodels$p.val.coef.bon
#' 371 results considered significant after adjustment.
#' goback.coxmodels <- goback.coxmodels[1:371, ]

#' pat1 <- 'conganomalies.'; pat2 <- '.any'; pat3 <- '.other'; pat4 <- 'chromosomalanomalies'
#' pat5 <- 'rvot.'; pat6 <- 'lvot.'; pat7 <- 'septal.defects' 

#' goback.coxmodels <- goback.coxmodels[!grepl(pat1, goback.coxmodels$defect), ]
#' goback.coxmodels <- goback.coxmodels[!grepl(pat4, goback.coxmodels$defect), ]

#' goback.coxmodels <- goback.coxmodels[!grepl(pat2, goback.coxmodels$cancer), ]
#' goback.coxmodels <- goback.coxmodels[!grepl(pat3, goback.coxmodels$cancer), ]

#' goback.coxmodels <- goback.coxmodels[!grepl(pat5, goback.coxmodels$defect), ]
#' goback.coxmodels <- goback.coxmodels[!grepl(pat6, goback.coxmodels$defect), ]
#' goback.coxmodels <- goback.coxmodels[!grepl(pat7, goback.coxmodels$defect), ]

#' top.hits <- goback.coxmodels; rm(goback.coxmodels, pat1, pat2, pat3, pat4, pat5, pat6, pat7)

#' save(top.hits, file = 'goback.cox.ph.top.hits.v20180124.1.rdata')

#' write.csv(top.hits, file = 'Z:/Jeremy/GOBACK/R outputs/goback.cox.ph.top.hits.v20180124.1.csv', row.names = FALSE)



# Filter down to specific BD-CC pairs (revised 20180131) ------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' THIS SECTION IS DEPRECATED.  USE THE VERSION REVISED 02/23/2018.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.01.31.
#' 
#' Some changes based on the meeting with Philip.
#' 
#' DO exclude septal defects.  
#' DO NOT exclude RVOT or LVOT defects.
#' DO NOT exclude non-RMS soft tissues sarcomas (i.e., 'soft.other').
#' 
#' Two sensitivity analyses:
#'  1. Adjust for BW in models where congenital hip dislocation is the
#'     exposure.
#'  2. Exclude children diagnosed before 1st birthday in models where
#'     hydrocephalus.wo.spina.bifida is the exposure.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

#' setwd('Z:/Jeremy/GOBACK/Datasets/') 

#' load('goback.cox.ph.results.v20180124.1.rdata')

#' Compute Bonferroni-Holm adjusted p-values.
#' goback.coxmodels <- arrange(goback.coxmodels, p.val.coef)
#' goback.coxmodels$j.rev <- 606:1
#' goback.coxmodels$p.val.coef.bon <- 0.05/goback.coxmodels$j.rev
#' goback.coxmodels$p.val.coef.bon.delta <- goback.coxmodels$p.val.coef - goback.coxmodels$p.val.coef.bon
#' 371 results considered significant after adjustment.
#' goback.coxmodels <- goback.coxmodels[1:371, 1:13]

#' pat1 <- 'conganomalies.'; pat2 <- 'chromosomalanomalies'; pat3 <- '.any'; pat4 <- '.other'; 
#' pat5 <- 'septal.defects'; pat6 <- 'soft.other'

#' goback.coxmodels <- goback.coxmodels[!grepl(pat1, goback.coxmodels$defect), ]
#' goback.coxmodels <- goback.coxmodels[!grepl(pat2, goback.coxmodels$defect), ]
#' goback.coxmodels <- goback.coxmodels[!grepl(pat3, goback.coxmodels$cancer), ]

#' Grab soft.other models, set them aside, them bind them back in.
#' tmp <- goback.coxmodels[grepl(pat6, goback.coxmodels$cancer), ]
#' goback.coxmodels <- goback.coxmodels[!grepl(pat4, goback.coxmodels$cancer), ]
#' goback.coxmodels <- rbind(goback.coxmodels, tmp)

#' goback.coxmodels <- goback.coxmodels[!grepl(pat5, goback.coxmodels$defect), ]

#' top.hits <- arrange(goback.coxmodels, defect) 
#' save(top.hits, file = 'goback.cox.ph.top.hits.v20180131.1.rdata')
#' write.csv(top.hits, file = 'Z:/Jeremy/GOBACK/R outputs/goback.cox.ph.top.hits.v20180131.1.csv', row.names = FALSE)

#' rm(top.hits, goback.coxmodels, pat1, pat2, pat3, pat4, pat5, pat6, tmp); gc()



# Filter down to specific BD-CC pairs (revised 20180223) ------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.02.23.
#' 
#' Philip wants to change the algorithm for selecting top results.
#' 1. Filter out non-specific models.
#' 2. Perform B-H adjustment only on the set of specific models.
#' 
#' DO exclude septal defects.  
#' DO NOT exclude RVOT or LVOT defects.
#' DO NOT exclude non-RMS soft tissues sarcomas (i.e., 'soft.other').
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

setwd('Z:/Jeremy/GOBACK/Datasets/') 

load('goback.cox.ph.results.v20180124.1.rdata')

pat1 <- 'conganomalies.'; pat2 <- 'chromosomalanomalies'; pat3 <- '.any'; pat4 <- '.other'; 
pat5 <- 'septal.defects'; pat6 <- 'soft.other'

goback.coxmodels <- goback.coxmodels[!grepl(pat1, goback.coxmodels$defect), ]
goback.coxmodels <- goback.coxmodels[!grepl(pat2, goback.coxmodels$defect), ]
goback.coxmodels <- goback.coxmodels[!grepl(pat3, goback.coxmodels$cancer), ]

#' Grab soft.other models, set them aside, them bind them back in.
tmp <- goback.coxmodels[grepl(pat6, goback.coxmodels$cancer), ]
goback.coxmodels <- goback.coxmodels[!grepl(pat4, goback.coxmodels$cancer), ]
goback.coxmodels <- rbind(goback.coxmodels, tmp)

goback.coxmodels <- goback.coxmodels[!grepl(pat5, goback.coxmodels$defect), ]

#' Compute Bonferroni-Holm adjusted p-values.
goback.coxmodels <- arrange(goback.coxmodels, p.val.coef)
goback.coxmodels$rank <- 1:74
goback.coxmodels$p.val.coef.bon <- 0.05/(74-goback.coxmodels$rank+1)
goback.coxmodels$p.val.coef.bon.delta <- goback.coxmodels$p.val.coef - goback.coxmodels$p.val.coef.bon
goback.coxmodels <- filter(goback.coxmodels, p.val.coef.bon.delta < 0)


top.hits <- arrange(goback.coxmodels, defect) 
save(top.hits, file = 'goback.cox.ph.top.hits.v20180223.1.rdata')
write.csv(top.hits, file = 'Z:/Jeremy/GOBACK/R outputs/goback.cox.ph.top.hits.v20180214.1.csv', row.names = FALSE)

rm(goback.coxmodels, top.hits, pat1, pat2, pat3, pat4, pat5, pat6, tmp); gc()



# Flag new BD-CC associations ---------------------------------------------

setwd('Z:/Jeremy/GOBACK/Datasets/') 

load("goback.cox.ph.top.hits.v20180223.1.rdata")
tmp <- top.hits

load('./Old Datasets/goback.cox.ph.top.hits.v20180131.1.rdata')

tmp2 <- left_join(select(tmp, defect, cancer, HR),
                  select(top.hits, defect, cancer, HR),
                  by = c('defect','cancer'))
tmp2$new <- ifelse(is.na(tmp2$HR.y), 1, 0)
tmp2 <- c(tmp2$new)

load("goback.cox.ph.top.hits.v20180223.1.rdata")

top.hits$is.new <- tmp2

rm(tmp, tmp2)

save(top.hits, file = 'goback.cox.ph.top.hits.v20180223.2.rdata')



# Sensitivity analyses: hip dislocation -----------------------------------

setwd('Z:/Jeremy/GOBACK/Datasets/') 

load('goback.no.chrom.v20180122.1.rdata')

goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                          cancer = goback.nochrom[,'gct.extra'],
                          defect = goback.nochrom[,'congenital.hip.dislocation'],
                          sex = factor(goback.nochrom$sex,
                                       levels = c(1,2),
                                       labels = c('Male','Female')),
                          m.age = goback.nochrom$m.age,
                          state = goback.nochrom$state.num,
                          birth.wt = goback.nochrom$birth.wt)      


cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state + birth.wt, data = goback.surv)
cox.coef <- summary(cox)$coefficients
test.ph <- cox.zph(cox)
test.ph <- test.ph$table['defect','p']

rm(goback.nochrom, cox, cox.coef, test.ph); gc()

# Sensitivity analyses: hydrocephalus -------------------------------------

setwd('Z:/Jeremy/GOBACK/Datasets/') 

load('goback.no.chrom.v20180122.1.rdata')

for (i in c('astro','ependymoma','epithe','nephro')){
  
  goback.nochrom[,i] <- ifelse(is.na(goback.nochrom[,i]), 0, goback.nochrom[,i])
  
  print(i); print(table(goback.nochrom[,i], useNA = 'always'))
  
  goback.nochrom <- filter(goback.nochrom, goback.nochrom[,i] == 0 | (goback.nochrom[,i] == 1 & person.yrs >= 1))
  
  goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                            cancer = goback.nochrom[,i],
                            defect = goback.nochrom[,'hydrocephalus.wo.spinabifida'],
                            sex = factor(goback.nochrom$sex,
                                         levels = c(1,2),
                                         labels = c('Male','Female')),
                            m.age = goback.nochrom$m.age,
                            state = goback.nochrom$state.num)   
  
  cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
  print(summary(cox))
  
  rm(cox, goback.surv); gc()
  
}













goback.nochrom <- filter(goback.nochrom, goback.nochrom[,i] == 0 | (goback.nochrom[,i] == 1 & person.yrs >= 1))

for (i in c('astro','ependymoma','epithe','nephro')){
  
  goback.nochrom[,i] <- ifelse(is.na(goback.nochrom[,i]), 0, goback.nochrom[,i])
  
}

for (i in c('astro','ependymoma','epithe','nephro')){
  
goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                          cancer = goback.nochrom[,i],
                          defect = goback.nochrom[,'hydrocephalus.wo.spinabifida'],
                          sex = factor(goback.nochrom$sex,
                                       levels = c(1,2),
                                       labels = c('Male','Female')),
                          m.age = goback.nochrom$m.age,
                          state = goback.nochrom$state.num)      


cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state + birth.wt, data = goback.surv)
cox.coef <- summary(cox)$coefficients
test.ph <- cox.zph(cox)
test.ph <- test.ph$table['defect','p']

}

# Sensitivity analyses: risk of any cancer given any defect ---------------

setwd('Z:/Jeremy/GOBACK/Datasets/')

load('goback.no.chrom.v20180122.1.rdata')

goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                          cancer = goback.nochrom$cancer,
                          defect = goback.nochrom$any.birthdefect,
                          sex = factor(goback.nochrom$sex,
                                       levels = c(1,2),
                                       labels = c('Male','Female')),
                          m.age = goback.nochrom$m.age,
                          state = goback.nochrom$state.num)

tab <- as.numeric(table(goback.nochrom$any.birthdefect, goback.nochrom$cancer)[2,2])

cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
cox.coef <- summary(cox)$coefficients
test.ph <- cox.zph(cox)
test.ph <- test.ph$table['defect','p']

estimates <- data.frame(defect = 'any.nonchromosomal.defect', 
                        cancer = 'any.cancer', 
                        HR = exp(cox.coef[1,1]), 
                        ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                        ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                        p.value.coef = cox.coef[1,5],
                        p.value.zph = test.ph,
                        num.comorbid = tab)

write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.coxph.models.csv', sep=',', append = TRUE, 
            row.names = FALSE, col.names = FALSE)

rm(cox.coef, estimates, goback.nochrom, goback.surv, cox, tab, test.ph)

load('goback.chrom.v20180122.1.rdata')

goback.surv <- data.frame(time = goback.chrom$person.yrs,
                          cancer = goback.chrom$cancer,
                          defect = goback.chrom$any.birthdefect,
                          sex = factor(goback.chrom$sex,
                                       levels = c(1,2),
                                       labels = c('Male','Female')),
                          m.age = goback.chrom$m.age,
                          state = goback.chrom$state.num)

tab <- as.numeric(table(goback.chrom$any.birthdefect, goback.chrom$cancer)[2,2])

cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
cox.coef <- summary(cox)$coefficients
test.ph <- cox.zph(cox)
test.ph <- test.ph$table['defect','p']

estimates <- data.frame(defect = 'any.nonchromosomal.defect', 
                        cancer = 'any.cancer', 
                        HR = exp(cox.coef[1,1]), 
                        ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                        ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                        p.value.coef = cox.coef[1,5],
                        p.value.zph = test.ph,
                        num.comorbid = tab)

write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.coxph.models.csv', sep=',', append = TRUE, 
            row.names = FALSE, col.names = FALSE)

rm(cox.coef, estimates, goback.chrom, goback.surv, cox, tab, test.ph)




# Sensitivity analyses: recreate table 3 for TX  --------------------------

setwd('Z:/Jeremy/GOBACK/Datasets/')

load('goback.no.chrom.v20180122.1.rdata')

goback.nochrom <- filter(goback.nochrom, state == 'TX'); gc()

tmp <- table(goback.nochrom[,16], goback.nochrom$cancer1)
tmp <- which(tmp[2, ] >= 5)
tmp <- tmp + 107

bw.adj <- c(108, 122, 127) # Column indices for ALL, hepatoblastoma and nephroblastoma.

for (j in tmp){
  
  tab <- as.numeric(table(goback.nochrom$any.birthdefect, goback.nochrom[,j])[2,2])

  if (j %in% bw.adj){

  goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                            cancer = goback.nochrom[,j],
                            defect = goback.nochrom$any.birthdefect,
                            sex = factor(goback.nochrom$sex,
                                         levels = c(1,2),
                                         labels = c('Male','Female')),
                            m.age = goback.nochrom$m.age,
                            birth.wt = goback.nochrom$birth.wt)

  cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + birth.wt, data = goback.surv)
  cox.coef <- summary(cox)$coefficients
  test.ph <- cox.zph(cox)
  test.ph <- test.ph$table['defect','p']

  }
  
  else{
    
    goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                              cancer = goback.nochrom[,j],
                              defect = goback.nochrom$any.birthdefect,
                              sex = factor(goback.nochrom$sex,
                                           levels = c(1,2),
                                           labels = c('Male','Female')),
                              m.age = goback.nochrom$m.age)
    
    cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex, data = goback.surv)
    cox.coef <- summary(cox)$coefficients
    test.ph <- cox.zph(cox)
    test.ph <- test.ph$table['defect','p']

  }
  
  estimates <- data.frame(defect = 'any.nonchromosomal.defect', 
                          cancer = names(goback.nochrom[j]), 
                          HR = exp(cox.coef[1,1]), 
                          ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                          ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                          p.value.coef = cox.coef[1,5],
                          p.value.zph = test.ph,
                          num.comorbid = tab,
                          state = 'TX')
  
  write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/statewise.coxph.models.csv', sep=',', append = TRUE, 
              row.names = FALSE, col.names = FALSE)

}

rm(bw.adj, cox.coef, estimates, goback.surv, goback.nochrom, cox, j, tab, test.ph, tmp); gc()

load('goback.chrom.v20180122.1.rdata')

goback.chrom <- filter(goback.chrom, state == 'TX'); gc()

tmp <- table(goback.chrom[,16], goback.chrom$cancer1)
tmp <- which(tmp[2, ] >= 5)
tmp <- tmp + 107

bw.adj <- c(108, 122, 127) # Column indices for ALL, hepatoblastoma and nephroblastoma.

for (j in tmp){
  
  tab <- as.numeric(table(goback.chrom$any.birthdefect, goback.chrom[,j])[2,2])
  
  if (j %in% bw.adj){
    
    goback.surv <- data.frame(time = goback.chrom$person.yrs,
                              cancer = goback.chrom[,j],
                              defect = goback.chrom$any.birthdefect,
                              sex = factor(goback.chrom$sex,
                                           levels = c(1,2),
                                           labels = c('Male','Female')),
                              m.age = goback.chrom$m.age,
                              birth.wt = goback.chrom$birth.wt)
    
    cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + birth.wt, data = goback.surv)
    cox.coef <- summary(cox)$coefficients
    test.ph <- cox.zph(cox)
    test.ph <- test.ph$table['defect','p']
    
  }
  
  else{
    
    goback.surv <- data.frame(time = goback.chrom$person.yrs,
                              cancer = goback.chrom[,j],
                              defect = goback.chrom$any.birthdefect,
                              sex = factor(goback.chrom$sex,
                                           levels = c(1,2),
                                           labels = c('Male','Female')),
                              m.age = goback.chrom$m.age)
    
    cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex, data = goback.surv)
    cox.coef <- summary(cox)$coefficients
    test.ph <- cox.zph(cox)
    test.ph <- test.ph$table['defect','p']
    
  }
  
  estimates <- data.frame(defect = 'any.chromosomal.defect', 
                          cancer = names(goback.chrom[j]), 
                          HR = exp(cox.coef[1,1]), 
                          ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                          ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                          p.value.coef = cox.coef[1,5],
                          p.value.zph = test.ph,
                          num.comorbid = tab,
                          state = 'TX')
  
  write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/statewise.coxph.models.csv', sep=',', append = TRUE, 
              row.names = FALSE, col.names = FALSE)
  
}

rm(cox.coef, estimates, goback.surv, goback.chrom, cox, j, tab, test.ph, tmp); gc()



# Sensitivity analyses: recreate table 3 for MI ---------------------------

setwd('Z:/Jeremy/GOBACK/Datasets/')

load('goback.no.chrom.v20180122.1.rdata')

goback.nochrom <- filter(goback.nochrom, state == 'MI'); gc()

tmp <- table(goback.nochrom[,16], goback.nochrom$cancer1)
tmp <- which(tmp[2, ] >= 5)
tmp <- tmp + 107

bw.adj <- c(108, 122, 127) # Column indices for ALL, hepatoblastoma and nephroblastoma.

for (j in tmp){
  
  tab <- as.numeric(table(goback.nochrom$any.birthdefect, goback.nochrom[,j])[2,2])
  
  if (j %in% bw.adj){
    
    goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                              cancer = goback.nochrom[,j],
                              defect = goback.nochrom$any.birthdefect,
                              sex = factor(goback.nochrom$sex,
                                           levels = c(1,2),
                                           labels = c('Male','Female')),
                              m.age = goback.nochrom$m.age,
                              birth.wt = goback.nochrom$birth.wt)
    
    cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + birth.wt, data = goback.surv)
    cox.coef <- summary(cox)$coefficients
    test.ph <- cox.zph(cox)
    test.ph <- test.ph$table['defect','p']
    
  }
  
  else{
    
    goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                              cancer = goback.nochrom[,j],
                              defect = goback.nochrom$any.birthdefect,
                              sex = factor(goback.nochrom$sex,
                                           levels = c(1,2),
                                           labels = c('Male','Female')),
                              m.age = goback.nochrom$m.age)
    
    cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex, data = goback.surv)
    cox.coef <- summary(cox)$coefficients
    test.ph <- cox.zph(cox)
    test.ph <- test.ph$table['defect','p']
    
  }
  
  estimates <- data.frame(defect = 'any.nonchromosomal.defect', 
                          cancer = names(goback.nochrom[j]), 
                          HR = exp(cox.coef[1,1]), 
                          ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                          ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                          p.value.coef = cox.coef[1,5],
                          p.value.zph = test.ph,
                          num.comorbid = tab,
                          state = 'MI')
  
  write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/statewise.coxph.models.csv', sep=',', append = TRUE, 
              row.names = FALSE, col.names = FALSE)
  
}

rm(bw.adj, cox.coef, estimates, goback.surv, goback.nochrom, cox, j, tab, test.ph, tmp); gc()

load('goback.chrom.v20180122.1.rdata')

goback.chrom <- filter(goback.chrom, state == 'MI'); gc()

tmp <- table(goback.chrom[,16], goback.chrom$cancer1)
tmp <- which(tmp[2, ] >= 5)
tmp <- tmp + 107

bw.adj <- c(108, 122, 127) # Column indices for ALL, hepatoblastoma and nephroblastoma.

for (j in tmp){
  
  tab <- as.numeric(table(goback.chrom$any.birthdefect, goback.chrom[,j])[2,2])
  
  if (j %in% bw.adj){
    
    goback.surv <- data.frame(time = goback.chrom$person.yrs,
                              cancer = goback.chrom[,j],
                              defect = goback.chrom$any.birthdefect,
                              sex = factor(goback.chrom$sex,
                                           levels = c(1,2),
                                           labels = c('Male','Female')),
                              m.age = goback.chrom$m.age,
                              birth.wt = goback.chrom$birth.wt)
    
    cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + birth.wt, data = goback.surv)
    cox.coef <- summary(cox)$coefficients
    test.ph <- cox.zph(cox)
    test.ph <- test.ph$table['defect','p']
    
  }
  
  else{
    
    goback.surv <- data.frame(time = goback.chrom$person.yrs,
                              cancer = goback.chrom[,j],
                              defect = goback.chrom$any.birthdefect,
                              sex = factor(goback.chrom$sex,
                                           levels = c(1,2),
                                           labels = c('Male','Female')),
                              m.age = goback.chrom$m.age)
    
    cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex, data = goback.surv)
    cox.coef <- summary(cox)$coefficients
    test.ph <- cox.zph(cox)
    test.ph <- test.ph$table['defect','p']
    
  }
  
  estimates <- data.frame(defect = 'any.chromosomal.defect', 
                          cancer = names(goback.chrom[j]), 
                          HR = exp(cox.coef[1,1]), 
                          ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                          ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                          p.value.coef = cox.coef[1,5],
                          p.value.zph = test.ph,
                          num.comorbid = tab,
                          state = 'MI')
  
  write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/statewise.coxph.models.csv', sep=',', append = TRUE, 
              row.names = FALSE, col.names = FALSE)
  
}

rm(cox.coef, estimates, goback.surv, goback.chrom, cox, j, tab, test.ph, tmp); gc()

# Sensitivity analyses: AnyBD-AnyCC by state, total data ------------------

setwd('Z:/Jeremy/GOBACK/Datasets/')
load('goback.v20180216.1.rdata')

for (j in unique(goback$state)){
  
  tmp <- filter(goback, state == j)
  
  tab <- as.numeric(table(tmp$any.birthdefect, tmp$cancer)[2,2])
  
  goback.surv <- data.frame(time = tmp$person.yrs,
                            cancer = tmp$cancer,
                            defect = tmp$any.birthdefect,
                            sex = factor(tmp$sex,
                                         levels = c(1,2),
                                         labels = c('Male','Female')),
                            m.age = tmp$m.age)
  
  cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex, data = goback.surv)
  cox.coef <- summary(cox)$coefficients
  test.ph <- cox.zph(cox)
  test.ph <- test.ph$table['defect','p']
  
  estimates <- data.frame(defect = 'any.defect', 
                          cancer = 'any.cancer', 
                          HR = exp(cox.coef[1,1]), 
                          ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                          ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                          p.value.coef = cox.coef[1,5],
                          p.value.zph = test.ph,
                          num.comorbid = tab,
                          state = j)
  
  write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/any.defect.any.cancer.by.state.and.chromstatus.csv', sep=',', append = TRUE, 
              row.names = FALSE, col.names = FALSE)
  
  rm(goback.surv, cox, cox.coef, test.ph, estimates, tmp); gc()
}

rm(j, tab, goback); gc()

# Sensitivity analyses: AnyBD-Any CC by state, non-chrom ------------------

setwd('Z:/Jeremy/GOBACK/Datasets/')
load('goback.no.chrom.v20180122.1.rdata')

for (j in unique(goback.nochrom$state)){
  
  tmp <- filter(goback.nochrom, state == j)
  
  tab <- as.numeric(table(tmp$any.birthdefect, tmp$cancer)[2,2])
  
  goback.surv <- data.frame(time = tmp$person.yrs,
                            cancer = tmp$cancer,
                            defect = tmp$any.birthdefect,
                            sex = factor(tmp$sex,
                                         levels = c(1,2),
                                         labels = c('Male','Female')),
                            m.age = tmp$m.age)
  
  cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex, data = goback.surv)
  cox.coef <- summary(cox)$coefficients
  test.ph <- cox.zph(cox)
  test.ph <- test.ph$table['defect','p']

  estimates <- data.frame(defect = 'any.nonchromosomal.defect', 
                          cancer = 'any.cancer', 
                          HR = exp(cox.coef[1,1]), 
                          ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                          ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                          p.value.coef = cox.coef[1,5],
                          p.value.zph = test.ph,
                          num.comorbid = tab,
                          state = j)
  
  write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/any.defect.any.cancer.by.state.and.chromstatus.csv', sep=',', append = TRUE, 
              row.names = FALSE, col.names = FALSE)
  
  rm(goback.surv, cox, cox.coef, test.ph, estimates, tmp); gc()
}

rm(goback.surv, tmp, estimates, cox.coef, cox, j, tab, test.ph, goback.nochrom); gc()



# Sensitivity analyses: AnyBD-Any CC by state, chrom ------------------

setwd('Z:/Jeremy/GOBACK/Datasets/')
load('goback.chrom.v20180122.1.rdata')

for (j in unique(goback.chrom$state)){
  
  tmp <- filter(goback.chrom, state == j)
  
  tab <- as.numeric(table(tmp$any.birthdefect, tmp$cancer)[2,2])
  
  goback.surv <- data.frame(time = tmp$person.yrs,
                            cancer = tmp$cancer,
                            defect = tmp$any.birthdefect,
                            sex = factor(tmp$sex,
                                         levels = c(1,2),
                                         labels = c('Male','Female')),
                            m.age = tmp$m.age)
  
  cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex, data = goback.surv)
  cox.coef <- summary(cox)$coefficients
  test.ph <- cox.zph(cox)
  test.ph <- test.ph$table['defect','p']
  
  estimates <- data.frame(defect = 'any.chromosomal.defect', 
                          cancer = 'any.cancer', 
                          HR = exp(cox.coef[1,1]), 
                          ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                          ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                          p.value.coef = cox.coef[1,5],
                          p.value.zph = test.ph,
                          num.comorbid = tab,
                          state = j)
  
  write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/any.defect.any.cancer.by.state.and.chromstatus.csv', sep=',', append = TRUE, 
              row.names = FALSE, col.names = FALSE)
  
  rm(goback.surv, cox, cox.coef, test.ph, estimates, tmp); gc()
}

rm(j, tab, goback.chrom); gc()
# As above, but with logistic regression ----------------------------------

setwd('Z:/Jeremy/GOBACK/Datasets/')
load('goback.no.chrom.v20180122.1.rdata')

for (j in unique(goback.nochrom$state)){
  
  tmp <- filter(goback.nochrom, state == j)
  
  tab <- as.numeric(table(tmp$any.birthdefect, tmp$cancer)[2,2])
  
  mod <- glm(cancer ~ any.birthdefect + sex + m.age, data = tmp, family = binomial(link = 'logit'))
  coef <- summary(mod)$coefficients
  
  estimates <- data.frame(defect = 'any.nonchromosomal.defect', 
                          cancer = 'any.cancer', 
                          HR = exp(coef[2,1]), 
                          ci.lower = exp(coef[2,1]-(1.96*coef[2,2])), 
                          ci.upper = exp(coef[2,1]+(1.96*coef[2,2])),
                          p.value.coef = coef[2,4],
                          num.comorbid = tab,
                          state = j)
  
  write.table(estimates, file = 'C:/Users/schraw/Desktop/any.defect.any.cancer.by.state.and.chromstatus.logreg.csv', sep=',', append = TRUE, 
              row.names = FALSE, col.names = FALSE)
  
  rm(coef, estimates, tmp, mod); gc()
}

rm(goback.nochrom, j, tab); gc()

load('goback.chrom.v20180122.1.rdata')

for (j in unique(goback.chrom$state)){
  
  tmp <- filter(goback.chrom, state == j)
  
  tab <- as.numeric(table(tmp$any.birthdefect, tmp$cancer)[2,2])
  
  mod <- glm(cancer ~ any.birthdefect + sex + m.age, data = tmp, family = binomial(link = 'logit'))
  coef <- summary(mod)$coefficients
  
  estimates <- data.frame(defect = 'any.chromosomal.defect', 
                          cancer = 'any.cancer', 
                          HR = exp(coef[2,1]), 
                          ci.lower = exp(coef[2,1]-(1.96*coef[2,2])), 
                          ci.upper = exp(coef[2,1]+(1.96*coef[2,2])),
                          p.value.coef = coef[2,4],
                          num.comorbid = tab,
                          state = j)
  
  write.table(estimates, file = 'C:/Users/schraw/Desktop/any.defect.any.cancer.by.state.and.chromstatus.logreg.csv', sep=',', append = TRUE, 
              row.names = FALSE, col.names = FALSE)
  
  rm(coef, estimates, tmp, mod); gc()
}

rm(goback.chrom, j, tab); gc()

load('goback.v20180216.1.rdata')

for (j in unique(goback$state)){
  
  tmp <- filter(goback, state == j)
  
  tab <- as.numeric(table(tmp$any.birthdefect, tmp$cancer)[2,2])
  
  mod <- glm(cancer ~ any.birthdefect + sex + m.age, data = tmp, family = binomial(link = 'logit'))
  coef <- summary(mod)$coefficients
  
  estimates <- data.frame(defect = 'any.defect', 
                          cancer = 'any.cancer', 
                          HR = exp(coef[2,1]), 
                          ci.lower = exp(coef[2,1]-(1.96*coef[2,2])), 
                          ci.upper = exp(coef[2,1]+(1.96*coef[2,2])),
                          p.value.coef = coef[2,4],
                          num.comorbid = tab,
                          state = j)
  
  write.table(estimates, file = 'C:/Users/schraw/Desktop/any.defect.any.cancer.by.state.and.chromstatus.logreg.csv', sep=',', append = TRUE, 
              row.names = FALSE, col.names = FALSE)
  
  rm(coef, estimates, tmp, mod); gc()
}

rm(goback.chrom, j, tab); gc()



