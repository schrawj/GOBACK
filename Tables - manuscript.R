#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.06.11.
#' 
#' Generate Tables 1-4 and Figure 1 for GOBACK manuscript.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

# Table 1 -----------------------------------------------------------------

require(gmodels)

setwd('Z:/Jeremy/GOBACK/Datasets/')
load('goback.v20180611.rdata')

table(goback$state, goback$birth.yr, useNA = 'ifany')
CrossTable(goback$state, prop.chisq = FALSE)
CrossTable(goback$state, goback$any.birthdefect, prop.chisq = FALSE)
CrossTable(goback$state, goback$cancer, prop.chisq = FALSE)
CrossTable(goback$any.birthdefect, goback$cancer, prop.chisq = FALSE)

for (i in unique(goback$state)){
  tmp <- filter(goback, state == i)
  print(i)
  CrossTable(tmp$cancer, tmp$any.birthdefect, prop.chisq = FALSE)
  rm(i, tmp)
}



# Table 2 -----------------------------------------------------------------

setwd('Z:/Jeremy/GOBACK/Datasets/')
load('goback.v20180611.rdata')

#' Pare down the datasaet to help with performance.
goback <- goback[, c(1:16,107)]

#' Collapse plurality to singleton vs multiple.
goback$plu.cat <- factor(ifelse(goback$plu > 1, 1, 0),
                         levels = c(0,1),
                         labels = c('singleton','multiple'))

#' The simplified version of the table will compare children only by birth defects status.
for (i in c(3,4,17,18)){
  print(names(goback[i]))
  print(gmodels::CrossTable(goback[,i], goback$any.birthdefect, prop.t = FALSE, prop.chisq = FALSE, prop.r = FALSE, chisq = TRUE))
}

for (i in c(7,9,6)){
  print(names(goback[i]))
  print(aggregate(goback[,i] ~ goback$any.birthdefect, data = goback, mean))
  print(aggregate(goback[,i] ~ goback$any.birthdefect, data = goback, sd))
  print(t.test(goback[,i] ~ goback$any.birthdefect, data = goback, na.rm = TRUE))
}

rm(list = ls()); gc()



# Table 4: chromosomal and genetic conditions -----------------------------

require(survival)

setwd('Z:/Jeremy/GOBACK/Datasets/')
load('goback.chrom.v20180611.rdata')

goback.chrom <- goback.chrom[, c(1:21,95:156)]

#' Any anomaly, any cancer.
goback.surv <- data.frame(time = goback.chrom$person.yrs,
                          cancer = goback.chrom$cancer,
                          defect = goback.chrom$any.chromosomal.anomaly,
                          sex = factor(goback.chrom$sex,
                                       levels = c(1,2),
                                       labels = c('Male','Female')),
                          m.age = goback.chrom$m.age,
                          state = goback.chrom$state.num)      

cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
cox.coef <- summary(cox)$coefficients

rm(cox, goback.surv); gc()

estimates <- data.frame(defect = 'any.chromosomal.anomaly', 
                        cancer = 'any.cancer', 
                        HR = exp(cox.coef[1,1]), 
                        ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                        ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                        p.value.coef = cox.coef[1,5],
                        num.comorbid = table(goback.chrom$any.chromosomal.anomaly, goback.chrom$cancer)[2,2])

write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.specific.cancer.by.any.chromosomal.defect.v20180611.csv', sep=',', append = TRUE, 
            row.names = FALSE, col.names = FALSE)

#' Models for cancers except ALL, Wilms, hepatoblastoma in children with 
#' chromosomal anomalies or single-gene syndromes.
for (j in c(40:52,54:57,59:78)){
  
  for (i in 22:31){
  
    tab <- table(goback.chrom[,i], goback.chrom[,j])[2,2]
    
    if (tab >= 5){
      
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
      test.ph <- cox.zph(cox)
      test.ph <- test.ph$table['defect','p']
      
      rm(cox, goback.surv); gc()
      
      estimates <- data.frame(defect = names(goback.chrom[i]), 
                              cancer = names(goback.chrom[j]), 
                              HR = exp(cox.coef[1,1]), 
                              ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                              ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                              p.value.coef = cox.coef[1,5],
                              p.value.zph = test.ph,
                              num.comorbid = tab)
      
      write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.specific.cancer.by.any.chromosomal.defect.v20180611.csv', sep=',', append = TRUE, 
                  row.names = FALSE, col.names = FALSE)
    }
    
    else{
      
      next
      
    }
    
  }
  
}

rm(cox.coef, estimates, test.ph, tab, i, j); gc()

#' Models for ALL, Wilms, hepatoblastoma, in children with 
#' chromosomal anomalies or single gene syndromes.  
for (j in c(39,53,58)){
  
  for (i in 22:31){
  
    tab <- table(goback.chrom[,i], goback.chrom[,j])[2,2]
    
    if (tab >= 5){
      
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
      test.ph <- cox.zph(cox)
      test.ph <- test.ph$table['defect','p']
      
      rm(cox, goback.surv); gc()
      
      estimates <- data.frame(defect = names(goback.chrom[i]), 
                              cancer = names(goback.chrom[j]), 
                              HR = exp(cox.coef[1,1]), 
                              ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                              ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                              p.value.coef = cox.coef[1,5],
                              p.value.zph = test.ph,
                              num.comorbid = tab)
      
      write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.specific.cancer.by.any.chromosomal.defect.v20180611.csv', sep=',', append = TRUE, 
                  row.names = FALSE, col.names = FALSE)
    }
    
    else{
      
      next
      
    }
    
  }
  
}

rm(list = ls()); gc()



# Table 4: Non-chromosomal defects ----------------------------------------

require(survival)

setwd('Z:/Jeremy/GOBACK/Datasets/')
load('goback.nochrom.v20180611.rdata')

#' Any anomaly, any cancer.
goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                          cancer = goback.nochrom$cancer,
                          defect = goback.nochrom$any.birthdefect,
                          sex = factor(goback.nochrom$sex,
                                       levels = c(1,2),
                                       labels = c('Male','Female')),
                          m.age = goback.nochrom$m.age,
                          state = goback.nochrom$state.num)      

cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
cox.coef <- summary(cox)$coefficients

test.ph <- cox.zph(cox)
test.ph <- test.ph$table['defect','p']

rm(cox, goback.surv); gc()

estimates <- data.frame(defect = 'any.nonchromosomal.anomaly', 
                        cancer = 'any.cancer', 
                        HR = exp(cox.coef[1,1]), 
                        ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                        ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                        p.value.coef = cox.coef[1,5],
                        p.value.zph = test.ph,
                        num.comorbid = table(goback.nochrom$any.birthdefect, goback.nochrom$cancer)[2,2])

write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.specific.cancer.by.any.nonchromosomal.defect.v20180611.csv', sep=',', append = TRUE, 
            row.names = FALSE, col.names = FALSE)

goback.nochrom <- goback.nochrom[,c(1,156,3,6,7,9,15,16,112:151)]; gc()

#' Models for cancers except ALL, Wilms, hepatoblastoma in children with non-chromosomal structural birth defects.
for (j in c(10:22,24:27,29:48)){
  
  tab <- table(goback.nochrom[,'any.birthdefect'], goback.nochrom[,j])[2,2]
  
  if (tab >= 5){
    
    goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                              cancer = goback.nochrom[,j],
                              defect = goback.nochrom$any.birthdefect,
                              sex = factor(goback.nochrom$sex,
                                           levels = c(1,2),
                                           labels = c('Male','Female')),
                              m.age = goback.nochrom$m.age,
                              state = goback.nochrom$state.num)      
    
    cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
    cox.coef <- summary(cox)$coefficients
    test.ph <- cox.zph(cox)
    test.ph <- test.ph$table['defect','p']
    
    rm(cox, goback.surv); gc()
    
    estimates <- data.frame(defect = 'any.nonchromosomal.defect', 
                            cancer = names(goback.nochrom[j]), 
                            HR = exp(cox.coef[1,1]), 
                            ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                            ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                            p.value.coef = cox.coef[1,5],
                            p.value.zph = test.ph,
                            num.comorbid = tab)
    
    write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.specific.cancer.by.any.nonchromosomal.defect.v20180611.csv', sep=',', append = TRUE, 
                row.names = FALSE, col.names = FALSE)
  }
  
  else{
    
    next
    
  }
  
}

rm(cox.coef, estimates, tab, test.ph, j); gc()

#' Models for ALL, Wilms and hepatoblastoma in children with non-chromosomal structural birth defects.
for (j in c(9,23,28)){
  
  tab <- table(goback.nochrom[,'any.birthdefect'], goback.nochrom[,j])[2,2]
  
  if (tab >= 5){
    
    goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                              cancer = goback.nochrom[,j],
                              defect = goback.nochrom$any.birthdefect,
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
    
    rm(cox, goback.surv); gc()
    
    estimates <- data.frame(defect = 'any.nonchromosomal.defect', 
                            cancer = names(goback.nochrom[j]), 
                            HR = exp(cox.coef[1,1]), 
                            ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                            ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                            p.value.coef = cox.coef[1,5],
                            p.value.zph = test.ph,
                            num.comorbid = tab)
    
    write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.specific.cancer.by.any.nonchromosomal.defect.v20180611.csv', sep=',', append = TRUE, 
                row.names = FALSE, col.names = FALSE)
  }
  
  else{
    
    next
    
  }
  
}

rm(list = ls()); gc()



# Table 3: Chromosomal defects --------------------------------------------

require(survival)

setwd('Z:/Jeremy/GOBACK/Datasets/')
load('goback.chrom.v20180611.rdata')

for (i in 95:104){
  
    tab <- table(goback.chrom[,i], goback.chrom$cancer)[2,2]

    if (tab >= 5){
      
      goback.surv <- data.frame(time = goback.chrom$person.yrs,
                                cancer = goback.chrom$cancer,
                                defect = goback.chrom[,i],
                                sex = factor(goback.chrom$sex,
                                             levels = c(1,2),
                                             labels = c('Male','Female')),
                                m.age = goback.chrom$m.age,
                                state = goback.chrom$state.num)      
      
      cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
      cox.coef <- summary(cox)$coefficients
      test.ph <- cox.zph(cox)
      test.ph <- test.ph$table['defect','p']
      
      rm(cox, goback.surv); gc()
      
      estimates <- data.frame(defect = names(goback.chrom[i]), 
                              cancer = 'any.cancer', 
                              HR = exp(cox.coef[1,1]), 
                              ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                              ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                              p.value.coef = cox.coef[1,5],
                              p.value.zph = test.ph,
                              num.comorbid = tab)
      
      write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.any.cancer.by.chromosomal.defect.v20180611.csv', sep=',', append = TRUE, 
                  row.names = FALSE, col.names = FALSE)
    }
    
    else{
      
          next
      
        }

}

rm(list = ls()); gc()



# Table 3: Non-chromosomal defects ----------------------------------------

require(survival)

setwd('Z:/Jeremy/GOBACK/Datasets/')
load('goback.nochrom.v20180611.rdata')

goback.nochrom <- goback.nochrom[,c(1,156,3,6,7,9,15,22:94,107)]; gc()

for (i in 8:80){
  
  tab <- table(goback.nochrom[,i], goback.nochrom$cancer)[2,2]
  
  if (tab >= 5){
    
    goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                              cancer = goback.nochrom$cancer,
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
    
    rm(cox, goback.surv); gc()
    
    estimates <- data.frame(defect = names(goback.nochrom[i]), 
                            cancer = 'any.cancer', 
                            HR = exp(cox.coef[1,1]), 
                            ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                            ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                            p.value.coef = cox.coef[1,5],
                            p.value.zph = test.ph,
                            num.comorbid = tab)
    
    write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.any.cancer.by.nonchromosomal.defect.v20180611.csv', sep=',', append = TRUE, 
                row.names = FALSE, col.names = FALSE)
  }
  
  else{
    
    next
    
  }
  
}

rm(list = ls()); gc()



# Reformat output for convenience -----------------------------------------

#' Reformat tables 3, 4 and 5 to allow easy copy-paste of HRs and CIs into a Word document.
tablenames <- c('defect','cancer','hr','ci.lower','ci.upper','p.val.coef','p.val.zph','n.comorbid')

table3part1 <- read.csv(file='Z:/Jeremy/GOBACK/R outputs/goback.any.cancer.by.chromosomal.defect.v20180611.csv', 
                        header = FALSE, stringsAsFactors = FALSE)
table3part2 <- read.csv(file='Z:/Jeremy/GOBACK/R outputs/goback.any.cancer.by.nonchromosomal.defect.v20180611.csv', 
                        header = FALSE, stringsAsFactors = FALSE)
table3 <- rbind(table3part1, table3part2)

names(table3) <- tablenames
for (i in 3:5){
  table3[,i] <- round(table3[,i], 1)
}
table3$hr.ci <- paste0(table3$hr, ' (',table3$ci.lower,'-',table3$ci.upper,')')
table3 <- table3[,c(1,8,9)]

write.csv(table3, file = 'Z:/Jeremy/GOBACK/Tables/table3.raw.v20180611.csv', row.names = FALSE)

table4part1 <- read.csv(file='Z:/Jeremy/GOBACK/R outputs/goback.specific.cancer.by.any.chromosomal.defect.v20180611.csv', 
                        header = FALSE, stringsAsFactors = FALSE)
table4part2 <- read.csv(file='Z:/Jeremy/GOBACK/R outputs/goback.specific.cancer.by.any.nonchromosomal.defect.v20180611.csv', 
                        header = FALSE, stringsAsFactors = FALSE)
table4 <- rbind(table4part1, table4part2)

names(table4) <- tablenames
for (i in 3:5){
  table4[,i] <- round(table4[,i], 1)
}
table4$hr.ci <- paste0(table4$hr, ' (',table4$ci.lower,'-',table4$ci.upper,')')
table4 <- table4[,c(1,2,8,9)]

write.csv(table4, file = 'Z:/Jeremy/GOBACK/Tables/table4.raw.v20180611.csv', row.names = FALSE)

load(file = 'Z:/Jeremy/GOBACK/Datasets/Expanded datasets/goback.coxph.top.hits.v20180612.rdata')

table5 <- data.frame(defect = top.hits$defect, cancer = top.hits$cancer, n.comorbid = top.hits$n.comorbid, 
                     hr.ci = paste0(round(top.hits$hr, 1),' (',round(top.hits$ci.lower, 1),'-',round(top.hits$ci.upper, 1),')'))
write.csv(table5, file = 'Z:/Jeremy/GOBACK/Tables/table5.raw.v20180611.csv', row.names = FALSE)

rm(list = ls()); gc()



# Figure 1: Body system defects and body system cancers -------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Philip wants the heatmap to be at the level of the body system defects
#' variables (e.g., 'any CNS anomaly') and body system level cancer (e.g.,
#' 'any CNS cancer').  This is a subset of the results already generated 
#' by the Cox models (see the script 'Cox regression models - GOBACK data) 
#' and can be astracted from them.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

load("Z:/Jeremy/GOBACK/Datasets/Expanded datasets/goback.coxph.results.v20180612.rdata")
#' TODO: patterns may need to change; or alternatively, rename variables to comply with patterns.
pat1 <- 'conganomalies.'; pat2 <- 'chromosomalanomalies'
pat3 <- '.any'; pat4 <- '.other'

tmp <- goback.coxmodels[grepl(pat1, goback.coxmodels$defect) & grepl(pat3, goback.coxmodels$cancer), ]
tmp <- tmp[!grepl(pat4, tmp$defect), ]
tmp2 <- goback.coxmodels[grepl(pat2, goback.coxmodels$defect) & grepl(pat3, goback.coxmodels$cancer), ]
tmp2 <- tmp2[!grepl(pat4, tmp2$defect), ]
tmp <- rbind(tmp, tmp2); rm(tmp2,pat1, pat2, pat3, pat4); gc()

heatmap <- tmp
save(heatmap, file = 'Z:/Jeremy/GOBACK/Datasets/figure1.v20180126.1.rdata')



# Etable 5 ----------------------------------------------------------------

require(survival); require(dplyr)

#' Risk of any cancer among children with specific BDs.  
#' Only Texas. Only cases DX'd 1 year of age or greater.
load('Z:/Jeremy/GOBACK/Datasets/goback.nochrom.v20180611.rdata')

goback.nochrom$exclude <- ifelse(goback.nochrom$cancer == 1 & goback.nochrom$person.yrs < 1, 1, 0)
goback.nochrom <- filter(filter(goback.nochrom, state == 'TX'),
                         exclude == 0)

goback.nochrom <- goback.nochrom[,c(1,156,3,6,7,9,15,16,112:151)]; gc()

for (i in 8:80){
  
  tab <- table(goback.nochrom[,i], goback.nochrom$cancer)[2,2]
  
  if (tab >= 1){
    
    goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                              cancer = goback.nochrom$cancer,
                              defect = goback.nochrom[,i],
                              sex = factor(goback.nochrom$sex,
                                           levels = c(1,2),
                                           labels = c('Male','Female')),
                              m.age = goback.nochrom$m.age)      
    
    cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex, data = goback.surv)
    cox.coef <- summary(cox)$coefficients

    rm(cox, goback.surv); gc()
    
    estimates <- data.frame(defect = names(goback.nochrom[i]), 
                            cancer = 'any.cancer', 
                            HR = exp(cox.coef[1,1]), 
                            ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                            ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                            p.value.coef = cox.coef[1,5],
                            num.comorbid = tab)
    
    write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/etable5.v20180611.csv', sep=',', append = TRUE, 
                row.names = FALSE, col.names = FALSE)
  }
  
  else{
    
    next
    
  }
  
}

rm(list = ls()); gc()

load('goback.chrom.v20180611.rdata')

goback.chrom$exclude <- ifelse(goback.chrom$cancer == 1 & goback.chrom$person.yrs < 1, 1, 0)
goback.chrom <- filter(filter(goback.chrom, state == 'TX'),
                       exclude == 0)

for (i in 95:101){
  
  tab <- table(goback.chrom[,i], goback.chrom$cancer)[2,2]
  
  if (tab >= 1){
    
    goback.surv <- data.frame(time = goback.chrom$person.yrs,
                              cancer = goback.chrom$cancer,
                              defect = goback.chrom[,i],
                              sex = factor(goback.chrom$sex,
                                           levels = c(1,2),
                                           labels = c('Male','Female')),
                              m.age = goback.chrom$m.age)      
    
    cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex, data = goback.surv)
    cox.coef <- summary(cox)$coefficients
    
    rm(cox, goback.surv); gc()
    
    estimates <- data.frame(defect = names(goback.chrom[i]), 
                            cancer = 'any.cancer', 
                            HR = exp(cox.coef[1,1]), 
                            ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                            ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                            p.value.coef = cox.coef[1,5],
                            num.comorbid = tab)
    
    write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/etable5.csv', sep=',', append = TRUE, 
                row.names = FALSE, col.names = FALSE)
  }
  
  else{
    
    next
    
  }
  
}

rm(list = ls()); gc()

#' Let's load the output back in and cbind it to the original table 4.
setwd('Z:/Jeremy/GOBACK/R outputs/')

#' Have rounded HR and CI columns to 2 decimal places in Excel for convenience.
tab.sens <- read.csv(file = './eTable5.csv', header = TRUE, stringsAsFactors = FALSE)
tab.og.nochrom <- read.csv(file = './table 4 nochrom.csv', header = TRUE, stringsAsFactors = FALSE)
tab.og.chrom <- read.csv(file = './table 4 chrom.csv', header = TRUE, stringsAsFactors = FALSE)

tab.og <- select(rbind(tab.og.chrom, tab.og.nochrom), -p.val.zph); rm(tab.og.chrom, tab.og.nochrom)

etable5 <- rename(left_join(tab.og, tab.sens, by = c('defect','cancer')),
                            num.comorbid.sensitivity = n.comorbid)

etable5$hr <- with(etable5, paste0(hr.x,' (',ci.lower.x,'-',ci.upper.x,')'))
etable5$hr.sensitivity <- with(etable5, paste0(hr.y,' (',ci.lower.y,'-',ci.upper.y,')'))
etable5$change.flag <- ifelse((etable5$hr.x/etable5$hr.y <= 0.5 | etable5$hr.x/etable5$hr.y >= 1.5) | 
                              (etable5$ci.lower.x > 1 & etable5$ci.lower.y < 1), 1, 0)

tmp <- select(etable5, defect, cancer, hr, hr.sensitivity, change.flag, num.comorbid, num.comorbid.sensitivity)
write.csv(tmp, file = 'Z:/Jeremy/GOBACK/R outputs/eTable5.sidebyside.csv', row.names = FALSE)
