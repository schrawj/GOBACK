#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.01.26.
#' 
#' Generate Tables 1-4 and Figure 1 for GOBACK manuscript.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

# Prep environment --------------------------------------------------------
require(dplyr); require(gmodels); require(ggplot2); require(survival)

setwd('Z:/Jeremy/GOBACK/Datasets/')



# Table 1 -----------------------------------------------------------------

load('goback.v20180125.1.rdata')

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

#' CMH test for sex, cancer, defect status.
tmp <- xtabs(cancer ~ sex + any.birthdefect, data = goback)
ftable(tmp)
mantelhaen.test(goback$cancer, goback$sex, goback$any.birthdefect)

tmp <- data.frame(sex = rep(c(1,2),2), cancer = rep(c(0,1),each = 2), birth.defect = rep(c(0,1), each = 4),
                  count = c(4880855,4747666,6872,6153,319649,218115,1189,938))
tmp.xtab <- xtabs(count ~ sex + cancer + birth.defect, data = tmp)
mantelhaen.test(tmp.xtab)

# Table 2 -----------------------------------------------------------------

load('goback.v20180125.1.rdata')

#' Pare down the datasaet to help with performance.
goback <- goback[, c(1:16,103)]

with(subset(goback, any.birthdefect == 0), CrossTable(sex, cancer, prop.chisq = FALSE, chisq = TRUE))
with(subset(goback, any.birthdefect == 1), CrossTable(sex, cancer, prop.chisq = FALSE, chisq = TRUE))

goback$multi.birth <- ifelse(goback$plu %in% c(2:8), 1, 
                             ifelse(goback$plu == 1, 0, 9))

with(subset(goback, any.birthdefect == 0 & multi.birth != 9), CrossTable(multi.birth, cancer, prop.chisq = FALSE, chisq = TRUE))
with(subset(goback, any.birthdefect == 1 & multi.birth != 9), CrossTable(multi.birth, cancer, prop.chisq = FALSE, chisq = TRUE))

aggregate(birth.wt ~ any.birthdefect + cancer, data = goback, mean) 
aggregate(birth.wt ~ any.birthdefect + cancer, data = goback, sd)
p <- ggplot(data = goback) + geom_histogram(aes(x = birth.wt), color = 'blue', fill = 'white', binwidth = 100) + facet_wrap(~any.birthdefect)
print(p)
t.test( birth.wt ~ cancer, data = subset(goback, any.birthdefect == 0), alternative = 'two.sided') 
t.test( birth.wt ~ cancer, data = subset(goback, any.birthdefect == 1), alternative = 'two.sided')

aggregate(gest.age ~ any.birthdefect + cancer, data = goback, mean); aggregate(gest.age ~ any.birthdefect + cancer, data = goback, sd)
p <- ggplot(data = goback) + geom_histogram(aes(x = gest.age), color = 'blue', fill = 'white', binwidth = 1) + facet_wrap(~any.birthdefect)
print(p)
t.test( gest.age ~ cancer, data = subset(goback, any.birthdefect == 0), alternative = 'two.sided') 
t.test( gest.age ~ cancer, data = subset(goback, any.birthdefect == 1), alternative = 'two.sided')

aggregate(m.age ~ any.birthdefect + cancer, data = goback, mean); aggregate(m.age ~ any.birthdefect + cancer, data = goback, sd)
p <- ggplot(data = goback) + geom_histogram(aes(x = m.age), color = 'blue', fill = 'white', binwidth = 1) + facet_wrap(~any.birthdefect)
print(p)
t.test( m.age ~ cancer, data = subset(goback, any.birthdefect == 0), alternative = 'two.sided') 
t.test( m.age ~ cancer, data = subset(goback, any.birthdefect == 1), alternative = 'two.sided')

with(subset(goback, any.birthdefect == 0), CrossTable(m.race, cancer, prop.chisq = FALSE, chisq = TRUE))
with(subset(goback, any.birthdefect == 1), CrossTable(m.race, cancer, prop.chisq = FALSE, chisq = TRUE))

with(subset(goback, any.birthdefect == 0), CrossTable(m.edu2, cancer, prop.chisq = FALSE, chisq = TRUE))
with(subset(goback, any.birthdefect == 1), CrossTable(m.edu2, cancer, prop.chisq = FALSE, chisq = TRUE))
# Table 3: Chromosomal defects --------------------------------------------

load('goback.chrom.v20180122.1.rdata')

goback.chrom <- goback.chrom[, c(1:21,95:152)]

#' Generate models for all cancers except ALL, Wilm's, hepatoblastoma.
#' These will additionally need to be adjusted for birthweight.
for (j in c(36:48,50:53,55:74)){
  
  tab <- table(goback.chrom[,22], goback.chrom[,j])[2,2]
  
  if (tab >= 5){
    
    goback.surv <- data.frame(time = goback.chrom$person.yrs,
                              cancer = goback.chrom[,j],
                              defect = goback.chrom$chromosomalanomalies,
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
    
    estimates <- data.frame(defect = 'chromosomalanomalies', 
                            cancer = names(goback.chrom[j]), 
                            HR = exp(cox.coef[1,1]), 
                            ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                            ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                            p.value.coef = cox.coef[1,5],
                            p.value.zph = test.ph,
                            num.comorbid = tab)
    
    write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.specific.cancer.by.any.chromosomal.defect.csv', sep=',', append = TRUE, 
                row.names = FALSE, col.names = FALSE)
  }
  
  else{
    
    next
    
  }
  
}

rm(cox.coef, estimates, test.ph, tab, i, j); gc()

#' Models for ALL, Wilm's, hepatoblastoma.  
for (j in c(35,49,54)){
  
  tab <- table(goback.chrom[,22], goback.chrom[,j])[2,2]
  
  if (tab >= 5){
    
    goback.surv <- data.frame(time = goback.chrom$person.yrs,
                              cancer = goback.chrom[,j],
                              defect = goback.chrom$chromosomalanomalies,
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
    
    estimates <- data.frame(defect = 'chromosomalanomalies', 
                            cancer = names(goback.chrom[j]), 
                            HR = exp(cox.coef[1,1]), 
                            ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                            ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                            p.value.coef = cox.coef[1,5],
                            p.value.zph = test.ph,
                            num.comorbid = tab)
    
    write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.specific.cancer.by.any.chromosomal.defect.csv', sep=',', append = TRUE, 
                row.names = FALSE, col.names = FALSE)
  }
  
  else{
    
    next
    
  }
  
}

rm(goback.chrom, cox.coef, estimates, test.ph, tab, i, j); gc()



# Table 3: Non-chromosomal defects ----------------------------------------

load('goback.nochrom.v20180419.rdata')

goback.nochrom <- goback.nochrom[,c(1:3,6,7,9,15,16,107:147,152)]; gc()

for (j in c(11:23,25:28,30:49)){
  
  tab <- table(goback.nochrom[,8], goback.nochrom[,j])[2,2]
  
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
    
    write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.specific.cancer.by.any.nonchromosomal.defect.csv', sep=',', append = TRUE, 
                row.names = FALSE, col.names = FALSE)
  }
  
  else{
    
    next
    
  }
  
}

rm(cox.coef, estimates, tab, test.ph, j); gc()

#' ALL, Wilm's, hepatoblastoma.
for (j in c(10,24,29)){
  
  tab <- table(goback.nochrom[,8], goback.nochrom[,j])[2,2]
  
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
    
    write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.specific.cancer.by.any.nonchromosomal.defect.csv', sep=',', append = TRUE, 
                row.names = FALSE, col.names = FALSE)
  }
  
  else{
    
    next
    
  }
  
}

rm(cox.coef, estimates, goback.nochrom, j, tab, test.ph); gc()



# Table 4: Chromosomal defects --------------------------------------------

load('goback.chrom.v20180122.1.rdata')

for (i in 95:101){
  
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
      
      write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.any.cancer.by.chromosomal.defect.csv', sep=',', append = TRUE, 
                  row.names = FALSE, col.names = FALSE)
    }
    
    else{
      
          next
      
        }

}

rm(cox.coef, estimates, i, tab, test.ph, goback.chrom); gc()



# Table 4: Non-chromosomal defects ----------------------------------------

load('goback.no.chrom.v20180122.1.rdata')

goback.nochrom <- goback.nochrom[,c(1:3,6,7,9,15,22:94,152,103)]; gc()

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
    
    write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.any.cancer.by.nonchromosomal.defect.csv', sep=',', append = TRUE, 
                row.names = FALSE, col.names = FALSE)
  }
  
  else{
    
    next
    
  }
  
}

rm(cox.coef, estimates, i, tab, test.ph, goback.nochrom); gc()

#' This is a long list.  Let's convert it into a format that requires less typing to transpose to a word doc.
results <- read.csv(file='Z:/Jeremy/GOBACK/R outputs/goback.any.cancer.by.nonchromosomal.defect.csv', 
                    header = TRUE, stringsAsFactors = FALSE)
results$hr.ci <- paste0(results$HR, ' (',results$ci.lower,'-',results$ci.upper,')')
results <- results[,c(1,8,9)]
write.csv(results, file = 'C:/Users/schraw/desktop/results.csv', row.names = FALSE)

# Figure 1: Body system defects and body system cancers -------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Philip wants the heatmap to be at the level of the body system defects
#' variables (e.g., 'any CNS anomaly') and body system level cancer (e.g.,
#' 'any CNS cancer').  This is a subset of the results already generated 
#' by the Cox models and can be astracted from them.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

load("Z:/Jeremy/GOBACK/Datasets/goback.cox.ph.results.v20180124.1.rdata")

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

#' Risk of any cancer among children with specific BDs.  Only in TX, only 
#' cases DX'd 1 year of age or greater.
load('goback.nochrom.v20180419.rdata')

goback.nochrom$exclude <- ifelse(goback.nochrom$cancer == 1 & goback.nochrom$person.yrs < 1, 1, 0)
goback.nochrom <- filter(goback.nochrom, state == 'TX')
goback.nochrom <- filter(goback.nochrom, exclude == 0)
goback.nochrom <- goback.nochrom[,c(1:3,6,7,9,15,22:94,152,103)]; gc()

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
    
    write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/etable5.csv', sep=',', append = TRUE, 
                row.names = FALSE, col.names = FALSE)
  }
  
  else{
    
    next
    
  }
  
}

rm(list = ls()); gc()

load('goback.chrom.v20180419.rdata')

goback.chrom$exclude <- ifelse(goback.chrom$cancer == 1 & goback.chrom$person.yrs < 1, 1, 0)
goback.chrom <- filter(goback.chrom, state == 'TX')
goback.chrom <- filter(goback.chrom, exclude == 0)

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


# Scratch paper -----------------------------------------------------------

#' This is a long list.  Let's convert it into a format that requires less typing to transpose to a word doc.
results <- read.csv(file='Z:/Jeremy/GOBACK/R outputs/goback.any.cancer.by.nonchromosomal.defect.csv', 
                    header = TRUE, stringsAsFactors = FALSE)
results$hr.ci <- paste0(results$HR, ' (',results$ci.lower,'-',results$ci.upper,')')
results <- results[,c(1,8,9)]
write.csv(results, file = 'C:/Users/schraw/desktop/results.csv', row.names = FALSE)