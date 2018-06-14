#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' This is a graveyard for various obsolete models and sensitivity analyses
#' that were performed on earlier versions of the data.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------



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

require(dplyr); require(survival)

load('Z:/Jeremy/GOBACK/Datasets/goback.nochrom.v20180611.rdata')

tmp <- filter(goback.nochrom, cancer == 0)
tmp2 <- filter(goback.nochrom, cancer == 1 & (cancer1 == 'all' | cancer1 == 'hepato' | cancer1 == 'nephro'))

rm(goback.nochrom); gc()
goback.filter <- rbind(tmp, tmp2); rm(tmp, tmp2); gc()

#' Remove unnecessary columns. 
#' This task is causing repeated issues with RAM usage.
goback.filter <- goback.filter[,c(1,156,3,6,7,15,22:94,112,126,131)]

for (i in 7:79){
  
  for (j in 80:82){
    
    tab <- table(goback.filter[,i], goback.filter[,j])[2,2]
    
    if (tab >= 5){
      
      goback.surv <- data.frame(time = goback.filter$person.yrs,
                                cancer = goback.filter[,j],
                                defect = goback.filter[,i],
                                sex = factor(goback.filter$sex,
                                             levels = c(1,2),
                                             labels = c('Male','Female')),
                                m.age = goback.filter$m.age,
                                state = goback.filter$state.num,
                                birth.wt = goback.filter$birth.wt)
      
      cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state + birth.wt, data = goback.surv)
      cox.coef <- summary(cox)$coefficients
      #'    test.ph <- cox.zph(cox)
      #'    test.ph <- test.ph$table['defect','p']
      rm(cox, goback.surv); gc()
      
      estimates <- data.frame(defect = names(goback.filter[i]), 
                              cancer = names(goback.filter[j]), 
                              HR = exp(cox.coef[1,1]), 
                              ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                              ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                              p.value.coef = cox.coef[1,5],
                              #' p.value.zph = test.ph,
                              num.comorbid = tab)
      
      write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.coxph.models.v20180612.csv', sep=',', append = TRUE, 
                  row.names = FALSE, col.names = FALSE)
      
    }
    
    else{
      
      next
      
    }
  }
}

rm(list = ls()); gc()



# BW-adjusted models for Wilms/ALL/hepatoblastoma: chromosomal ------------

require(dplyr); require(survival)

load('Z:/Jeremy/GOBACK/Datasets/goback.chrom.v20180611.rdata')

tmp <- filter(goback.chrom, cancer == 0)
tmp2 <- filter(goback.chrom, cancer == 1 & (cancer1 == 'all' | cancer1 == 'hepato' | cancer1 == 'nephro'))

rm(goback.chrom); gc()

goback.filter <- rbind(tmp, tmp2); rm(tmp, tmp2); gc()

goback.filter <- goback.filter[,c(1,156,3,6,7,15,97:104,112,126,131)]

for (i in 7:13){
  
  for (j in 14:16){
    
    tab <- table(goback.filter[,i], goback.filter[,j])[2,2]
    
    if (tab >= 5){
      
      goback.surv <- data.frame(time = goback.filter$person.yrs,
                                cancer = goback.filter[,j],
                                defect = goback.filter[,i],
                                sex = factor(goback.filter$sex,
                                             levels = c(1,2),
                                             labels = c('Male','Female')),
                                m.age = goback.filter$m.age,
                                state = goback.filter$state.num,
                                birth.wt = goback.filter$birth.wt)
      
      cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state + birth.wt, data = goback.surv)
      cox.coef <- summary(cox)$coefficients
      test.ph <- cox.zph(cox)
      test.ph <- test.ph$table['defect','p']
      rm(cox, goback.surv); gc()
      
      estimates <- data.frame(defect = names(goback.filter[i]), 
                              cancer = names(goback.filter[j]), 
                              HR = exp(cox.coef[1,1]), 
                              ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                              ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                              p.value.coef = cox.coef[1,5],
                              p.value.zph = test.ph,
                              num.comorbid = tab)
      
      write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/goback.coxph.models.v20180612.csv', sep=',', append = TRUE, 
                  row.names = FALSE, col.names = FALSE)
      
    }
    
    else{
      
      next
      
    }
  }
}

rm(list = ls()); gc()

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
#' 
#' IT IS RETAINED HERE FOR RECORD-KEEPING PURPOSES.
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
#' 
#' IT IS RETAINED HERE FOR RECORD-KEEPING PURPOSES.
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

setwd('Z:/Jeremy/GOBACK/Datasets/') 

load('goback.cox.ph.results.v20180124.1.rdata')

#' Compute Bonferroni-Holm adjusted p-values.
goback.coxmodels <- arrange(goback.coxmodels, p.val.coef)
goback.coxmodels$j.rev <- 606:1
goback.coxmodels$p.val.coef.bon <- 0.05/goback.coxmodels$j.rev
goback.coxmodels$p.val.coef.bon.delta <- goback.coxmodels$p.val.coef - goback.coxmodels$p.val.coef.bon
#' 371 results considered significant after adjustment.
goback.coxmodels <- goback.coxmodels[1:371, 1:13]

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

top.hits <- arrange(goback.coxmodels, defect) 
save(top.hits, file = 'goback.cox.ph.top.hits.v20180131.1.rdata')
write.csv(top.hits, file = 'Z:/Jeremy/GOBACK/R outputs/goback.cox.ph.top.hits.v20180131.1.csv', row.names = FALSE)

rm(top.hits, goback.coxmodels, pat1, pat2, pat3, pat4, pat5, pat6, tmp); gc()




# Flag new BD-CC associations ---------------------------------------------

setwd('Z:/Jeremy/GOBACK/Datasets/') 

load("goback.cox.ph.top.hits.v20180612.rdata")
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




# eTable 2 ----------------------------------------------------------------

#' Need to re-run model for risk of hepatoblastoma among kids with ASD in
#' non-chromosomal set, adjusting for BW.
load('./goback.no.chrom.v20180122.1.rdata')


goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                          cancer = goback.nochrom$hepato,
                          defect = goback.nochrom$atrialseptaldefect,
                          sex = factor(goback.nochrom$sex,
                                       levels = c(1,2),
                                       labels = c('Male','Female')),
                          m.age = goback.nochrom$m.age,
                          state = goback.nochrom$state,
                          birth.wt = goback.nochrom$birth.wt/100)

cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state + birth.wt, data = goback.surv)
summary(cox)
