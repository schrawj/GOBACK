#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.08.23.
#' 
#' Contains the code used to generate the R Markdown documents describing 
#' QC measures in response to co-author comments, plus some additional 
#' code to output data in other formats.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

# Review BD codes for cancer cases w/22q or trisomy 13 --------------------

#' Sonja Rasmussen asked to see birth defects codes for kids with 22q11.2 deletions and cancer.
#' We also need to go back in time to retrieve their birth defects registry ID for possible
#' review by Angela Scheurle, the TBDR medical geneticist.
require(dplyr); require(kableExtra); require(xlsx)

load('Z:/Jeremy/GOBACK/Datasets/goback.v20180711.rdata')
load('Z:/Jeremy/GOBACK/Datasets/Old Datasets/Texas/tx.raw.data.rdata')
load('Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.v20180606.rdata')
load('Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.v20180712.rdata')


#' Export TX kids only to an Excel file, and include the CASE_ID variable from the raw data.
tx.raw <- select(tx.raw, birthID, CASE_ID)
tx.raw$studyid <- paste0('tx',tx.raw$birthID)

#' Pull IDs for kids with 22q or 13q and cancer, then export their BPA codes for Sonja.
comorbid.22q <- filter(goback, del.22q == 1 & cancer == 1 & state == 'TX')
comorbid.13q <- filter(goback, del.13q == 1 & cancer == 1 & state == 'TX')

comorbid.13q <- left_join(filter(bd.codes.txnc, studyid %in% comorbid.13q$studyid), tx.raw, 'studyid')
comorbid.13q <- comorbid.13q[, c(1,69,2:67)]
comorbid.22q <- left_join(filter(bd.codes.txnc, studyid %in% comorbid.22q$studyid), tx.raw, 'studyid')
comorbid.22q <- comorbid.22q[, c(1,69,2:67)]

write.xlsx(comorbid.13q, file = 'W:/Old_genepi2/Jeremy/GOBACK/R outputs/22q.and.13q.comorbid.cases.xlsx', sheetName = '13q-cancer cases', showNA = FALSE, row.names = FALSE)
write.xlsx(comorbid.22q, file = 'W:/Old_genepi2/Jeremy/GOBACK/R outputs/22q.and.13q.comorbid.cases.xlsx', sheetName = '22q-cancer cases', showNA = FALSE, row.names = FALSE,
           append = TRUE)

#' Philip had originally requested I re-do the speciifc cancer association analyses for kids with chromosomally confirmed diagnoses.
#' That actually won't be productive, as we have so few.

rm(list = ls()); gc()

#check <- select(rbind(filter(goback.chrom, di.george.syndrome == 1 & cancer == 1 & state != 'AR'),
#                      filter(goback.chrom, trisomy13 == 1 & cancer == 1 & state != 'AR')),
#                studyid, cancer1, di.george.syndrome, trisomy13)

#check.bpa <- select(filter(bd.codes.txnc, studyid %in% check$studyid), 1:32)
#check.icd <- filter(bd.codes.mi, studyid %in% check$studyid)

#check.icd[26:32] <- as.character(NA)

#names <- c('studyid',rep(paste0('code',1:31)))
#names(check.bpa) <- names 
#names(check.icd) <- names

#check <- arrange(left_join(check, rbind(check.bpa, check.icd), by = 'studyid'), studyid)

#kable(check) %>% kable_styling() %>% scroll_box(width = '100%', height = '400px')



# Review specific BD codes for kids w/obstructive GU defects --------------

require(dplyr); require(xlsx)

load('Z:/Jeremy/GOBACK/Datasets/goback.v20180711.rdata')
load('Z:/Jeremy/GOBACK/Datasets/Old Datasets/Texas/tx.raw.data.rdata')
load('Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.transpose.v20180614.rdata')
load('Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.v20180606.rdata')

ob <- filter(goback, obstructive.genitourinary.defects == 1 & state %in% c('TX','NC'))
ob.bpa <- filter(bd.codes.txnc.transpose, studyid %in% ob$studyid)

cols <- c(1, which(grepl('753.2', names(ob.bpa))))

ob.bpa <- select(ob.bpa, cols)
ob.bpa$`753.200` <- ifelse(rowSums(ob.bpa[2:7], na.rm = TRUE) >= 1, 1, 0)
ob.bpa$`753.210` <- ifelse(rowSums(ob.bpa[8:12], na.rm = TRUE) >= 1, 1, 0)
ob.bpa$`753.220` <- ifelse(rowSums(ob.bpa[13:18], na.rm = TRUE) >= 1, 1, 0)
ob.bpa$`753.290` <- ifelse(rowSums(ob.bpa[19:24], na.rm = TRUE) >= 1, 1, 0)
ob.bpa <- ob.bpa[, c(1,2,8,13,20)]

gu.defect.counts <- data.frame(code = names(ob.bpa[2:5]),
                               count = colSums(ob.bpa[2:5]))

gu.names <- read.xlsx('Z:/Jeremy/GOBACK/Data dictionaries and data definitions/dxmap_current.xlsx', sheetName = 'Current')

gu.defect.counts <- left_join(gu.defect.counts, 
                              select(gu.names, BPA_number, BPANAME), 
                              by = c('code' = 'BPA_number'))

kable(gu.defect.counts) %>% kable_styling() %>% scroll_box(width = '100%')

ob <- filter(goback, obstructive.genitourinary.defects == 1 & cancer == 1 & state == 'TX')

tab <- table(ob$cancer1)

gu.cancers <- data.frame(cancer = names(tab), count = tab)
gu.cancers <- rename(select(gu.cancers, -count.Var1), cases = count.Freq)

ob <- left_join(select(ob, studyid, cancer1), bd.codes.txnc, by = 'studyid')

write.csv(ob, file = 'Z:/Jeremy/GOBACK/R outputs/bd.codes.for.tx.obstructive.gu.defects.cancer.cases.csv')

rm(list = ls()); gc()



# Confirm that no non-chromosomal CHD cases have Alagille -----------------

#' Do any of the children with pulmonary valve atresia/stenosis and hepatoblastoma have Alagille?

require(dplyr)

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.transpose.v20180614.rdata')
load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.v20180606.rdata')
load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.nochrom.v20180711.rdata')

which(grepl('759.87', names(bd.codes.txnc.transpose)))

bd.codes.txnc.transpose <- bd.codes.txnc.transpose[, c(1,2318,2319)]

possible.alagille.ids <- filter(bd.codes.txnc.transpose, `759.870` == 1 | `759.878`== 1)
possible.alagille.ids <- c(possible.alagille.ids$studyid)

#' I remember that we elected not to remove kids with this code from the non-chromosomal dataset
#' because it was unclear what condition they might have based on their codes, the conditions in general are 
#' rare and they are not all cancer-predisposing.
#' 1 of 5 pulmonary valve atresia-hepatoblastoma cases has a code that MAY indicate Alagille.
tmp <- filter(goback.nochrom, studyid %in% possible.alagille.ids)

tmp2 <- filter(goback.nochrom, hepato == 1 & pulmvalveatresiaandstenosis == 1)
tmp2$studyid %in% possible.alagille.ids

#' The BD codes in that particular child are suggestive of Alagille: 
#' multiple forms of CHD, hepatosplenomegaly, craniofacial features.
#' What would happen if we dropped that child?
require(survival)

rm(bd.codes.txnc, bd.codes.txnc.transpose, tmp, tmp2, possible.alagille.ids)

goback.nochrom <- goback.nochrom[c(1:237109,237111:10155647), ]

goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                          cancer = goback.nochrom$hepato,
                          defect = goback.nochrom$pulmvalveatresiaandstenosis,
                          sex = factor(goback.nochrom$sex,
                                       levels = c(1,2),
                                       labels = c('Male','Female')),
                          m.age = goback.nochrom$m.age,
                          state = goback.nochrom$state.num,
                          birth.wt = goback.nochrom$birth.wt)

cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state + birth.wt, data = goback.surv)
summary(cox)$coefficients


# Pick out some MBD kids at random ----------------------------------------

require(dplyr); require(kableExtra)

load('Z:/Jeremy/GOBACK/Datasets/goback.nochrom.v20180711.rdata')
load('Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.v20180606.rdata')
load('Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.v20180712.rdata')

tmp <- arrange(filter(goback.nochrom, defect.total >= 4 & state %in% c('TX','NC')), runif)
tmp <- tmp[30000:30099,]
mbd.codes <- filter(bd.codes.txnc, studyid %in% tmp$studyid)
kable(mbd.codes) %>% kable_styling() %>% scroll_box(width = '100%', height = '400px')

tmp <- arrange(filter(goback.nochrom, defect.total >= 4 & state == 'MI', runif))
tmp <- tmp[10900:10999,]
mbd.codes <- filter(bd.codes.mi, studyid %in% tmp$studyid)
kable(mbd.codes) %>% kable_styling() %>% scroll_box(width = '100%', height = '400px')



# Re-run PDA models adjusting for BW and GA -------------------------------

require(dplyr); require(survival)

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.nochrom.v20180829.rdata')

associations <- read.csv(file = 'W:/Old_genepi2/Jeremy/GOBACK/R outputs/goback.cox.ph.top.hits.v20180612.csv',
                         header = TRUE, stringsAsFactors = FALSE)

associations <- filter(associations, defect %in% c('atrialseptaldefect','patentductusarteriosis','ventricularseptaldefect'))

for (i in 1:nrow(associations)){
  
  index.defect <- associations$defect[i]
  index.cancer <- associations$cancer[i]
  
  print(paste('Computing BW- and GA-adjusted model for risk of',index.cancer, 'in kids with', index.defect))
  
  goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                            cancer = goback.nochrom[,index.cancer],
                            defect = goback.nochrom[,index.defect],
                            sex = factor(goback.nochrom$sex,
                                         levels = c(1,2),
                                         labels = c('Male','Female')),
                            m.age = goback.nochrom$m.age,
                            state = goback.nochrom$state.num,
                            birth.wt = goback.nochrom$birth.wt,
                            gest.age = goback.nochrom$gest.age)
  
  cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state + birth.wt + gest.age, data = goback.surv)
  cox.coef <- summary(cox)$coefficients

  estimates <- data.frame(defect = index.defect, 
                          cancer = index.cancer, 
                          HR = exp(cox.coef[1,1]), 
                          ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                          ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                          p.value.coef = cox.coef[1,5],
                          num.comorbid = as.numeric(table(goback.nochrom[,index.defect], goback.nochrom[,index.cancer])[2,2]))
  
  write.table(estimates, file = 'W:/Old_genepi2/Jeremy/GOBACK/R Outputs/goback.coxph.models.for.selected.chd.bw.and.ga.adjusted.v20180829.csv', 
              sep=',', append = TRUE, row.names = FALSE, col.names = FALSE)
}

write.table(associations[,1:7], file = 'W:/Old_genepi2/Jeremy/GOBACK/R Outputs/goback.coxph.models.for.selected.chd.v20180829.csv', 
            sep=',', append = TRUE, row.names = FALSE)

# Re-do all models for specific associations ------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' This time, we will restrict all analyses to Texas children older than 1
#' at cancer DX.
#' 
#' Wel will also adjust all models for PDA, ASD, and VSD for birthweight
#' and gestational age.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(dplyr); require(survival)

chrom.syndromes <- c('di.george.syndrome','down.syndrome','nf','trisomy18')

associations <- read.csv(file = 'W:/Old_genepi2/Jeremy/GOBACK/R outputs/goback.cox.ph.top.hits.v20180612.csv',
                         header = TRUE, stringsAsFactors = FALSE)
nonchrom.associations <- filter(associations, !(defect %in% chrom.syndromes))
chrom.associations <- filter(filter(associations, defect %in% chrom.syndromes), 
                             defect != 'di.george.syndrome') #' Old DiGeorge variable is moot.

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.nochrom.v20180829.rdata')
goback.nochrom <- filter(goback.nochrom, state %in% c('TX','NC') & person.yrs > 1)

for (i in 1:nrow(nonchrom.associations)){
  
  index.defect <- nonchrom.associations$defect[i]
  index.cancer <- nonchrom.associations$cancer[i]
  
  if (index.defect %in% c('atrialseptaldefect','patentductusarteriosis','ventricularseptaldefect')){
  
    print(paste('Computing BW- and GA-adjusted model for risk of',index.cancer, 'in kids with', index.defect))
  
    goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                              cancer = goback.nochrom[,index.cancer],
                              defect = goback.nochrom[,index.defect],
                              sex = factor(goback.nochrom$sex,
                                           levels = c(1,2),
                                           labels = c('Male','Female')),
                              m.age = goback.nochrom$m.age,
                              birth.wt = goback.nochrom$birth.wt,
                              gest.age = goback.nochrom$gest.age,
                              state = goback.nochrom$state.num)
    
    cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + birth.wt + gest.age + state, data = goback.surv)
    cox.coef <- summary(cox)$coefficients
    
    estimates <- data.frame(defect = index.defect, 
                            cancer = index.cancer, 
                            HR = exp(cox.coef[1,1]), 
                            ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                            ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                            p.value.coef = cox.coef[1,5],
                            num.comorbid = as.numeric(table(goback.nochrom[,index.defect], goback.nochrom[,index.cancer])[2,2]))
    
    write.table(estimates, file = 'W:/Old_genepi2/Jeremy/GOBACK/R Outputs/goback.coxph.top.hits.over.one.yr.v20180904.csv', 
                sep=',', append = TRUE, row.names = FALSE, col.names = FALSE)
    
  }
  
  else if (index.cancer %in% c('all','hepato','nephro')){
    
    print(paste('Computing BW-adjusted model for risk of',index.cancer, 'in kids with', index.defect))
    
    goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                              cancer = goback.nochrom[,index.cancer],
                              defect = goback.nochrom[,index.defect],
                              sex = factor(goback.nochrom$sex,
                                           levels = c(1,2),
                                           labels = c('Male','Female')),
                              m.age = goback.nochrom$m.age,
                              birth.wt = goback.nochrom$birth.wt,
                              state = goback.nochrom$state.num)
    
    cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + birth.wt + state, data = goback.surv)
    cox.coef <- summary(cox)$coefficients
    
    estimates <- data.frame(defect = index.defect, 
                            cancer = index.cancer, 
                            HR = exp(cox.coef[1,1]), 
                            ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                            ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                            p.value.coef = cox.coef[1,5],
                            num.comorbid = as.numeric(table(goback.nochrom[,index.defect], goback.nochrom[,index.cancer])[2,2]))
    
    write.table(estimates, file = 'W:/Old_genepi2/Jeremy/GOBACK/R Outputs/goback.coxph.top.hits.over.one.yr.v20180904.csv', 
                sep=',', append = TRUE, row.names = FALSE, col.names = FALSE)
    
  }
  
  else {
    
    print(paste('Computing model for risk of',index.cancer, 'in kids with', index.defect))
    
    goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                              cancer = goback.nochrom[,index.cancer],
                              defect = goback.nochrom[,index.defect],
                              sex = factor(goback.nochrom$sex,
                                           levels = c(1,2),
                                           labels = c('Male','Female')),
                              m.age = goback.nochrom$m.age,
                              state = goback.nochrom$state.num)
    
    cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
    cox.coef <- summary(cox)$coefficients
    
    estimates <- data.frame(defect = index.defect, 
                            cancer = index.cancer, 
                            HR = exp(cox.coef[1,1]), 
                            ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                            ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                            p.value.coef = cox.coef[1,5],
                            num.comorbid = as.numeric(table(goback.nochrom[,index.defect], goback.nochrom[,index.cancer])[2,2]))
    
    write.table(estimates, file = 'W:/Old_genepi2/Jeremy/GOBACK/R Outputs/goback.coxph.top.hits.over.one.yr.v20180904.csv', 
                sep=',', append = TRUE, row.names = FALSE, col.names = FALSE)
    
  }
  
}

rm(goback.nochrom, cox, cox.coef, goback.surv); gc()

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.chrom.v20180829.rdata')
goback.chrom <- filter(goback.chrom, state %in% c('TX','NC') & person.yrs > 1)

for (i in 1:nrow(chrom.associations)){
  
  index.defect <- chrom.associations$defect[i]
  index.cancer <- chrom.associations$cancer[i]
  
  if (index.cancer %in% c('all','hepato','nephro')){
    
    print(paste('Computing BW-adjusted model for risk of',index.cancer, 'in kids with', index.defect))
    
    goback.surv <- data.frame(time = goback.chrom$person.yrs,
                              cancer = goback.chrom[,index.cancer],
                              defect = goback.chrom[,index.defect],
                              sex = factor(goback.chrom$sex,
                                           levels = c(1,2),
                                           labels = c('Male','Female')),
                              m.age = goback.chrom$m.age,
                              birth.wt = goback.chrom$birth.wt,
                              state = goback.chrom$state.num)
    
    cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + birth.wt + state, data = goback.surv)
    cox.coef <- summary(cox)$coefficients
    
    estimates <- data.frame(defect = index.defect, 
                            cancer = index.cancer, 
                            HR = exp(cox.coef[1,1]), 
                            ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                            ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                            p.value.coef = cox.coef[1,5],
                            num.comorbid = as.numeric(table(goback.chrom[,index.defect], goback.chrom[,index.cancer])[2,2]))
    
    write.table(estimates, file = 'W:/Old_genepi2/Jeremy/GOBACK/R Outputs/goback.coxph.top.hits.over.one.yr.v20180904.csv', 
                sep=',', append = TRUE, row.names = FALSE, col.names = FALSE)
    
  }
  
  else {
    
    print(paste('Computing model for risk of',index.cancer, 'in kids with', index.defect))
    
    goback.surv <- data.frame(time = goback.chrom$person.yrs,
                              cancer = goback.chrom[,index.cancer],
                              defect = goback.chrom[,index.defect],
                              sex = factor(goback.chrom$sex,
                                           levels = c(1,2),
                                           labels = c('Male','Female')),
                              m.age = goback.chrom$m.age,
                              state = goback.chrom$state.num)
    
    cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
    cox.coef <- summary(cox)$coefficients
    
    estimates <- data.frame(defect = index.defect, 
                            cancer = index.cancer, 
                            HR = exp(cox.coef[1,1]), 
                            ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                            ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                            p.value.coef = cox.coef[1,5],
                            num.comorbid = as.numeric(table(goback.chrom[,index.defect], goback.chrom[,index.cancer])[2,2]))
    
    write.table(estimates, file = 'W:/Old_genepi2/Jeremy/GOBACK/R Outputs/goback.coxph.top.hits.over.one.yr.v20180904.csv', 
                sep=',', append = TRUE, row.names = FALSE, col.names = FALSE)
    
  }
  
}

rm(list = ls()); gc()



# Print side-by-side comparisons ------------------------------------------

require(dplyr); require(xlsx)

overall <- read.csv(file = 'W:/Old_genepi2/Jeremy/GOBACK/R outputs/goback.cox.ph.top.hits.v20180612.csv',
                    stringsAsFactors = FALSE)
overall <- overall[, 1:7]

sensitivity <- read.csv(file = 'W:/Old_genepi2/Jeremy/GOBACK/R outputs/goback.coxph.top.hits.over.one.yr.v20180904.csv',
                        stringsAsFactors = FALSE, header = FALSE)

names(sensitivity) <- c('defect','cancer',paste0('sensitivity.',names(overall)[3:7]))
names(overall)[3:7] <- paste0('overall.',names(overall)[3:7])

compare <- left_join(overall, sensitivity, by = c('defect','cancer'))

write.xlsx(compare, file = 'W:/Old_genepi2/Jeremy/GOBACK/R outputs/goback.cox.ph.top.hits.over.1.v20180904.xlsx',
           row.names = FALSE, sheetName = 'All Associations')

compare <- filter(compare, sensitivity.n.comorbid >= 5)

write.xlsx(compare, file = 'W:/Old_genepi2/Jeremy/GOBACK/R outputs/goback.cox.ph.top.hits.over.1.v20180904.xlsx',
           row.names = FALSE, sheetName = 'Five Comorbid Cases', append = TRUE)

rm(list = ls()); gc()