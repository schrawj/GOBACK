#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.06.07.
#' 
#' Analyses for BTEC abstract.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------


# Any cancer and any CBT risk by any defect and any CNS defect ------------

require(dplyr); require(survival)

#' There are no striking associations of chromosomal defects with CBT, 
#' so we will only perform this analysis in the non-chromosomal set.
#' We'll also want to exclude children with neurofibromatosis or TSC 
#' from these analyses.
setwd('Z:/Jeremy/GOBACK/Datasets/')
load('./Expanded datasets/list.of.syndromic.kids.in.tx.mi.nc.rdata')
load('./goback.nochrom.v20180530.2.rdata')

exclusions <- c(syndrome.ids$tuberous.sclerosis, syndrome.ids$neurofibromatosis)

goback.nochrom <- filter(goback.nochrom, !(studyid %in% exclusions))

cbts <- c('astro','medullo','ependymoma','pnet')

#' Compute a variable 'any.cbt', i.e. any childhood brain tumor.
#' Set to 1 if child DX'd with astro, medullo, ependymoma or PNET.
goback.nochrom$any.cbt <- 0

for (i in cbts){
  goback.nochrom[,'any.cbt'] <- ifelse(goback.nochrom[,i] == 1, 1, goback.nochrom[,'any.cbt'])
}

defects <- c('any.birthdefect', 'conganomalies.cns')

estimates <- as.data.frame(matrix(nrow = 1, ncol = 8))
names(estimates) <- c('defect', 'cancer', 'hr', 'ci.lower', 'ci.upper', 'p.val.coef', 'num.comorbid', 'set')

#' Cox PH model for any defect - any cancer.
tab <- as.numeric(table(goback.nochrom$any.birthdefect, goback.nochrom$cancer)[2,2])

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

new.estimates <- data.frame(defect = 'any.birthdefect', 
                            cancer = 'any.cancer', 
                            hr = exp(cox.coef[1,1]), 
                            ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                            ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                            p.val.coef = cox.coef[1,5],
                            num.comorbid = tab, 
                            set = 'goback.nochrom.nosyn')

estimates <- rbind(estimates, new.estimates)

#' Cox PH models for any CBT in children with any defect and any CNS defect.
for (i in defects){
  
  tab <- as.numeric(table(goback.nochrom[,i], goback.nochrom$any.cbt)[2,2])
  
  goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                            cancer = goback.nochrom$any.cbt,
                            defect = goback.nochrom[,i],
                            sex = factor(goback.nochrom$sex,
                                         levels = c(1,2),
                                         labels = c('Male','Female')),
                            m.age = goback.nochrom$m.age,
                            state = goback.nochrom$state.num)
  
  cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
  cox.coef <- summary(cox)$coefficients

  new.estimates <- data.frame(defect = as.character(i), 
                          cancer = 'any.cbt', 
                          hr = exp(cox.coef[1,1]), 
                          ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                          ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                          p.val.coef = cox.coef[1,5],
                          num.comorbid = tab, 
                          set = 'goback.nochrom')
  
  estimates <- rbind(estimates, new.estimates)
  
}

estimates <- estimates[2:4, ]

write.table(estimates, file = 'Z:/Jeremy/GOBACK/R outputs/btec.models.csv', sep = ',', row.names = FALSE)

rm(defects, tab, goback.surv, cox, cox.coef, new.estimates, i, exclusions, cbts); gc()



# Organ system defects and CBTs -------------------------------------------

require(dplyr); require(survival)

setwd('Z:/Jeremy/GOBACK/Datasets/')
load('./goback.nochrom.v20180530.2.rdata')
load('./Expanded datasets/list.of.syndromic.kids.in.tx.mi.nc.rdata')

#' Exclude children with neurofibromatosis or TSC.
exclusions <- c(syndrome.ids$tuberous.sclerosis, syndrome.ids$neurofibromatosis)

goback.nochrom <- filter(goback.nochrom, !(studyid %in% exclusions))

#' Compute a variable 'any.cbt', i.e. any childhood brain tumor.
#' Set to 1 if child DX'd with astro, medullo, ependymoma or PNET.
goback.nochrom$any.cbt <- 0

cbts <- c('astro','medullo','ependymoma','pnet')

for (i in cbts){
  goback.nochrom[,'any.cbt'] <- ifelse(goback.nochrom[,i] == 1, 1, goback.nochrom[,'any.cbt'])
}

#' Initialize empty data frame to hold results for Cox models below.
estimates <- as.data.frame(matrix(nrow = 1, ncol = 8))
names(estimates) <- c('defect', 'cancer', 'hr', 'ci.lower', 'ci.upper', 'p.val.coef', 'num.comorbid', 'set')

#' Values in i represent organ system defect variables.
#' Values in j represent CBT variables.
for (i in c(22,30,35,38,61,65,68,76,83,94)){                       
  
  for (j in c(153,111,126,114,132)){ 
    
    comorbid.cases <- table(goback.nochrom[,i], goback.nochrom[,j])[2,2]
    
    if (comorbid.cases >= 5){
      
      print(paste('computing model for', names(goback.nochrom[i]),'and',names(goback.nochrom[j])))
      
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

      new.estimates <- data.frame(defect = names(goback.nochrom[i]), 
                              cancer = names(goback.nochrom[j]), 
                              hr = exp(cox.coef[1,1]), 
                              ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                              ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                              p.val.coef = cox.coef[1,5],
                              num.comorbid = comorbid.cases,
                              set = 'goback.nochrom.nosyn')
      
      estimates <- rbind(estimates, new.estimates)
      
    }
    
    else{
      
      next
      
    }
  }
}

rm(i, j, exclusions, cbts, comorbid.cases, goback.surv, cox, cox.coef, new.estimates, goback.nochrom); gc()
              
#' Generate results for any chromosomal anomaly.
load('./goback.chrom.v20180530.2.rdata')
load('./Expanded datasets/list.of.syndromic.kids.in.tx.mi.nc.rdata')

#' Exclude children with neurofibromatosis or TSC.
exclusions <- c(syndrome.ids$tuberous.sclerosis, syndrome.ids$neurofibromatosis)

goback.chrom <- filter(goback.chrom, !(studyid %in% exclusions))

#' Compute a variable 'any.cbt', i.e. any childhood brain tumor.
#' Set to 1 if child DX'd with astro, medullo, ependymoma or PNET.
goback.chrom$any.cbt <- 0

cbts <- c('astro','medullo','ependymoma','pnet')

for (i in cbts){
  goback.chrom[,'any.cbt'] <- ifelse(goback.chrom[,i] == 1, 1, goback.chrom[,'any.cbt'])
}

#' It turns out there are not 5 co-occuring cases of chromosomal anomalies and CBTs.
for (j in c(153,111,126,114,132)){ 
  
  comorbid.cases <- table(goback.chrom$chromosomalanomalies, goback.chrom[,j])[2,2]
  
  if (comorbid.cases >= 5){
    
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

    new.estimates <- data.frame(defect = 'chromosomal.anomalies', 
                                cancer = names(goback.chrom[j]), 
                                hr = exp(cox.coef[1,1]), 
                                ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                                ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                                p.val.coef = cox.coef[1,5],
                                num.comorbid = comorbid.cases,
                                set = 'goback.chrom.nosyn')
    
    estimates <- rbind(estimates, new.estimates)
    
  }
  
  else{
    
    next
    
  }
}

estimates <- estimates[2:32, ]

write.table(estimates, file = 'Z:/Jeremy/GOBACK/R outputs/btec.models.csv', sep = ',', append = TRUE, row.names = FALSE, col.names = FALSE)

rm(list = ls()); gc()



# CBT heatmaps ------------------------------------------------------------

require(dplyr); require(ggplot2)

cbt.results <- read.csv('Z:/Jeremy/GOBACK/R outputs/btec.models.csv', header = TRUE)
cbt.results <- cbt.results[4:37, ]

#' Character vectors housing the defects and cancers of interest.
#defects <- names(goback.nochrom[c(22,30,35,38,61,65,68,76,83,94)])
defects <- c( "conganomalies.cns",               "conganomalies.eye",               "conganomalies.ear.face.neck",     "conganomalies.heart.circsys",    
              "conganomalies.respsys",           "oral.clefts",                     "conganomalies.digestivesystem",   "conganomalies.genitalandurinary",
              "conganomalies.musculoskelsys",    "conganomalies.integument")
cancers <- c('any.cbt','astro','medullo','ependymoma','pnet')

#' Generate a data frame housing these associations.
bd.cc.associations <- data.frame(defect = rep(defects, each = 5), cancer = rep(cancers, 10))
bd.cc.associations <- left_join(bd.cc.associations, cbt.results, by = c('defect','cancer'))
colnames(bd.cc.associations) <- tolower(colnames(bd.cc.associations))

#' Generate an indicator variable for strength of association.
bd.cc.associations$signif.cat.int <- 9999
bd.cc.associations$signif.cat.int <- factor(
  
  ifelse(bd.cc.associations$signif.cat.int == 9999 & is.na(bd.cc.associations$hr), 7,
  ifelse(bd.cc.associations$ci.lower < 1 & bd.cc.associations$ci.upper > 1, 2,
  ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$hr > 1 & bd.cc.associations$hr <= 2.5, 3,
  ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$hr > 2.5 & bd.cc.associations$hr <= 5, 4,
  ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$hr > 5 & bd.cc.associations$hr <= 10, 5,
  ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$hr > 10, 6,
  ifelse(bd.cc.associations$ci.lower < 1 & bd.cc.associations$ci.upper < 1, 1, NA))))))),
  
  levels = c(1:7),
  labels = c('< 1','Null','1.01-2.50', '2.51-5.00','5.01-10.00','>10.00','Not Tested'))

#' Convert axis variables to factors for better display.
bd.cc.associations$cancer <- factor(bd.cc.associations$cancer, 
                                    levels = c('pnet', 'ependymoma', 'medullo', 'astro', 'any.cbt'),
                                    labels = c('PNET', 'Ependymoma', 'Medulloblastoma', 'Astrocytoma', 'Any CBT'))

bd.cc.associations$defect <- factor(bd.cc.associations$defect, 
                                    levels = c( "conganomalies.cns",               "conganomalies.eye",               "conganomalies.ear.face.neck",     "conganomalies.heart.circsys",    
                                                "conganomalies.respsys",           "oral.clefts",                     "conganomalies.digestivesystem",   "conganomalies.genitalandurinary",
                                                "conganomalies.musculoskelsys",    "conganomalies.integument"),
                                    labels = c('CNS Anomaly', 'Eye Anomaly', 'Ear, Face or Neck Anomaly', 'Heart or Circulatory System Anomaly',
                                               'Respiratory System Anomaly', 'Oral Clefts', 'Digestive System Anomaly', 'Genitourinary Anomaly',
                                               'Musculoskeletal Anomaly', 'Integument Anomaly'))

#' For landscape display.
plot <-  ggplot(data = bd.cc.associations, aes(x = defect, y = cancer)) +
  geom_tile(aes(fill = signif.cat.int), color = 'black') + 
  scale_fill_manual(values = c('grey75','indianred1','red','darkred','white')) +
  guides(fill = guide_legend(title='Hazard Ratio')) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_blank()) 
print(plot)

rm(list = ls()); gc()



# CBT K-M curves by number of defects -------------------------------------

require(survival); require(survminer); require(tictoc)

setwd('Z:/Jeremy/GOBACK/Datasets/')
load('./goback.nochrom.v20180530.2.rdata')
load('./Expanded datasets/list.of.syndromic.kids.in.tx.mi.nc.rdata')

#' Exclude children with neurofibromatosis or TSC.
exclusions <- c(syndrome.ids$tuberous.sclerosis, syndrome.ids$neurofibromatosis)

goback.nochrom <- filter(goback.nochrom, !(studyid %in% exclusions))

#' Generate a categorical variable for number of birth defects.
goback.nochrom$majordefect.cat <- factor(ifelse(goback.nochrom$majordefect.total == 0, 0,
                                         ifelse(goback.nochrom$majordefect.total == 1, 1,
                                         ifelse(goback.nochrom$majordefect.total == 2, 2, 
                                         ifelse(goback.nochrom$majordefect.total == 3, 3,
                                         ifelse(goback.nochrom$majordefect.total >= 4, 4, NA))))),
                                         levels = c(0:4),
                                         labels = c('0', '1', '2', '3', '4 or more'))

#' Compute a variable 'any.cbt', i.e. any childhood brain tumor.
#' Set to 1 if child DX'd with astro, medullo, ependymoma or PNET.
goback.nochrom$any.cbt <- 0

cbts <- c('astro','medullo','ependymoma','pnet')

for (i in cbts){
  goback.nochrom[,'any.cbt'] <- ifelse(goback.nochrom[,i] == 1, 1, goback.nochrom[,'any.cbt'])
}

#' Initialize an empty list to house plots.
l <- list()
j <- 1

tic()

for (i in c(154,111,126,114,132)){
  
  print(paste('starting iteration', j))
  
  #' Generate survival object.
  goback.surv <- data.frame(studyid = goback.nochrom$studyid,
                            time = goback.nochrom$person.yrs, 
                            cancer = goback.nochrom[,i], 
                            defect = goback.nochrom$majordefect.cat,
                            sex = factor(goback.nochrom$sex, 
                                         levels = c(1,2),
                                         labels = c('Male','Female')),
                            m.age = goback.nochrom$m.age,
                            state = goback.nochrom$state)
  
  #' Generate K-M curves from survival object.
  fit <- survfit(Surv(time, cancer) ~ defect, data = goback.surv)
  
  plot <- ggsurvplot(fit,
             ylim = c(0.995, 1),
             xlim = c(0,18),
             legend.labs = c('No anomaly', '1 anomaly', '2 anomalies', '3 anomalies', '4 or more anomalies'))
  
  l[j] <- plot
  
  rm(goback.surv, fit, plot)
  
  j <- j + 1
}

toc()

save(l, file = 'Z:/Jeremy/GOBACK/R Outputs/CBT KM curves/cbt.km.curves.wo.syndromic.cases.rdata')

rm(list = ls()); gc()



# CBT Cox PH models by number of defects ----------------------------------

require(survival); require(tictoc)

setwd('Z:/Jeremy/GOBACK/Datasets/')
load('./goback.nochrom.v20180530.2.rdata')
load('./Expanded datasets/list.of.syndromic.kids.in.tx.mi.nc.rdata')

#' Exclude children with neurofibromatosis or TSC.
exclusions <- c(syndrome.ids$tuberous.sclerosis, syndrome.ids$neurofibromatosis)

goback.nochrom <- filter(goback.nochrom, !(studyid %in% exclusions))

#' Generate a categorical variable for number of birth defects.
goback.nochrom$majordefect.cat <- factor(ifelse(goback.nochrom$majordefect.total == 0, 0,
                                                ifelse(goback.nochrom$majordefect.total == 1, 1,
                                                       ifelse(goback.nochrom$majordefect.total == 2, 2, 
                                                              ifelse(goback.nochrom$majordefect.total == 3, 3,
                                                                     ifelse(goback.nochrom$majordefect.total >= 4, 4, NA))))),
                                         levels = c(0:4),
                                         labels = c('0', '1', '2', '3', '4 or more'))

#' Compute a variable 'any.cbt', i.e. any childhood brain tumor.
#' Set to 1 if child DX'd with astro, medullo, ependymoma or PNET.
goback.nochrom$any.cbt <- 0

cbts <- c('astro','medullo','ependymoma','pnet')

for (i in cbts){
  goback.nochrom[,'any.cbt'] <- ifelse(goback.nochrom[,i] == 1, 1, goback.nochrom[,'any.cbt'])
}

#' Initialize an empty list to house plots.
cbt.models <- list()
j <- 1

for (i in c(154,111,126,114,132)){
  
  print(paste('starting iteration', j))
  
  goback.surv <- data.frame(studyid = goback.nochrom$studyid,
                            time = goback.nochrom$person.yrs, 
                            cancer = goback.nochrom[,i], 
                            defect = goback.nochrom$majordefect.cat,
                            sex = factor(goback.nochrom$sex, 
                                         levels = c(1,2),
                                         labels = c('Male','Female')),
                            m.age = goback.nochrom$m.age,
                            state = goback.nochrom$state)
  
  cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
  
  cox.coef <- as.data.frame(summary(cox)$coefficients)
  cox.coef$ci.lower <- exp(cox.coef$coef - (1.96*cox.coef$`se(coef)`))
  cox.coef$ci.upper <- exp(cox.coef$coef + (1.96*cox.coef$`se(coef)`))
  
  cbt.models[[j]] <- cox.coef
  names(cbt.models)[j] <- names(goback.nochrom[i])
  
  j <- j + 1
  
}

save(cbt.models, file = 'Z:/Jeremy/GOBACK/Datasets/Expanded Datasets/cbt.risk.by.num.defects.no.syndromic.cases.rdata')
