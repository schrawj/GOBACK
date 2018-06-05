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

defects <- c('any.birthdefect', 'conganomalies.cns')



# Generate Cox model: non-chromosomal -------------------------------------

#' There are no striking associations of chromosomal defects with CBT, 
#' so we will only perform this analysis in the non-chromosomal set.
load('./goback.no.chrom.v20180122.1.rdata')

#' Compute a variable 'any.cbt', i.e. any childhood brain tumor.
#' Set to 1 if child DX'd with astro, medullo, ependymoma or PNET.
goback.nochrom$any.cbt <- 0

for (i in cbts){
  goback.nochrom[,'any.cbt'] <- ifelse(goback.nochrom[,i] == 1, 1, goback.nochrom[,'any.cbt'])
}

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
  test.ph <- cox.zph(cox)
  test.ph <- test.ph$table['defect','p']
  
  estimates <- data.frame(defect = as.character(i), 
                          cancer = 'any.cbt', 
                          HR = exp(cox.coef[1,1]), 
                          ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                          ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                          p.value.coef = cox.coef[1,5],
                          p.value.zph = test.ph,
                          num.comorbid = tab, 
                          set = 'goback.nochrom')
  
  write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/btec.models.csv', sep=',', append = TRUE, 
              row.names = FALSE, col.names = FALSE)
  
}

rm(cox, cox.coef, estimates, goback.nochrom, goback.surv, tab, test.ph, i, cbts, defects)



# CBT heatmaps ------------------------------------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.06.01.
#' 
#' Analyses for BTEC presentation.
#' 
#' Heatmaps with X axis restricted to CBTs.
#' 
#' K-M curves and HRs for CBT risk by number of birth defects.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(dplyr); require(ggplot2); require(ggthemes)

setwd('Z:/Jeremy/GOBACK/Datasets/')
load('./goback.nochrom.v20180530.2.rdata')
load('./Expanded datasets/goback.cox.ph.results.v20180124.1.rdata')

#' Character vectors housing the defects and cancers of interest.
defects <- names(goback.nochrom[c(22,30,35,38,61,65,68,76,83,94,95)])
cancers <- c('cns.any','astro','medullo','ependymoma','pnet')

#' Generate a data frame housing these associations.
bd.cc.associations <- data.frame(defect = rep(defects, each = 5), cancer = rep(cancers, 11))
bd.cc.associations <- left_join(bd.cc.associations, goback.coxmodels, by = c('defect','cancer'))
colnames(bd.cc.associations) <- tolower(colnames(bd.cc.associations))

#' Generate an indicator variable for strength of association.
bd.cc.associations$signif.cat.int <- 9999
bd.cc.associations$signif.cat.int <- factor(
  
  ifelse(bd.cc.associations$signif.cat.int == 9999 & is.na(bd.cc.associations$hr), 7,
  ifelse(bd.cc.associations$ci.lower < 1 & bd.cc.associations$ci.upper > 1, 2,
  ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$ci.upper > 1 & bd.cc.associations$hr > 1 & bd.cc.associations$hr <= 5, 3,
  ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$ci.upper > 1 & bd.cc.associations$hr > 5 & bd.cc.associations$hr <= 10, 4,
  ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$ci.upper > 1 & bd.cc.associations$hr > 10 & bd.cc.associations$hr <= 20, 5,
  ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$ci.upper > 1 & bd.cc.associations$hr > 20, 6,
  ifelse(bd.cc.associations$ci.lower < 1 & bd.cc.associations$ci.upper < 1, 1, NA))))))),
  
  levels = c(1:7),
  labels = c('< 1','Null','1.01-5', '5.01-10','10.01-20','>20','Not Tested'))

#' Convert axis variables to factor form for prettier display.
bd.cc.associations$cancer <- factor(bd.cc.associations$cancer, 
                                    levels = c('pnet', 'ependymoma', 'medullo', 'astro', 'cns.any'),
                                    labels = c('PNET', 'Ependymoma', 'Medulloblastoma', 'Astrocytoma', 'Any CNS Tumor'))

bd.cc.associations$defect <- factor(bd.cc.associations$defect, 
                                    levels = c( "conganomalies.cns",               "conganomalies.eye",               "conganomalies.ear.face.neck",     "conganomalies.heart.circsys",    
                                                "conganomalies.respsys",           "oral.clefts",                     "conganomalies.digestivesystem",   "conganomalies.genitalandurinary",
                                                "conganomalies.musculoskelsys",    "conganomalies.integument",        "chromosomalanomalies"),
                                    labels = c('CNS Anomaly', 'Eye Anomaly', 'Ear, Face or Neck Anomaly', 'Heart or Circulatory System Anomaly',
                                               'Respiratory System Anomaly', 'Oral Clefts', 'Digestive System Anomaly', 'Genitourinary Anomaly',
                                               'Musculoskeletal Anomaly', 'Integument Anomaly', 'Chromosomal Anomaly'))

#' For landscape display.
plot <-  ggplot(data = bd.cc.associations, aes(x = defect, y = cancer)) +
  geom_tile(aes(fill = signif.cat.int), color = 'black') + 
  scale_fill_manual(values = c('grey75','indianred1','red','darkred','white')) +
  guides(fill = guide_legend(title='Hazard Ratio')) +
#labs(x = "Anomaly", y = "Cancer") +
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

#' Generate a categorical variable for number of birth defects.
goback.nochrom$majordefect.cat <- factor(ifelse(goback.nochrom$majordefect.total == 0, 0,
                                         ifelse(goback.nochrom$majordefect.total == 1, 1,
                                         ifelse(goback.nochrom$majordefect.total == 2, 2, 
                                         ifelse(goback.nochrom$majordefect.total == 3, 3,
                                         ifelse(goback.nochrom$majordefect.total >= 4, 4, NA))))),
                                         levels = c(0:4),
                                         labels = c('0', '1', '2', '3', '4 or more'))

#' Initialize an empty list to house plots.
l <- list()
j <- 1

tic()

for (i in c(140,111,126,114,132)){
  
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
             legend.labs = c('No birth defect', '1 defect', '2 defects', '3 defects', '4 or more defects'))
  
  l[j] <- plot
  
  rm(goback.surv, fit, plot)
  
  j <- j + 1
}

toc()

save(l, file = 'Z:/Jeremy/GOBACK/R Outputs/CBT KM curves/cbt.km.curves.rdata')

rm(list = ls()); gc()



# Cox PH model for any defect-any CBT -------------------------------------

require(survival)

setwd('Z:/Jeremy/GOBACK/Datasets/')
load('./goback.nochrom.v20180530.2.rdata')

table(goback.nochrom$any.birthdefect, goback.nochrom$cns.any, useNA = 'ifany')

goback.surv <- data.frame(studyid = goback.nochrom$studyid,
                          time = goback.nochrom$person.yrs, 
                          cancer = goback.nochrom$cns.any, 
                          defect = goback.nochrom$any.birthdefect,
                          sex = factor(goback.nochrom$sex, 
                                       levels = c(1,2),
                                       labels = c('Male','Female')),
                          m.age = goback.nochrom$m.age,
                          state = goback.nochrom$state)

cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)

cox.coef <- as.data.frame(summary(cox)$coefficients)
cox.coef$ci.lower <- exp(cox.coef$coef - (1.96*cox.coef$`se(coef)`))
cox.coef$ci.upper <- exp(cox.coef$coef + (1.96*cox.coef$`se(coef)`))
print(cox.coef)

rm(list = ls()); gc()



# Cox PH model for any CNS anomaly-any CBT --------------------------------

require(survival)

setwd('Z:/Jeremy/GOBACK/Datasets/')
load('./goback.nochrom.v20180530.2.rdata')

table(goback.nochrom$conganomalies.cns, goback.nochrom$cns.any, useNA = 'ifany')

goback.surv <- data.frame(studyid = goback.nochrom$studyid,
                          time = goback.nochrom$person.yrs, 
                          cancer = goback.nochrom$cns.any, 
                          defect = goback.nochrom$conganomalies.cns,
                          sex = factor(goback.nochrom$sex, 
                                       levels = c(1,2),
                                       labels = c('Male','Female')),
                          m.age = goback.nochrom$m.age,
                          state = goback.nochrom$state)

cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)

cox.coef <- as.data.frame(summary(cox)$coefficients)
cox.coef$ci.lower <- exp(cox.coef$coef - (1.96*cox.coef$`se(coef)`))
cox.coef$ci.upper <- exp(cox.coef$coef + (1.96*cox.coef$`se(coef)`))
print(cox.coef)

rm(list = ls()); gc()



# CBT Cox PH models by number of defects ----------------------------------

require(survival); require(tictoc)

setwd('Z:/Jeremy/GOBACK/Datasets/')
load('./goback.nochrom.v20180530.2.rdata')

#' Generate a categorical variable for number of birth defects.
goback.nochrom$majordefect.cat <- factor(ifelse(goback.nochrom$majordefect.total == 0, 0,
                                                ifelse(goback.nochrom$majordefect.total == 1, 1,
                                                       ifelse(goback.nochrom$majordefect.total == 2, 2, 
                                                              ifelse(goback.nochrom$majordefect.total == 3, 3,
                                                                     ifelse(goback.nochrom$majordefect.total >= 4, 4, NA))))),
                                         levels = c(0:4),
                                         labels = c('0', '1', '2', '3', '4 or more'))

#' Initialize an empty list to house plots.
cbt.models <- list()
j <- 1

for (i in c(140,111,126,114,132)){
  
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

save(cbt.models, file = 'Z:/Jeremy/GOBACK/Datasets/Expanded Datasets/cbt.risk.by.num.defects.rdata')