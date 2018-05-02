
# Prep environment --------------------------------------------------------

setwd('Z:/Jeremy/GOBACK/')

require(dplyr); require(survival)


# Load in data ------------------------------------------------------------

load("./Datasets/goback.cox.ph.top.hits.v20180223.1.rdata")

#' Create vectors of defects and cancers to loop through.
nonchrom.defects <- c(top.hits$defect[-c(15:17,42)])
chrom.defects <- c(top.hits$defect[c(15:17,42)])
nonchrom.cancers <- c(top.hits$cancer[-c(15:17,42)])
chrom.cancers <- c(top.hits$cancer[c(15:17,42)])

load('./Datasets/goback.chrom.v20180122.1.rdata')
load('./Datasets/cancer.codes.v20180226.1.rdata')



# Generate Cox models using only cases > 1 yr (chromosomal) ---------------

cancer.ids <- list()

bw.adj <- c('hepato','all','nephro') 

for (i in 1:4){
  
  index.defect <- chrom.defects[i]
  index.cancer <- chrom.cancers[i]

  surv.data <- data.frame(studyid = goback.chrom$studyid, 
                    state = goback.chrom$state, 
                    sex = factor(goback.chrom$sex, 
                                 levels = c(1,2),
                                 labels = c('Male','Female')),
                    m.age = goback.chrom$m.age,
                    birth.wt = goback.chrom$birth.wt,
                    time = goback.chrom$person.yrs,
                    defect = goback.chrom[,index.defect],
                    cancer = goback.chrom[,index.cancer])
  surv.data$studyid <- as.character(surv.data$studyid)
  surv.data$state <- as.numeric(surv.data$state)
  
  tmp <- filter(surv.data, surv.data$cancer == 1 & surv.data$defect == 1)
  tmp <- c(tmp$studyid)
  
  cancer.ids[[i]] <- data.frame(select(filter(cancer.codes, studyid %in% tmp), studyid, morph31, site_code1, behavior1))
  names(cancer.ids)[i] <- paste0(index.defect,'-',index.cancer)

  surv.data <- filter(surv.data, cancer == 0 | (cancer == 1 & time > 1))
  
  tab <- as.numeric(table(surv.data[,'defect'], surv.data[,'cancer'])[2,2])
  
  if (tab > 0){
    
    if (index.cancer %in% bw.adj){
      
      cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + birth.wt + state, data = surv.data)
      cox.coef <- summary(cox)$coefficients
      test.ph <- cox.zph(cox)
      test.ph <- test.ph$table['defect','p']
      
    }
    
    else{

      cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = surv.data)
      cox.coef <- summary(cox)$coefficients
      test.ph <- cox.zph(cox)
      test.ph <- test.ph$table['defect','p']
      
    }
    
    estimates <- data.frame(defect = index.defect, 
                            cancer = index.cancer, 
                            hr = exp(cox.coef[1,1]), 
                            ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                            ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                            p.value.coef = cox.coef[1,5],
                            p.value.zph = test.ph,
                            num.comorbid = tab,
                            set = 'GOBACK.chrom')
    
    write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/top.hits.excluding.infants.csv', sep=',', append = TRUE, 
                row.names = FALSE, col.names = FALSE)
    
  }
  
  else{
    
    next
  }
}
  
rm(cox, cox.coef, test.ph, tab, estimates, tmp, surv.data, index.defect, index.cancer, i, goback.chrom, chrom.cancer, chrom.defects); gc()


  

# Generate Cox models using only cases >1 yr (non-chromosomal) ------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' This code will fail if you run it without running the section above, 
#' because it refers to the cancer.ids and bw.adj objects.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

load('./Datasets/goback.no.chrom.v20180122.1.rdata')

start.time <- Sys.time()

for (i in 1:39){
  
  index.defect <- nonchrom.defects[i]
  index.cancer <- nonchrom.cancers[i]
  
  surv.data <- data.frame(studyid = goback.nochrom$studyid, 
                          state = goback.nochrom$state, 
                          sex = factor(goback.nochrom$sex, 
                                       levels = c(1,2),
                                       labels = c('Male','Female')),
                          m.age = goback.nochrom$m.age,
                          birth.wt = goback.nochrom$birth.wt,
                          time = goback.nochrom$person.yrs,
                          defect = goback.nochrom[,index.defect],
                          cancer = goback.nochrom[,index.cancer])
  surv.data$studyid <- as.character(surv.data$studyid)
  surv.data$state <- as.numeric(surv.data$state)
  
  tmp <- filter(surv.data, surv.data$cancer == 1 & surv.data$defect == 1)
  tmp <- c(tmp$studyid)
  
  cancer.ids[[i]] <- data.frame(select(filter(cancer.codes, studyid %in% tmp), studyid, morph31, site_code1, behavior1))
  names(cancer.ids)[i] <- paste0(index.defect,'-',index.cancer)
  
  surv.data <- filter(surv.data, cancer == 0 | (cancer == 1 & time > 1))
  
  tab <- as.numeric(table(surv.data[,'defect'], surv.data[,'cancer'])[2,2])
  
  if (tab > 0){
    
    if (index.cancer %in% bw.adj){
      
      cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + birth.wt + state, data = surv.data)
      cox.coef <- summary(cox)$coefficients
      test.ph <- cox.zph(cox)
      test.ph <- test.ph$table['defect','p']
      
    }
    
    else{
      
      cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = surv.data)
      cox.coef <- summary(cox)$coefficients
      test.ph <- cox.zph(cox)
      test.ph <- test.ph$table['defect','p']
      
    }
    
    estimates <- data.frame(defect = index.defect, 
                            cancer = index.cancer, 
                            hr = exp(cox.coef[1,1]), 
                            ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                            ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                            p.value.coef = cox.coef[1,5],
                            p.value.zph = test.ph,
                            num.comorbid = tab,
                            set = 'GOBACK.NON.CHROM')
    
    write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/top.hits.excluding.infants.csv', sep=',', append = TRUE, 
                row.names = FALSE, col.names = FALSE)
    
  }
  
  else{
    
    next
  }
}

end.time <- Sys.time()

print(end.time - start.time)

rm(goback.nochrom, cox, cox.coef, test.ph, surv.data, estimates, tab, i, nonchrom.cancers, nonchrom.defects, index.cancer, index.defect, 
   cancer.codes, bw.adj, start.time, end.time, tmp); gc()



# Generate Cox models for any CBT in kids >1 yr in TX ---------------------

require(dplyr); require(survival)

setwd('Z:/Jeremy/GOBACK/')
load('./Datasets/goback.nochrom.v20180419.rdata')

goback.nochrom <- filter(goback.nochrom, state == 'TX')
goback.nochrom <- filter(goback.nochrom, person.yrs >= 1)

tab <- table(goback.nochrom$hydrocephalus.wo.spinabifida, goback.nochrom$cns.any)[2,2]

goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                          cancer = goback.nochrom$cns.any,
                          defect = goback.nochrom$hydrocephalus.wo.spinabifida,
                          sex = factor(goback.nochrom$sex,
                                       levels = c(1,2),
                                       labels = c('Male','Female')),
                          m.age = goback.nochrom$m.age)      

cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex, data = goback.surv)
cox.coef <- summary(cox)$coefficients

rm(cox, goback.surv); gc()

estimates <- data.frame(defect = 'hydrocephalus.wo.spinabifida', 
                        cancer = 'any.cancer', 
                        HR = exp(cox.coef[1,1]), 
                        ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                        ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                        p.value.coef = cox.coef[1,5],
                        num.comorbid = tab)
print(estimates)

#' Have just entered directly into 'GOBACK Jama Peds Tables v2018502.docx'.
#' Too lazy to write one row to a csv file.
rm(list = ls()); gc()



# Side-by-side comparison of effect estimates overall and w/o infa --------

load('./Datasets/goback.cox.ph.top.hits.v20180223.2.rdata')

top.hits.sub <- read.csv(file = 'Z:/Jeremy/GOBACK/R Outputs/top.hits.excluding.infants.csv', 
                         stringsAsFactors = FALSE, sep = ',', 
                         col.names = c('defect','cancer','over1.hr','over1.ci.lower','over1.ci.upper','over1.p.val.coef','over1.p.val.zph','over1.num.comorbid','set'))

top.hits$estimate <- paste0(top.hits$HR, ' (',top.hits$CI.lower,'-',top.hits$CI.upper,')')
top.hits.sub$over1.estimate <- paste0(top.hits.sub$over1.hr, ' (',top.hits.sub$over1.ci.lower,'-',top.hits.sub$over1.ci.upper,')')

compare <- left_join(select(top.hits, defect, cancer, HR, CI.lower, CI.upper, n.comorbid, is.new),
                     select(top.hits.sub, defect, cancer, over1.hr, over1.ci.lower, over1.ci.upper, over1.num.comorbid, set),
                     by = c('defect','cancer'))
compare <- compare[,c(1:6,8:12,7)]
write.csv(compare, file = './R Outputs/top.hits.with.and.wo.infants.v20180226.1.csv', row.names = FALSE)



# Reformat HR and CI for easy import to MS Word ---------------------------

#' I have manually rounded the HR and CI values down to 2 decimal places in Excel.
tmp <- read.csv('./R Outputs/top.hits.excluding.infants.csv', header = TRUE)

#' A column with a nicely formatted HR and 95% CI.
tmp$estimate <- paste0(tmp$hr,' (',tmp$ci.lower,'-',tmp$ci.upper,')')
tmp <- select(tmp, defect, cancer, estimate, n.comorbid)

write.csv(tmp, file = './R Outputs/top.hits.excluding.infants.v20180227.1.csv', row.names = FALSE)



# Export cancer codes for comorbid cases ----------------------------------

cancer.ids <- lapply(cancer.ids, function(x) {arrange(x, by = morph31)})
#' TODO: bind in narrative info about diagnostic codes.

load("./Datasets/goback.cox.ph.top.hits.v20180223.1.rdata")

#' Create vectors of defects and cancers to loop through.
nonchrom.defects <- c(top.hits$defect[-c(15:17,42)])
nonchrom.cancers <- c(top.hits$cancer[-c(15:17,42)])

for( i in 1:length(cancer.ids)){
  xlsx::write.xlsx(cancer.ids[i], file = paste0('Z:/Jeremy/GOBACK/R Outputs/Cancer codes in comorbid cases/',nonchrom.cancers[i],'.xlsx'), 
             append = TRUE, row.names = FALSE, sheetName = nonchrom.defects[i])
}



# A version that will only pull the IDs of comorbid cases -----------------

load('./Datasets/goback.no.chrom.v20180122.1.rdata')

start.time <- Sys.time()

cancer.ids <- list()

for (i in 1:39){
  
  index.defect <- nonchrom.defects[i]
  index.cancer <- nonchrom.cancers[i]
  
  surv.data <- data.frame(studyid = goback.nochrom$studyid, 
                          defect = goback.nochrom[,index.defect],
                          cancer = goback.nochrom[,index.cancer])
  surv.data$studyid <- as.character(surv.data$studyid)

  tmp <- filter(surv.data, surv.data$cancer == 1 & surv.data$defect == 1)
  tmp <- c(tmp$studyid)
  
  cancer.ids[[i]] <- data.frame(select(filter(cancer.codes, studyid %in% tmp), studyid, morph31, site_code1, behavior1))
  names(cancer.ids)[i] <- paste0(index.defect,'-',index.cancer)
  
  print(i)
  
 }

end.time <- Sys.time()

print(end.time - start.time)

cancer.ids <- lapply(cancer.ids, function(x) {arrange(x, by = morph31)})

for( i in 1:length(cancer.ids)){
  xlsx::write.xlsx(cancer.ids[i], file = paste0('Z:/Jeremy/GOBACK/R Outputs/Cancer codes in comorbid cases/',nonchrom.cancers[i],'.xlsx'), 
                   append = TRUE, row.names = FALSE, sheetName = nonchrom.defects[i])
}

rm(goback.nochrom, nonchrom.cancers, nonchrom.defects, surv.data, tmp, i, start.time, end.time, index.defect, index.cancer); gc()

cancer.ids <- list()

load('./Datasets/goback.chrom.v20180122.1.rdata')

for (i in 1:4){
  
  index.defect <- chrom.defects[i]
  index.cancer <- chrom.cancers[i]
  
  surv.data <- data.frame(studyid = goback.chrom$studyid, 
                          defect = goback.chrom[,index.defect],
                          cancer = goback.chrom[,index.cancer])
  surv.data$studyid <- as.character(surv.data$studyid)

  tmp <- filter(surv.data, surv.data$cancer == 1 & surv.data$defect == 1)
  tmp <- c(tmp$studyid)
  
  cancer.ids[[i]] <- data.frame(select(filter(cancer.codes, studyid %in% tmp), studyid, morph31, site_code1, behavior1))
  names(cancer.ids)[i] <- paste0(index.defect,'-',index.cancer)
  
  print(i)
  
}

cancer.ids <- lapply(cancer.ids, function(x) {arrange(x, by = morph31)})

for( i in 1:length(cancer.ids)){
  xlsx::write.xlsx(cancer.ids[i], file = paste0('Z:/Jeremy/GOBACK/R Outputs/Cancer codes in comorbid cases/',chrom.cancers[i],'.xlsx'), 
                   append = TRUE, row.names = FALSE, sheetName = chrom.defects[i])
  }
