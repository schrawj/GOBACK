#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2019.02.26.
#' 
#' One of the reviewers at JAMA Oncology requested that we stratify 
#' by race-ethnicity where possible. This is probably only feasible for 
#' broader analyses. We will do this for the major groups in table 2.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(survival); require(dplyr)

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20180829.rdata')

#' Maternal race-ethnicity will need to be collapsed for this analysis.
#' NHW, NHB, Hispanic. Reclassify NC AIAIN children first.
goback$m.race <- as.character(goback$m.race)
goback$m.race <- ifelse(goback$m.race == 'AIAN' & goback$state == 'NC', 'Hispanic', goback$m.race)

#' A variable for any major non-chromosomal defect would be helpful.
#' 75,203 NAs are 49777 kids with isolated minor defects and 25427 kids with genetic condiitons.
goback$major.nonchromosomal.defect <- ifelse((goback$any.birthdefect == 1 & goback$majordefect.total >= 1 & is.na(goback$any.genetic.anomaly)), 1, 
                                             ifelse(goback$any.genetic.anomaly == 1, NA, 0))

#' Optionally, filter out unnecessary rows and some unnecessary columns. May help performance.
goback <- filter(goback, m.race %in% c('NHW','NHB','Hispanic'))
goback <- goback[, c(1:95, 107,156)]; gc()

defects <- subset(names(goback), grepl('conganomalies.', names(goback)))
defects <- subset(defects, !grepl('.other', defects))
defects <- c('major.nonchromosomal.defect', defects, 'oral.clefts','any.genetic.anomaly')
defect.cols <- which(names(goback) %in% defects)

groups <- unique(goback$m.race)

#' Compute Cox regression models for select categories of defects.
for (i in defect.cols){
  
  for (j in groups){
    
    tmp <- subset(goback, goback$m.race == j)
  
    tab <- table(tmp[,i], tmp$cancer)[2,2]
    
    if (tab >= 5){
      
      tmp.surv <- data.frame(time = tmp$person.yrs,
                                cancer = tmp$cancer,
                                defect = tmp[,i],
                                sex = factor(tmp$sex,
                                             levels = c(1,2),
                                             labels = c('Male','Female')),
                                m.age = tmp$m.age,
                                state = tmp$state.num)      
      
      cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = tmp.surv)
      cox.coef <- summary(cox)$coefficients
      test.ph <- cox.zph(cox)
      test.ph <- test.ph$table['defect','p']
      
      rm(cox, tmp.surv); gc()
      
      estimates <- data.frame(defect = names(tmp[i]), 
                              cancer = 'any.cancer', 
                              race.eth = j,
                              HR = exp(cox.coef[1,1]), 
                              ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                              ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                              p.value.coef = cox.coef[1,5],
                              p.value.zph = test.ph,
                              num.comorbid = tab)
      
      write.table(estimates, file = 'W:/Old_genepi2/Jeremy/GOBACK/R outputs/risk.of.any.cancer.by.birth.defects.and.race.ethnicity.v20190226.csv', 
                  sep = ',', append = TRUE, row.names = FALSE, col.names = FALSE)
    }
    
    else{
      
      next
      
    }
  
  }
  
}

#' Round the HR and CI values to 2 decimal places using Excel.
models <- read.csv(file = 'W:/Old_genepi2/Jeremy/GOBACK/R outputs/risk.of.any.cancer.by.birth.defects.and.race.ethnicity.v20190226.csv',
                   stringsAsFactors = F, header = F)
names(models) <- c('defect','cancer','race.eth','hr','ci.lower','ci.upper','pval.coef','pval.zph','n.comorbid')
models$estimate <- paste0(models$hr,' (',models$ci.lower,'-',models$ci.upper,')')

write.csv(models, file = 'W:/Old_genepi2/Jeremy/GOBACK/R outputs/risk.of.any.cancer.by.birth.defects.and.race.ethnicity.v20190227.csv', row.names = F)

rm(list = ls()); gc()
