#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.08.23.
#' 
#' Cumulative incidence of any cancer, any heme cancer, any cns cancer, and
#' any non-cns solid tissue cancer by number of major defects.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(survival); require(survminer)

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.nochrom.v20180829.rdata')

goback.nochrom$majordefect.cat <- factor(ifelse(goback.nochrom$majordefect.total == 0, 0,
                                                ifelse(goback.nochrom$majordefect.total == 1, 1,
                                                       ifelse(goback.nochrom$majordefect.total == 2, 2, 
                                                              ifelse(goback.nochrom$majordefect.total == 3, 3,
                                                                     ifelse(goback.nochrom$majordefect.total >= 4, 4, NA))))),
                                         levels = c(0:4),
                                         labels = c('0', '1', '2', '3', '4 or more'))

#' Need new variables for any heme cancer (any leukemia or any lymphoma) and forany non-CNS solid tumor.
goback.nochrom$any.heme.cancer <- ifelse(goback.nochrom$leu.any == 1 | goback.nochrom$lym.any == 1, 1, 0)

solid.tumors <- unique(goback.nochrom$cancer1)
solid.tumors <- subset(solid.tumors, !(solid.tumors %in% c(NA, 'all','leu.other','aml','hl','nhl','cns.other','medullo','pnet','lym.other',
                                                           'gct.intra','astro','ependymoma')))

goback.nochrom$any.non.cns.solid.tumor <- ifelse(goback.nochrom$cancer1 %in% solid.tumors, 1, 0)

outcomes <- c('cancer','any.heme.cancer','cns.any','any.non.cns.solid.tumor')
outcomes.fancy <- c('Cancer','Hematologic Cancer', "Central Nervous Sytem (CNS) Tumor", 'Non-CNS Solid Tumor')

the.plots.thicken <- list()

for (i in 1:length(outcomes)){
  
  goback.surv <- data.frame(studyid = goback.nochrom$studyid,
                            time = goback.nochrom$person.yrs, 
                            cancer = goback.nochrom[,outcomes[i]], 
                            defect = goback.nochrom$majordefect.cat,
                            sex = factor(goback.nochrom$sex, 
                                         levels = c(1,2),
                                         labels = c('Male','Female')),
                            m.age = goback.nochrom$m.age,
                            state = goback.nochrom$state)
  
  fit <- survfit(Surv(time, cancer) ~ defect, data = goback.surv)
  
  if (i == 1){ # y-axis limits are different for "any cancer" than for other panels.
    
    new.plot <- ggsurvplot(fit, 
                           fun = 'event',
                           xlim = c(0,18),
                           ylim = c(0,0.015),
                           ylab = paste("Cumulative Incidence of any",outcomes.fancy[i]),
                           xlab = "Time in Years",
                           font.x = c(15, 'bold', 'black'),
                           font.y = c(15, 'bold', 'black'),
                           font.legend = c(15, 'bold', 'black'),
                           font.tickslab = c(15, 'bold', 'black'),
                           conf.int = FALSE,
                           linetype = 'strata', 
                           legend.labs = c('No defect', '1 defect', '2 defects', '3 defects', '4 or more defects'))
    
    the.plots.thicken[i] <- new.plot
    
  }
  
  else{
    
    new.plot <- ggsurvplot(fit, 
                           fun = 'event',
                           xlim = c(0,18),
                           ylim = c(0,0.01),
                           ylab = paste("Cumulative Incidence of any",outcomes.fancy[i]),
                           xlab = "Time in Years",
                           font.x = c(15, 'bold', 'black'),
                           font.y = c(15, 'bold', 'black'),
                           font.legend = c(15, 'bold', 'black'),
                           font.tickslab = c(15, 'bold', 'black'),
                           conf.int = FALSE,
                           linetype = 'strata', 
                           legend.labs = c('No defect', '1 defect', '2 defects', '3 defects', '4 or more defects'))
    
    the.plots.thicken[i] <- new.plot
  }
  
}

names(the.plots.thicken) <- outcomes




# Scratch paper -----------------------------------------------------------

#' If needed, this line will perform a log-rank test for trend on the data.
surv_pvalue(fit)