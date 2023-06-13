#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.08.23.
#' 
#' Cumulative incidence of any cancer, any heme cancer, any cns cancer, and
#' any non-cns solid tissue cancer by number of major defects.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(survival); require(survminer); require(ggsci) 
require(extrafont); require(grid); require(gridExtra)

loadfonts(device = 'win')

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.nochrom.v20190429.rdata')

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
outcomes.fancy <- c('Cancer','Hematologic Malignancies', "CNS Tumors", 'Non-CNS Solid Tumors')
panels <- c('A.', 'B.','C.','D.')

the.plots.thicken <- list()

goback.surv <- data.frame(studyid = goback.nochrom$studyid,
                          time = goback.nochrom$person.yrs, 
                          defect = goback.nochrom$majordefect.cat,
                          sex = factor(goback.nochrom$sex, 
                                       levels = c(1,2),
                                       labels = c('Male','Female')),
                          m.age = goback.nochrom$m.age,
                          state = goback.nochrom$state)

for (i in 1:length(outcomes)){
  
  goback.surv$cancer = goback.nochrom[,outcomes[i]]
  
  fit <- survfit(Surv(time, cancer) ~ defect, data = goback.surv)
  
  if (i == 1){ # y-axis limits are different for "any cancer" than for other panels.
    
    new.plot <- ggsurvplot(fit, 
                           fun = 'event',
                           title = paste(panels[i],'Cumulative Incidence of',outcomes.fancy[i]),
                           font.title = c(20, 'bold','black'),
                           legend.labs = c('None', 'One', 'Two', 'Three', 'Four or more'),
                           xlim = c(0,18),
                           ylim = c(0,0.015),
                           ylab = '',
                           xlab = "Time in Years",
                           font.x = c(20, 'bold', 'black'),
                           font.y = c(20, 'bold', 'black'),
                           font.tickslab = c(20, 'bold', 'black'),
                           conf.int = FALSE,
                           linetype = 'strata',
                           risk.table = TRUE)
    
  }
  
  else{
    
    new.plot <- ggsurvplot(fit, 
                           fun = 'event',
                           title = paste(panels[i],'Cumulative Incidence of',outcomes.fancy[i]),
                           font.title = c(20, 'bold','black'),
                           legend.labs = c('None', 'One', 'Two', 'Three', 'Four or more'),
                           xlim = c(0,18),
                           ylim = c(0,0.01),
                           ylab = '',
                           xlab = "Time in Years",
                           font.x = c(20, 'bold', 'black'),
                           font.y = c(20, 'bold', 'black'),
                           font.tickslab = c(20, 'bold', 'black'),
                           conf.int = FALSE,
                           linetype = 'strata',
                           risk.table = TRUE)

  }

  new.plot$plot <- new.plot$plot + theme(legend.position = 'none')
  
  new.plot$table <- new.plot$table +
    ggtitle('Number of Individuals at Risk, According to Number of Major Birth Defects') + 
    ylab('No. of Birth Defects') +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.text = element_blank(),
          legend.position = 'none',
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 15, face = 'bold'),
          plot.title = element_text(size = 15, face = 'bold', hjust = 0.5))
  
  the.plots.thicken[[i]] <- new.plot
  
}



# Scratch paper -----------------------------------------------------------

#' This version prints the legend.
new.plot <- ggsurvplot(fit, 
                       fun = 'event',
                       title = paste(panels[i],'Cumulative Incidence of',outcomes.fancy[i]),
                       font.title = c(25, 'bold','black'),
                       xlim = c(0,18),
                       ylim = c(0,0.015),
                       ylab = '',
                       xlab = "Time in Years",
                       font.x = c(20, 'bold', 'black'),
                       font.y = c(20, 'bold', 'black'),
                       font.legend = c(20, 'bold', 'black'),
                       font.tickslab = c(20, 'bold', 'black'),
                       conf.int = FALSE,
                       linetype = 'strata', 
                       legend.title = 'Number of birth defects',
                       legend.labs = c('None', '1', '2', '3', '4 or more'))

#' Performs a log-rank test for trend on the data.
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
  p.trend <- surv_pvalue(fit)
  print(paste('p-value, test for trend, number of birth defects and risk of', outcomes[i],'is',p.trend))
  
}

#' Solid instead of dashed lines, and smaller font for events/censoring.
the.plots.thicken <- list()

goback.surv <- data.frame(studyid = goback.nochrom$studyid,
                          time = goback.nochrom$person.yrs, 
                          defect = goback.nochrom$majordefect.cat,
                          sex = factor(goback.nochrom$sex, 
                                       levels = c(1,2),
                                       labels = c('Male','Female')),
                          m.age = goback.nochrom$m.age,
                          state = goback.nochrom$state)

for (i in 1:length(outcomes)){
  
  goback.surv$cancer = goback.nochrom[,outcomes[i]]
  
  fit <- survfit(Surv(time, cancer) ~ defect, data = goback.surv)
  
  if (i == 1){ # y-axis limits are different for "any cancer" than for other panels.
    
    new.plot <- ggsurvplot(fit, 
                           fun = 'event',
                           title = paste(panels[i],'Cumulative Incidence of',outcomes.fancy[i]),
                           legend.labs = c('None', 'One', 'Two', 'Three', 'Four or more'),
                           xlim = c(0,18),
                           ylim = c(0,0.015),
                           ylab = '',
                           xlab = "Time, y",
                           conf.int = FALSE,
                           risk.table = TRUE,
                           censor.shape= 124, # Set to empty character string to suppress. Set to 124 for vertical line.
                           censor.size = 0.75) # Default is 4.5.

  }
  
  else{
    
    new.plot <- ggsurvplot(fit, 
                           fun = 'event',
                           title = paste(panels[i],'Cumulative Incidence of',outcomes.fancy[i]),
                           legend.labs = c('None', 'One', 'Two', 'Three', 'Four or more'),
                           xlim = c(0,18),
                           ylim = c(0,0.01),
                           ylab = '',
                           xlab = "Time, y",
                           conf.int = FALSE,
                           risk.table = TRUE,
                           censor.shape= 124, # Set to empty character string to suppress. Set to 124 for vertical line.
                           censor.size = 0.75) # Default is 4.5.
                                                   
    
  }
  
  new.plot$plot <- new.plot$plot + 
                   theme(legend.position = 'none',
                         panel.grid.major.y = element_line(color = 'grey75'),
                         text = element_text(family = 'Arial')) +
                   scale_color_jama() 

  new.plot$table <- new.plot$table +
    ggtitle('No. at Risk') + 
    theme(legend.text = element_blank(),
          legend.position = 'none',
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(color = 'black'),
          plot.title = element_text(hjust = 0),
          text = element_text(color = 'black', family = 'Arial'))

  the.plots.thicken[[i]] <- new.plot

}

grid.arrange(the.plots.thicken[[1]]$plot,
             the.plots.thicken[[2]]$plot,
             the.plots.thicken[[3]]$plot,
             the.plots.thicken[[4]]$plot, 
             the.plots.thicken[[1]]$table,
             layout_matrix = (rbind(c(1,2),
                                    c(1,2),
                                    c(3,4),
                                    c(3,4),
                                    c(5,5)))
             )
