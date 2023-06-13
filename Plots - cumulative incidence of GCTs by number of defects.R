
require(dplyr)
require(survival); require(survminer); require(ggsci) 
require(extrafont); require(grid); require(gridExtra)

load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/goback.nochrom.v20191114.rdata")

goback.nochrom$majordefect.cat <- factor(ifelse(goback.nochrom$majordefect.total == 0, 0,
                                                ifelse(goback.nochrom$majordefect.total == 1, 1,
                                                       ifelse(goback.nochrom$majordefect.total == 2, 2, 
                                                              ifelse(goback.nochrom$majordefect.total == 3, 3,
                                                                     ifelse(goback.nochrom$majordefect.total >= 4, 4, NA))))),
                                         
                                         levels = c(0:4),
                                         labels = c('0', '1', '2', '3', '4 or more'))

goback.surv <- data.frame(studyid = goback.nochrom$studyid,
                          time = goback.nochrom$person.yrs, 
                          cancer = goback.nochrom$gct.any,
                          defect = goback.nochrom$majordefect.cat,
                          sex = factor(goback.nochrom$sex, 
                                       levels = c(1,2),
                                       labels = c('Male','Female')),
                          m.age = goback.nochrom$m.age,
                          state = goback.nochrom$state)

fit <- survfit(Surv(time, cancer) ~ defect, data = goback.surv)

#' shape will be used for events/censoring. Empty character string to suppress. 124 for vertical line.
shape <- ''

new.plot <- ggsurvplot(fit, 
                       fun = 'event',
                       #title = paste(panels[i],'Cumulative Incidence of',outcomes.fancy[i]),
                       #legend.labs = c('None', 'One', 'Two', 'Three', 'Four or more'),
                       xlim = c(0,18),
                       ylim = c(0,0.01),
                       ylab = '',
                       xlab = "Time, y",
                       surv.scale = 'percent',
                       conf.int = FALSE,
                       risk.table = TRUE,
                       censor.shape = shape, # Set to empty character string to suppress. Set to 124 for vertical line.
                       censor.size = 0.75) # Default is 4.5.
print(new.plot)
