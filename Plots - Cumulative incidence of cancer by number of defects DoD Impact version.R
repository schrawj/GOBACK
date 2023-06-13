#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2019.08.13.
#' 
#' Minimal production version of the cumulative incidence plot for all 
#' cancers as provided in the Fall 2019 DoD Impact Award grant submission.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(survival); require(survminer); require(ggsci) 

load('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/goback.nochrom.v20190606.rdata')

goback.nochrom$majordefect.cat <- factor(ifelse(goback.nochrom$majordefect.total == 0, 0,
                                                ifelse(goback.nochrom$majordefect.total == 1, 1,
                                                       ifelse(goback.nochrom$majordefect.total == 2, 2, 
                                                              ifelse(goback.nochrom$majordefect.total == 3, 3,
                                                                     ifelse(goback.nochrom$majordefect.total >= 4, 4, NA))))),
                                         
                                         levels = c(0:4),
                                         labels = c('0', '1', '2', '3', '4 or more'))


goback.surv <- data.frame(studyid = goback.nochrom$studyid,
                          cancer = goback.nochrom$cancer, 
                          time = goback.nochrom$person.yrs, 
                          defect = goback.nochrom$majordefect.cat,
                          sex = factor(goback.nochrom$sex, 
                                       levels = c(1,2),
                                       labels = c('Male','Female')),
                          m.age = goback.nochrom$m.age,
                          state = goback.nochrom$state)

fit <- survfit(Surv(time, cancer) ~ defect, data = goback.surv)

#' shape is used for events/censoring. Empty character string to suppress. 124 for vertical line.
shape <- ''

new.plot <- ggsurvplot(fit, 
                       fun = 'event',
                       font.title = c(20, 'bold','black'),
                       xlim = c(0,18),
                       ylab = 'Cumulative Incidence of Cancer (%)',
                       xlab = "Time (Years)",
                       font.x = c(20, 'bold', 'black'),
                       font.y = c(20, 'bold', 'black'),
                       font.tickslab = c(15, 'bold', 'black'),
                       conf.int = FALSE,
                       censor.shape = shape,
                       risk.table = F)

new.plot$plot <- new.plot$plot + 
                 theme(legend.position = 'none',
                       panel.grid.major.y = element_line(color = 'grey75')) +
                 scale_color_jama() +
                 scale_y_continuous(limits= c(0, 0.015),
                                    breaks = c(0, 0.005, 0.01,0.015),
                                    labels = c(0, 0.5, 1, 1.5))

print(new.plot)
