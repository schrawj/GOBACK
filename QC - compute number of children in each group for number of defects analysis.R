#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.10.24.
#' 
#' Reviewers asked to know the number of children with 1, 2, 3, and 4+ 
#' defects in the number of birth defects analysis.
#' 
#' This analysis is performed in the 
#' 'Plots - Cumulative incidence of cancer by number of defects' script. 
#' 
#' This is basically a companion piece to that script.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.nochrom.v20180829.rdata')

goback.nochrom$majordefect.cat <- factor(ifelse(goback.nochrom$majordefect.total == 0, 0,
                                                ifelse(goback.nochrom$majordefect.total == 1, 1,
                                                       ifelse(goback.nochrom$majordefect.total == 2, 2, 
                                                              ifelse(goback.nochrom$majordefect.total == 3, 3,
                                                                     ifelse(goback.nochrom$majordefect.total >= 4, 4, NA))))),
                                         levels = c(0:4),
                                         labels = c('0', '1', '2', '3', '4 or more'))

goback.nochrom$any.heme.cancer <- ifelse(goback.nochrom$leu.any == 1 | goback.nochrom$lym.any == 1, 1, 0)

solid.tumors <- unique(goback.nochrom$cancer1)
solid.tumors <- subset(solid.tumors, !(solid.tumors %in% c(NA, 'all','leu.other','aml','hl','nhl','cns.other','medullo','pnet','lym.other',
                                                           'gct.intra','astro','ependymoma')))
goback.nochrom$any.non.cns.solid.tumor <- ifelse(goback.nochrom$cancer1 %in% solid.tumors, 1, 0)

outcomes <- c('cancer','any.heme.cancer','cns.any','any.non.cns.solid.tumor')

for (i in outcomes){
  
  tab <- table(goback.nochrom$majordefect.cat, goback.nochrom[,i], useNA = 'ifany', dnn = c('Number of major defects',i))
  
  write.xlsx(tab, file = 'W:/Old_genepi2/Jeremy/GOBACK/R outputs/Cancer risk by number of birth defects/n.by.cancer.and.number.of.birth.defects.xlsx',
             sheetName = i, row.names = FALSE, append = TRUE)
  
}