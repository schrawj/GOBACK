#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2018.04.26.
#' 
#' Counts of the number of CNS tumors diagnosed at each age overall and in Texas.
#' 
#' Philip suggested we do a sensitivity analysis looking at these cancers only in
#' children diagnosed after 4 years and only in Texas.  I'm skeptical that this will be
#' very useful.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------


load("Z:/Jeremy/GOBACK/Datasets/goback.v20180419.rdata")

require(dplyr)

for (i in c(111,114,126,132)){
  
  tmp <- names(goback.nochrom[i])
  
  tab <- table(floor(goback.nochrom$person.yrs), goback.nochrom[, i], deparse.level = 0)
  
  out <- data.frame(age = rownames(tab),
                    dx.count = tab[,2])
  
  write.csv(tab, file = paste0('Z:/Jeremy/GOBACK/R outputs/CBT person years/',tmp,'.personyrs.csv'))
  
}

goback.nochrom <- filter(goback.nochrom, state == 'TX')

for (i in c(111,114,126,132)){
  
  tmp <- names(goback.nochrom[i])
  
  tab <- table(floor(goback.nochrom$person.yrs), goback.nochrom[, i], deparse.level = 0)
  
  out <- data.frame(age = rownames(tab),
                    dx.count = tab[,2])
  
  write.csv(tab, file = paste0('Z:/Jeremy/GOBACK/R outputs/CBT person years/',tmp,'.personyrs.tx.csv'))
  
}

rm(list = ls()); gc()

