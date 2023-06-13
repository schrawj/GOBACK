require(descr); require(magrittr); require(tidyverse)

setwd('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/')

goback <- readRDS('Datasets/goback.v20220215.rds')

goback %<>% mutate(state = as.factor(state))

bd.report <- data.frame()

cancer.report <- data.frame()

bd.vars <- c(15, 19:95)

cancer.vars <- c(96, 100:139)

bd.pop.tab <- table(goback$state)/10000

cancer.pop.tab <- table(goback$state)/100000

tabulate.percentages <- function(my.vector, denominator, output.df){
  
  for (i in my.vector){
    
    print(paste('tabulating prevalence for',names(goback)[i]))
    
    tab <- crosstab(goback[,i], goback$state, prop.c = T, drop.levels = F, plot = F)
    
    new.report <- data.frame(variable = names(goback)[i],
                             AR = round(tab$tab[2,1]/denominator[1], 1),
                             MA = round(tab$tab[2,2]/denominator[2], 1),
                             MI = round(tab$tab[2,3]/denominator[3], 1),
                             NC = round(tab$tab[2,4]/denominator[4], 1),
                             OK = round(tab$tab[2,5]/denominator[5], 1),
                             TX = round(tab$tab[2,6]/denominator[6], 1)
                             )

    output.df <- rbind(output.df, new.report)
    
  }
  
  return(output.df)

}

cancer.report <- tabulate.percentages(cancer.vars, cancer.pop.tab, cancer.report)

bd.report <- tabulate.percentages(bd.vars, bd.pop.tab, bd.report)

saveRDS(cancer.report, 'Datasets/Expanded datasets/cancer.frequencies.by.state.v20220215.rds')
write_csv(cancer.report, 'R outputs/cancer.frequencies.by.state.v20220215.csv')

saveRDS(bd.report, 'Datasets/Expanded datasets/birth.defect.frequencies.by.state.v20220215.rds')
write_csv(bd.report, 'R outputs/birth.defect.frequencies.by.state.v20220215.csv')
