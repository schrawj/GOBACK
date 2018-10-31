#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.10.25.
#' 
#' One reviewer asked for mean/median and total person time contributed by
#' children in each group.
#' 
#' Simple enough.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(tictoc); require(dplyr)

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20180829.rdata')

tab <- aggregate(person.yrs ~ any.birthdefect + cancer, data = goback, median)
tab <- rename(tab, median.person.yrs = person.yrs)

tab2 <- aggregate(person.yrs ~ any.birthdefect + cancer, data = goback, mean)
tab2 <- rename(tab2, mean.person.yrs = person.yrs)

tab3 <- aggregate(person.yrs ~ any.birthdefect + cancer, data = goback, sum)
tab3 <- rename(tab3, total.person.yrs = person.yrs)

tab <- left_join(left_join(tab, tab2, by = c('any.birthdefect','cancer')), 
                           tab3, 
                           by = c('any.birthdefect','cancer'))

for (i in 3:5){
  tab[,i] <- floor(tab[,i])
}

write.csv(tab, file = 'W:/Old_genepi2/Jeremy/GOBACK/R outputs/goback.person.years.by.birth.defect.and.cancer.status.csv', row.names = FALSE)
