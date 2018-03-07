#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' GOBACK logistic regression modeling
#' 
#' Liftover of code used to generate preliminary data in AR, MI and TX.
#' 
#' Will generate logistic regression models for all potential 
#' cancer-birth defect associations with at least 5 cormorbid cases.
#' 
#' Two sets of tables: one for kids with chromosomal birth defects, one 
#' for kids with non-chromosomal birth defects.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------



# prep environment --------------------------------------------------------
require(dplyr)

setwd('Z:/Jeremy/GOBACK/Datasets/')



# Logistic regression: cancer in children w/o chromosomal defects ---------
load('goback.no.chrom.v20171211.1.rdata')

for (i in 22:102){
  
  tmp <- table(goback.nochrom[,i], goback.nochrom$cancer1)
  tmp <- which(tmp[2, ] >= 5)
  tmp <- tmp + 107
  
  if (length(tmp) > 0){
    
    for (j in tmp){
      
      z <- names(goback.nochrom[i])
      y <- names(goback.nochrom[j])
      
      x <- glm(goback.nochrom[,j] ~ goback.nochrom[,i], data = goback.nochrom, family = binomial(link = 'logit'))
      x.summary <- summary(x)$coefficients
      tab <- as.numeric(table(goback.nochrom[,i], goback.nochrom[,j])[2,2])
      estimates <- data.frame(defect = z, cancer = y, OR = exp(x.summary[2,1]), 
                              ci.lower = exp(x.summary[2,1]-(1.96*x.summary[2,2])), 
                              ci.upper = exp(x.summary[2,1]+(1.96*x.summary[2,2])),
                              p.value = x.summary[2,4],
                              num.pos.events = tab)
      
      write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/bd.cc.models.csv', sep=',', append = TRUE, 
                  row.names = FALSE, col.names = FALSE)
    }
    
  }
  
  else{
    
    sink(file = 'Z:/Jeremy/GOBACK/R Outputs/list of unmodeled defects.txt', append = TRUE)
    
    print(paste('There were no cancers with 5 or more comorbid cases for',names(goback.nochrom[i])))
    
    sink()
  }
  
}

rm(estimates, x.summary, i, j, tab, tmp, x, y, z)
gc()

#' Models for individual non-chromosomal birth defects and [cancer].any variables
for (i in 22:102){
  for (j in 138:147){
    
    z <- names(goback.nochrom[i])
    y <- names(goback.nochrom[j])
    comorbid.cases <- table(goback.nochrom[,i], goback.nochrom[,j])[2,2]
    
    if (comorbid.cases >= 5){
      x <- glm(goback.nochrom[,j] ~ goback.nochrom[,i], data = goback.nochrom, family = binomial(link='logit'))
      x.summary <- summary(x)$coefficients
      estimates <- data.frame(defect = z, cancer = y, OR = exp(x.summary[2,1]), 
                              ci.lower = exp(x.summary[2,1]-(1.96*x.summary[2,2])), 
                              ci.upper = exp(x.summary[2,1]+(1.96*x.summary[2,2])),
                              p.value = x.summary[2,4],
                              num.pos.events = as.numeric(comorbid.cases))
      write.table(estimates, file = 'C:/Users/schraw/Desktop/goback models/bd.cc.models.csv', sep=',', append = TRUE, 
                  row.names = FALSE, col.names = FALSE)
    }
    
    else{
      
      sink(file = 'Z:/Jeremy/GOBACK/R Outputs/list of unmodeled defects.txt', append = TRUE)
      
      print(paste('There were less than five comorbid instances of',z,'and',y))
      
      sink()
      
    }
  }
}

rm(goback.nochrom, i, j, estimates, x, x.summary, y, z, comorbid.cases)
gc()



# Logistic regression: cancer in children w/chromosomal defects -----------
load('./goback.chrom.v20171211.1.rdata')

#' Models for individual chromosomal birth defects and individual cancers
for (i in 96:101){
  
  tmp <- table(goback.chrom[,i], goback.chrom$cancer1)
  tmp <- which(tmp[2, ] >= 5)
  tmp <- tmp + 107
  
  if (length(tmp) > 0){
    
    for (j in tmp){
      
      z <- names(goback.chrom[i])
      y <- names(goback.chrom[j])
      
      x <- glm(goback.chrom[,j] ~ goback.chrom[,i], data = goback.chrom, family = binomial(link = 'logit'))
      x.summary <- summary(x)$coefficients
      tab <- as.numeric(table(goback.chrom[,i], goback.chrom[,j])[2,2])
      estimates <- data.frame(defect = z, cancer = y, OR = exp(x.summary[2,1]), 
                              ci.lower = exp(x.summary[2,1]-(1.96*x.summary[2,2])), 
                              ci.upper = exp(x.summary[2,1]+(1.96*x.summary[2,2])),
                              p.value = x.summary[2,4],
                              num.pos.events = tab)
      
      write.table(estimates, file = 'Z:/Jeremy/GOBACK/R Outputs/bd.cc.models.csv', sep=',', append = TRUE, 
                  row.names = FALSE, col.names = FALSE)
    }
    
  }
  
  else{
    
    sink(file = 'Z:/Jeremy/GOBACK/R Outputs/list of unmodeled defects.txt', append = TRUE)
    
    print(paste('There were no cancers with 5 or more comorbid cases for',names(goback.chrom[i])))
    
    sink()
  }
  
}

rm(estimates, x.summary, i, j, tab, tmp, x, y, z)
gc()

#' Models for individual chromosomal birth defects and [cancer].any variables
for (i in 96:101){
  for (j in 138:147){
    
    z <- names(goback.chrom[i])
    y <- names(goback.chrom[j])
    comorbid.cases <- table(goback.chrom[,i], goback.chrom[,j])[2,2]
    
    if (comorbid.cases >= 5){
      x <- glm(goback.chrom[,j] ~ goback.chrom[,i], data = goback.chrom, family = binomial(link='logit'))
      x.summary <- summary(x)$coefficients
      estimates <- data.frame(defect = z, cancer = y, OR = exp(x.summary[2,1]), 
                              ci.lower = exp(x.summary[2,1]-(1.96*x.summary[2,2])), 
                              ci.upper = exp(x.summary[2,1]+(1.96*x.summary[2,2])),
                              p.value = x.summary[2,4],
                              num.pos.events = as.numeric(comorbid.cases))
      write.table(estimates, file = 'C:/Users/schraw/Desktop/goback models/BD-CC associations.csv', sep=',', append = TRUE, 
                  row.names = FALSE, col.names = FALSE)
    }
    
    else{
      
      sink(file = 'Z:/Jeremy/GOBACK/R Outputs/list of unmodeled defects.txt', append = TRUE)
      
      print(paste('There were less than five comorbid instances of',z,'and',y))
      
      sink()
      
    }
  }
}

rm(goback.chrom, estimates, x, x.summary, i, j, y, z, comorbid.cases)
gc()