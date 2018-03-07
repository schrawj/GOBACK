#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' GOBACK logistic regression modeling
#' 
#' Discussed initial modeling approach at meeting on 10/23/2017.
#' 
#' Will generate a table of logistic regression models for all potential 
#' cancer x birth defect associations with at least 5 cormorbid cases, 
#' and heatmaps based on the one in the WA state paper.
#' 
#'  Two sets of tables: one for kids with chromosomal birth defects, one 
#'  for kids with non-chromosomal birth defects.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------



# prep environment --------------------------------------------------------
require(dplyr)

#' For desktop
setwd('Z:/Jeremy/GOBACK/Datasets/Combined Arkansas Michigan and Texas/')



# Logistic regression: cancer in children w/o chromosomal defects ---------
load('./ar.mi.tx.nochromosomaldefects.v20171025.1.rdata')

for (i in 22:103){
  
  tmp <- table(armitx.nochrom[,i], armitx.nochrom$cancer1)
  tmp <- which(tmp[2, ] >= 5)
  tmp <- tmp + 117
  
  if (length(tmp) > 0){
    
    for (j in tmp){
      
      z <- names(armitx.nochrom[i])
      y <- names(armitx.nochrom[j])
      
      x <- glm(armitx.nochrom[,j] ~ armitx.nochrom[,i], data = armitx.nochrom, family = binomial(link = 'logit'))
      x.summary <- summary(x)$coefficients
      tab <- as.numeric(table(armitx.nochrom[,i], armitx.nochrom[,j])[2,2])
      estimates <- data.frame(defect = z, cancer = y, OR = exp(x.summary[2,1]), 
                              ci.lower = exp(x.summary[2,1]-(1.96*x.summary[2,2])), 
                              ci.upper = exp(x.summary[2,1]+(1.96*x.summary[2,2])),
                              num.pos.events = tab)
      
      write.table(estimates, file = 'C:/Users/schraw/Desktop/goback models/BD-CC associations.csv', sep=',', append = TRUE, 
                  row.names = FALSE, col.names = FALSE)
    }
    
  }
  
  else{
    
    sink(file = 'C:/Users/schraw/Desktop/goback models/list of defects with no models.txt', append = TRUE)
    
    print(paste('There were no cancers with 5 or more comorbid cases for',names(armitx.nochrom[i])))
    
    sink()
  }
  
}

rm(estimates, x.summary, i, j, tab, tmp, x, y, z)
gc()

#' Models for individual non-chromosomal birth defects and [cancer].any variables
for (i in 22:103){
  for (j in 148:157){
    
    z <- names(armitx.nochrom[i])
    y <- names(armitx.nochrom[j])
    comorbid.cases <- table(armitx.nochrom[,i], armitx.nochrom[,j])[2,2]
    
    if (comorbid.cases > 5){
      x <- glm(armitx.nochrom[,j] ~ armitx.nochrom[,i], data = armitx.nochrom, family = binomial(link='logit'))
      x.summary <- summary(x)$coefficients
      estimates <- data.frame(defect = z, cancer = y, OR = exp(x.summary[2,1]), 
                              ci.lower = exp(x.summary[2,1]-(1.96*x.summary[2,2])), 
                              ci.upper = exp(x.summary[2,1]+(1.96*x.summary[2,2])),
                              num.pos.events = as.numeric(comorbid.cases))
      write.table(estimates, file = 'C:/Users/schraw/Desktop/goback models/BD-CC associations.csv', sep=',', append = TRUE, 
                  row.names = FALSE, col.names = FALSE)
    }
    
    else{
      
      sink(file = 'C:/Users/schraw/Desktop/goback models/list of defects with no models.txt', append = TRUE)
      
      print(paste('There were less than five comorbid instances of',z,'and',y))
      
      sink()
      
    }
  }
}

rm(armitx.nochrom, i, j, estimates, x, x.summary, y, z, comorbid.cases)
gc()



# Logistic regression: cancer in children w/chromosomal defects -----------
load('./ar.mi.tx.chromosomaldefects.v20171025.1.rdata')

#' Models for individual chromosomal birth defects and individual cancers
for (i in 104:111){
  
  tmp <- table(armitx.chrom[,i], armitx.chrom$cancer1)
  tmp <- which(tmp[2, ] >= 5)
  tmp <- tmp + 117
  
  if (length(tmp) > 0){
    
    for (j in tmp){
      
      z <- names(armitx.chrom[i])
      y <- names(armitx.chrom[j])
      
      x <- glm(armitx.chrom[,j] ~ armitx.chrom[,i], data = armitx.chrom, family = binomial(link = 'logit'))
      x.summary <- summary(x)$coefficients
      tab <- as.numeric(table(armitx.chrom[,i], armitx.chrom[,j])[2,2])
      estimates <- data.frame(defect = z, cancer = y, OR = exp(x.summary[2,1]), 
                              ci.lower = exp(x.summary[2,1]-(1.96*x.summary[2,2])), 
                              ci.upper = exp(x.summary[2,1]+(1.96*x.summary[2,2])),
                              num.pos.events = tab)
      
      write.table(estimates, file = 'C:/Users/schraw/Desktop/goback models/BD-CC associations.csv', sep=',', append = TRUE, 
                  row.names = FALSE, col.names = FALSE)
    }
    
  }
  
  else{
    
    sink(file = 'C:/Users/schraw/Desktop/goback models/list of defects with no models.txt', append = TRUE)
    
    print(paste('There were no cancers with 5 or more comorbid cases for',names(armitx.chrom[i])))
    
    sink()
  }
  
}

rm(estimates, x.summary, i, j, tab, tmp, x, y, z)
gc()

#' Models for individual chromosomal birth defects and [cancer].any variables
for (i in 104:111){
  for (j in 148:157){
    
    z <- names(armitx.chrom[i])
    y <- names(armitx.chrom[j])
    comorbid.cases <- table(armitx.chrom[,i], armitx.chrom[,j])[2,2]
    
    if (comorbid.cases > 5){
      x <- glm(armitx.chrom[,j] ~ armitx.chrom[,i], data = armitx.chrom, family = binomial(link='logit'))
      x.summary <- summary(x)$coefficients
      estimates <- data.frame(defect = z, cancer = y, OR = exp(x.summary[2,1]), 
                              ci.lower = exp(x.summary[2,1]-(1.96*x.summary[2,2])), 
                              ci.upper = exp(x.summary[2,1]+(1.96*x.summary[2,2])),
                              num.pos.events = as.numeric(comorbid.cases))
      write.table(estimates, file = 'C:/Users/schraw/Desktop/goback models/BD-CC associations.csv', sep=',', append = TRUE, 
                  row.names = FALSE, col.names = FALSE)
    }
    
    else{
      
      sink(file = 'C:/Users/schraw/Desktop/goback models/list of defects with no models.txt', append = TRUE)
      
      print(paste('There were less than five comorbid instances of',z,'and',y))
      
      sink()
      
    }
  }
}

rm(armitx.chrom, estimates, x, x.summary, y, z, comorbid.cases)
gc()

# Model QC: Re-run a few models manually ----------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Some of these ORs are quite dramatic.
#' 
#' Hopefully that reflects the biology of these associations.
#' 
#' Check the diagnostic codes for some cancer and birth defects diagnoses
#' to make sure there are no errors in our variables.  Just want to rule
#' out that there is some error in the input to the data.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

model <- glm(armitx.nochrom$hepato ~ armitx.nochrom$atrialseptaldefect, data = armitx.nochrom, family = binomial(link = 'logit'))
summary(model)

model <- glm(armitx.nochrom$gct.any ~ armitx.nochrom$digestivesystem.other.major, data = armitx.nochrom, family = binomial(link = 'logit'))
summary(model)

model <- glm(armitx.nochrom$pns.any ~ armitx.nochrom$musculoskelsys.other.major, data = armitx.nochrom, family = binomial(link = 'logit'))
summary(model)

model <- glm(armitx.nochrom$all ~ armitx.nochrom$microcephalus, data = armitx.nochrom, family = binomial(link = 'logit'))
summary(model)

rm(model)



# Model QC: Verifying birth defects diagnoses -----------------------------

for (i in 22:112){
  tmp <- table(armitx.nochrom[,i], useNA = 'always')
  tmp <- data.frame(defect = names(armitx.nochrom[i]),
                    negative.for.def = tmp[1],
                    positive.for.def = tmp[2],
                    num.actually.na = tmp[3],
                    num.should.be.na = 479467-(tmp[2]))
  write.table(tmp, file = 'C:/Users/schraw/Desktop/goback models/number of missing observations by defect.csv', sep= ',', 
              row.names = FALSE, col.names = FALSE, append = TRUE)
}

rm(i, tmp)

ids <- select(armitx.nochrom, studyid)

#' Look through some of the original ICD codes in MI data and verify they match the number of 
#' children DX'd with that anomaly.
load("Z:/Jeremy/GOBACK/Datasets/Michigan/mi.birthdefects.codes.rdata")

mi.bd$ebstein.code <- as.numeric(NA)

for (i in mi.bd){
  for (j in 110:133){
    mi.bd$ebstein.code <- ifelse(is.na(mi.bd$ebstein.code) & round(mi.bd[,j], digits = 1) == 746.2, 1, mi.bd$ebstein.code)
  }
}

table(mi.bd$ebstein.code, useNA = 'ifany')
table(mi.bd$EbsteinAnomaly, useNA = 'ifany')

tmp$sb.code <- as.numeric(NA)
sb.codes <- c(741.0, 741.9)

for (i in tmp){
  for (j in 110:133){
    tmp$sb.code <- ifelse(is.na(tmp$sb.code) & round(tmp[,j], digits = 1) %in% sb.codes, 1, tmp$sb.code)
    
  }
}

table(tmp$sb.code)

rm(mi.bd, tmp, sb.codes, i, j)



# Model QC: verifying some cancer diagnoses -------------------------------
hepato <- filter(filter(armitx.nochrom, cancer1 == 'hepato'), state == 'TX')
hepato <- c(hepato$studyid)

nhl <- filter(filter(armitx.nochrom, cancer1 == 'nhl'), state == 'TX')
nhl <- c(nhl$studyid)

all <- filter(filter(armitx.nochrom, cancer1 == 'all'), state == 'TX')
all <- c(all$studyid)

load('Z:/Jeremy/GOBACK/Datasets/Texas/tx.cancer1.codes.rdata')

tx.can$birthid <- paste0('tx',tx.can$birthid)

tx.hepato <- tx.can[tx.can$birthid %in% hepato, ]
unique(tx.hepato$morph31)
table(tx.hepato$morph31, useNA = 'ifany')

tx.nhl <- tx.can[tx.can$birthid %in% nhl, ]
tx.nhl <- arrange(tx.nhl, morph31)
unique(tx.nhl$morph31)
print(tx.nhl[,2:3])

tx.all <- tx.can[tx.can$birthid %in% all, ]
tx.all <- arrange(tx.all, morph31)
unique(tx.all$morph31)

tmp <- filter(tx.all, morph31 == 9811)
unique(tmp$site.code1)

tmp <- filter(tx.all, morph31 == 9823)
unique(tmp$site.code1)

rm(tx.can, tmp, tx.all, tx.nhl, tx.hepato, all, hepato, nhl)


