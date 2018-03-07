#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#'                                       
#'                            MICHIGAN + TEXAS
#'                                         
#' Certain issues exist in both the TX and MI data.  For efficiency, these
#' are addressed once after joining the files rather than separately before
#' joining the files.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------



# Prep environment --------------------------------------------------------
require(dplyr)

setwd('Z:/Jeremy/GOBACK/Datasets/Combined Michigan and Texas/')



# User-defined functions --------------------------------------------------
gen.cancer.var <- function(x){
  x <- 0
}
populate.cancer.var <- function(x,y,z){
  if(missing(z)){
    ifelse((x | y) == 1, 1, 0)
  }
  else{
    ifelse((x | y | z) == 1, 1, 0)
  }
}



# Load in and merge data --------------------------------------------------
load("Z:/Jeremy/GOBACK/Datasets/Michigan/mi.v20171006.4.rdata")
load('Z:/Jeremy/GOBACK/Datasets/Texas/tx.v20171019.3.rdata')

mitx <- rbind(tx, mi)

#' Appears that the structures of the data frames were unaffected by the merge.
flag <- as.logical(identical(str(mitx), str(mi)))
flag <- as.logical(identical(str(mitx), str(tx)))

rm(flag, mi, tx)

save(mitx, file = './mi.tx.v20171019.1.rdata')



# Demographics variables --------------------------------------------------
load('./mi.tx.v20171019.1.rdata')

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Implementing rules we established previously for demographics variables.
#' - Maternal age is set to NA if < 13 or > 50.
#' - Gestational age is set to NA if < 22 or > 44.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
mitx$gest.age <- ifelse(mitx$gest.age < 22 | mitx$gest.age > 44, NA, mitx$gest.age)
range(mitx$gest.age, na.rm = TRUE)
table(mitx$gest.age, useNA = 'ifany')

mitx$m.age <- ifelse(mitx$m.age < 13 | mitx$m.age > 50, NA, mitx$m.age)
range(mitx$m.age, na.rm = TRUE)
table(mitx$m.age, useNA = 'ifany')

save(mitx, file = './mi.tx.v20171020.1.rdata')



# Compute individual cancer variables -------------------------------------
load("./mi.tx.v20171020.1.rdata")

l <- c(colnames(mitx[118:157]))
mitx <- mitx[,-c(118:157)]

for (i in l){
  x <- (paste0('mitx$',i))
  mitx[i] <- gen.cancer.var(x)
}

tmp <- filter(mitx, cancer == 0)
tmp2 <- filter(mitx, cancer == 1)

rm(mitx)

for (i in 118:157){
  can.name <- (names(tmp2)[i])
  for (j in tmp2){
    tmp2[,i] <- ifelse(tmp2$cancer1 == can.name, 1, 0)
  }
}

rm(can.name, i, j, l, x)

mitx <- rbind(tmp, tmp2)

rm(tmp, tmp2)

save(mitx, file = './mi.tx.v20171023.1.rdata')



# Compute [cancer].any variables ------------------------------------------
load('./mi.tx.v20171023.1.rdata')

mitx$leu.any <- populate.cancer.var(mitx$all, mitx$aml, mitx$leu.other)
mitx$lym.any <- populate.cancer.var(mitx$hl, mitx$nhl, mitx$lym.other)
mitx$cns.any <- populate.cancer.var(mitx$astro, mitx$medullo, mitx$cns.other)
mitx$pns.any <- populate.cancer.var(mitx$neuro, mitx$pns.other)
mitx$renal.any <- populate.cancer.var(mitx$nephro, mitx$renal.other)
mitx$hepatic.any <- populate.cancer.var(mitx$hepato, mitx$hepatic.other)
mitx$bone.any <- populate.cancer.var(mitx$osteo, mitx$ewing, mitx$bone.other)
mitx$rms.any <- populate.cancer.var(mitx$erms, mitx$arms, mitx$rms.other)
mitx$soft.any <- populate.cancer.var(mitx$rms.any, mitx$soft.other)
mitx$gct.any <- populate.cancer.var(mitx$gct.extra, mitx$gct.gonad, mitx$gct.intra)

save(mitx, file = './mi.tx.v.20171023.2.rdata')



# Generate a variable for repeatable random sampling ----------------------
load('./mi.tx.v20171023.2.rdata')

mitx$runif <- runif(8312882, 0, 1)

save(mitx, file = './mi.tx.v20171023.3.rdata')



# Set birth defects variables to zero for children w/o defects ------------
load('./mi.tx.v20171023.3.rdata')

#' Fix an issue where non-BD children in TX have any.birthdefect set to NA.
mitx$any.birthdefect <- ifelse(mitx$state == 'TX' & is.na(mitx$any.birthdefect), 0, mitx$any.birthdefect)

#' Update any.birthdefect variable.  Set to 1 if any of rows 22:112 are 1.
mitx$defect.sum <- rowSums(mitx[,22:112], na.rm = TRUE)
mitx$any.birthdefect <- ifelse(mitx$defect.sum >= 1 & mitx$any.birthdefect == 0, 1, mitx$any.birthdefect)
mitx <- select(mitx, -defect.sum)

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Fill in zero values for children with no defects.
#' 
#' These variables are coded as:
#'  1 = child has that defect
#'  0 = child has no defects
#'  NA = child does not have that defect, but has at least one other
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
for (i in 22:112){
  mitx[, i] <- ifelse(mitx$any.birthdefect == 0, 0, mitx[, i])
}

save(mitx, file = './mi.tx.v20171023.4.rdata')








