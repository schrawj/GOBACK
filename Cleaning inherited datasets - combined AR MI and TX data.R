#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#'                                       
#'                      ARKANSAS + MICHIGAN + TEXAS
#'                                         
#' Having previously resolved issues that were specific to MI or TX, then
#' issues which affected both MI and TX in the same way, I am ready to 
#' join AR, MI and TX data together into a single data frame.
#' 
#' These data will then be used for preliminary analysis for Philip's 
#' 11/1/2017 deadline.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------



# Prep environment --------------------------------------------------------
require(dplyr)
require(stringr)

#' For desktop:
setwd('Z:/Jeremy/GOBACK/Datasets/Combined Arkansas Michigan and Texas/')

#' For laptop:
setwd("//discovery2/bcm-pedi-genepi2/Jeremy/GOBACK/Datasets/Combined Arkansas Michigan and Texas/")



# User-defined functions --------------------------------------------------
get.affected.ids <- function(codes, data.frame,first.defect.col.index, other.defect.col.indices, id.var.index){
  
  tmp2 <- data.frame[first.defect.col.index %in% codes, ]
  
  for (i in other.defect.col.indices){
    tmp3 <- filter(data.frame, i %in% codes)
    tmp2 <- rbind(tmp2, tmp3) 
  }
  
  tmp2 <- tmp2[!duplicated(tmp2[,id.var.index]), ]
  
  affected.ids <- c(paste0('tx', paste0(tmp2[,id.var.index])))
  
  rm(tmp2, tmp3)
  
  return(affected.ids)
}



# Load in and merge state-level data --------------------------------------
load("Z:/Jeremy/GOBACK/Datasets/Combined Michigan and Texas/mi.tx.v20171023.4.rdata")
load("Z:/Jeremy/GOBACK/Datasets/Arkansas/arkansas.v20170913.2.rdata")

#' Add runif variable to AR data.
ar$runif <- runif(629119, 0, 1)

armitx <- rbind(ar, mitx)

flag <- identical(str(ar),str(armitx))
flag2 <- identical(str(mitx), str(armitx))

rm(ar, mitx, flag, flag2)
gc()

save(armitx, file = './ar.mi.tx.v20171023.1.rdata')



# Clean demographics variables --------------------------------------------
load('./ar.mi.tx.v20171023.1.rdata')

#' Rules for filtering extreme maternal and gestational ages were not 
#' applied in AR data because they were agreed upon later.
armitx$m.age <- ifelse(armitx$m.age < 13 | armitx$m.age > 50, NA, armitx$m.age)
armitx$gest.age <- ifelse(armitx$gest.age < 22 | armitx$gest.age > 44, NA, armitx$gest.age)

save(armitx, file = './ar.mi.tx.v20171023.2.rdata')



# Person-years ------------------------------------------------------------
load('./ar.mi.tx.v20171023.2.rdata')

#' Decision was made 10/23/17 to exclude children with negative person-years
#' values from time-to-event analysis completely, as there are so few.
armitx$person.yrs <- ifelse(armitx$person.yrs < 0, NA, armitx$person.yrs)

save(armitx, file = './ar.mi.tx.v20171023.3.rdata')



# Clean overall.chrom -----------------------------------------------------
load('./ar.mi.tx.v20171023.3.rdata')

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' The overall.chrom variable in the TX data is a mess.
#' 
#' It is an empty string for kids with no birth defect.
#' 
#' It is 'Has chromosomal abnormality' for kids WITHOUT a chromosomal 
#' abnormality and it is 'No chromosomal abnormality' for kids WITH 
#' a chromosomal abnormality.
#' 
#' Should be 0 if no chromosomal defect, 1 otherwise.
#' 
#' Appears to be parameterized correctly in MI and AR.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

armitx$overall.chrom <- as.numeric(
                                    ifelse(armitx$overall.chrom == 1, 1,
                                             ifelse(armitx$overall.chrom == 0, 0,
                                                    ifelse(armitx$overall.chrom == 'Has Chromosomal Abnormality', 0, 
                                                            ifelse(armitx$overall.chrom == 'No Chromosomal Abnormality', 1, 
                                                                   ifelse(armitx$overall.chrom == '', 0, 
                                                                          ifelse(is.na(armitx$overall.chrom), 0, 
                                                                                 ifelse(armitx$overall.chrom == 1, 1,
                                                                                        ifelse(armitx$overall.chrom == 0, 0, 9999)))))))))
armitx$overall.chrom <- ifelse(is.na(armitx$overall.chrom), 0, armitx$overall.chrom)

save(armitx, file = './ar.mi.tx.v20171023.4.rdata')



# Fix errors in select TX birth defects -----------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' I don't know why, but there are no cases of PDA, coarctation, pulmonary
#' artery anomalies or chromosomal anomalies other minor in the TX data.
#' 
#' This is clearly an error as there ARE codes for these defects in TX 
#' kids.
#' 
#' Search for kids with these codes and update these variables.
#' 
#' This issue also affects NC, but (oddly) not MI.  It is addressed in the 
#' NC data cleaning script.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
load("Z:/Jeremy/GOBACK/Datasets/Texas/tx.raw.data.bpa.codes.v20171106.rdata")
load('./ar.mi.tx.v20171023.4.rdata')

#' Load in a file that has the birthID variable so we can link to main TX data frame.
tx.bd.raw <- read.csv(file = 'C:/Users/schraw/Desktop/TX BPA codes.csv', header = TRUE, stringsAsFactors = FALSE)
tx.bd.raw <- select(tx.bd.raw, X, CASE_ID, birthID)

tx.bd <- left_join(tx.bd, tx.bd.raw, by = c('X','CASE_ID'))
tx.bd <- tx.bd[!is.na(tx.bd$birthID), ]

rm(tx.bd.raw)

coarc.codes <- c(747.100, 747.110, 747.190)
pulm.codes <- c(747.300, 747.310, 747.320, 747.325, 747.330, 747.340, 747.380, 747.390)
pda.codes <- c(747.000, 747.008, 747.010, 747.060)

tmp <- tx.bd[!is.na(tx.bd$bpa1),]

#' Coarctation.
get.affected.ids <- function(codes, data.frame,first.defect.col.index, other.defect.col.indices, id.var.index){
  
  tmp2 <- data.frame[first.defect.col.index %in% codes, ]

  for (i in other.defect.col.indices){
    tmp3 <- filter(data.frame, i %in% codes)
    tmp2 <- rbind(tmp2, tmp3) 
  }
  
  tmp2 <- tmp2[!duplicated(tmp2[,id.var.index]), ]
  
  affected.ids <- c(paste0('tx', paste0(tmp2[,id.var.index])))
  
  rm(tmp2, tmp3)
  
  return(affected.ids)
}

coarc.ids <- get.affected.ids(coarc.codes, tmp, tmp[,3], tmp[,4:68], 69)
pulm.ids <- get.affected.ids(pulm.codes, tmp, tmp[,3], tmp[,4:68], 69)
pda.ids <- get.affected.ids(pda.codes, tmp, tmp[,3], tmp[,4:68], 69)

#' Update the main data set.
armitx$coarctationofaorta <- ifelse(armitx$state == 'TX' & armitx$studyid %in% coarc.ids, 1, armitx$coarctationofaorta)
armitx$patentductusarteriosis <- ifelse(armitx$state == 'TX' & armitx$studyid %in% pda.ids, 1, armitx$patentductusarteriosis)
armitx$pulmonaryarteryanomalies <- ifelse(armitx$state == 'TX' & armitx$studyid %in% pulm.ids, 1, armitx$pulmonaryarteryanomalies)

save(armitx, file = './ar.mi.tx.v20171106.1.rdata')

rm(tx.bd, tx.raw, coarc.codes, pulm.codes, pda.codes, coarc.ids, pda.ids, pulm.ids)



# Split dataset by chromosomal defect status ------------------------------
load('./ar.mi.tx.v20171106.1.rdata')

chrom <- armitx[armitx$any.birthdefect == 1 & armitx$overall.chrom == 1, ]
non.chrom <- armitx[armitx$any.birthdefect == 1 & armitx$overall.chrom == 0, ]
no.defect <- filter(armitx, any.birthdefect == 0)

rm(armitx)

armitx.chrom <- rbind(chrom, no.defect)

rm(chrom)
gc()

save(armitx.chrom, file = './ar.mi.tx.chromosomaldefects.v20171106.1.rdata')

rm(armitx.chrom)

armitx.nochrom <- rbind(non.chrom, no.defect)

rm(non.chrom, no.defect)

save(armitx.nochrom, file = './ar.mi.tx.nochromosomaldefects.v20171106.1.rdata')

rm(armitx.nochrom)
gc()



# Re-order cancer variables alphabetically --------------------------------
load('./ar.mi.tx.chromosomaldefects.v20171106.1.rdata')
load('./ar.mi.tx.nochromosomaldefects.v20171106.1.rdata')

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Once I started running logistic regression models for specific 
#' birth defect-cancer associations I realized I was wasting a lot of time
#' by computing regression models for situations where there were no 
#' comorbid cases.
#' 
#' To avoid this, I am going to edit the loop that generates these models
#' to only build models for associations where there are >= 5 comorbid 
#' cases.
#' 
#' The easiest way of doing this seems to be to generate a crosstab of 
#' cancer-birth defect counts, select the column indices that enough
#' observations, and run regression models for cancers in that list.
#' 
#' The crosstab will display cancers alphabetically rather than by their
#' actual column index in the data frame, hence I need to reorder cancers
#' alphabetically in the data frame to create an easy 1:1 column index
#' mapping.
#' 
#' Reorder the cancers that show up in cancer1 alphabetically.  
#' Move the .any variables to the end of the section.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

armitx.nochrom <- armitx.nochrom[,c(1:119,
                                    145, 126, 143, 129, 156, 154, 144, 141, 151, 152, 150, 139, 137, 122, 121, 125, 127, 134, 130, 123, 140,
                                    155, 157, 132, 136, 133, 147, 149, 
                                    120, 124, 128, 131, 135, 138, 142, 146, 148, 153, 
                                    158)]

save(armitx.nochrom, file = './ar.mi.tx.nochromosomaldefects.v20171025.1.rdata')

armitx.chrom <- armitx.chrom[,  c(1:119,
                                  145, 126, 143, 129, 156, 154, 144, 141, 151, 152, 150, 139, 137, 122, 121, 125, 127, 134, 130, 123, 140,
                                  155, 157, 132, 136, 133, 147, 149, 
                                  120, 124, 128, 131, 135, 138, 142, 146, 148, 153, 
                                  158)]

save(armitx.chrom, file = './ar.mi.tx.chromosomaldefects.v20171025.1.rdata')

