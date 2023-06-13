#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Authored: 2018.06.13.
#' 
#' Last updated: 2019.12.05.
#' 
#' It occurred to me to generate an alternative form of the data frame 
#' holding birth defects codes: one with a row for every subject and a 
#' column for every defect.  This would be better suited to looking up all
#' children who have a given birth defect, whereas the existing data frame
#' is well suited to looking up all codes in a given child.
#' 
#' Updated 2019.12.05 to combine the OK data with TX and NC.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

# Texas, North Carolina, Oklahoma -----------------------------------------

require(stringr); require(dplyr); require(tictoc)

#load("//discovery2/bcm-pedi-genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.v20180606.rdata")
load('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txncok.v20191205.rdata')

#' Generate a list of all the unique codes.
codes <- as.character()

#' Slow, as you can imagine.
for (i in 1:nrow(bd.codes.ok.nc.tx)){
  
  tmp <- as.character(bd.codes.ok.nc.tx[i,2:ncol(bd.codes.ok.nc.tx)])
  codes <- c(codes, subset(tmp, !is.na(tmp) & tmp != 'NA'))
  codes <- subset(codes, !duplicated(codes))
  
}

codes <- sort(codes)

#' Initialize an empty data frame with the proper dimensions and names.
bd.codes.ok.nc.tx.transpose <- as.data.frame(matrix(nrow = nrow(bd.codes.ok.nc.tx), ncol = 1 + length(codes)))

names(bd.codes.ok.nc.tx.transpose) <- c('studyid', codes)

bd.codes.ok.nc.tx.transpose[, 1] <- bd.codes.ok.nc.tx$studyid

#' Ensure both data frames have rows sorted in the same order.
bd.codes.ok.nc.tx <- arrange(bd.codes.ok.nc.tx, studyid)

bd.codes.ok.nc.tx.transpose <- arrange(bd.codes.ok.nc.tx.transpose, studyid)

#' Initialize an empty list.
#' Extract non-missing, unique codes from each row of the data frame.  
#' Map them to the column index they correspond to in the transposed data frame.
#' Save these indices as elements of a list.
l <- list()

#' Ran in ~3.5 minutes for TX, NC, and OK.
tic()
for (i in 1:nrow(bd.codes.ok.nc.tx)){
  
  tmp <- as.character(bd.codes.ok.nc.tx[i,2:ncol(bd.codes.ok.nc.tx)])
  tmp <- subset(tmp, !duplicated(tmp) & tmp != 'NA')

  for (j in 1:length(tmp)){
   
     tmp[j] <- which(codes == tmp[j])
     
  }
  
  tmp <- as.numeric(tmp) + 1
  
  l[[i]] <- tmp
  
}
toc()

#' For every row in the transposed data frame, set all columns representing that child's diagnoses to 1.
#' Slow. Ran in ~10.5 minutes for Tx, NC, and OK.
tic()
for (i in 1:nrow(bd.codes.ok.nc.tx.transpose)){
  
  for (j in l[[i]]){
    
    bd.codes.ok.nc.tx.transpose[i, j] <- 1
    
  }
  
}
toc()

save(bd.codes.ok.nc.tx.transpose, 
     file = '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txncok.transpose.v20191205.rdata')

rm(list = ls()); gc()



# Michigan ----------------------------------------------------------------

require(stringr); require(dplyr); require(tictoc)

load("//discovery2/bcm-pedi-genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.v20180712.rdata")

#' Generate a list of all the unique codes.
codes <- as.character()

for (i in 1:nrow(bd.codes.mi)){
  tmp <- as.character(bd.codes.mi[i,2:25])
  codes <- c(codes, subset(tmp, !is.na(tmp)))
  codes <- subset(codes, !duplicated(codes))
}

codes <- subset(codes, codes != 'NA')
codes <- sort(codes)

#' Initialize an empty data frame with the proper dimensions and names.
bd.codes.mi.transpose <- as.data.frame(matrix(nrow = nrow(bd.codes.mi), ncol = 1 + length(codes)))
names(bd.codes.mi.transpose) <- c('studyid', codes)
bd.codes.mi.transpose[, 1] <- bd.codes.mi$studyid

#' Ensure both data frames have rows sorted in the same order.
bd.codes.mi <- arrange(bd.codes.mi, studyid)
bd.codes.mi.transpose <- arrange(bd.codes.mi.transpose, studyid)

l <- list()

tic()
for (i in 1:nrow(bd.codes.mi)){
  
  tmp <- as.character(bd.codes.mi[i,2:25])
  tmp <- subset(tmp, !duplicated(tmp))
  tmp <- subset(tmp, tmp != 'NA')
  
  for (j in 1:length(tmp)){
    tmp[j] <- which(codes == tmp[j])
  }
  
  tmp <- as.numeric(tmp) + 1
  
  l[[i]] <- tmp
  
}
toc()

tic()
for (i in 1:nrow(bd.codes.mi.transpose)){
  
  for (j in l[[i]]){
    bd.codes.mi.transpose[i, j] <- 1
  }
}

toc()

save(bd.codes.mi.transpose, file = 'Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.transpose.v20180712.rdata')

rm(list = ls()); gc()