#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.06.13.
#' 
#' It occurred to me to generate an alternative form of the data frame 
#' holding birth defects codes: one with a row for every subject and a 
#' column for every defect.  This would be better suited to looking up all
#' children who have a given birth defect, whereas the existing data frame
#' is well suited to looking up all codes in a given child.
#' 
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------


# Texas and North Carolina ------------------------------------------------

require(stringr); require(dplyr); require(tictoc)

load("//discovery2/bcm-pedi-genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.v20180606.rdata")

#' Generate a list of all the unique codes.
codes <- as.character()

for (i in 1:nrow(bd.codes.txnc)){
  tmp <- as.character(bd.codes.txnc[i,2:67])
  codes <- c(codes, subset(tmp, !is.na(tmp)))
  codes <- subset(codes, !duplicated(codes))
}

codes <- sort(codes)

#' Initialize an empty data frame with the proper dimensions and names.
bd.codes.txnc.transpose <- as.data.frame(matrix(nrow = nrow(bd.codes.txnc), ncol = 1 + length(codes)))
names(bd.codes.txnc.transpose) <- c('studyid', codes)
bd.codes.txnc.transpose[, 1] <- bd.codes.txnc$studyid

#' Ensure both data frames have rows sorted in the same order.
bd.codes.txnc <- arrange(bd.codes.txnc, studyid)
bd.codes.txnc.transpose <- arrange(bd.codes.txnc.transpose, studyid)

#' Initialize an empty list.
#' Extract non-missing, unique codes from each row of the data frame.  
#' Map them to the column index they correspond to in the transposed data frame.
#' Save these indices as elements of a list.
l <- list()

tic()
for (i in 1:nrow(bd.codes.txnc)){
  
  tmp <- as.character(bd.codes.txnc[i,2:67])
  tmp <- subset(tmp, !duplicated(tmp))
  tmp <- subset(tmp, !is.na(tmp))
  
  for (j in 1:length(tmp)){
    tmp[j] <- which(codes == tmp[j])
  }
  
  tmp <- as.numeric(tmp) + 1
  
  l[[i]] <- tmp
  
}
toc()

#' For every row in the transposed data frame, set any column representing one of that child's diagnoses to 1.
tic()
for (i in 1:nrow(bd.codes.txnc.transpose)){
  
  for (j in l[[i]]){
    bd.codes.txnc.transpose[i, j] <- 1
  }
}
toc()

save(bd.codes.txnc.transpose, file = 'Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.transpose.v20180614.rdata')

rm(list = ls); gc()



# Michigan ----------------------------------------------------------------

require(stringr); require(dplyr); require(tictoc)

load("//discovery2/bcm-pedi-genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.v20180606.rdata")

#' Generate a list of all the unique codes.
codes <- as.character()

for (i in 1:nrow(bd.codes.mi)){
  tmp <- as.character(bd.codes.mi[i,2:25])
  codes <- c(codes, subset(tmp, !is.na(tmp)))
  codes <- subset(codes, !duplicated(codes))
}

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
  tmp <- subset(tmp, !is.na(tmp))
  
  for (j in 1:length(tmp)){
    tmp[j] <- which(codes == tmp[j])
  }
  
  tmp <- as.numeric(tmp) + 1
  
  l[[i]] <- tmp
  
}
toc()




