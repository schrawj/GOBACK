#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2019.08.14.
#' 
#' Pull BD codes and names for children wtih RMS and any birth defect.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(dplyr); require(xlsx)

load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20190606.1.rdata")
load('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.v20180606.rdata')
load('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.v20180712.rdata')
load('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.bpa.to.icd.mappings.rdata')

ids <- c(filter(goback, rms.any == 1 & any.birthdefect == 1)$studyid)

#' Define a function that finds the first column in which no children have a BD recorded.
#' This is one more than the number of columns we need to select.
find.column <- function(df){
  
  for (i in 1:ncol(df)){
    
    sum <- sum(is.na(df[,i]))
    
    if (sum == nrow(df)){ return(i) }
    
    else { next }
    
  }
  
}

#' Define a function that finds the longer data frame, adds an appropriate number of columns to the shorter one, 
#' harmonizes the names, rbinds them, and spits if to the global environment with the name 'codes.'
consolidate.dfs <- function(df1, df2){
  
  condition <- ncol(tmp) > ncol(tmp2)
  
  if (condition == T){
    
    for (i in (ncol(df2)+1):ncol(df1)){
      
      df2[,i] <- as.character()
      
    }
    
  }
  
  else {
    
    for (i in (ncol(df1)+1):ncol(df2)){
      
      df1[,i] <- as.character()
      
    }
    
  }
  
  names(df2) <- c('id',paste0(rep('code',ncol(df2)-1), 1:(ncol(df2)-1)))  
  names(df1) <- c('id',paste0(rep('code',ncol(df2)-1), 1:(ncol(df2)-1)))  
  
  codes <<- rbind(df1,df2)
}

tmp <- filter(bd.codes.txnc, studyid %in% ids)
c <- find.column(tmp)
tmp <- tmp[, 1:(c-1)]

tmp2 <- filter(bd.codes.mi, studyid %in% ids)
c <- find.column(tmp2)
tmp2 <- tmp2[ , 1:(c-1)]

consolidate.dfs(tmp, tmp2)

out <- data.frame(studyid = codes[,1])
out <- left_join(out, select(goback, studyid, cancer1), by = 'studyid')

#' Append defect names for children diagnosed using CDC/BPA codes.
for (i in 2:ncol(codes)){
  
  tmp <- select(codes, 1, i)
  names(tmp)[2] <- 'code'
  tmp <- left_join(tmp, 
                   select(map, bpa.number, bpaname), 
                   by = c('code' = 'bpa.number'))
  names(tmp)[2:3] <- paste0(names(tmp)[2:3],i-1)
  
  out <- left_join(out, tmp, c('studyid' = 'id'))
  
}

write.xlsx(out, 
          file = '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/R outputs/birth.defect.codes.for.rms.cases.v20190914.xlsx',
          row.names = F, showNA = F)
