#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2019.01.31.
#' 
#' Sharon, Philip, and I have picked out several additional specific BD-CC
#' associations we want to follow up in the family-based cohort.
#' 
#' Need to find IDs and print list of comorbid defects for children with 
#' these conditions.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(xlsx); require(dplyr)

#' W:/ points to //smb-main.ad.bcm.edu/
load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20180829.rdata')
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/linked.registry.ids.v20180908.rdata")
load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.v20180712.rdata')
load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.v20180606.rdata')
load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.bpa.to.icd.mappings.rdata')

goback.ids <- goback.ids[!duplicated(goback.ids$studyid), ]

#' Update MI dimensions to match TX and NC.
for (i in 25:67){
  bd.codes.mi[ ,i] <- as.numeric()
}

for (i in 2:67){
  colnames(bd.codes.mi)[i] <- paste0('bd.code',as.character(i))
  colnames(bd.codes.txnc)[i] <- paste0('bd.code',as.character(i))
}

bds <- c('choanal.atresia','spinabifida.wo.anencephaly','hydrocephalus.wo.spinabifida','pyloric.stenosis','lvot.defects',rep('pulmvalveatresiaandstenosis',2),'ventricularseptaldefect',rep('craniosynostosis',4))
cancers <- c('leu.any','soft.other','nephro','medullo','neuro','neuro','hepato',rep('hepato',2),'medullo','nephro','neuro')

#' Generates a list. Each element is a data frame with the IDs and BD codes for children with the index BD-CC event.
l <- list()

for (i in 1:length(bds)){
  
  bd.col <- which(colnames(goback) == bds[i])
  ca.col <- which(colnames(goback) == cancers[i])
  
  tag <- paste0(bds[i],'-',cancers[i])
  
  new.comorbid <- select(subset(goback, goback[,ca.col] == 1 & goback[,bd.col] == 1), studyid)

  new.comorbid <- left_join(new.comorbid, 
                            select(goback.ids, studyid, bd.registry.id, cancer.registry.id),
                            by = 'studyid')
  
  new.comorbid <- rbind(left_join(subset(new.comorbid, substr(new.comorbid$studyid,1,2) %in% c('tx','nc')),
                               bd.codes.txnc,
                               by = 'studyid'),
                         left_join(subset(new.comorbid, substr(new.comorbid$studyid,1,2) == 'mi'),
                               bd.codes.mi,
                               by = 'studyid'))

  l[[tag]] <- new.comorbid

}

rm(goback, bd.codes.mi, bd.codes.txnc, bd.col, ca.col, i, goback.ids); gc()

#' Removes birth defects columns with pure NA values.
for (i in 1:length(l)){
  
  last.col <-c()
  missing.count <- as.numeric()
  
  for (j in 4:ncol(l[[i]])){
    
    missing.count <- as.numeric(sum(is.na(l[[i]][,j])))
    
    if (missing.count == nrow(l[[i]]) & length(last.col) == 0) {
      
      last.col <- j-1
      
    }
    
    else if (missing.count < nrow(l[[i]]) & length(last.col) == 0) {
      
      next
      
    }
    
    else {
      
      break
      
    }
    
  }
  
  l[[i]] <- l[[i]][,1:last.col]
  
}

rm(last.col, missing.count, i, j)

#' Add text descriptions for BD codes from TX DSHS. This step is optional and comparatively slow.
map.bpa <- map[,1:2]
map.icd <- map[,3:4]

for (i in 1:length(l)){
  
  out <- data.frame(l[[i]])
  
  out <- data.frame(studyid = out[,1])
  
  for (j in 4:ncol(l[[i]])){
    
    new.col.names <- c(paste0('code',as.character(j-3)), paste0('name',as.character(j-3)))
    
    new.out.txnc <- select(filter(l[[i]], substr(l[[i]]$studyid,1,2) %in% c('tx','nc')), 1, j)
    names(new.out.txnc)[2] <- 'code'
    new.out.txnc <- left_join(new.out.txnc, map.bpa, by = c('code' = 'bpa.number'))
    names(new.out.txnc)[2:3] <- new.col.names
    
    new.out.mi <- select(filter(l[[i]], substr(l[[i]]$studyid,1,2) == 'mi'), 1, j)
    names(new.out.mi)[2] <- 'code'
    new.out.mi <- left_join(new.out.mi, map.icd, by = c('code' = 'icd9.code'))
    new.out.mi <- new.out.mi[!duplicated(new.out.mi$studyid), ]
    names(new.out.mi)[2:3] <- new.col.names
    
    new.out <- rbind(new.out.txnc, new.out.mi)
    
    out <<- left_join(out, new.out, 'studyid')
    
    rm(new.out.mi, new.out.txnc, new.col.names); gc()
    
    print(paste('done with column',j,'of data frame', i))
    
  }  
  
  l[[i]] <- out
  
}

#' Write to Excel.
for (i in 1:length(l)){
  
  if (nrow(l[[i]]) != 0){
    
    write.xlsx(l[[i]], file = 'W:/Old_genepi2/Jeremy/GOBACK/Family-based cohort/children.w.target.bd.cc.associations.v20190206.xlsx',
               sheetName = names(l[i]), 
               row.names = FALSE, 
               showNA = FALSE,
               append = ifelse(i == 1, FALSE, TRUE))
    
  }
  
  else {
    
    next
  }
  
}

rm(list = ls()); gc()
