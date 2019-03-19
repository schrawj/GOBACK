#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2017.11.09.
#' 
#' Pull together the cancer diagnostic codes for easy reference.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

# Load in state-level data ------------------------------------------------

load("Z:/Jeremy/GOBACK/Datasets/Arkansas/arkansas.v20170913.1.rdata")
ar.trim <- select(ar, studyid, morphology, site.code, behavior3, laterality)
ar.trim <- ar.trim[!is.na(ar.trim$morphology), ]
for (i in 6:21){
  ar.trim[,i] <- as.numeric(NA)
}
ar.trim <- ar.trim[, c(1,2,6:9,3,10:13,4,14:17,5,18:21)]
ar.trim$laterality <- as.numeric(ar.trim$laterality)
ar.trim$behavior3 <- as.numeric(ar.trim$behavior3)

load("Z:/Jeremy/GOBACK/Datasets/North Carolina/nc.cancer.data.rdata")
nc.trim <- select(nc.cancer, ncid, histologic.type.icdo3.1:histologic.type.icdo3.5, primary.site.1:primary.site.5, behavior.code.icdo3.1:behavior.code.icdo3.5)
nc.trim$ncid <- paste0('nc',nc.trim$ncid)
nc.trim <- rename(nc.trim, studyid = ncid)
#' Removes the leading "C" from site codes.
for (i in 7:11){
  nc.trim[,i] <- substr(nc.trim[,i], 2, 4)
}
for (i in 2:16){
  nc.trim[,i] <- ifelse(nc.trim[,i] == "", NA, nc.trim[,i])
  nc.trim[,i] <- as.numeric(nc.trim[,i])
}
#' Add empty variables for laterality, which I don't think we got from NC.
for (i in 17:21){
  nc.trim[,i] <- as.numeric(NA)
}
nc.trim <- nc.trim[!is.na(nc.trim$histologic.type.icdo3.1), ]

load("Z:/Jeremy/GOBACK/Datasets/Michigan/mi.raw.data.v20170818.rdata")
mi.trim <- select(mi, studyid, morph31:morph34,site_code1:site_code4, behave31:behave34, laterality1:laterality4)
mi.trim <- mi.trim[!is.na(mi.trim$morph31), ]
mi.trim$studyid <- paste0('mi',mi.trim$studyid)
mi.trim$morph35 <- as.numeric(NA)
mi.trim$site_code5 <- as.numeric(NA)
mi.trim$behave35 <- as.numeric(NA)
mi.trim$laterality5 <- as.numeric(NA)
mi.trim <- mi.trim[,c(1:5,18,6:9,19,10:13,20,14:17,21)]

tx <- as.data.frame(haven::read_dta('Z:/Birth Defects-Childhood Cancer Projects/GOBACK Project (Tiffany)/Texas/TX_Cancer_added can variables.dta'))
tx.can <- as.data.frame(haven::read_sas('Z:/Jeremy/GOBACK/Datasets/Texas/tcr_17_007_july20.sas7bdat'))
tx$birthID <- as.character(tx$birthID)
tx.can$birthID <- as.character(tx.can$birthID)
tx.trim <- left_join(select(tx, birthID, morph31, morph32, morph33, morph34,site_code1, site_code2, site_code3, site_code4),
                     select(tx.can, birthID, BEHAVIORICDO3, LATERAL),
                     by = 'birthID')
tx.trim$birthID <- as.character(paste0('tx',tx.trim$birthID))
tx.trim$LATERAL <- as.numeric(tx.trim$LATERAL)
tx.trim <- rename(tx.trim, studyid = birthID, behavior1 = BEHAVIORICDO3, laterality1 = LATERAL)
tx.trim <- tx.trim[!duplicated(tx.trim$studyid), ]
for (i in 2:10){
  tx.trim[,i] <- ifelse(tx.trim[,i] == "", NA, tx.trim[,i])
  tx.trim[,i] <- as.numeric(tx.trim[,i])
}
tx.trim$morph35 <- as.numeric(NA)
tx.trim$site_code5 <- as.numeric(NA)
tx.trim$behavior2 <- as.numeric(NA)
tx.trim$behavior3 <- as.numeric(NA)
tx.trim$behavior4 <- as.numeric(NA)
tx.trim$behavior5 <- as.numeric(NA)
tx.trim$laterality2 <- as.numeric(NA)
tx.trim$laterality3 <- as.numeric(NA)
tx.trim$laterality4 <- as.numeric(NA)
tx.trim$laterality5 <- as.numeric(NA)
tx.trim <- tx.trim[,c(1:5,12,6:9,13,10,14:17,11,18:21)]

names <- c(colnames(tx.trim))

colnames(ar.trim) <- names
colnames(mi.trim) <- names
colnames(nc.trim) <- names

cancer.codes <- rbind(ar.trim, mi.trim, nc.trim, tx.trim)

rm(ar, mi, tx, tx.can, nc.cancer, ar.trim, mi.trim, nc.trim, tx.trim, i, names); gc()

save(cancer.codes, file = 'Z:/Jeremy/GOBACK/Datasets/cancer.codes.v20180227.1.rdata')

