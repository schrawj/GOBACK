#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.09.05.
#' 
#' Sharon asked for a list of all birth defects codes in kids with 
#' choanal atresia and acute leukemia. That's easy.
#' 
#' Philip asked for all available IDs for each of these kids. That's more 
#' involved. Dani says she can't use their StudyID variable to get their 
#' info.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

# Find kids, get list of defects. -----------------------------------------

load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20180829.rdata")
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.v20180712.rdata")
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.v20180606.rdata")
load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/linked.registry.ids.v20180908.rdata')

require(dplyr); require(xlsx); require(stringr)

tmp <- filter(goback, choanal.atresia == 1 & leu.any == 1)
ids <- c(tmp$studyid)
 
tmp <- filter(bd.codes.txnc, studyid %in% ids)
tmp <- tmp[, 1:19]
for (i in 20:22){tmp[, i] <- as.character()}
names(tmp) <- c('studyid',rep(paste0('code',1:21)))

tmp2 <- filter(bd.codes.mi, studyid %in% ids)
tmp2 <- tmp2[, 1:22]
names(tmp2) <- c('studyid',rep(paste0('code',1:21)))

choanal <- rbind(tmp, tmp2)

map <- read.xlsx(file = 'W:/Old_genepi2/Jeremy/GOBACK/Data dictionaries and data definitions/dxmap_current.xlsx',
                 sheetName = 'Current', colIndex = c(1,2,13,14), stringsAsFactors = FALSE)
names(map) <- str_replace_all(tolower(names(map)), '_', '.')

#' For the love of God, save this file. I will use it again.
save(map, file = 'W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.bpa.to.icd.mappings.rdata')

out <- data.frame(studyid = choanal[,1])

for (i in 2:22){
  
  tmp <- select(choanal, 1, i)
  names(tmp)[2] <- 'code'
  tmp <- left_join(tmp, 
                   select(map, bpa.number, bpaname), 
                   by = c('code' = 'bpa.number'))
  names(tmp)[2:3] <- paste0(names(tmp)[2:3],i-1)
  
  out <- left_join(out, tmp, 'studyid')
  
}

write.xlsx(out, file = 'W:/Old_genepi2/Jeremy/GOBACK/R outputs/bd.codes.for.choanal.atresia.leukemia.cases.xlsx',
           row.names = FALSE, showNA = FALSE, sheetName = 'BDCodes')

choanal.ids <- filter(goback.ids, studyid %in% choanal$studyid)

write.csv(choanal.ids, file = 'W:/Old_genepi2/Jeremy/GOBACK/R outputs/ids.for.choanal.atresia.leukemia.cases.csv',
          row.names = FALSE)



# Retrieve their birth defects and cancer registry IDs --------------------

load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20180829.rdata")
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.v20180712.rdata")
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.v20180606.rdata")

require(dplyr); require(readstata13)

goback.ids <- select(goback, studyid, state); rm(goback)

#' This file has TX BD registry IDs: CASE_ID.
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/Texas/tx.raw.data.rdata")

tx.raw$birthID <- paste0('tx',tx.raw$birthID)
tx.raw$CASE_ID <- ifelse(tx.raw$CASE_ID == "", NA, tx.raw$CASE_ID)

goback.ids <- left_join(goback.ids,
                        select(tx.raw, birthID, CASE_ID),
                        by = c('studyid' = 'birthID'))

goback.ids <- rename(goback.ids, bd.registry.id = CASE_ID); rm(tx.raw)

#' This file has birthID and the TCR ID variables, patientid.
tx.can <- read.dta13("W:/Old_genepi2/Birth Defects-Childhood Cancer Projects/GOBACK Project (Tiffany)/Texas/tcr_17_007_july20 (STATA).dta")
tx.can$studyid <- paste0('tx',tx.can$birthid)

goback.ids <- left_join(goback.ids, 
                        select(tx.can, studyid, patientid),
                        by = 'studyid')
goback.ids <- rename(goback.ids, cancer.registry.id = patientid); rm(tx.can)

save(goback.ids, file = 'W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/linked.registry.ids.v20180908.rdata')
