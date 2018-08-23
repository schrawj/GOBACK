#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.08.23.
#' 
#' Contains the code used to generate the R Markdown documents describing 
#' QC measures in response to co-author comments, plus some additional 
#' code to output data in other formats.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

# Review BD codes for cancer cases w/22q or trisomy 13 --------------------

#' Sonja Rasmussen asked to see birth defects codes for kids with 22q11.2 deletions and cancer.
#' We also need to go back in time to retrieve their birth defects registry ID for possible
#' review by Angela Scheurle, the TBDR medical geneticist.
require(dplyr); require(kableExtra)

load('Z:/Jeremy/GOBACK/Datasets/goback.chrom.v20180711.rdata')
load('Z:/Jeremy/GOBACK/Datasets/Old Datasets/Texas/tx.raw.data.rdata')
load('Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.v20180606.rdata')
load('Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.v20180712.rdata')

check <- select(rbind(filter(goback.chrom, di.george.syndrome == 1 & cancer == 1 & state != 'AR'),
                      filter(goback.chrom, trisomy13 == 1 & cancer == 1 & state != 'AR')),
                studyid, cancer1, di.george.syndrome, trisomy13)

check.bpa <- select(filter(bd.codes.txnc, studyid %in% check$studyid), 1:32)
check.icd <- filter(bd.codes.mi, studyid %in% check$studyid)

check.icd[26:32] <- as.character(NA)

names <- c('studyid',rep(paste0('code',1:31)))
names(check.bpa) <- names 
names(check.icd) <- names

check <- arrange(left_join(check, rbind(check.bpa, check.icd), by = 'studyid'), studyid)

kable(check) %>% kable_styling() %>% scroll_box(width = '100%', height = '400px')

#' Export TX kids only to a csv file, and include the CASE_ID variable from the raw data.
tx.raw <- select(tx.raw, birthID, CASE_ID)
tx.raw$studyid <- paste0('tx',tx.raw$birthID)

check <- select(filter(goback.chrom, di.george.syndrome == 1 & cancer == 1 & state == 'TX'),
                studyid, cancer1, di.george.syndrome)
check <- left_join(check, select(tx.raw, CASE_ID, studyid), by = 'studyid')
check.bpa <- filter(bd.codes.txnc, studyid %in% check$studyid)
check <- left_join(select(check, studyid, cancer1, CASE_ID), check.bpa, by = 'studyid')

write.csv(check, file = 'Z:/Jeremy/GOBACK/R outputs/bd.codes.for.tx.22q.cancer.cases.csv', row.names = FALSE)

rm(list = ls()); gc()



# Review specific BD codes for kids w/obstructive GU defects --------------

require(dplyr); require(xlsx)

load('Z:/Jeremy/GOBACK/Datasets/goback.v20180711.rdata')
load('Z:/Jeremy/GOBACK/Datasets/Old Datasets/Texas/tx.raw.data.rdata')
load('Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.transpose.v20180614.rdata')
load('Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.v20180606.rdata')

ob <- filter(goback, obstructive.genitourinary.defects == 1 & state %in% c('TX','NC'))
ob.bpa <- filter(bd.codes.txnc.transpose, studyid %in% ob$studyid)

cols <- c(1, which(grepl('753.2', names(ob.bpa))))

ob.bpa <- select(ob.bpa, cols)
ob.bpa$`753.200` <- ifelse(rowSums(ob.bpa[2:7], na.rm = TRUE) >= 1, 1, 0)
ob.bpa$`753.210` <- ifelse(rowSums(ob.bpa[8:12], na.rm = TRUE) >= 1, 1, 0)
ob.bpa$`753.220` <- ifelse(rowSums(ob.bpa[13:18], na.rm = TRUE) >= 1, 1, 0)
ob.bpa$`753.290` <- ifelse(rowSums(ob.bpa[19:24], na.rm = TRUE) >= 1, 1, 0)
ob.bpa <- ob.bpa[, c(1,2,8,13,20)]

gu.defect.counts <- data.frame(code = names(ob.bpa[2:5]),
                               count = colSums(ob.bpa[2:5]))

gu.names <- read.xlsx('Z:/Jeremy/GOBACK/Data dictionaries and data definitions/dxmap_current.xlsx', sheetName = 'Current')

gu.defect.counts <- left_join(gu.defect.counts, 
                              select(gu.names, BPA_number, BPANAME), 
                              by = c('code' = 'BPA_number'))

kable(gu.defect.counts) %>% kable_styling() %>% scroll_box(width = '100%')

ob <- filter(goback, obstructive.genitourinary.defects == 1 & cancer == 1 & state == 'TX')

tab <- table(ob$cancer1)

gu.cancers <- data.frame(cancer = names(tab), count = tab)
gu.cancers <- rename(select(gu.cancers, -count.Var1), cases = count.Freq)

ob <- left_join(select(ob, studyid, cancer1), bd.codes.txnc, by = 'studyid')

write.csv(ob, file = 'Z:/Jeremy/GOBACK/R outputs/bd.codes.for.tx.obstructive.gu.defects.cancer.cases.csv')

rm(list = ls()); gc()



# Pick out some MBD kids at random ----------------------------------------

require(dplyr); require(kableExtra)

load('Z:/Jeremy/GOBACK/Datasets/goback.nochrom.v20180711.rdata')
load('Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.v20180606.rdata')
load('Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.v20180712.rdata')

tmp <- arrange(filter(goback.nochrom, defect.total >= 4 & state %in% c('TX','NC')), runif)
tmp <- tmp[30000:30099,]
mbd.codes <- filter(bd.codes.txnc, studyid %in% tmp$studyid)
kable(mbd.codes) %>% kable_styling() %>% scroll_box(width = '100%', height = '400px')

tmp <- arrange(filter(goback.nochrom, defect.total >= 4 & state == 'MI', runif))
tmp <- tmp[10900:10999,]
mbd.codes <- filter(bd.codes.mi, studyid %in% tmp$studyid)
kable(mbd.codes) %>% kable_styling() %>% scroll_box(width = '100%', height = '400px')
