
require(dplyr)

load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txncok.transpose.v20191205.rdata")

col.names <- subset(names(bd.codes.ok.nc.tx.transpose), grepl('748.3', names(bd.codes.ok.nc.tx.transpose)))

tmp <- bd.codes.ok.nc.tx.transpose %>% select(studyid, all_of(col.names))
tmp$case <- ifelse(rowSums(tmp[2:ncol(tmp)], na.rm = T) >=1, 1, 0)
tmp <- filter(tmp, case == 1)
tmp <- c(tmp$studyid)

rm(bd.codes.ok.nc.tx.transpose); gc()

load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20200123.rdata")

goback <- goback %>% filter(state %in% c('NC','TX','OK'))
goback$larynx.anom <- ifelse(goback$studyid %in% tmp, 1, 0)

table(goback$larynx.anom, goback$cancer)
table(goback$larynx.anom, goback$cancer1)
