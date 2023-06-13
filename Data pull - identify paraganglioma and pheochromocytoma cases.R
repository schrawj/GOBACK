
require(dplyr)

load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/cancer.codes.v20191008.1.rdata")
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20200123.rdata")
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/goback.coxph.results.v20180612.rdata")

paragang <- filter(cancer.codes, histtypeicdo3.1 %in% c(8680,8700,8693))
paragang <- paragang$studyid
paragang <- filter(goback, studyid %in% paragang)

print(filter(goback.coxmodels, cancer == 'pns.other'))
