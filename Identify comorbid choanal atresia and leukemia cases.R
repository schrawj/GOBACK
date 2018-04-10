require(dplyr)

goback.nochrom <- filter(goback.nochrom, any.birthdefect == 1)
cho.leu <- filter(goback.nochrom, leu.any == 1 & choanal.atresia == 1)
cho.leu <- c(cho.leu$studyid)

load("//discovery2/bcm-pedi-genepi2/Jeremy/GOBACK/Datasets/cancer.codes.v20180227.1.rdata")

cho.leu <- filter(cancer.codes, studyid %in% cho.leu)
cho.leu <- select(cho.leu, studyid, morph31, site_code1, behavior1)

write.csv(cho.leu, file = "//discovery2/bcm-pedi-genepi2/Jeremy/GOBACK/R outputs/Cancer codes in comorbid cases/leu.any.csv",
          row.names = FALSE)
