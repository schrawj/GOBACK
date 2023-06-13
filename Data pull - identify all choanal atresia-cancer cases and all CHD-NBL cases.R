
require(dplyr)

load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20190606.1.rdata")
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/linked.registry.ids.v20190327.rdata")

choanal <- filter(goback, choanal.atresia == 1 & cancer == 1)

choanal <- left_join(select(choanal, studyid, cancer1),
                     select(goback.ids, studyid, recruitment.id),
                     by = 'studyid')

write.csv(choanal, file = '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Family-based cohort/choanal.atresia.cancer.comorbid.cases.v20190709.csv',
          row.names = F)

chd <- filter(goback, neuro == 1 & conganomalies.heart.circsys == 1)

chd <- right_join(select(goback.ids, studyid, recruitment.id),
                 select(chd, studyid, cancer1, conganomalies.heart.circsys:conganomalies.heart.circsys.other),
                 by = 'studyid')

write.csv(chd, file = '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Family-based cohort/chd.neuroblastoma.comorbid.cases.v20190709.csv',
          row.names = F)
