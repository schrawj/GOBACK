#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.02.14.
#' 
#' A few of the top hits in GOBACK were for specific defects in association
#' with either 'non-RMS soft tissue sarcomas' or 'extracranial germ cell
#' tumors.'
#' 
#' These definitions are both somewhat heterogeneous.  Find IDs for these 
#' comorbid cases and use them to pull their site and morphology codes to 
#' get more granular information on their diagnosis.
#' 
#' As it happens, all the defects that popped out for these two cancers
#' are non-chromosomal.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------



# prep environment --------------------------------------------------------

require(dplyr); require(xlsx); require(haven)



# Generate data frames with cancer info for non-RMS soft tissue -----------

setwd('Z:/Jeremy/GOBACK/Datasets'); load('goback.no.chrom.v20180122.1.rdata'); load("./Old Datasets/cancer.codes.rdata")

soft.defects <- c('atrialseptaldefect','craniosynostosis','obstructive.genitourinary.defects','spinabifida.wo.anencephaly')

sarcoma <- list()

for (i in 1:length(soft.defects)){
  tmp <- filter(goback.nochrom, goback.nochrom$cancer1 == 'soft.other' & goback.nochrom[,as.character(soft.defects[i])] == 1)
  tmp <- c(tmp$studyid)
  sarcoma[[i]] <- data.frame(select(filter(cancer.codes, studyid %in% tmp), studyid, morph31, site_code1))
}

sarcoma <- lapply(sarcoma, function(x) {arrange(x, by = morph31)})
names(sarcoma) <- soft.defects

tumor.types <- data.frame(morph31 = c(8800,8806,8814,8832,8850,8851,9540,9560, 9581, 9121, 9131, 9133),
                          tumor.type = c('Unspecified soft tissue sarcomas',
                                         'miscellaneous soft tissue sarcomas',
                                         'fibroblastic and myofibroblastic tumors', #Not always the case, but true for all site codes in these sets.
                                         'fibrohistiocytic tumors','liposarcoma','liposarcoma','nerve sheath tumor','nerve sheath tumor', 'alveolar soft parts sarcoma',
                                         'blood vessel tumors','blood vessel tumors','blood vessel tumors'))
sites <- data.frame(site_code1 = c(23,220,445,481,
                                   490,492,496,
                                   710,715,717,719,
                                   720,721,724,729),
                    site = c(' Anterior 2/3 of tongue, NOS','Liver','Skin of trunk','Specified parts of peritoneum',
                           'Connective, Subcutaneous and other soft tissues of head, face, and neck','Connective, Subcutaneous and other soft tissues of lower limb and hip','Connective, Subcutaneous and other soft tissues of trunk, NOS',
                           'Cerebrum','Brain, Ventricle, NOS','Brain, Brain stem','Brain, NOS',
                           'Spinal cord','Cauda equina','Acoustic nerve','Nervous System, NOS'))

sarcoma <- lapply(sarcoma, function(x) {left_join(x, tumor.types, by = 'morph31')})
sarcoma <- lapply(sarcoma, function(x) {left_join(x, sites, by = 'site_code1')})

for( i in 1:length(sarcoma)){
  write.xlsx(sarcoma[i], file = 'C:/Users/schraw/Desktop/non.rms.sarcomas.xlsx', append = TRUE, row.names = FALSE, sheetName = names(sarcoma[i]))
}



# Generate data frames with cancer info for extracranial GCTs -------------

setwd('Z:/Jeremy/GOBACK/Datasets'); load('goback.no.chrom.v20180122.1.rdata'); load("./Old Datasets/cancer.codes.rdata")

ect.defects <- c('congenital.hip.dislocation','obstructive.genitourinary.defects','patentductusarteriosis')

ect <- list()

for (i in 1:length(ect.defects)){
  tmp <- filter(goback.nochrom, goback.nochrom$cancer1 == 'gct.extra' & goback.nochrom[,as.character(ect.defects[i])] == 1)
  tmp <- c(tmp$studyid)
  ect[[i]] <- data.frame(select(filter(cancer.codes, studyid %in% tmp), studyid, morph31, site_code1))
}

ect <- lapply(ect, function(x) {arrange(x, by = morph31)})
names(ect) <- ect.defects

tumor.types <- data.frame(morph31 = c(rep(9080, 3),rep(9085,2)),
                          site_code1 = c(381,383,495,495,763),
                          tumor.type = c(rep('malignant teratomas of extracranial and extragonadal sites',3),'Intracranial and intraspinal tumors of mixed form',
                                         'Malignant gonadal tumors of mixed forms'))

ect <- lapply(ect, function(x) {left_join(x, tumor.types, by = c('morph31','site_code1'))})

for( i in 1:length(ect)){
  write.xlsx(ect[i], file = 'C:/Users/schraw/Desktop/extracranial.gcts.xlsx', append = TRUE, row.names = FALSE, sheetName = names(ect[i]))
}

# Compute mean # defects for kids w/top defects ---------------------------

setwd('Z:/Jeremy/GOBACK/Datasets/')
load("goback.cox.ph.top.hits.v20180131.1.rdata")
load('goback.no.chrom.v20180122.1.rdata')

def <- c(top.hits$defect[c(1:8,12:32,34)]) #' Remove chromosomal anomalies for now.
can <- c(top.hits$cancer[c(1:8,12:32,34)])

tmp <- filter(goback.nochrom, cancer == 1)

num.def <- data.frame(defect = as.character(),
                      cancer = as.character(),
                      av.num.defects = as.numeric())

for (i in 1:30){
  tmp2 <- filter(tmp, cancer1 == can[i])
  av.num.defects <- aggregate(defect.total ~ tmp2[,def[i]], data = tmp2, mean)[2,2]
  tmp2 <- data.frame(defect = def[i],
                     cancer = can[i],
                     av.num.defects = av.num.defects)
  num.def <- rbind(num.def, tmp2)
}

rm(goback.nochrom); gc()

load('goback.chrom.v20180122.1.rdata')

tmp <- filter(goback.chrom, cancer == 1)

def <- c(top.hits$defect[c(9:11,33)]) #' Only chromosomal anomalies.
can <- c(top.hits$cancer[c(9:11,33)])

for (i in 1:4){
  tmp2 <- filter(tmp, cancer1 == can[i])
  av.num.defects <- aggregate(defect.total ~ tmp2[,def[i]], data = tmp2, mean)[2,2]
  tmp2 <- data.frame(defect = def[i],
                     cancer = can[i],
                     av.num.defects = av.num.defects)
  num.def <- rbind(num.def, tmp2)
}

top.hits <- left_join(top.hits, num.def, by = c('defect','cancer'))

save(top.hits, file = 'goback.cox.ph.top.hits.v20180214.1.rdata')
write.csv(top.hits, file = 'Z:/Jeremy/GOBACK/R outputs/goback.cox.ph.top.hits.v20180214.1.csv', row.names = FALSE)



# Do all non-RMS soft tissue sarcomas have malignant codes? ---------------

setwd('Z:/Jeremy/GOBACK/Datasets/')
load('./Old Datasets/cancer.codes.rdata')

tmp2 <- read_sas('./Texas/tcr_17_007_july20.sas7bdat')
tmp2$studyid <- paste0('tx',tmp2$birthID)
tmp2 <- arrange(tmp2, studyid, DXDATE) #' Arrange duplicated rows in order of DX.
tmp2 <- tmp2[!duplicated(tmp2$studyid), ] # Remove duplicates (will remove cancers DX'd after first primary)

tmp <- left_join(cancer.codes[grepl('tx', cancer.codes$studyid), ],
                 select(tmp2, studyid, BEHAVIORICDO3),
                 by = 'studyid')
table(is.na(tmp$BEHAVIORICDO3)) # 0 NA values implies correct matching on ID.

cancer.codes.update <- rename(tmp, behavior1 = BEHAVIORICDO3); rm(tmp, tmp2)

soft.defects <- c('atrialseptaldefect','craniosynostosis','obstructive.genitourinary.defects','spinabifida.wo.anencephaly')

sarcoma <- list()

goback.nochrom <- filter(goback.nochrom, state == 'TX')

for (i in 1:length(soft.defects)){
  tmp <- filter(goback.nochrom, goback.nochrom$cancer1 == 'soft.other' & goback.nochrom[,as.character(soft.defects[i])] == 1)
  tmp <- c(tmp$studyid)
  sarcoma[[i]] <- data.frame(select(filter(cancer.codes.update, studyid %in% tmp), studyid, morph31, site_code1, behavior1))
}

sarcoma <- lapply(sarcoma, function(x) {arrange(x, by = morph31)})
names(sarcoma) <- soft.defects

tumor.types <- data.frame(morph31 = c(8800,8806,8814,8832,8850,8851,9540,9560, 9581, 9121, 9131, 9133),
                          tumor.type = c('Unspecified soft tissue sarcomas',
                                         'miscellaneous soft tissue sarcomas',
                                         'fibroblastic and myofibroblastic tumors', #Not always the case, but true for all site codes in these sets.
                                         'fibrohistiocytic tumors','liposarcoma','liposarcoma','nerve sheath tumor','nerve sheath tumor', 'alveolar soft parts sarcoma',
                                         'blood vessel tumors','blood vessel tumors','blood vessel tumors'))
sites <- data.frame(site_code1 = c(23,220,445,481,
                                   490,492,496,
                                   710,715,717,719,
                                   720,721,724,729),
                    site = c(' Anterior 2/3 of tongue, NOS','Liver','Skin of trunk','Specified parts of peritoneum',
                             'Connective, Subcutaneous and other soft tissues of head, face, and neck','Connective, Subcutaneous and other soft tissues of lower limb and hip','Connective, Subcutaneous and other soft tissues of trunk, NOS',
                             'Cerebrum','Brain, Ventricle, NOS','Brain, Brain stem','Brain, NOS',
                             'Spinal cord','Cauda equina','Acoustic nerve','Nervous System, NOS'))

sarcoma <- lapply(sarcoma, function(x) {left_join(x, tumor.types, by = 'morph31')})
sarcoma <- lapply(sarcoma, function(x) {left_join(x, sites, by = 'site_code1')})

for( i in 1:length(sarcoma)){
  write.xlsx(sarcoma[i], file = 'C:/Users/schraw/Desktop/non.rms.sarcomas.xlsx', append = TRUE, row.names = FALSE, sheetName = names(sarcoma[i]))
}




