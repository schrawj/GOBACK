#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Generate a data frame of every unique morphology/behavior code.
#' 
#' Then, map these to the cancer diagnoses we've defined.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Prep environment.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
require(dplyr)
require(stringr)
require(readstata13)

pull.hist <- function(x){
  select(x, histology.behavior)
}

restring.columns <- function(x){
  str_replace_all(tolower(colnames(x)),'_','.')
}

#' Searches through a vector of ICD-O-3 morph/behavior codes and updates the master list.
update.dx <- function(x, y){
  ifelse(is.na(histology.behavior.codes$cancer1) & histology.behavior.codes$histology.behavior %in% x, y, histology.behavior.codes$cancer1)
}

load("Y:/Jeremy Schraw/GOBACK project/Datasets/icdo3.codes.descriptions.from.seer.rdata")
load("Y:/Jeremy Schraw/GOBACK project/Datasets/tx.nc.cancer.rdata")
load("Y:/Jeremy Schraw/GOBACK project/Datasets/Michigan/michigan.rdata")
load("Y:/Jeremy Schraw/GOBACK project/Datasets/North Carolina/nc.cancer.data.v20170821.1.rdata")
ar.new <- read.dta13('C:/Users/schraw/Downloads/bth9511vs_update.dta', convert.underscore = TRUE)



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Generate the list of unique codes.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

#' Begin by removing non-cancer children from cancer datasets.
ar.tmp <- ar.new[ar.new$primarysit != "", ]
mi.tmp <- mi[mi$any_cancer == 1, ]
nc.tmp <- nc.cancer[nc.cancer$cancer == 1, ]
tx.tmp <- tx.nc.ca[tx.nc.ca$state == 'TX', ]
tx.tmp <- tx.tmp[tx.tmp$cancer == 1, ]

#' Concatenate morphology and behavior codes.
ar.tmp$histology.behavior <- as.character(paste0(ar.tmp$HISTOLOGY3, '/', ar.tmp$BEHAVIOR3))
mi.tmp$histology.behavior <- as.character(paste0(mi.tmp$morph31, '/', mi.tmp$behave31))
nc.tmp$histology.behavior <- as.character(paste0(nc.tmp$morph31, '/', nc.tmp$behavior.code.icdo3.1))
tx.tmp$histology.behavior <- as.character(paste0(tx.tmp$morph31, '/', tx.tmp$behavior1))

#' Keep only the new variable.
ar.codes <- pull.hist(ar.tmp)
mi.codes <- pull.hist(mi.tmp)
nc.codes <- pull.hist(nc.tmp)
tx.codes <- pull.hist(tx.tmp)

histology.behavior.codes <- rbind(tx.codes, nc.codes, mi.codes, ar.codes)
histology.behavior.codes <- data.frame(histology.behavior = 
                                         histology.behavior.codes[!duplicated(histology.behavior.codes$histology.behavior), ])
histology.behavior.codes$histology.behavior <- as.character(histology.behavior.codes$histology.behavior)
histology.behavior.codes <- arrange(histology.behavior.codes, histology.behavior)

save(histology.behavior.codes, file = 'Y:/Jeremy Schraw/GOBACK project/Datasets/icdo3.codes.all.unique.codes.v20170823.1.rdata')

histology.behavior.codes[1:100, ]
histology.behavior.codes[101:200, ]
histology.behavior.codes[201:300, ]
histology.behavior.codes[301:346, ]



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Pair these with the SEER descriptions.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
histology.behavior.codes <- left_join(histology.behavior.codes, 
                                      
                                      select(seer.codes, histology.behavior, histology.behavior.description),
                                      
                                      by = 'histology.behavior')

histology.behavior.codes <- arrange(histology.behavior.codes, histology.behavior)

histology.behavior.codes$cancer1 <- NA

histology.behavior.codes[1:100, ]
histology.behavior.codes[101:200, ]
histology.behavior.codes[201:300, ]
histology.behavior.codes[301:346, ]



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Start matching these up with the cancer variables we defined. 
#' 
#' Generate a data frame of all mapped and unmapped ICD-O-3 codes in the TX cancer data. 
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

#' Prep TX cancer data for join.
tx.tmp <- tx.nc.ca[tx.nc.ca$cancer == 1, ]
tx.tmp <- tx.tmp[tx.tmp$state == 'TX', ]
tx.tmp <- subset(tx.tmp, !is.na(tx.tmp$morph31) & is.na(tx.tmp$morph32))
tx.tmp$histology.behavior <- as.character(paste0(tx.tmp$morph31, '/', tx.tmp$behavior1))

tx.tmp$cancer1 <- apply(tx.tmp[, c(60,61,63:65,67:69,71:72,74:76,78,79,81,82,83,85,86,87,89,91:94,96:99)], 1, function(x){
  names(which(x == 1))
})  

tx.tmp$cancer1 <- as.character(tx.tmp$cancer1)
tx.tmp$cancer1 <- ifelse(tx.tmp$cancer1 == 'character(0)', NA, tx.tmp$cancer1)

unique(tx.tmp$cancer1)

tx.icd03.codes <- tx.tmp[, c(204,205)]
tx.icd03.codes <- tx.icd03.codes[!duplicated(tx.icd03.codes$histology.behavior), ]

#' Join tx cancer data to the subset of codes missing cancer diagnoses.
histology.behavior.codes <- left_join(histology.behavior.codes, tx.icd03.codes, by = 'histology.behavior')

save(histology.behavior.codes, file = 'Y:/Jeremy Schraw/GOBACK project/Datasets/icdo3.codes.all.unique.codes.v20170823.2.rdata')



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' A decent head start I suppose.  We have 118 of these 294 codes mapped to our list of 
#' diagnoses.
#' 
#' Start concatenating vectors of codes that map to each DX, and update those rows with 
#' the DX.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
load("Y:/Jeremy Schraw/GOBACK project/Datasets/icdo3.codes.all.unique.codes.v20170823.2.rdata")

#' Per IARC, certain codes re-direct to others.  These are essentially duplicates and 
#' can be deleted.
#' Use 9751/3 for 9751/1.  Both codes for LCH.
histology.behavior.codes <- filter(histology.behavior.codes, histology.behavior != '9751/1')

#' Sort and view the unmapped codes.
histology.behavior.codes <- arrange(histology.behavior.codes, cancer1)

#' Generate vectors of all the codes that map to a DX.
#' CAUTION: Some codes are not unique to a single site.  These will need to be paired with site codes for a definitive mapping.

#' Vetted.
leu.other.codes <- c('9989/3','9987/3','9985/3','9983/3','9982/3','9980/3','9964/3','9962/3','9961/3','9950/3','9946/3', '9945/3',
                     '9930/3','9863/3','9860/3','9809/3','9805/3','9971/3', '9831/3', '9801/3','9800/3')

#' Vetted.
aml.codes <- c('9920/3','9865/3')

#' Vetted.
lym.other.codes <- c('9754/3','9751/3','9750/3','9741/3', '9687/3','9755/3','9590/3')

#' Vetted.
nhl.codes <- c('9719/3','9709/3','9708/3','9698/3','9684/3','9679/3','9671/3')

#' Vetted.
hl.codes <- c('9667/3')

#' Vetted.
soft.other.codes <- c('9540/3','9540/0','9540/1','9550/0','9560/0','9560/1','9560/3','9561/3','9570/0','9252/3','9161/1','9133/3',
                      '9131/0','9130/3','9121/0','9120/0','9120/3','9044/3','9043/3','9041/3','9041/1','9040/3','8990/3','8890/0','8850/0','8840/3',
                      '8833/3','8832/3','8806/3','8711/3','8710/3', '8814/3')

#' Vetted.
cns.other.codes <- c('9539/1','9538/3','9537/0','9530/0','9508/3','9505/1','9505/3','9506/1','9492/0','9451/3','9450/3','9444/1',
                     '9430/3','9413/0','9413/3','9412/1','9382/3','9362/3','9361/1','9530/3','9538/1','9351/1','9350/3','9350/1','9350/0',
                     '9390/0','9390/1','9390/3')
#' Vetted.
astro.codes <- c('9442/1')

#' Vetted.
neuro.codes <- c('9500/2')

#' Vetted.
bone.other.codes <- c('9370/3','9261/3', '9330/3', '9310/3')

#' Vetted.
osteo.codes <- c('9193/3')

#' Vetted.
renal.other.codes <- c('8964/3','8312/3')

#' Vetted.
other.any.codes <- c('8936/3','8930/3','8921/3')

#' Vetted.
pns.other.codes <- c('8700/3','9501/3')

#' Vetted.
epithe.codes <- c('8500/3')

#' Vetted.
rms.other.codes <- c('8912/3','8902/3','8901/3','8900/3')



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Update the cancer1 variable.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
histology.behavior.codes$cancer1 <- update.dx(leu.other.codes, 'leu.other')
histology.behavior.codes$cancer1 <- update.dx(aml.codes, 'aml')
histology.behavior.codes$cancer1 <- update.dx(lym.other.codes, 'lym.other')
histology.behavior.codes$cancer1 <- update.dx(nhl.codes, 'nhl')
histology.behavior.codes$cancer1 <- update.dx(hl.codes, 'hl')
histology.behavior.codes$cancer1 <- update.dx(soft.other.codes, 'soft.other')
histology.behavior.codes$cancer1 <- update.dx(cns.other.codes, 'cns.other')
histology.behavior.codes$cancer1 <- update.dx(astro.codes, 'astro')
histology.behavior.codes$cancer1 <- update.dx(neuro.codes, 'neuro')
histology.behavior.codes$cancer1 <- update.dx(bone.other.codes, 'bone.other')
histology.behavior.codes$cancer1 <- update.dx(osteo.codes, 'osteo')
histology.behavior.codes$cancer1 <- update.dx(renal.other.codes, 'renal.other')
histology.behavior.codes$cancer1 <- update.dx(other.any.codes, 'other.any')
histology.behavior.codes$cancer1 <- update.dx(pns.other.codes, 'pns.other')
histology.behavior.codes$cancer1 <- update.dx(epithe.codes, 'epithe')
histology.behavior.codes$cancer1 <- update.dx(rms.other.codes, 'rms.other')

histology.behavior.codes <- arrange(histology.behavior.codes, cancer1)

save(histology.behavior.codes, file = 'Y:/Jeremy Schraw/GOBACK project/Datasets/icdo3.codes.all.unique.codes.v20170823.3.rdata')

histology.behavior.codes[101:200, ]
histology.behavior.codes[201:300, ]
histology.behavior.codes[301:346, ]


#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Generate the second data frame: for unmapped morphology codes, one instance of every
#' unique combination of morphology and site code.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
unmapped.codes <- histology.behavior.codes[is.na(histology.behavior.codes$cancer1), ]
unmapped.codes <- c(unmapped.codes$histology.behavior)
unmapped.codes <- str_replace_all(unmapped.codes, '/[0123]', '')
unmapped.codes <- as.numeric(unmapped.codes)



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Generate the list of unique codes.
#' 
#' Remove non-cancer cases and select site and morphology code variables.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

#' AR.
ar.tmp <- ar.new[ar.new$primarysit != "", ]
ar.tmp <- select(ar.tmp, primarysit, HISTOLOGY3)
ar.tmp$primarysit <- as.numeric(str_replace(ar.tmp$primarysit, 'C', ''))
ar.tmp <- rename(ar.tmp, morphology = HISTOLOGY3, site.code = primarysit)
ar.tmp$morphology <- as.numeric(ar.tmp$morphology)

#' MI.
mi.tmp <- mi[mi$any_cancer == 1, ]
mi.tmp <- select(mi.tmp, site_code1, morph31)
mi.tmp <- rename(mi.tmp, site.code = site_code1, morphology = morph31)

#' NC.
nc.tmp <- nc.cancer[nc.cancer$cancer == 1, ]
nc.tmp <- select(nc.tmp, site.code1, morph31)
nc.tmp <- rename(nc.tmp, site.code = site.code1, morphology = morph31)

#' TX.
tx.tmp <- tx.nc.ca[tx.nc.ca$state == 'TX', ]
tx.tmp <- tx.tmp[tx.tmp$cancer == 1, ]
tx.tmp <- select(tx.tmp, site.code1, morph31)
tx.tmp <- rename(tx.tmp, site.code = site.code1, morphology = morph31)

#' Bind together; select only unmapped morph codes; remove rows with duplicated site codes; 
#' initialize a cancer variable.
site.histology.codes <- rbind(ar.tmp, mi.tmp, nc.tmp, tx.tmp)
site.histology.codes <- site.histology.codes[site.histology.codes$morphology %in% unmapped.codes, ]
site.histology.codes <- site.histology.codes[!duplicated(site.histology.codes$site.code), ]
site.histology.codes <- arrange(site.histology.codes, morphology)
site.histology.codes$cancer1 <- NA

print(site.histology.codes)

save(site.histology.codes, file = 'Y:/Jeremy Schraw/GOBACK project/Datasets/icdo3.codes.all.unique.sitemorph.codes.v20170824.1.rdata')



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' An unfortunate number of ifelse statements for mapping these codes.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
load("Y:/Jeremy Schraw/GOBACK project/Datasets/icdo3.codes.all.unique.sitemorph.codes.v20170824.1.rdata")

site.histology.codes$cancer1 <- with(site.histology.codes,
                                     ifelse(is.na(cancer1) & morphology %in% 9501:9504 & site.code %in% c(0:699,739:768,809),
                                       'pns.other', cancer1))

site.histology.codes$cancer1 <- with(site.histology.codes,
                                     ifelse(is.na(cancer1) & morphology == 9380 & site.code == 723, 'astro', 
                                       ifelse(is.na(cancer1) & morphology == 9380 & site.code %in% c(700:722,724:729,751,753),
                                              'cns.other', cancer1)))

site.histology.codes$cancer1 <- with(site.histology.codes,
                                     ifelse(is.na(cancer1) & morphology == 9364 & site.code == 649, 'renal.other',
                                       ifelse(is.na(cancer1) & morphology == 9364 & site.code %in% 400:419, 'bone.other',
                                              ifelse(is.na(cancer1) & morphology == 9364 & site.code %in% c(000:399,470:639,659:699,739:768,809),
                                                     'soft.other', cancer1))))

site.histology.codes$cancer1 <- with(site.histology.codes,
                                     ifelse(is.na(cancer1) & morphology %in% c(9070,9072) & site.code %in% c(700:729,751:753),'gct.intra',
                                            ifelse(is.na(cancer1) & morphology %in% c(9070,9072) & site.code %in% c(0:559,570:619,630:699,739:750,754:768,809),'gct.extra',
                                                   ifelse(is.na(cancer1) & morphology %in% c(9070,9072) & site.code %in% c(569,620:629), 'gct.gonad',cancer1))))

site.histology.codes$cancer1 <- with(site.histology.codes,
                                     ifelse(is.na(cancer1) & morphology == 8963 & site.code == 649, 'renal.other',
                                            ifelse(is.na(cancer1) & morphology == 8963 & site.code %in% c(0:639,659:699,739:768,809),
                                                   'soft.other', cancer1)))

site.histology.codes$cancer1 <- with(site.histology.codes,
                                     ifelse(is.na(cancer1) & morphology == 8830 & site.code %in% c(0:399,440:768,809), 'soft.other',
                                            ifelse(is.na(cancer1) & morphology == 8830 & site.code %in% 400:419, 'bone.other', cancer1)))

site.histology.codes$cancer1 <- with(site.histology.codes,
                                     ifelse(is.na(cancer1) & morphology %in% 8810:8811 & site.code %in% 400:419, 'bone.other',
                                            ifelse(is.na(cancer1) & morphology %in% 8810:8811 & site.code %in% c(0:399,440:768,809), 'soft.other',cancer1)))

site.histology.codes$cancer1 <- with(site.histology.codes,
                                     ifelse(is.na(cancer1) & morphology %in% c(8800,8801,8803:8805) & site.code %in% 400:419, 'bone.other',
                                            ifelse(is.na(cancer1) & morphology %in% 8800:8805 & site.code %in% c(0:399,440:768,809),'soft.other',cancer1)))

site.histology.codes$cancer1 <- with(site.histology.codes,
                                     ifelse(is.na(cancer1) & morphology == 8620 & site.code %in% 0:809, 'gct.gonad',cancer1))

site.histology.codes$cancer1 <- with(site.histology.codes,
                                     ifelse(is.na(cancer1) & morphology %in% 8480:8490 & site.code == 649, 'renal.other',
                                            ifelse(is.na(cancer1) & morphology %in% 8480:8490 & site.code %in% 220:221, 'hepatic.other',
                                                   ifelse(is.na(cancer1) & morphology %in% 8480:8490 & site.code %in% c(569,620:629), 'gct.gonad', 
                                                          ifelse(is.na(cancer1) & morphology %in% 8480:8490 & site.code %in% c(180,182:189,199,209,210:218),'epithe',cancer1)))))

site.histology.codes$cancer1 <- with(site.histology.codes,
                                     ifelse(is.na(cancer1) & morphology %in% 8000:8005 & site.code %in% c(700:729,751:753), 'cns.other',
                                            ifelse(is.na(cancer1) & morphology %in% 8000:8005 & site.code == 649, 'renal.other',
                                                   ifelse(is.na(cancer1) & morphology %in% 8000:8005 & site.code %in% 220:221, 'hepatic.other',
                                                          ifelse(is.na(cancer1) & morphology %in% 8000:8005 & site.code %in% 400:419, 'bone.other',
                                                                 ifelse(is.na(cancer1) & morphology %in% 8000:8005 & site.code %in% c(569,620:629),'gct.gonad',
                                                                       ifelse(is.na(cancer1) & morphology %in% 8000:8005 & site.code %in% c(0:218,239:399,420:559,570:619,630:639,659:699,739:750,754:809),'other.any', cancer1)))))))

site.histology.codes$cancer1 <- with(site.histology.codes,
                                     ifelse(is.na(cancer1) & morphology %in% c(8010:8041,8082) & site.code == 649, 'renal.other',
                                            ifelse(is.na(cancer1) & morphology %in% c(8010:8041,8082) & site.code %in% 220:221, 'hepatic.other',
                                                   ifelse(is.na(cancer1) & morphology %in% c(8010:8041,8082) & site.code %in% c(569,620:629), 'gct.gonad',
                                                          ifelse(is.na(cancer1) & morphology %in% c(8010:8041,8082) & site.code %in% c(110:199,440:449,739), 'epithe', cancer1)))))
site.histology.codes$cancer1 <- with(site.histology.codes,
                                     ifelse(is.na(cancer1) & morphology %in% c(8010:8041,8130:8141) & site.code == 649, 'renal.other',
                                            ifelse(is.na(cancer1) & morphology %in% c(8010:8041,8130:8141) & site.code %in% c(569,620:629), 'gct.gonad',
                                                   ifelse(is.na(cancer1) & morphology %in% c(8010:8041,8130:8141) & site.code %in% c(739,110:119) ,'epithe',cancer1))))

site.histology.codes$cancer1 <- with(site.histology.codes, ifelse(is.na(cancer1) & morphology == 8130 & site.code %in% 670:679, 'epithe', cancer1))
site.histology.codes$cancer1 <- with(site.histology.codes, ifelse(is.na(cancer1) & morphology == 8077 & site.code %in% c(500:539), 'epithe', cancer1))
site.histology.codes$cancer1 <- with(site.histology.codes, ifelse(is.na(cancer1) & morphology %in% 8010:8084 & site.code %in% c(90:109,239:339), 'epithe', cancer1))                                                
site.histology.codes$cancer1 <- with(site.histology.codes, ifelse(is.na(cancer1) & morphology == 8171 & site.code %in% c(0:809), 'hepatic.other', cancer1))
site.histology.codes$cancer1 <- with(site.histology.codes, ifelse(is.na(cancer1) & morphology == 8290 & site.code %in% c(739), 'epithe', cancer1))
site.histology.codes$cancer1 <- with(site.histology.codes, ifelse(is.na(cancer1) & morphology %in% 8120:8157 & site.code %in% c(0:69,90:109,129:179,239:339,380:399,480:488), 'epithe', cancer1))
site.histology.codes$cancer1 <- with(site.histology.codes, ifelse(is.na(cancer1) & morphology == 8272 & site.code %in% c(0:809), 'cns.other', cancer1))
site.histology.codes$cancer1 <- with(site.histology.codes, ifelse(is.na(cancer1) & morphology == 8310 & site.code %in% c(530:539), 'epithe', cancer1))
site.histology.codes$cancer1 <- with(site.histology.codes, ifelse(is.na(cancer1) & morphology == 8440 & site.code %in% c(739,110:119,79:89,180:189,199,209,210:218,340:349,379,500:509,670:679,690:699), 'epithe', cancer1))
site.histology.codes$cancer1 <- with(site.histology.codes, ifelse(is.na(cancer1) & morphology == 8200 & site.code %in% c(79:89), 'epithe', cancer1))
site.histology.codes$cancer1 <- with(site.histology.codes, ifelse(is.na(cancer1) & morphology %in% 8390:8420 & site.code %in% c(440:449), 'epithe', cancer1))

#' Philip clarified: pPNET outside the CNS are not included in the PNET variable.
site.histology.codes$cancer1 <- with(site.histology.codes, ifelse(is.na(cancer1) & morphology == 9364 & site.code %in% c(0:399,470:639,659:699,739:768,809), 'soft.other', cancer1))

print(site.histology.codes)

save(site.histology.codes, file = 'Y:/Jeremy Schraw/GOBACK project/Datasets/icdo3.codes.all.unique.sitemorph.codes.v20170824.2.rdata')

