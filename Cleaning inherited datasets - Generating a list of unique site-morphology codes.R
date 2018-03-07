#' Searches through a vector of ICD-O-3 morph/behavior codes and updates the master list.
update.dx <- function(x, y){
  ifelse(is.na(site.morph.codes$cancer1) & site.morph.codes$morphology %in% x, y, site.morph.codes$cancer1)
}








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
ar.codes <- rename(
                    select(ar.tmp, primarysit, HISTOLOGY3 ),
                    
                    site.code = primarysit, morphology = HISTOLOGY3)
ar.codes$site.code <- str_replace_all(ar.codes$site.code, '^C','')
ar.codes <- transmute(ar.codes, site.code = as.numeric(site.code), morphology = as.numeric(morphology))

mi.codes <- rename(
                    select(mi.tmp, site_code1, morph31),
                    
                    site.code = site_code1, morphology = morph31)

nc.codes<- rename(
                    select(nc.tmp, site.code1, morph31),
                    
                    site.code = site.code1, morphology = morph31)

tx.codes <- rename(
                    select(tx.tmp, site.code1, morph31),
                    
                    site.code = site.code1, morphology = morph31)


site.morph.codes <- rbind(tx.codes, nc.codes, mi.codes, ar.codes)
site.morph.codes <- site.morph.codes[!duplicated(paste0(site.morph.codes$site.code, site.morph.codes$morphology)), ]
site.morph.codes <- arrange(site.morph.codes, morphology)



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

tx.icdo3.codes <- tx.tmp[, c(36,31,205)]
tx.icdo3.codes <- rename(tx.icdo3.codes, site.code = site.code1, morphology = morph31)
tx.icdo3.codes <- tx.icdo3.codes[!duplicated(paste0(tx.icdo3.codes$site.code,tx.icdo3.codes$morphology)), ]

#' Join tx cancer data to the subset of codes missing cancer diagnoses.
site.morph.codes <- left_join(site.morph.codes, tx.icdo3.codes, by = c('site.code','morphology'))

site.morph.codes <- arrange(site.morph.codes, morphology)



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' A decent head start I suppose.  
#' 
#' Generate vectors of all the codes that map to a DX.
#' CAUTION: Some codes are not unique to a single site.  
#' These will need to be paired with site codes for a definitive mapping.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
load("Y:/Jeremy Schraw/GOBACK project/Datasets/icdo3.codes.all.unique.codes.v20170823.2.rdata")

leu.other.codes <- c(9989,9987,9985,9983,9982,9980,9964,9962,9961,9950,9946, 9945,
                     9930,9898,9863,9860,9809,9805,9971, 9831, 9801,9800, 9986, 9975, 9960:9964)

all.codes <- c(9811:9818,9826,9835:9837)

aml.codes <- c('9920','9865')

lym.other.codes <- c(9754,9751,9750,9740:9742, 9687,9755,9596,9590)

nhl.codes <- c(9727:9729,9716:9719,9714,9709,9708,9700:9702,9699,9698,9689:9691,9684,9678:9680,9671,9591)

hl.codes <- c(9667,9661:9665, 9659, 9650:9655)

soft.other.codes <- c(9581,9540,9550,9560,9560,9560,9561,9570,9252,9161,9133,
                      9131,9130,9121,9120,9120,9044,9043,9041,9041,9040,8990,8890,8850,8840,
                      8833,8832,8806,8711,8710, 8814)

cns.other.codes <- c(9530:9539,9508,9505,9506,9493,9492,9451,9450,9444,9430,9413,9412,9382,9362,9361,9538,9351,9350,9390)

retino.codes <- 9510:9514
#' Vetted.
astro.codes <- c(9440:9442, 9420:9424, 9400:9411)

#' Vetted.
neuro.codes <- c(9490,9500)

#' Vetted.
bone.other.codes <- c('9370','9261', '9330', '9310')

#' Vetted.
osteo.codes <- c('9193')

#' Vetted.
renal.other.codes <- c('8964','8312')

#' Vetted.
other.any.codes <- c('8936','8930','8921')

#' Vetted.
pns.other.codes <- c('8700','9501')

#' Vetted.
epithe.codes <- c('8500')

#' Vetted.
rms.other.codes <- c('8912','8902','8901','8900')

pnet.codes <- 9473

ependymoma.codes <- c(9383, 9391:9394)

medullo.codes <- c(9470:9472, 9474, 9480)




site.morph.codes$cancer1 <- update.dx(leu.other.codes, 'leu.other')
site.morph.codes$cancer1 <- update.dx(all.codes, 'all')
site.morph.codes$cancer1 <- update.dx(aml.codes, 'aml')
site.morph.codes$cancer1 <- update.dx(lym.other.codes, 'lym.other')
site.morph.codes$cancer1 <- update.dx(nhl.codes, 'nhl')
site.morph.codes$cancer1 <- update.dx(hl.codes, 'hl')
site.morph.codes$cancer1 <- update.dx(soft.other.codes, 'soft.other')
site.morph.codes$cancer1 <- update.dx(cns.other.codes, 'cns.other')
site.morph.codes$cancer1 <- update.dx(retino.codes, 'retino')
site.morph.codes$cancer1 <- update.dx(astro.codes, 'astro')
site.morph.codes$cancer1 <- update.dx(neuro.codes, 'neuro')
site.morph.codes$cancer1 <- update.dx(bone.other.codes, 'bone.other')
site.morph.codes$cancer1 <- update.dx(osteo.codes, 'osteo')
site.morph.codes$cancer1 <- update.dx(renal.other.codes, 'renal.other')
site.morph.codes$cancer1 <- update.dx(other.any.codes, 'other.any')
site.morph.codes$cancer1 <- update.dx(pns.other.codes, 'pns.other')
site.morph.codes$cancer1 <- update.dx(epithe.codes, 'epithe')
site.morph.codes$cancer1 <- update.dx(rms.other.codes, 'rms.other')
site.morph.codes$cancer1 <- update.dx(pnet.codes, 'pnet')
site.morph.codes$cancer1 <- update.dx(medullo.codes, 'medullo')
site.morph.codes$cancer1 <- update.dx(ependymoma.codes, 'ependymoma')

site.morph.codes$cancer1 <- with(site.morph.codes,
                                     ifelse(is.na(cancer1) & morphology %in% 9501:9504 & site.code %in% c(0:699,739:768,809), 'pns.other',
                                            ifelse(is.na(cancer1) & morphology %in% 9501:9504 & site.code %in% 700:729,'cns.other', cancer1)))

site.morph.codes$cancer1 <- with(site.morph.codes,
                                     ifelse(is.na(cancer1) & morphology == 9380 & site.code == 723, 'astro', 
                                            ifelse(is.na(cancer1) & morphology == 9380 & site.code %in% c(700:722,724:729,751,753),
                                                   'cns.other', cancer1)))

site.morph.codes$cancer1 <- with(site.morph.codes,
                                     ifelse(is.na(cancer1) & morphology == 9364 & site.code == 649, 'renal.other',
                                            ifelse(is.na(cancer1) & morphology == 9364 & site.code %in% 400:419, 'bone.other',
                                                   ifelse(is.na(cancer1) & morphology == 9364 & site.code %in% c(000:399,470:639,659:699,739:768,809),
                                                          'soft.other', cancer1))))

site.morph.codes$cancer1 <- with(site.morph.codes,
                                     ifelse(is.na(cancer1) & morphology %in% c(9070,9072) & site.code %in% c(700:729,751:753),'gct.intra',
                                            ifelse(is.na(cancer1) & morphology %in% c(9070,9072) & site.code %in% c(0:559,570:619,630:699,739:750,754:768,809),'gct.extra',
                                                   ifelse(is.na(cancer1) & morphology %in% c(9070,9072) & site.code %in% c(569,620:629), 'gct.gonad',cancer1))))

site.morph.codes$cancer1 <- with(site.morph.codes,
                                     ifelse(is.na(cancer1) & morphology == 8963 & site.code == 649, 'renal.other',
                                            ifelse(is.na(cancer1) & morphology == 8963 & site.code %in% c(0:639,659:699,739:768,809),
                                                   'soft.other', cancer1)))

site.morph.codes$cancer1 <- with(site.morph.codes,
                                     ifelse(is.na(cancer1) & morphology == 8830 & site.code %in% c(0:399,440:768,809), 'soft.other',
                                            ifelse(is.na(cancer1) & morphology == 8830 & site.code %in% 400:419, 'bone.other', cancer1)))

site.morph.codes$cancer1 <- with(site.morph.codes,
                                     ifelse(is.na(cancer1) & morphology %in% 8810:8811 & site.code %in% 400:419, 'bone.other',
                                            ifelse(is.na(cancer1) & morphology %in% 8810:8811 & site.code %in% c(0:399,440:768,809), 'soft.other',cancer1)))

site.morph.codes$cancer1 <- with(site.morph.codes,
                                     ifelse(is.na(cancer1) & morphology %in% c(8800,8801,8803:8805) & site.code %in% 400:419, 'bone.other',
                                            ifelse(is.na(cancer1) & morphology %in% 8800:8805 & site.code %in% c(0:399,440:768,809),'soft.other',cancer1)))

site.morph.codes$cancer1 <- with(site.morph.codes,
                                     ifelse(is.na(cancer1) & morphology == 8620 & site.code %in% 0:809, 'gct.gonad',cancer1))

site.morph.codes$cancer1 <- with(site.morph.codes,
                                     ifelse(is.na(cancer1) & morphology %in% 8480:8490 & site.code == 649, 'renal.other',
                                            ifelse(is.na(cancer1) & morphology %in% 8480:8490 & site.code %in% 220:221, 'hepatic.other',
                                                   ifelse(is.na(cancer1) & morphology %in% 8480:8490 & site.code %in% c(569,620:629), 'gct.gonad', 
                                                          ifelse(is.na(cancer1) & morphology %in% 8480:8490 & site.code %in% c(180,182:189,199,209,210:218),'epithe',cancer1)))))

site.morph.codes$cancer1 <- with(site.morph.codes,
                                     ifelse(is.na(cancer1) & morphology %in% 8000:8005 & site.code %in% c(700:729,751:753), 'cns.other',
                                            ifelse(is.na(cancer1) & morphology %in% 8000:8005 & site.code == 649, 'renal.other',
                                                   ifelse(is.na(cancer1) & morphology %in% 8000:8005 & site.code %in% 220:221, 'hepatic.other',
                                                          ifelse(is.na(cancer1) & morphology %in% 8000:8005 & site.code %in% 400:419, 'bone.other',
                                                                 ifelse(is.na(cancer1) & morphology %in% 8000:8005 & site.code %in% c(569,620:629),'gct.gonad',
                                                                        ifelse(is.na(cancer1) & morphology %in% 8000:8005 & site.code %in% c(0:218,239:399,420:559,570:619,630:639,659:699,739:750,754:809),'other.any', cancer1)))))))

site.morph.codes$cancer1 <- with(site.morph.codes,
                                     ifelse(is.na(cancer1) & morphology %in% c(8010:8041,8082) & site.code == 649, 'renal.other',
                                            ifelse(is.na(cancer1) & morphology %in% c(8010:8041,8082) & site.code %in% 220:221, 'hepatic.other',
                                                   ifelse(is.na(cancer1) & morphology %in% c(8010:8041,8082) & site.code %in% c(569,620:629), 'gct.gonad',
                                                          ifelse(is.na(cancer1) & morphology %in% c(8010:8041,8082) & site.code %in% c(110:199,440:449,739), 'epithe', cancer1)))))
site.morph.codes$cancer1 <- with(site.morph.codes,
                                     ifelse(is.na(cancer1) & morphology %in% c(8010:8041,8130:8141) & site.code == 649, 'renal.other',
                                            ifelse(is.na(cancer1) & morphology %in% c(8010:8041,8130:8141) & site.code %in% c(569,620:629), 'gct.gonad',
                                                   ifelse(is.na(cancer1) & morphology %in% c(8010:8041,8130:8141) & site.code %in% c(739,110:119) ,'epithe',cancer1))))

site.morph.codes$cancer1 <- with(site.morph.codes, ifelse(is.na(cancer1) & morphology == 8130 & site.code %in% 670:679, 'epithe', cancer1))
site.morph.codes$cancer1 <- with(site.morph.codes, ifelse(is.na(cancer1) & morphology == 8077 & site.code %in% c(500:539), 'epithe', cancer1))
site.morph.codes$cancer1 <- with(site.morph.codes, ifelse(is.na(cancer1) & morphology %in% 8010:8084 & site.code %in% c(90:109,239:339), 'epithe', cancer1))                                                
site.morph.codes$cancer1 <- with(site.morph.codes, ifelse(is.na(cancer1) & morphology == 8171 & site.code %in% c(0:809), 'hepatic.other', cancer1))
site.morph.codes$cancer1 <- with(site.morph.codes, ifelse(is.na(cancer1) & morphology == 8290 & site.code %in% c(739), 'epithe', cancer1))
site.morph.codes$cancer1 <- with(site.morph.codes, ifelse(is.na(cancer1) & morphology %in% 8120:8157 & site.code %in% c(0:69,90:109,129:179,239:339,380:399,480:488), 'epithe', cancer1))
site.morph.codes$cancer1 <- with(site.morph.codes, ifelse(is.na(cancer1) & morphology == 8272 & site.code %in% c(0:809), 'cns.other', cancer1))
site.morph.codes$cancer1 <- with(site.morph.codes, ifelse(is.na(cancer1) & morphology == 8310 & site.code %in% c(530:539), 'epithe', cancer1))
site.morph.codes$cancer1 <- with(site.morph.codes, ifelse(is.na(cancer1) & morphology == 8440 & site.code %in% c(739,110:119,79:89,180:189,199,209,210:218,340:349,379,500:509,670:679,690:699), 'epithe', cancer1))
site.morph.codes$cancer1 <- with(site.morph.codes, ifelse(is.na(cancer1) & morphology == 8200 & site.code %in% c(79:89), 'epithe', cancer1))
site.morph.codes$cancer1 <- with(site.morph.codes, ifelse(is.na(cancer1) & morphology %in% 8390:8420 & site.code %in% c(440:449), 'epithe', cancer1))

#' Philip clarified: pPNET outside the CNS are not included in the PNET variable.
site.morph.codes$cancer1 <- with(site.morph.codes, ifelse(is.na(cancer1) & morphology == 9364 & site.code %in% c(0:399,470:639,659:699,739:768,809), 'soft.other', cancer1))








site.morph.codes <- arrange(site.morph.codes, cancer1)
site.morph.codes <- arrange(site.morph.codes, morphology)


site.morph.codes[1:250, ]
site.morph.codes[251:500, ]
site.morph.codes[501:750, ]
site.morph.codes[751:1000, ]
site.morph.codes[1001:1250, ]
site.morph.codes[1250:1397, ]

table(is.na(site.morph.codes))

save(site.morph.codes, file = 'Y:/Jeremy Schraw/GOBACK project/Datasets/icdo3.to.dx.mappings.v20170824.1.rdata')



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' The first time I assigned cancer1 in the AR dataset, 5 site/morphology codes did not 
#' map to a DX.
#' 
#' update our site.morph codes with these AR cancer codes that didn't match a DX.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
load("Z:/GOBACK/Jeremy/Datasets/Arkansas/unmatched.arkansas.cancers.v20170908.rdata")
print(ar.unmatched.cancers[,c(14,15,122)])

str(site.morph.codes)

#' Build a new data frame to match the codes and assign them a DX.
site.morph.codes2 <- data.frame(morphology = as.numeric(c(9591,9381,9380,9133,9084)),
                                site.code = as.numeric(c(380,718,700,172,710)),
                                cancer1 = as.character(c('nhl', 'cns.other','cns.other','soft.other','gct.intra')))

site.morph.codes <- rbind(site.morph.codes, site.morph.codes2)
site.morph.codes$cancer1 <- str_replace_all(site.morph.codes$cancer1, '_','.')


save(site.morph.codes, file= 'Z:/GOBACK/Jeremy/Datasets/iccc.codes.to.dx.mappings.v20170908.rdata')

