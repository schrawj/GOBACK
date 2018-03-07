#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#'                                       
#'                              NORTH CAROLINA
#'                                    
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------



# Prep environment --------------------------------------------------------
setwd('Z:/Jeremy/GOBACK/Datasets/North Carolina/')

require(dplyr)
require(stringr)



# User-defined functions --------------------------------------------------
restring.columns <- function(x){
  stringr::str_replace_all(tolower(colnames(x)),'_','.')
}

rename.tiff.vars <- function(data.frame){
  rename(data.frame, 
         
         conganomalies.ear.face.neck.other.major = ear.face.neck.other.major,
         conganomalies.ear.face.neck.other.minor = ear.face.neck.other.minor,
         total.anomalous.pulmonary.venous.return = totalanompulmvenousreturn,
         rvot.defects = rvot, lvot.defects = lvot, 
         interrupted.aortic.arch.type.a.or.c = iaa.aorc,
         lung.agenesis.hypoplasia = lungagenesis.hypoplasia,
         choanal.atresia = choanalatresia,
         cleft.palate.wo.cleft.lip = cleftpalate,
         esophageal.atre.tracheofist = esophagealatre.tracheofist,
         pyloric.stenosis = pyloricstenosis,
         hirshsprung.disease = hirshsprungdisease,
         biliary.atresia = biliaryatresia,
         small.intestinal.atresia = smallintestatresiastenosis,
         conganomalies.genitalandurinary = conganomgenitalandurinary,
         renal.agenesis.hypoplasia = renalagenesis.hypoplasia,
         bladder.exstrophy = bladderexstrophy,
         obstructive.genitourinary.defects = obstructivegenitourinarydefects,
         upper.limb.reduction.deformities = upperlimbreductiondeformities,
         lower.limb.reduction.deformities = lowerlimbreductiondeformities,
         limb.deformities.unspecified = limbdeformities.unspecified,
         congenital.hip.dislocation = congenitalhipdislocation,
         diagphragmatic.hernia = diaphragmatichernia,
         conganomalies.integument = conganomaliesintegument,
         di.george.syndrome = georgesyndrome,
         chromosomalanomalies.other.major = chromanom.other.major,
         chromosomalanomalies.other.minor = chromanom.other.minor,
         common.truncus = commontruncus,
         aortic.valve.stenosis = aorticvalvestenosis,
         oral.clefts = cleftpalateandcleftlip,
         cleft.lip.w.and.wo.cleft.palate = cleftlip.w.and.wo.cleftpalate,
         conganomalies.digestivesystem = conganomaliesdigestivesystem,
         conganomalies.musculoskelsys = conganomaliesmusculoskelsys,
         down.syndrome = downsyndrome,
         turner.syndrome = turnersyndrome,
         transposition.of.greatvessels = transpositionofgreatvessels,
         tetralogy.of.fallot = tetralogyoffallot,
         
         leu.any = leu, lym.any = lym, pns.any = pns, rms.any = rms, 
         soft.any = soft, gct.any = gct)
}

compute.cancer1 <- function(hist, site){
ifelse(hist %in% cancer.codes$all.non.contingent.codes, 'all',
ifelse(hist %in% cancer.codes$all.contingent.codes & site %in% cancer.codes$all.contingent.codes.sites, 'all', 
ifelse(hist %in% cancer.codes$aml.codes, 'aml', 
ifelse(hist %in% cancer.codes$other.hem.codes, 'leu.other', 
ifelse(hist %in% cancer.codes$hl.codes, 'hl',
ifelse(hist %in% cancer.codes$nhl.non.contingent.codes, 'nhl',
ifelse(hist %in% cancer.codes$nhl.contingent.codes & site %in% cancer.codes$nhl.contingent.codes.sites, 'nhl',
ifelse(hist %in% cancer.codes$lym.other.codes, 'lym.other',
ifelse(hist %in% cancer.codes$ependymoma.codes, 'ependymoma',
ifelse(hist %in% cancer.codes$astro.non.contingent.codes, 'astro',
ifelse(hist %in% cancer.codes$astro.contingent.codes & site %in% cancer.codes$astro.contingent.codes.sites, 'astro',
ifelse(hist %in% cancer.codes$medullo.codes, 'medullo',
ifelse(hist %in% cancer.codes$pnet.codes, 'pnet',
ifelse(hist %in% cancer.codes$cns.other.non.contingent.codes, 'cns.other',
ifelse(hist %in% cancer.codes$cns.other.contingent.codes1 & site %in% cancer.codes$cns.other.contingent.codes1.sites, 'cns.other',
ifelse(hist %in% cancer.codes$cns.other.contingent.codes2 & site %in% cancer.codes$cns.other.contingent.codes2.sites, 'cns.other',
ifelse(hist %in% cancer.codes$cns.other.contingent.codes3 & site %in% cancer.codes$cns.other.contingent.codes3.sites, 'cns.other',
ifelse(hist %in% cancer.codes$neuroblast.codes, 'neuro',
ifelse(hist %in% cancer.codes$pns.other.non.contingent.codes, 'pns.other',
ifelse(hist %in% cancer.codes$pns.other.contingent.codes & site %in% cancer.codes$pns.other.contingent.codes.sites, 'pns.other',      
ifelse(hist %in% cancer.codes$retino.codes, 'retino',
ifelse(hist %in% cancer.codes$nephro.codes, 'nephro',
ifelse(hist %in% cancer.codes$renal.other.non.contingent.codes, 'renal.other',
ifelse(hist %in% cancer.codes$renal.other.contingent.codes & site %in% cancer.codes$renal.other.contingent.codes.sites, 'renal.other',
ifelse(hist %in% cancer.codes$hepatoblast.codes, 'hepato',
ifelse(hist %in% cancer.codes$hepatic.other.non.contingent.codes, 'hepatic.other',
ifelse(hist %in% cancer.codes$hepatic.other.contingent.codes & site %in% cancer.codes$hepatic.other.contingent.codes.sites, 'hepatic.other',
ifelse(hist %in% cancer.codes$osteo.codes & site %in% cancer.codes$osteo.codes.sites, 'osteo',
ifelse(hist %in% cancer.codes$ewing.contingent.codes1 & site %in% cancer.codes$ewing.contingent.codes1.sites, 'ewing',
ifelse(hist %in% cancer.codes$ewing.contingent.codes2 & site %in% cancer.codes$ewing.contingent.codes2.sites, 'ewing', 'other.cancer'))))))))))))))))))))))))))))))
}
compute.cancer2 <- function(hist, site, newvar){
ifelse(hist %in% cancer.codes$bone.other.non.contingent.codes, 'bone.other',
ifelse(hist %in% cancer.codes$bone.other.contingent.codes1 & site %in% cancer.codes$bone.other.contingent.codes1.sites, 'bone.other',
ifelse(hist %in% cancer.codes$bone.other.contingent.codes2 & site %in% cancer.codes$bone.other.contingent.codes2.sites, 'bone.other',  
ifelse(hist %in% cancer.codes$erms.codes, 'erms',
ifelse(hist %in% cancer.codes$arms.codes, 'arms',      
ifelse(hist %in% cancer.codes$rms.other.codes, 'rms.other',
ifelse(hist %in% cancer.codes$soft.other.non.contingent.codes, 'soft.other',
ifelse(hist %in% cancer.codes$soft.other.contingent.codes1 & site %in% cancer.codes$soft.other.contingent.codes1.sites, 'soft.other',
ifelse(hist %in% cancer.codes$soft.other.contingent.codes2 & site %in% cancer.codes$soft.other.contingent.codes2.sites, 'soft.other',
ifelse(hist %in% cancer.codes$soft.other.contingent.codes3 & site %in% cancer.codes$soft.other.contingent.codes3.sites, 'soft.other',
ifelse(hist %in% cancer.codes$soft.other.contingent.codes4 & site %in% cancer.codes$soft.other.contingent.codes4.sites, 'soft.other',
ifelse(hist %in% cancer.codes$soft.other.contingent.codes5 & site %in% cancer.codes$soft.other.contingent.codes5.sites, 'soft.other',
ifelse(hist %in% cancer.codes$soft.other.contingent.codes6 & site %in% cancer.codes$soft.other.contingent.codes6.sites, 'soft.other',
ifelse(hist %in% cancer.codes$soft.other.contingent.codes7 & site %in% cancer.codes$soft.other.contingent.codes7.sites, 'soft.other',
ifelse(hist %in% cancer.codes$soft.other.contingent.codes8 & site %in% cancer.codes$soft.other.contingent.codes8.sites, 'soft.other',
ifelse(hist %in% cancer.codes$intra.gct.codes & site %in% cancer.codes$intra.gct.codes.sites, 'gct.intra',
ifelse(hist %in% cancer.codes$extra.gct.codes & site %in% cancer.codes$extra.gct.codes.sites, 'gct.extra',
ifelse(hist %in% cancer.codes$gonad.gct.codes & site %in% cancer.codes$gonad.gct.codes.sites, 'gct.gonad',
ifelse(hist %in% cancer.codes$other.unspec.non.contingent.codes, 'other.any',
ifelse(hist %in% cancer.codes$other.unspec.contingent.codes1 & site %in% cancer.codes$other.unspec.contingent.codes1.sites, 'other.any',
ifelse(hist %in% cancer.codes$other.unspec.contingent.codes2 & site %in% cancer.codes$other.unspec.contingent.codes2.sites, 'other.any',
ifelse(hist %in% cancer.codes$other.unspec.contingent.codes3 & site %in% cancer.codes$other.unspec.contingent.codes3.sites, 'other.any',
ifelse(hist %in% cancer.codes$other.unspec.contingent.codes4 & site %in% cancer.codes$other.unspec.contingent.codes4.sites, 'other.any',
ifelse(hist %in% cancer.codes$epithe.non.contingent.codes, 'epithe',
ifelse(hist %in% cancer.codes$epithe.contingent.codes1 & site %in% cancer.codes$epithe.contingent.codes1.sites, 'epithe',
ifelse(hist %in% cancer.codes$epithe.contingent.codes2 & site %in% cancer.codes$epithe.contingent.codes2.sites, 'epithe',
ifelse(hist %in% cancer.codes$epithe.contingent.codes3 & site %in% cancer.codes$epithe.contingent.codes3.sites, 'epithe',
ifelse(hist %in% cancer.codes$epithe.contingent.codes4 & site %in% cancer.codes$epithe.contingent.codes4.sites, 'epithe',
ifelse(hist %in% cancer.codes$epithe.contingent.codes5 & site %in% cancer.codes$epithe.contingent.codes5.sites, 'epithe',
ifelse(hist %in% cancer.codes$epithe.contingent.codes6 & site %in% cancer.codes$epithe.contingent.codes6.sites, 'epithe',
ifelse(is.na(hist), NA, newvar)))))))))))))))))))))))))))))))
}

get.standard.dem.vars <- function(x){
  require(dplyr)
  require(stringr)
  colnames(x) <- str_replace_all(colnames(x),'_','.')
  colnames(x) <- tolower(colnames(x))
  select(x, studyid, state, sex, m.race, m.edu2, m.age, birth.wt, birth.wt.cat, gest.age, birth.yr, plu, f.race, f.edu2, f.age, person.yrs)
}

populate.cancer.var <- function(x,y,z){
  if(missing(z)){
    ifelse((x | y) == 1, 1, 0)
  }
  else{
    ifelse((x | y | z) == 1, 1, 0)
  }
}

get.standard.dem.vars <- function(x){
  require(dplyr)
  require(stringr)
  colnames(x) <- str_replace_all(colnames(x),'_','.')
  colnames(x) <- tolower(colnames(x))
  select(x, studyid, state, sex, m.race, m.edu2, m.age, birth.wt, birth.wt.cat, gest.age, birth.yr, plu, f.race, f.edu2, f.age, person.yrs)
}
get.standard.def.vars <- function(x){
  require(dplyr)
  require(stringr)
  colnames(x) <- str_replace_all(colnames(x),'_','.')
  colnames(x) <- tolower(colnames(x))
  select(x,   
         studyid, any.birthdefect,                                   minor.status,                                     
         overall.chrom,                                     defect.total,                                     
         majordefect.total,                                 minordefect.total,                                
         conganomalies.cns,                                 anencephalus,                                     
         spinabifida.wo.anencephaly,                        hydrocephalus.wo.spinabifida,                     
         encephalocele,                                     microcephalus,                                    
         holoprosencephaly,                                 conganomalies.cns.other.major,                    
         conganomalies.cns.other.minor,                     conganomalies.eye,                                
         anopthalmos.micropthalmos,                         congenitalcataract,                              
         aniridia,                                          conganomalies.eye.other.major,                    
         conganomalies.eye.other.minor,                     conganomalies.ear.face.neck,                      
         anotia.microtia,                                   conganomalies.ear.face.neck.other.major,          
         conganomalies.ear.face.neck.other.minor,           conganomalies.heart.circsys,
         conotruncal.defects,                               common.truncus,
         tetralogy.of.fallot,                               transposition.of.greatvessels,
         avsd,                                              endocardialcushiondefect,
         apvr,                                              total.anomalous.pulmonary.venous.return,
         lvot.defects,                                      hypoplasticleftheartsyndrome,
         interrupted.aortic.arch.type.a.or.c,               coarctationofaorta,
         aortic.valve.stenosis,                             rvot.defects,
         pulmvalveatresiaandstenosis,                       ebsteinanomaly,
         septal.defects,                                    ventricularseptaldefect,
         atrialseptaldefect,                                singleventricle,
         trivalveatresiaandstenosis,                        patentductusarteriosis,
         pulmonaryarteryanomalies,                          heart.circsys.other.major,
         heart.circsys.other.minor,
         conganomalies.respsys,                            
         lung.agenesis.hypoplasia,                          choanal.atresia,                                  
         respsys.other.major,                               respsys.other.minor,                              
         oral.clefts,                                       cleft.palate.wo.cleft.lip,                        
         cleft.lip.w.and.wo.cleft.palate,                                      
         conganomalies.digestivesystem,                     esophageal.atre.tracheofist,                      
         rectal.largeintestatresia.sten,                    pyloric.stenosis,                                 
         hirshsprung.disease,                               biliary.atresia,                                  
         small.intestinal.atresia,                          digestivesystem.other.major,                      
         digestivesystem.other.minor,                       conganomalies.genitalandurinary,                  
         renal.agenesis.hypoplasia,                         bladder.exstrophy,                                
         obstructive.genitourinary.defects,                 hypospadias,              
         epispadias,               genitalandurinary.other.major,                    
         genitalandurinary.other.minor,                     conganomalies.musculoskelsys,                     
         limb.deformities.unspecified,                      upper.limb.reduction.deformities,                 
         lower.limb.reduction.deformities,                  gastroschisis.omphalocele,                        
         congenital.hip.dislocation,                        diagphragmatic.hernia,                            
         clubfoot,                                          craniosynostosis,                                 
         musculoskelsys.other.major,                        musculoskelsys.other.minor,                       
         conganomalies.integument,                          chromosomalanomalies,                             
         trisomy13,                                         down.syndrome,                                    
         trisomy18,                                         turner.syndrome,                                  
         di.george.syndrome,                                chromosomalanomalies.other.major,                 
         chromosomalanomalies.other.minor,                  other.unspeccongenitalanomalies) 
}
get.standard.cancer.vars <- function(x){
  require(dplyr)
  require(stringr)
  colnames(x) <- str_replace_all(colnames(x),'_','.')
  colnames(x) <- tolower(colnames(x))
  select(x, studyid,  cancer,        num.diagnoses, cancer1,       cancertime1,   laterality1,   all,           aml,           leu.any,       leu.other,    
         hl,            nhl,           lym.any,       lym.other,     astro,         medullo,       cns.any,       cns.other,     neuro,        
         pns.any,       pns.other,     retino,        nephro,        renal.any,     renal.other,   hepato,        hepatic.any,   hepatic.other,
         osteo,         ewing,         bone.any,      bone.other,    erms,          arms,          rms.any,       rms.other,     soft.any,     
         soft.other,    gct.intra,     gct.extra,     gct.gonad,     gct.any,       epithe,        other.any,     ependymoma,    pnet)  
}



# Load and save raw data --------------------------------------------------
nc.raw <- read.csv(file = './nc.raw.data.csv', header = TRUE, stringsAsFactors = FALSE)

colnames(nc.raw) <- restring.columns(nc.raw)
nc.raw <- rename.tiff.vars(nc.raw)

save(nc.raw, file = './nc.raw.data.rdata')

nc <- nc.raw
rm(nc.raw)
save(nc, file = './nc.v20171102.1.rdata')



# Demographics variables --------------------------------------------------
load('./nc.v20171102.1.rdata')

#' Lay of the land.

sink('C:/Users/schraw/Desktop/nc.var.report.txt')

for (i in 136:156){
  
  if ((is.character(nc[,i])) | is.logical(nc[,i])){
    print(names(nc[i]))
    print(table(nc[,i], useNA = 'always'))
  }
  
  else{
    print(names(nc[i]))
    print(summary(nc[,i]))
  }
}

sink()

#' Set NA sex values to 9.
nc$sex <- ifelse(is.na(nc$sex), 9, nc$sex)

#' Filter extreme maternal ages.
nc$m.age <- ifelse(nc$m.age < 13 | nc$m.age > 50, NA, nc$m.age)

#' Set NA plu values to 9.
nc$plu <- ifelse(is.na(nc$plu), 9, nc$plu)

#' Filter extreme paternal ages.  LIkely value 99 has to be a missing value code.
nc$f.age <- ifelse(nc$f.age == 99, NA, nc$f.age)

#' append state to studyid.
nc$studyid <- paste0('nc',nc$studyid)

save(nc, file = './nc.v20171102.2.rdata')



# Person-years ------------------------------------------------------------
load('./nc.v20171102.2.rdata')
load('./nc.birth.certificate.data.rdata')

nc.nocancer.nodeath <- filter(nc, cancer == 2 & infant.death == 0)
nc.nocancer.death <- filter(nc, cancer == 2 & infant.death == 1)

nc.single.cancer <- left_join(filter(nc, cancer == 1 & num.diagnoses == 1), select(nc.sas, studyid, age.at.diagnosis.1), by = 'studyid')
nc.single.cancer$age.at.diagnosis.1 <- as.numeric(nc.single.cancer$age.at.diagnosis.1)

nc.mult.cancer <- left_join(filter(nc, cancer == 1 & num.diagnoses > 1), select(nc.sas, studyid, age.at.diagnosis.1, age.at.diagnosis.2, age.at.diagnosis.3, age.at.diagnosis.4, age.at.diagnosis.5), by = 'studyid')
for (i in 199:203){
  nc.mult.cancer[,i] <- ifelse(nc.mult.cancer[,i] == '',NA,nc.mult.cancer[,i])
  nc.mult.cancer[,i] <- as.numeric(nc.mult.cancer[,i])
}
nc.mult.cancer$person.years <- apply(nc.mult.cancer[,199:203], 1, min, na.rm = TRUE)
nc.mult.cancer$person.years <- ifelse(nc.mult.cancer$person.years == 0, 0.5, nc.mult.cancer$person.years)

nc.nocancer.nodeath$person.years <- 2012 - nc.nocancer.nodeath$yearbth
nc.nocancer.nodeath$person.years <- ifelse(nc.nocancer.nodeath$person.years == 0, 0.5, nc.nocancer.nodeath$person.years)

nc.nocancer.death$person.years <- (nc.nocancer.death$agedth/365)

nc.single.cancer$person.years <- nc.single.cancer$age.at.diagnosis.1
nc.single.cancer$person.years <- ifelse(nc.single.cancer$person.years == 0, 0.5, nc.single.cancer$person.years)

person.years <- rbind(select(nc.single.cancer, studyid, person.years),
                select(nc.mult.cancer, studyid, person.years),
                select(nc.nocancer.nodeath, studyid, person.years),
                select(nc.nocancer.death, studyid, person.years))

nc <- left_join(nc, person.years, by = 'studyid')

nc <- rename(select(nc, -person.years.x), person.yrs = person.years.y)

rm(nc.nocancer.nodeath, nc.nocancer.death, nc.single.cancer, nc.mult.cancer, person.years, nc.sas)

save(nc, file = './nc.v20171103.1.rdata')



# Birth defects variables -------------------------------------------------
load('./nc.v20171103.1.rdata')

#' Lay of the land.
sink('C:/Users/schraw/Desktop/nc.bd.var.report.txt')
for (i in 4:135){
  print(names(nc[i]))
  print(table(nc[,i], useNA = 'always'))
}
sink()

#' Overall.chrom is a disaster.
nc <- select(nc, -overall.chrom)

nc$any.birthdefect <- ifelse(is.na(nc$any.birthdefect), 0, nc$any.birthdefect)
nc$minor.status <- ifelse(nc$any.birthdefect == 0 & nc$minor.status == '', "No Minor Defects", nc$minor.status)
nc$overall.chrom <- apply(nc[,103:111], 1, sum, na.rm = TRUE)

#' Major and minor defect total.
for (i in 5:6){
    nc[,i] <- ifelse(is.na(nc[,i]) & nc$any.birthdefect == 0, 0, nc[,i])
}

nc$defect.total <- ifelse(is.na(nc$defect.total) & nc$any.birthdefect == 0, 0, nc$defect.total)

#' Set individual birth defects variables to 0 for kids w/o a birth defect.
for (i in 9:134){
  nc[,i] <- ifelse(is.na(nc[,i]) & nc$any.birthdefect == 0, 0, nc[,i])
}

save(nc, file = './nc.v20171103.2.rdata')



# Try to recover data for missing CHDs ------------------------------------
nc.bd <- haven::read_sas('C:/Users/schraw/Downloads/nc_linked_forbcm.sas7bdat')
nc.bd <- as.data.frame(nc.bd)
save(nc.bd, file = './nc.birth.defect.data.rdata')

tmp <- as.data.frame(nc.bd[nc.bd$DX1 != "",])

tmp$DX1 <- as.numeric(tmp$DX1)

#' Coarctation.
coarc.codes <- c(747100, 747110, 747190)

coarctations <- tmp[tmp$DX1 %in% coarc.codes, ]
coarctations <- filter(tmp, tmp[,21] == 747100 | DX1 == 747110 | DX1 == 747190)
for (i in 22:70){
  new.coarcs <- filter(tmp, tmp[,i] == 747100 | DX1 == 747110 | DX1 == 747190)
  coarctations <- rbind(coarctations, new.coarcs) 
}
coarctations <- as.data.frame(coarctations[!duplicated(coarctations$NCID), ])

coarc.ids <- c(paste0('nc',coarctations$NCID))

#' Pulmonary artery anomalies.
pulm.codes <- c(747300, 747310, 747320, 747325, 747330, 747340, 747380, 747390)

pulm <- tmp[tmp[,21] %in% pulm.codes, ]
for (i in 22:70){
  new.pulm <- tmp[tmp[,i] %in% pulm.codes, ]
  pulm <- rbind(pulm, new.pulm) 
}
pulm <- pulm[!duplicated(pulm$NCID), ]

pulm.ids <- c(paste0('nc',pulm$NCID))

#' Patent ductus arteriosus
pda.codes <- c(747000, 747008, 747010, 747060)

pda <- tmp[tmp[,21] %in% pda.codes, ]

for (i in 22:70){
  new.pda <- tmp[tmp[,i] %in% pda.codes, ]
  pda <- rbind(pda, new.pda) 
}

pda <- pda[!duplicated(pda$NCID), ]

pda.ids <- c(paste0('nc',pda$NCID))

assign.chd <- function(vector.of.ids, var.to.update){ifelse(nc$studyid %in% vector.of.ids, 1, var.to.update)}

nc$coarctationofaorta <- assign.chd(coarc.ids, nc$coarctationofaorta)
nc$pulmonaryarteryanomalies <- assign.chd(pulm.ids, nc$pulmonaryarteryanomalies)
nc$patentductusarteriosis <- assign.chd(pda.ids, nc$patentductusarteriosis)

#' Recalculate lvot.defects with new coarctation cases.
nc$lvot.defects <- ifelse(nc$any.birthdefect == 0, 0,
                          ifelse(nc$any.birthdefect == 1 & (nc$hypoplasticleftheartsyndrome == 1 | nc$interrupted.aortic.arch.type.a.or.c == 1 | nc$coarctationofaorta == 1 | nc$aortic.valve.stenosis == 1), 1, nc$lvot.defects))

rm(coarc.ids, pda.ids, pulm.ids, coarctations, new.coarcs, new.pda, new.pulm, pda, pulm, tmp, coarc.codes, pda.codes, pulm.codes, i, nc.bd)

save(nc, file = './nc.v2017.11.03.3.rdata')



# Compute cancer and cancer1 ----------------------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' I need to compute cancer1 but there is an extra layer of complexity for
#' NC: the age at diagnosis variables indicate that cancers may not be 
#' listed in the appropriate order for children with more than one 
#' diagnosis.  
#' 
#' Fortunately this is a small group, but I'll need to make sure to use
#' the right diagnosis.  In the event of a tie, we agreed to take the DX 
#' marked as having occurred sooner.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
load('./nc.v2017.11.03.3.rdata')
load("./nc.cancer.data.rdata")
load("Z:/Jeremy/GOBACK/Datasets/iccc.codes.to.dx.mappings.v20171018.rdata")

nc$cancer <- ifelse(nc$cancer == 2, 0, nc$cancer)

for (i in 30:34){
  nc.cancer[,i] <- ifelse(nc.cancer[,i] == "", NA, nc.cancer[,i])
  nc.cancer[,i] <- as.numeric(nc.cancer[,i])
}

for (i in 35:39){
  nc.cancer[,i] <- stringr::str_replace(nc.cancer[,i], '^C', "")
  nc.cancer[,i] <- as.numeric(nc.cancer[,i]) 
}

nc.cancer$studyid <- paste0('nc',nc.cancer$ncid)
nc.cancer <- left_join(nc.cancer, select(nc, studyid, num.diagnoses), by = 'studyid')

nc.nocancer <- filter(nc.cancer, cancer == 2)
nc.nocancer$cancer1 <- as.character(NA)

#' 3 'other cancers'.  They don't fit any definition under ICCC WHO recode.  Set to 'other.any'.
nc.singlecancer <- filter(nc.cancer, cancer == 1 & num.diagnoses == 1)
nc.singlecancer$cancer1 <- as.character(NA)
nc.singlecancer$cancer1 <- compute.cancer1(nc.singlecancer$histologic.type.icdo3.1, nc.singlecancer$primary.site.1)
nc.singlecancer$cancer1 <- compute.cancer2(nc.singlecancer$histologic.type.icdo3.1, nc.singlecancer$primary.site.1, newvar = nc.singlecancer$cancer1)
nc.singlecancer$cancer1 <- ifelse(nc.singlecancer$cancer1 == 'other.cancer','other.any', nc.singlecancer$cancer1)

nc.multicancer <- filter(nc.cancer, cancer == 1 & num.diagnoses > 1)

for (i in 10:14){
  nc.multicancer[,i] <- as.numeric(nc.multicancer[,i]) 
  nc.multicancer[,i] <- ifelse(nc.multicancer[,i] == 0, 0.5, nc.multicancer[,i])
}

#' There are more than 3 age at diagnosis columns, but this is the farthest we have to go to find the minimum for any given row.
nc.multicancer$earliest.diagnosis <- 20
nc.multicancer$earliest.diagnosis <- ifelse(!is.na(nc.multicancer$age.at.diagnosis.1) & nc.multicancer$age.at.diagnosis.1 < nc.multicancer$earliest.diagnosis, 10, nc.multicancer$iearliest.diagnosis)
nc.multicancer$earliest.diagnosis <- ifelse(!is.na(nc.multicancer$age.at.diagnosis.2) & nc.multicancer$age.at.diagnosis.2 < nc.multicancer$age.at.diagnosis.1, 11, nc.multicancer$earliest.diagnosis)
nc.multicancer$earliest.diagnosis <- ifelse(!is.na(nc.multicancer$age.at.diagnosis.3) & nc.multicancer$age.at.diagnosis.3 < nc.multicancer$age.at.diagnosis.2, 12, nc.multicancer$earliest.diagnosis)

nc.multicancer$morphology <- ifelse(nc.multicancer$earliest.diagnosis == 10, nc.multicancer$histologic.type.icdo3.1,
                                    ifelse(nc.multicancer$earliest.diagnosis == 11, nc.multicancer$histologic.type.icdo3.2,
                                           ifelse(nc.multicancer$earliest.diagnosis == 12, nc.multicancer$histologic.type.icdo3.3, NA)))

nc.multicancer$site <- ifelse(nc.multicancer$earliest.diagnosis == 10, nc.multicancer$primary.site.1,
                                    ifelse(nc.multicancer$earliest.diagnosis == 11, nc.multicancer$primary.site.2,
                                           ifelse(nc.multicancer$earliest.diagnosis == 12, nc.multicancer$primary.site.3, NA)))

nc.multicancer$cancer1 <- as.character(NA)
nc.multicancer$cancer1 <- compute.cancer1(nc.multicancer$morphology, nc.multicancer$site)
nc.multicancer$cancer1 <- compute.cancer2(nc.multicancer$morphology, nc.multicancer$site, nc.multicancer$cancer1)
nc.multicancer$cancer1 <- ifelse(nc.multicancer$cancer1 == 'other.cancer','other.any', nc.multicancer$cancer1)

print(nc.multicancer[,c(61:64)])

nc.cancer1 <- rbind(select(nc.nocancer, ncid, cancer1),
                    select(nc.singlecancer, ncid, cancer1),
                    select(nc.multicancer, ncid, cancer1))

nc.cancer1$studyid <- paste0('nc',nc.cancer1$ncid)

nc <- left_join(select(nc, -cancer1), select(nc.cancer1, studyid, cancer1), by = 'studyid')

rm(nc.cancer, nc.cancer1, nc.multicancer, nc.nocancer, nc.singlecancer, cancer.codes, i, j)

save(nc, file = './nc.v20171107.1.rdata')



# Compute [cancer].any variables ------------------------------------------
load('./nc.v20171107.1.rdata')

nc$leu.any <- populate.cancer.var(nc$all, nc$aml, nc$leu.other)
nc$lym.any <- populate.cancer.var(nc$hl, nc$nhl, nc$lym.other)
nc$cns.any <- populate.cancer.var(nc$astro, nc$medullo, nc$cns.other)
nc$pns.any <- populate.cancer.var(nc$neuro, nc$pns.other)
nc$renal.any <- populate.cancer.var(nc$nephro, nc$renal.other)
nc$hepatic.any <- populate.cancer.var(nc$hepato, nc$hepatic.other)
nc$bone.any <- populate.cancer.var(nc$osteo, nc$ewing, nc$bone.other)
nc$rms.any <- populate.cancer.var(nc$erms, nc$arms, nc$rms.other)
nc$soft.any <- populate.cancer.var(nc$rms.any, nc$soft.other)
nc$gct.any <- populate.cancer.var(nc$gct.extra, nc$gct.gonad, nc$gct.intra)

save(nc, file = './nc.v20171107.2.rdata')




# Extract standard variables ----------------------------------------------
load('./nc.v20171107.2.rdata')

nc.dem <- get.standard.dem.vars(nc)
nc.def <- get.standard.def.vars(nc)

nc$laterality1 <- as.numeric(NA)
nc.can <- get.standard.cancer.vars(nc)

tmp <- left_join(nc.dem, nc.def, by = 'studyid')
nc <- left_join(tmp, nc.can, by = 'studyid')

nc$runif <- runif(1239652, 0, 1)

save(nc, file = './nc.v20171107.3.rdata')

rm(nc.dem, nc.can, nc.def, tmp)

go <- 
