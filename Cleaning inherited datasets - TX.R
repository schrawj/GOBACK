#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#'                                       
#'                                TEXAS
#'                                         
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------



# Prep environment --------------------------------------------------------
require(dplyr)
require(stringr)

setwd('Z:/Jeremy/GOBACK/Datasets/Texas/')



# User-defined functions --------------------------------------------------
restring.columns <- function(x){
  stringr::str_replace_all(tolower(colnames(x)),'_','.')
}
compute.birth.wt.cat <- function(y){
  ifelse(is.na(y), NA, 
         ifelse(!is.na(y) & y > 2499, 'Normal or High Birthweight',
                ifelse(!is.na(y) & floor(y) %in% 1500:2499, 'Low Birthweight',
                       ifelse(!is.na(y) & floor(y) %in% 400:1500, '400-1499 grams',
                              ifelse(!is.na(y) & y < 400, '<400 grams', '99999')))))
}
rename.tiff.vars <- function(data.frame){
  rename(data.frame, 
         
         m.edu2 = m.edu,
         
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

get.standard.dem.vars <- function(x){
#  require(dplyr)
#  require(stringr)
#  colnames(x) <- str_replace_all(colnames(x),'_','.')
#  colnames(x) <- tolower(colnames(x))
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



# Load in and save raw data -----------------------------------------------
tx.raw <- read.csv('./tx.raw.data.v20171005.csv', header = TRUE, stringsAsFactors = FALSE)
save(tx.raw, file = './tx.raw.data.v20171005.rdata')



# Rename variables --------------------------------------------------------
load("./tx.raw.data.v20171005.rdata")
tx <- tx.raw
colnames(tx) <- restring.columns(tx)

tx <- rename(tx, birth.yr = birth.year, person.yrs = person.years, 
             double.outlet.right.ventricle = dououtletrightvent)
tx <- rename.tiff.vars(tx)



# Demographic variables ---------------------------------------------------
tx$studyid <- paste0(tx$state, tx$X)
tx$birth.wt.cat <- compute.birth.wt.cat(tx$birth.wt)
tx$f.edu2 <- as.numeric(NA)

save(tx, file = './tx.v20171006.1.rdata')



# Paternal ethnicity ------------------------------------------------------
tx.bc.9903 <- readstata13::read.dta13('Z:/Birth Defects-Childhood Cancer Projects/GOBACK Project (Tiffany)/Texas/Work on SAS datasets/TX99-03.dta',
                                      convert.factors = TRUE, convert.underscore = TRUE)
tx.bc.04 <- readstata13::read.dta13('Z:/Birth Defects-Childhood Cancer Projects/GOBACK Project (Tiffany)/Texas/Work on SAS datasets/TX04.dta',
                                    convert.factors = TRUE, convert.underscore = TRUE)
tx.bc.05 <- readstata13::read.dta13('Z:/Birth Defects-Childhood Cancer Projects/GOBACK Project (Tiffany)/Texas/Work on SAS datasets/TX05.dta',
                                    convert.factors = TRUE, convert.underscore = TRUE)
tx.bc.0610 <- readstata13::read.dta13('Z:/Birth Defects-Childhood Cancer Projects/GOBACK Project (Tiffany)/Texas/Work on SAS datasets/TX06-10.dta',
                                    convert.factors = TRUE, convert.underscore = TRUE)
tx.bc.1113 <- readstata13::read.dta13('Z:/Birth Defects-Childhood Cancer Projects/GOBACK Project (Tiffany)/Texas/Work on SAS datasets/TX11-13.dta',
                                      convert.factors = TRUE, convert.underscore = TRUE)

tx.bc.0613 <- rbind(tx.bc.0610, tx.bc.1113)

rm(tx.bc.0610, tx.bc.1113)

print(setdiff(colnames(tx.bc.05), colnames(tx.bc.0613)))

tmp <- c(paste(colnames(tx.bc.0613), ','))

tmp <- select(tx.bc.05,  
              C.DOB ,         birthmonth ,    birthday ,      birthyear ,     C.SEX ,         plurality ,     M.EDU ,         M.HISNOT ,     
              M.HISMEX ,      M.HISPR ,       M.HISCUB ,      M.HISOTH ,      M.HISDES ,      M.RWHITE ,      M.RBLACK ,      M.RAMIND ,     
              M.RINDES ,      M.RASNIN ,      M.RCHINA ,      M.RFILIP ,      M.RJAPAN ,      M.RKOREA ,      M.RVIET ,       M.ROTHAS ,     
              M.RASDES ,      M.RHAWAI ,      M.RGUAM ,       M.RSAMOA ,      M.ROTHPA ,      M.RPACIS ,      M.ROTHER ,      M.ROTDES ,     
              M.RACEETHN ,    M.BR.RAC ,      F.HISNOT ,      F.HISMEX ,      F.HISPR ,       F.HISCUB ,      F.HISOTH ,      F.HISDES ,     
              F.RWHITE ,      F.RBLACK ,      F.RAMIND ,      F.RINDES ,      F.RASNIN ,      F.RCHINA ,      F.RFILIP ,      F.RJAPAN ,     
              F.RKOREA ,      F.RVIET ,       F.ROTHAS ,      F.RASDES ,      F.RHAWAI ,      F.RGUAM ,       F.RSAMOA ,      F.ROTHPA ,     
              F.RPACIS ,      F.ROTHER ,      F.ROTDES ,      F.BR.RAC ,      C.ES.GES ,      F.AGE ,         M.AGE ,         BW.GROUP ,     
              BW.CALGRAMS ,   B.PREG.LENGTH , birthID ,       state ,         datayear )

tx.bc.0513 <- rbind(tmp, tx.bc.0613)

rm(tx.bc.05, tx.bc.0613, tmp)

print(setdiff(colnames(tx.bc.9903), colnames(tx.bc.04)))

table(tx.bc.04$M.HISP) #' Hispanic if == 1
table(tx.bc.04$M.HISPC)
table(tx.bc.04$M.HISP, tx.bc.04$M.HISPC)

tmp <- filter(tx, datayr == '2004')
table(tmp$m.race)

table(tx.bc.9903$M.HISP) #' Hispanic if == 1
table(tx.bc.9903$M.HISPC)
table(tx.bc.9903$M.HISP, tx.bc.9903$M.HISPC)

tx.bc.9904 <- rbind(
                 select(tx.bc.9903, birthID, birthmonth, birthday, birthyear, F.HISP),
                 select(tx.bc.04, birthID, birthmonth, birthday, birthyear, F.HISP))

rm(tx.bc.9903, tx.bc.04, tmp)

tx.bc.0513 <- select(tx.bc.0513, birthID, birthmonth, birthday, birthyear, 
                     M.HISNOT, M.HISMEX, M.HISPR, M.HISCUB, M.HISOTH, M.HISDES, M.RACEETHN,
                     F.HISNOT, F.HISMEX, F.HISPR, F.HISCUB, F.HISOTH, F.HISDES)

tx.bc.0513[1:100,5:11] # Think column 5 indicates non-Hispanic.  6:9 are subcategories of Hispanic ethnicity.

for (i in 5:9){
  tx.bc.0513[,i] <- as.numeric(tx.bc.0513[,i])
}

for (i in 13:16){
  tx.bc.0513[,i] <- as.numeric(tx.bc.0513[,i])
}

tx.bc.0513$M.HISP <- ifelse(rowSums(tx.bc.0513[,6:9]) >= 1, 1, 2)
tx.bc.0513$F.HISP <- ifelse(rowSums(tx.bc.0513[,13:16]) >= 1, 1, 2)

table(tx.bc.0513$M.HISP)
table(tx.bc.0513$M.RACEETHN)
table(tx.bc.0513$F.HISP)

tx.bc.9913 <- rbind(
                    select(tx.bc.9904, birthID, F.HISP),
                    select(tx.bc.0513, birthID, F.HISP))

colnames(tx.bc.9913) <- restring.columns(tx.bc.9913)

save(tx.bc.9913, file = './tx.raw.data.paternal.ethnicity.v20171018.rdata')

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Setting f.race to Hispanic if f.hisp equals Hispanic results in Hispanic
#' becoming the largest category of f.race, and draws subjects from Asian,
#' NHB and AIAN categories into Hispanic category.
#' 
#' Secondly, the number of father's missing race in tx is exactly equal to 
#' the number recorded as Hispanic in 99-04.  Suspicious.
#' 
#' Update f.race in situations where it is missing and father is Hispanic.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
tx <- left_join(tx, tx.bc.9913, by = 'birthid')

table(tx$f.race, useNA = 'always')
tx$f.race <- ifelse(is.na(tx$f.race) & tx$f.hisp == 1, 1, tx$f.race)
table(tx$f.race, useNA = 'always')

save(tx, file = './tx.v20171018.1.rdata')

rm(tx.bc.0513, tx.bc.9904, tx.bc.9913)



# Compute cancer1 ---------------------------------------------------------
load("./tx.v20171018.1.rdata")
load('./tx.cancer1.codes.rdata')
load("Z:/Jeremy/GOBACK/Datasets/iccc.codes.to.dx.mappings.v20171018.rdata")

tx.can <- rename(tx.can, morphology = morph31, site.code = site.code1)

tx.can$cancer1 <- compute.cancer1(tx.can$morphology, tx.can$site.code)
tx.can$cancer1 <- compute.cancer2(tx.can$morphology, tx.can$site.code, tx.can$cancer1)

#' 17 other.cancer instances.
#' I can't reconcile these with the ICCC WHO recode sheet.  
#' Maybe Philip has ideas.
tmp <- arrange(filter(tx.can, cancer1 == 'other.cancer'), morphology)
write.csv(tmp[,1:3], file = 'C:/Users/schraw/Desktop/unmapped.cancer.codes.csv', row.names = FALSE)

#' Set to 'other.any' for now.
tx.can$cancer1 <- ifelse(tx.can$cancer1 == 'other.cancer', 'other.any', tx.can$cancer1)

tx <- select(tx, -cancer1)
tx <- left_join(tx, select(tx.can, birthid, cancer1), by = 'birthid')

save(tx, file = './tx.v20171018.2.rdata')






# Person-years part 1: recover DOB ----------------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' The current person-years variable is not computed correctly.  It just 
#' gives the number of years between birth and end of the study period.
#' 
#' To fix it, I need to recover any available information on DOB, 
#' date of cancer DX and date of death.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
tx.bc.9903 <- readstata13::read.dta13('Z:/Birth Defects-Childhood Cancer Projects/GOBACK Project (Tiffany)/Texas/Work on SAS datasets/TX99-03.dta',
                                      convert.factors = TRUE, convert.underscore = TRUE)
tx.bc.04 <- readstata13::read.dta13('Z:/Birth Defects-Childhood Cancer Projects/GOBACK Project (Tiffany)/Texas/Work on SAS datasets/TX04.dta',
                                    convert.factors = TRUE, convert.underscore = TRUE)
tx.bc.05 <- readstata13::read.dta13('Z:/Birth Defects-Childhood Cancer Projects/GOBACK Project (Tiffany)/Texas/Work on SAS datasets/TX05.dta',
                                    convert.factors = TRUE, convert.underscore = TRUE)
tx.bc.0610 <- readstata13::read.dta13('Z:/Birth Defects-Childhood Cancer Projects/GOBACK Project (Tiffany)/Texas/Work on SAS datasets/TX06-10.dta',
                                      convert.factors = TRUE, convert.underscore = TRUE)
tx.bc.1113 <- readstata13::read.dta13('Z:/Birth Defects-Childhood Cancer Projects/GOBACK Project (Tiffany)/Texas/Work on SAS datasets/TX11-13.dta',
                                      convert.factors = TRUE, convert.underscore = TRUE)

tx.dob.0413 <- rbind(
                      select(tx.bc.04, birthID, C.DOB),
                      select(tx.bc.05, birthID, C.DOB),
                      select(tx.bc.0610, birthID, C.DOB),
                      select(tx.bc.1113, birthID, C.DOB))

tx.bc.9903$C.DOB <- as.character(paste0(tx.bc.9903$birthmonth, tx.bc.9903$birthday, tx.bc.9903$birthyear))

tx.dob.9913 <- rbind(tx.dob.0413, select(tx.bc.9903, birthID, C.DOB))

rm(tx.bc.9903, tx.bc.04, tx.bc.05, tx.bc.0610, tx.bc.1113, tx.dob.0413)

#' Have all DOBs.  Need to convert to date.  Verify all are 8 characters long.
tx.dob.9913$dob.str.len <- stringr::str_length(tx.dob.9913$C.DOB)
table(tx.dob.9913$dob.str.len, useNA = 'always')

tx.dob.9913$C.DOB <- paste0(substr(tx.dob.9913$C.DOB, 1, 2), '/', substr(tx.dob.9913$C.DOB, 3, 4), '/', substr(tx.dob.9913$C.DOB, 5, 8))

tx.dob.9913$C.DOB <- as.Date(tx.dob.9913$C.DOB, format = '%m/%d/%Y')

save(tx.dob.9913, file = './tx.raw.data.birth.dates.v20171018.rdata')



# Person-years part 2: recover date DX ------------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' A few cancer cases are missing the day and/or month of the year in 
#' which they were diagnosed.
#' 
#' For cases with month and year, arbitrarily set day to 15th.
#' This is approximately the middle of the month.
#' 
#' For cases with only year, arbitrarily set month and day to 07/2.
#' This is the middle of the year in a non-leap year.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
tx.datedx <- readstata13::read.dta13('Z:/Birth Defects-Childhood Cancer Projects/GOBACK Project (Tiffany)/Texas/tcr_17_007_july20 (STATA).dta',
                                      convert.factors = TRUE, convert.underscore = TRUE)

tx.datedx$dxdate.str.len <- stringr::str_length(tx.datedx$dxdate)

#' A record of the cases for which we're imputing dates.
tx.datedx.cases.with.imputed.dates <- filter(tx.datedx, dxdate.str.len < 8)
save(tx.datedx.cases.with.imputed.dates, file = './tx.cancer1.cases.with.imputed.dates.rdata')
rm(tx.datedx.cases.with.imputed.dates)

tx.datedx$date.was.imputed <- ifelse(tx.datedx$dxdate.str.len < 8, 1, 0)

tx.datedx$dxdate <- ifelse(tx.datedx$dxdate.str.len == 4, paste0(tx.datedx$dxdate, '0702'),
                           ifelse(tx.datedx$dxdate.str.len == 6, paste0(tx.datedx$dxdate,'15'), tx.datedx$dxdate))

tx.datedx$dxdate.str.len <- stringr::str_length(tx.datedx$dxdate)
table(tx.datedx$dxdate.str.len, useNA = 'always')

tx.datedx$dxdate <- paste0(substr(tx.datedx$dxdate, 1, 4), '/', substr(tx.datedx$dxdate, 5, 6), '/', substr(tx.datedx$dxdate, 7, 8))
tx.datedx$dxdate <- as.Date(tx.datedx$dxdate)

head(tx.datedx$dxdate, 250)

save(tx.datedx, file = './tx.raw.data.date.cancer.dx.v20171018.rdata')



# Person-years part 3: compute person.yrs ---------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Contrary to my initial recollection, DOD is not given in the TX data.
#' 
#' Person.yrs is therefore end of study period 
#' (date last cancer DX is 03/10/2017) minus date of birth if no cancer;
#' date of diagnosis - date of of birth if cancer.
#' 
#' I do not recall if we defined an end of the study period for TX.
#' Tiffany used 2014 as the endpoint in her calculations (I think).
#' 
#' For now I will use 03/10/2017 as end of study period, but I made a note
#' to confirm this at our meeting on 10/23.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
load("./tx.raw.data.date.cancer.dx.v20171018.rdata")
load("./tx.raw.data.birth.dates.v20171018.rdata")
load("./tx.v20171018.2.rdata")

tx <- left_join(tx, select(tx.datedx, birthid, dxdate), by = 'birthid')
tx.dob.9913 <- rename(tx.dob.9913, dob = C.DOB, birthid = birthID)
tx <- left_join(tx, select(tx.dob.9913, birthid, dob), by = 'birthid')

rm(tx.datedx, tx.dob.9913)

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Set a new rule: for kids missing elements of date of DX who were DX'd 
#' during the year of their birth, impute the midpoint of the difference
#' between their DOB and the end of year as the date DX.
#' 
#' Given the small number of kids for whom this is an issue and the fact
#' that I'll probably never have to solve this problem again, I'm just
#' going to fix them manually.
#' 
#' Applies to birthid's 20040119799 and 20130000126.
#' DOBs are 2004-08-16 and 2013-08-10 respectively.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
tx[tx$birthid == 20040119799, "dxdate"] <- as.Date('2004-10-23')
tx[tx$birthid == 20130000126, 'dxdate'] <- as.Date('2013-10-20')

#' Compute person.yrs.
tx$person.yrs <- ifelse(tx$cancer == 0 | (tx$cancer == 1 & tx$dxdate > as.Date('2015-12-31')), (as.Date('2015-12-31') - tx$dob)/365,
                        ifelse(tx$cancer == 1 & tx$dxdate <= as.Date('2015-12-31'), (tx$dxdate - tx$dob)/365, 9999))

save(tx, file = './tx.v20171018.3.rdata')



# Person-years part 4: validate computed variable -------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' This results in a small number of individuals with negative time at 
#' risk.
#' 
#' These may have to be filtered out.
#' 
#' It also seems likely that they might be taken from the population of 
#' kids for whom I imputed date or month values.
#' 
#' I need to go back and get a list of those IDs, then cross-check them
#' with the list of kids with negative person.yrs.
#' 
#' There are 6 kids with negative person.yrs.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
load("./tx.cancer1.cases.with.imputed.dates.rdata")
load("./tx.raw.data.birth.dates.v20171018.rdata")
load("./tx.raw.data.date.cancer.dx.v20171018.rdata")

colnames(tx.dob.9913) <- restring.columns(tx.dob.9913)

tmp <- filter(tx, person.yrs < 0)
tmp <- c(tmp$birthid)

tmp2 <- tx.datedx.cases.with.imputed.dates[tx.datedx.cases.with.imputed.dates$birthid %in% tmp, ]

dates <- left_join(
                      select(tx.datedx, birthid, dxdate),
                      select(tx.dob.9913, birthid, c.dob),
                    by = 'birthid')

dates <- dates[dates$birthid %in% tmp, ]
print(dates)

rm(dates, tmp, tmp2, tx.datedx, tx.datedx.cases.with.imputed.dates, tx.dob.9913)



# Cancer variables --------------------------------------------------------
load('./tx.v20171018.3.rdata')

#' Set num.diagnoses to 0 if no HX of cancer.
table(tx$num.diagnoses, useNA = 'ifany')
tx$num.diagnoses <- ifelse(tx$cancer == 0, 0, tx$num.diagnoses)
table(tx$num.diagnoses, useNA = 'ifany')

save(tx, file = './tx.v20171019.1.rdata')



# Standardize birthid variable --------------------------------------------
load('./tx.v20171019.1.rdata')

tx <- select(tx, -studyid)
tx$studyid <- as.character(paste0('tx',tx$birthid))
tx <- tx[!duplicated(tx$studyid), ]

save(tx, file = './tx.v20171019.2.rdata')

# Extract standard variables ----------------------------------------------
load('./tx.v20171019.2.rdata')

tx.dem <- get.standard.dem.vars(tx)
tx.def <- get.standard.def.vars(tx)
tx.can <- get.standard.cancer.vars(tx)

tx <- left_join(tx.dem, tx.def, by = 'studyid')
tx <- left_join(tx, tx.can, by = 'studyid')

rm(tx.can, tx.def, tx.dem)

save(tx, file = './tx.v20171019.3.rdata')



