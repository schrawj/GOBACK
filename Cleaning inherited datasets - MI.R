#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#'                                       
#'                              MICHIGAN
#' 
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------



# Prep environment --------------------------------------------------------
require(dplyr)
require(stringr)
require(ggplot2)
require(readstata13)

setwd('Z:/Jeremy/GOBACK/Datasets/Old Datasets/Michigan/')

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
compute.birth.wt.cat <- function(y){
  ifelse(is.na(y), NA, 
         ifelse(!is.na(y) & y > 2499, 'Normal or High Birthweight',
                ifelse(!is.na(y) & floor(y) %in% 1500:2499, 'Low Birthweight',
                       ifelse(!is.na(y) & floor(y) %in% 400:1500, '400-1499 grams',
                              ifelse(!is.na(y) & y < 400, '<400 grams', '99999')))))
}
cancerflag <- function(cancer.time.var){
  ifelse(!is.na(cancer.time.var), 1, 0)
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
gen.cancer.var <- function(x){
  x <- 0
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



# Load in and save raw data -----------------------------------------------
mi.raw <-read.csv('./mi.raw.data.v20171005.csv', header = TRUE, stringsAsFactors = FALSE)
save(mi.raw, file = './mi.raw.data.v20171005.rdata')



# Rename variables --------------------------------------------------------
load("./mi.raw.data.v20171005.rdata")
mi <- mi.raw
colnames(mi) <- restring.columns(mi)

#' Remove a duplicated column.
mi <- mi[,-232]
mi <- rename(mi, 
             conganomgenitalandurinary = conganomgenandurinary,
             double.outlet.right.ventricle = dououtletrightvent,
             person.yrs = person.years, 
             cancer = any.cancer)
mi <- rename.tiff.vars(mi)

save(mi, file = './mi.v20171005.1.rdata')



# Demographic variables ---------------------------------------------------
load('./mi.v20171005.1.rdata')

mi$m.edu2 <- ifelse(mi$m.edu2 == 99, NA, mi$m.edu2)
mi$m.edu.yrs <- mi$m.edu2
mi$m.edu2 <- with(mi,
                  ifelse(m.edu.yrs < 12, '< HS',
                         ifelse(m.edu.yrs == 12, 'HS',
                                ifelse(m.edu.yrs > 12, '> HS','This is not a valid category'))))

#' NOTE: birth weight and gestational age are missing for non-cancer children.
#' This is addressed by code near the end of this script.
mi$birth.wt.cat <- compute.birth.wt.cat(mi$birth.wt)

save(mi, file = './mi.v20171005.2.rdata')



# Person-years ------------------------------------------------------------
load('./mi.v20171005.2.rdata')

#' Need valid age.diag1 variable in order to generate person.yrs.
mi$calculated.person.yrs <- mi$yeardiag1 - mi$birth.yr
mi$person.yrs.diff <- abs(mi$calculated.person.yrs - mi$age.diag1)

#' 4 observations that appear to have unreliable age.diag1.
table(mi$person.yrs.diff)
tmp <- filter(mi, person.yrs.diff > 1)
tmp[ ,c(2,12,65,68,70:72,96,97,224,225)]

mi$age.diag1 <- ifelse(mi$person.yrs.diff > 1, mi$cancertime,  mi$age.diag1)

#' Verify problem solved.
list.of.bad.agediag1 <- c(tmp$studyid)
tmp <- mi[mi$studyid %in% list.of.bad.agediag1, ]
tmp[ ,c(2,12,65,68,96,97,224,225,239,240)]
mi$person.yrs.diff <- abs(mi$calculated.person.yrs - mi$age.diag1)
table(mi$person.yrs.diff)
rm(list.of.bad.agediag1)

#' Revise person.yrs.
mi$person.yrs <- ifelse(mi$cancer == 1, mi$age.diag1, 
                        ifelse(mi$cancer == 0, 2011 - mi$birth.yr, 9999))
table(mi$person.yrs, mi$birth.yr)

save(mi, file = './mi.v20171005.3.rdata')

# Compute birth defects variables -----------------------------------------
load('./mi.v20171005.3.rdata')

mi$overall.chrom <- ifelse(mi$chromosomalanomalies == 1, 1, 0)

save(mi, file = './mi.v20171005.4.rdata')



# Load in MI cancer codes -------------------------------------------------

require(readstata13); require(dplyr)

mi.can <- read.dta13('Z:/Birth Defects-Childhood Cancer Projects/GOBACK Project (Tiffany)/Michigan/MI_cancer_062017_rg.dta', 
                     convert.factors = TRUE, convert.underscore = TRUE) 
mi.can <- mi.can[,c(1,52,68)]
mi.can$morph31 <- as.numeric(mi.can$morph31)
mi.can$icdoii1 <- str_replace(mi.can$icdoii1, '^C','')
mi.can$icdoii1 <- as.numeric(mi.can$icdoii1)
mi.can <- rename(mi.can, morphology = morph31, site.code = icdoii1)

save(mi.can, file = './mi.cancer1.codes.rdata')

# Compute num.diagnoses ---------------------------------------------------
load('./mi.v20171005.4.rdata')

#' Need to compute num.diagnoses.  
#' At least one older dataset has info on serial cancer diagnoses.  
tmp <- mi
load('./mi.raw.data.v20170818.rdata')
mi.old <- mi
mi <- tmp
rm(tmp)

mi <- left_join(mi, select(mi.old, studyid, cancertime2, cancertime3, cancertime4), by = 'studyid')

mi$cancerflag1 <- cancerflag(mi$cancertime1)
mi$cancerflag2 <- cancerflag(mi$cancertime2)
mi$cancerflag3 <- cancerflag(mi$cancertime3)
mi$cancerflag4 <- cancerflag(mi$cancertime4)

mi$num.diagnoses <- rowSums(mi[,246:249], na.rm = TRUE)
mi <- mi[ ,-c(246:249)]

save(mi, file = './mi.v20171005.5.rdata')



# Compute cancer1 ---------------------------------------------------------
load('./mi.v20171005.5.rdata')
load('./mi.cancer1.codes.rdata')

mi <- left_join(mi, mi.can, by = 'studyid')

load("Z:/Jeremy/GOBACK/Datasets/iccc.codes.to.dx.mappings.v20171018.rdata")

mi$cancer1 <- compute.cancer1(mi$morphology, mi$site.code)
mi$cancer1 <- compute.cancer2(mi$morphology, mi$site.code, mi$cancer1)

#' There are 5 cancers with odd combinations of site and morphology codes.
#' It's unclear what these might be.
#' Move them to 'other.any'.
tmp <- filter(mi, cancer1 == 'other.cancer')
table(tmp$any.birthdefect, useNA = 'always')
mi$cancer1 <- ifelse(mi$cancer1 == 'other.cancer','other.any',mi$cancer1)

save(mi, file = './mi.v20171006.1.rdata')



# Compute individual cancers ----------------------------------------------
load("./mi.v20171006.1.rdata")

l <- c(colnames(mi[20:59]))
mi <- mi[,-c(20:59)]

for (i in l){
  x <- (paste0('mi$',i))
  mi[i] <- gen.cancer.var(x)
}

mi$cancer <- ifelse(!is.na(mi$cancer1), 1, 0)

tmp <- filter(mi, cancer == 0)
tmp2 <- filter(mi, cancer == 1)

for (i in 209:248){
  can.name <- (names(tmp2)[i])
  for (j in tmp2){
    tmp2[,i] <- ifelse(tmp2$cancer1 == can.name, 1, 0)
  }
}

mi <- arrange(rbind(tmp, tmp2), studyid)



# Compute [cancer].any variables ------------------------------------------
mi$leu.any <- populate.cancer.var(mi$all, mi$aml, mi$leu.other)
mi$lym.any <- populate.cancer.var(mi$hl, mi$nhl, mi$lym.other)
mi$cns.any <- populate.cancer.var(mi$astro, mi$medullo, mi$cns.other)
mi$pns.any <- populate.cancer.var(mi$neuro, mi$pns.other)
mi$renal.any <- populate.cancer.var(mi$nephro, mi$renal.other)
mi$hepatic.any <- populate.cancer.var(mi$hepato, mi$hepatic.other)
mi$bone.any <- populate.cancer.var(mi$osteo, mi$ewing, mi$bone.other)
mi$rms.any <- populate.cancer.var(mi$erms, mi$arms, mi$rms.other)
mi$soft.any <- populate.cancer.var(mi$rms.any, mi$soft.other)
mi$gct.any <- populate.cancer.var(mi$gct.extra, mi$gct.gonad, mi$gct.intra)

save(mi, file = './mi.v20171006.2.rdata')




# Inspect cancer variables ------------------------------------------------

#' Tiffany should have re-coded some gct.gonad cases as other.any.
#' Compare frequencies of these diagnoses between the two sets.
#' Her numbers don't make sense though.
table(mi.old$any_cancer)
table(mi.old$gct_gonad)
table(mi.old$other_any)

table(mi$cancer)
table(mi$gct.gonad)
table(mi$other.any)
table(mi$cancer1)

table(mi$cancer1)

#' 5 'other cancers'.
tmp <- filter(mi, cancer1 == 'other.cancer')
tmp <- arrange(tmp, morphology)
print(tmp[,c(2,247,248)])

for (i in 22:59){
  print(table(mi[,i]))
}

for (i in 209:248){
  print(names(mi)[i])
  print(table(mi[i]))
}

tmp <- filter(mi, gct.gonad == 1)
print(unique(tmp$morphology))

# Recover missing birthweight data ----------------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Due to some sort of error when Tiffany was joining the raw files, 
#' birthweight data was lost for all children without cancer.  She says she
#' fixed this on 10/12/2017.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

mi.new.bw <- read.csv(file = 'C:/Users/schraw/Desktop/mi.updated.birthweights.csv', header = TRUE, stringsAsFactors = FALSE)
colnames(mi.new.bw) <- restring.columns(mi.new.bw)
mi.new.bw <- mi.new.bw[,-230]
save(mi.new.bw, file = './mi.raw.data.new.birthweights.v20171017.rdata')



# Update birthweight and gest age variables -------------------------------
load("Z:/Jeremy/GOBACK/Datasets/Michigan/mi.v20171006.2.rdata")
load("Z:/Jeremy/GOBACK/Datasets/Michigan/mi.raw.data.new.birthweights.v20171017.rdata")

tmp <- select(mi.new.bw, studyid, birth.wt, gest.age)
tmp$birth.wt <- ifelse(tmp$birth.wt == 9999, NA, tmp$birth.wt)

mi <- select(mi, -birth.wt, -gest.age)
mi <- left_join(mi, tmp, by = 'studyid')

#' birth.wt.cat must now be recalculated.
mi$birth.wt.cat <- compute.birth.wt.cat(mi$birth.wt)

rm(tmp, mi.new.bw)

save(mi, file = './mi.v20171006.3.rdata')



# Extract standard variables ----------------------------------------------
load("./mi.v20171006.3.rdata")

mi.dem <- get.standard.dem.vars(mi)
mi.can <- get.standard.cancer.vars(mi)
mi.def <- get.standard.def.vars(mi)

mi <- left_join(mi.dem, mi.def, by = 'studyid')
mi <- left_join(mi, mi.can, by = 'studyid')

rm(mi.can, mi.def, mi.dem)

mi$studyid <- paste0('mi',mi$studyid)

save(mi, file = './mi.v20171006.4.rdata')

# Format the updated (May 2022) cancer registry data ----------------------

require(tidyverse); require(magrittr); require(haven)

setwd('//smb-main.ad.bcm.edu/GenEpi3/TiffanyChambers/GOBACK Data/Michigan Data 2021/Updated Data (June 2021)/Stata/')

mi <- read_dta('all years (harmonized).dta')

saveRDS(mi,
        '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old_Datasets/Michigan/mi.raw.data.v20220511.rds')

mi <- readRDS('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old_Datasets/Michigan/mi.raw.data.v20220511.rds')

#' Compute maternal race/ethnicity.
mi <- mi %>% 
  select(-M_RACE) %>% 
  filter(BIRTH_YR > 1991) %>% 
  mutate(race.all.fields = paste0(Birth_MomRace, Birth_MomRace2b2006, Birth_MomRace3b2006, Birth_MomRace4b2006, Birth_MomRace5b2006) %>% 
           str_remove_all('[\\Q.\\E09*]') %>% # remove values that do not correspond to valid categories per emails from MI registry staff
           str_trim(side = 'both') %>% # remove white spaces
           str_replace_all('[5678AGKOSV]', '4') %>% # consolidate different categories for Asians
           gsub('([[:digit:]])\\1+', '\\1', .), # 11 or 111 etc. will be converted to 1.
         m.race = ifelse(race.all.fields == '1', 2, 
                         ifelse(race.all.fields == '2', 3, 
                                ifelse(race.all.fields == '4', 4, 
                                       ifelse(race.all.fields == '', 7, 
                                              ifelse(race.all.fields == '3', 6, 5))))),
         Birth_MomHispb2006 = ifelse(is.na(Birth_MomHispb2006), 9, Birth_MomHispb2006))

mi <- mi %>% 
  mutate(m.race = ifelse(Birth_MomHispb2006 == 1, 1, m.race))

#' Some cancer cases were diagnosed after the age of 20 years.
#' They should be considered unaffected.
#mi %<>% mutate(cancer = ifelse(Cancer_AgeDiag1 > 18, 0, cancer))
mi %<>% mutate(cancer = ifelse(ageatdx > 216, 0, cancer))

#' Cancer should be set to zero for unaffected children.
#' Age at cancer diagnosis should be set to missing in unaffected children.
mi %<>% mutate(cancer = ifelse(is.na(cancer), 0, cancer),
               ageatdx = ifelse(cancer == 0, NA, ageatdx))

#' Strip tumor data for cancers diagnosed after 18 years of age.
mi %<>% 
  mutate(Cancer_Histology1 = ifelse(Cancer_AgeDiag1 > 18, NA, Cancer_Histology1),
         Cancer_BehaviorICDO31 = ifelse(Cancer_AgeDiag1 > 18, NA, Cancer_BehaviorICDO31),
         Cancer_PrimarySite1 = ifelse(Cancer_AgeDiag1 > 18, NA, Cancer_PrimarySite1),
         Cancer_Laterality1 = ifelse(Cancer_AgeDiag1 > 18, NA, Cancer_Laterality1),
         Cancer_SEERSumStg20001 = ifelse(Cancer_AgeDiag1 > 18, NA, Cancer_SEERSumStg20001),
         
         Cancer_Histology2 = ifelse(Cancer_AgeDiag2 > 18, NA, Cancer_Histology2),
         Cancer_BehaviorICDO32 = ifelse(Cancer_AgeDiag2 > 18, NA, Cancer_BehaviorICDO32),
         Cancer_PrimarySite2 = ifelse(Cancer_AgeDiag2 > 18, NA, Cancer_PrimarySite2),
         Cancer_Laterality2 = ifelse(Cancer_AgeDiag2 > 18, NA, Cancer_Laterality2),
         Cancer_SEERSumStg20002 = ifelse(Cancer_AgeDiag2 > 18, NA, Cancer_SEERSumStg20002),
         
         Cancer_Histology3 = ifelse(Cancer_AgeDiag3 > 18, NA, Cancer_Histology3),
         Cancer_BehaviorICDO33 = ifelse(Cancer_AgeDiag3 > 18, NA, Cancer_BehaviorICDO33),
         Cancer_PrimarySite3 = ifelse(Cancer_AgeDiag3 > 18, NA, Cancer_PrimarySite3),
         Cancer_Laterality3 = ifelse(Cancer_AgeDiag3 > 18, NA, Cancer_Laterality3),
         Cancer_SEERSumStg20003 = ifelse(Cancer_AgeDiag3 > 18, NA, Cancer_SEERSumStg20003),
         
         Cancer_Histology4 = ifelse(Cancer_AgeDiag4 > 18, NA, Cancer_Histology4),
         Cancer_BehaviorICDO34 = ifelse(Cancer_AgeDiag4 > 18, NA, Cancer_BehaviorICDO34),
         Cancer_PrimarySite4 = ifelse(Cancer_AgeDiag4 > 18, NA, Cancer_PrimarySite4),
         Cancer_Laterality4 = ifelse(Cancer_AgeDiag4 > 18, NA, Cancer_Laterality4),
         Cancer_SEERSumStg20004 = ifelse(Cancer_AgeDiag4 > 18, NA, Cancer_SEERSumStg20004),
         
         Cancer_Histology5 = ifelse(Cancer_AgeDiag5 > 18, NA, Cancer_Histology5),
         Cancer_BehaviorICDO35 = ifelse(Cancer_AgeDiag5 > 18, NA, Cancer_BehaviorICDO35),
         Cancer_PrimarySite5 = ifelse(Cancer_AgeDiag5 > 18, NA, Cancer_PrimarySite5),
         Cancer_Laterality5 = ifelse(Cancer_AgeDiag5 > 18, NA, Cancer_Laterality5),
         Cancer_SEERSumStg20005 = ifelse(Cancer_AgeDiag5 > 18, NA, Cancer_SEERSumStg20005),
         
         Cancer_Histology6 = ifelse(Cancer_AgeDiag6 > 18, NA, Cancer_Histology6),
         Cancer_BehaviorICDO36 = ifelse(Cancer_AgeDiag6 > 18, NA, Cancer_BehaviorICDO36),
         Cancer_PrimarySite6 = ifelse(Cancer_AgeDiag6 > 18, NA, Cancer_PrimarySite6),
         Cancer_Laterality6 = ifelse(Cancer_AgeDiag6 > 18, NA, Cancer_Laterality6),
         Cancer_SEERSumStg20006 = ifelse(Cancer_AgeDiag6 > 18, NA, Cancer_SEERSumStg20006)
         )

mi %<>% mutate(person.yrs = ifelse(cancer == 1, ageatdx/12, 
                            ifelse(cancer == 0 & timeatrisk_nocancer <= 216, timeatrisk_nocancer/12,
                            ifelse(timeatrisk_nocancer > 216, 18, NA))))

mi.cancer <- mi %>% 
  filter(cancer == 1) %>% 
  select(STUDYID, 
         starts_with('Cancer_Histology'), 
         starts_with('Cancer_Behavior'), 
         starts_with('Cancer_PrimarySite'), 
         starts_with('Cancer_Laterality'),
         starts_with('Cancer_SEERSumStg2000'))

saveRDS(mi,
        '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old_Datasets/Michigan/mi.v20220516.rds')

saveRDS(mi.cancer,
        '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded_datasets/mi.cancer.codes.v20220516.rds')

rm(list = ls()); gc()
