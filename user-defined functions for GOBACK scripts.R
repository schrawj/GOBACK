#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' This is a central list of all the user-defined functions I wrote during
#' data cleaning for GOBACK, as well as a place to hold a TO-DO list for unresolved 
#' data cleaning concerns.
#' 
#' TODO: Handle missing paternal variables.  
#' TODO: Update get.standard.vars function to fetch variables in the updated final data schema.
#' TODO: Include canceryr1 variable?  
#' TODO: Use the compute.hem.dx function and hem.codes list in the AR data?  Rename as 
#' compute.cancer.dx or something more apropos.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
pull.hist <- function(x){
  select(x, histology.behavior)
}

check.var.levels <- function(x, y){
  tmp <- x[, colnames(x) %in% y]
  apply(tmp[,-1], 2, function(x){
    str(x)
    unique(x)
  })
}

#' Takes argument X: a data frame, and converts its column names to lower case while replacing underscores with periods.
restring.columns <- function(x){
  stringr::str_replace_all(tolower(colnames(x)),'_','.')
}

#' Given a data frame with Tiffany's typical variable names, converts their names to those specified 
#' in the data dictionary.
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

compute.m.edu2 <- function(x, y){
  factor(
    ifelse(x < 12, 1,
           ifelse(x == 12, 2,
                  ifelse(x %in% 12:98, 3,
                         ifelse(x == y, 7, NA)))),
    levels = c(1:3,7),
    labels = c('< HS', 'HS', '> HS', 'Uknown'))
}

#' Computes the birth defects groupings variables.  Takes arguemnt X: a data frame
#' and Y: a vector of column indices to sum across.
compute.defect <- function(x, y){
  ifelse(rowSums(x[,y]) >= 1, 1, 0)
}

#' Computes birth weight category according to the ranges provided in the August AR data.
#' Takes argument Y: joint data frame-column name for numeric BW data.
compute.birth.wt.cat <- function(y){
  ifelse(is.na(y), NA, 
         ifelse(!is.na(y) & y > 2499, 'Normal or High Birthweight',
                ifelse(!is.na(y) & floor(y) %in% 1500:2499, 'Low Birthweight',
                       ifelse(!is.na(y) & floor(y) %in% 400:1500, '400-1499 grams',
                              ifelse(!is.na(y) & y < 400, '<400 grams', '99999')))))
}

#' Takes argument X: a variable in a data frame which should have certain values 
#' recoded as NA, and Y: a vector of numeric values which should be recoded as NA.
set.to.na <- function(x, y){
  ifelse(x %in% y, NA, x)
}

#' A trio of functions for abstracting and arranging studyid + the standard demographic, birth defect and cancer variables from 
#' a data frame containing them.
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

#' Searches through a vector of ICD-O-3 morph/behavior codes and updates the master list.
update.dx <- function(x, y){
  ifelse(is.na(site.morph.codes$cancer1) & site.morph.codes$morphology %in% x, y, site.morph.codes$cancer1)
}

#' Used to initialize empty cancer variables.
gen.cancer.var <- function(x){
  x <- 0
}

#' A function to compute the [cancertype].any variables.
#' Requires two arguments, X & Y and can accept a third, Z.  All 3 should be 
#' columns in a data frame, which, if equal to 1, indicate the variable being 
#' computed should also be set to 1.
populate.cancer.var <- function(x,y,z){
  if(missing(z)){
    ifelse((x | y) == 1, 1, 0)
  }
  else{
    ifelse((x | y | z) == 1, 1, 0)
  }
}

#' As above, but takes up to 4 arguments.
populate.chd.var <- function(w,x,y,z){
  if(missing(z) & missing(w)){
    ifelse((x | y) == 1, 1, 0)
  }
  if(missing(z)){
    ifelse((w | x | y) == 1, 1, 0)
  }
  else{
    ifelse((x | y | z | w) == 1, 1, 0)
  }
}

#' Function used to compute num.diagnoses in MI data.
#' Create a variable == 1 if information on the ith cancer diagnosis is not missing.
cancerflag <- function(cancer.time.var){
  ifelse(!is.na(cancer.time.var), 1, 0)
}

#' A two-part function to compute cancer diagnosis.
#' Had to be split into two parts because there are too many ifelse statements in my shitty code to run it in a single call.
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

#' Used to recalculate several incorrectly generated birth defects variables in the TX and NC datasets.
get.affected.ids <- function(codes, data.frame,first.defect.col.index, other.defect.col.indices, id.var.index){
  
  tmp2 <- data.frame[first.defect.col.index %in% codes, ]
  
  for (i in other.defect.col.indices){
    tmp3 <- filter(data.frame, i %in% codes)
    tmp2 <- rbind(tmp2, tmp3) 
  }
  
  tmp2 <- tmp2[!duplicated(tmp2[,id.var.index]), ]
  
  affected.ids <- c(paste0('tx', paste0(tmp2[,id.var.index])))
  
  rm(tmp2, tmp3)
  
  return(affected.ids)
}

sort.codes <- function(dataframe){
  dataframe <- data.frame(codes = dataframe)
  dataframe <- arrange(dataframe, codes)
  dataframe <- c(dataframe$codes)
}

compute.other.var <- function(var1, var2){
  ifelse(var1 == 1 | var2 == 1, 1,
         ifelse(goback$any.birthdefect == 1 & is.na(var1) & is.na(var2), NA, 0))
}



# Miscellaneous vectors for these functions -------------------------------


birth.cert.vars <- c('studyid', 'state','sex','m.race','m.edu','m.edu2','m.age','birth.wt',
                     'gest.age','birth.yr','plu','f.race','f.edu','f.edu2','f.age')

birth.cert.char.vars <- c('studyid','state')

birth.cert.numeric.vars <- c('sex','m.race','m.edu','m.edu2','m.age','birth.wt','gest.age',
                             'birth.yr','plu','f.race','f.edu','f.edu2','f.age')

newcol <- c('conganomalies.cns.other.major','conganomalies.cns.other.minor',
            'conganomalies.eye.other.major','conganomalies.eye.other.minor',
            'conganomalies.ear.face.neck.other.major', 'conganomalies.ear.face.neck.other.minor',
            'heart.circsys.other.major', 'heart.circsys.other.minor',
            'respsys.other.major', 'respsys.other.minor',
            'cleftpalateandlip.other.major','cleftpalateandlip.other.minor',
            'digestivesystem.other.major', 'digestivesystem.other.minor',
            'genitalandurinary.other.major','genitalandurinary.other.minor',
            'musculoskelsys.other.major','musculoskelsys.other.minor',
            'chromosomalanomalies.other.major','chromosomalanomalies.other.minor')

#' For loading in a vector of all the standard birth defects variables.
load("Z:/GOBACK/Jeremy/Datasets/list.of.standard.birth.defects.variables.rdata")




# Cancer codes ------------------------------------------------------------

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Vectors for the site and morphology codes that define each of our cancers.
#' 
#' Used to create the list object cancer.codes, which is required for the functions 
#' compute.cancer1 and compute.cancer2.
#' 
#' Updated on 10.05.2017 to include gonadal carcinomas and other and unspecified gonadal
#' tumors codes in the 'other.any' variable.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
all.non.contingent.codes <- c(9831,9834,9948,9835,9836,9826,9832:9833,9940,9820) 
all.contingent.codes <- c(9811:9818,9837,9823,9827)
all.contingent.codes.sites <- c(420,421,444) 
aml.codes <- c(9840,9861,9865:9867,9869:9874,9891,9895:9898,9910,9911,9920,9931)
hl.codes <- c(9650:9655,9659,9661:9665,9667)
nhl.non.contingent.codes <- c(9727:9729,9597,9670:9671,9673,9675,9678:9680,9684,9688:9691,9695,9698:9699,9712,9731:9735,9737:9738,9761:9762,9764:9766,9769,9970,9700:9702,9705,9708:9709,9714,9716:9719,9724:9726,9767:9768,9591,9760)
nhl.contingent.codes <- c(9811:9818,9837,9823,9827)
nhl.contingent.codes.sites <- c(0:419,422:423,425:809)
other.hem.codes <-  c(9863,9875,9876,9950,9960:9964,9945,9946,9975,9980,9982:9987,9989,9991:9992,9800,9801,9805:9809,9860,9930,9965:9967,9971)
lym.other.codes <- c(9687,9740:9742,9750:9759,9590,9596)
ependymoma.codes <- c(9383,9391:9394)
astro.non.contingent.codes <- c(9384,9400:9411,9420:9424,9440:9442)
astro.contingent.codes <- 9380
astro.contingent.codes.sites <- 723
medullo.codes <- c(9470:9472,9474,9480)
pnet.codes <- 9473
cns.other.non.contingent.codes <- c(9390,9508,9450,9451,9460,9382,9381,9430,9444,8270:8281,8300,9350:9352,9582,9360:9362,9412:9413,
                                    9492:9493,9505:9507,9530:9539)
cns.other.contingent.codes1 <- 9501:9504
cns.other.contingent.codes1.sites <- 700:729
cns.other.contingent.codes2 <- 9380
cns.other.contingent.codes2.sites <- c(700:722,724:729,751,753)
cns.other.contingent.codes3 <- 8000:8005
cns.other.contingent.codes3.sites <- c(700:729,751:753)
neuroblast.codes <- c(9490,9500)
pns.other.non.contingent.codes <- c(8680:8683,8690:8693,8700,9520:9523)
pns.other.contingent.codes <- 9501:9504
pns.other.contingent.codes.sites <- c(000:699,739:768,809)
retino.codes <- 9510:9514
nephro.codes <- 8959:8960
renal.other.non.contingent.codes <- c(8964:8967,8311:8312,8316:8319,8361)
renal.other.contingent.codes <- c(8963,9364,8010:8041,8050:8075,8082,8120:8122,8130:8141,8143,8155,8190:8201,8210:8211,8221:8231,8240:8241,8244:8246,8260:8263,8290,8310,8320,8323,8401,8430,8440,8480:8490,8504,8510,8550,8560:8576,
                                  8000:8005)
renal.other.contingent.codes.sites <- 649
hepatoblast.codes <- 8970
hepatic.other.non.contingent.codes <- 8160:8180
hepatic.other.contingent.codes <- c(8010:8041,8050:8075,8120:8122,8140:8141,8143,8155,8190:8201,8210:8211,8230,8231,8240:8241,8244:8246,8260:8264,8310,8320,8323,8401,8430,8440,8480:8490,8504,8510,8550,8560:8576,8000:8005)
hepatic.other.contingent.codes.sites <- 220:221
osteo.codes <- c(9180:9187,9191:9195,9200)
osteo.codes.sites <- c(400:419,760:768,809)
ewing.contingent.codes1 <- c(9260,9365)
ewing.contingent.codes1.sites <- c(400:419)
ewing.contingent.codes2 <- 9260
ewing.contingent.codes2.sites <- c(760:768,809)
bone.other.non.contingent.codes <- c(9221,9230,9241:9243,8812,9250,9261,9262,9370:9372,9270:9275,9280:9282,9290,9300:9302,9310:9312,9320:9322,9330,9340:9342)
bone.other.contingent.codes1 <- c(9210,9220,9240)
bone.other.contingent.codes1.sites <- c(400:419,760:768,809)
bone.other.contingent.codes2 <- c(9363:9364,8810:8811,8823,8830,8000:8005,8800:8801,8803:8805)
bone.other.contingent.codes2.sites <- 400:419
erms.codes <- 8910
arms.codes <- 8920
rms.other.codes <- c(8900:8905,8912,8991)
soft.other.non.contingent.codes <- c(8820,8822,8824:8827,9150,9160,9540:9571,9491,9580,9140,8850:8858,8860:8862,8870,8880:8881,8831:8833,8836,9251:9252,8890:8898,9040:9044,9120:9125,9130:9133,9135:9136,9141:9142,
                                     9161,9170:9175,9231,9581,8587,8710:8713,8806,8840:8842,8921,8982,8990,9373) 
soft.other.contingent.codes1.sites <- c(0:399,440:768,809)
soft.other.contingent.codes1 <- c(8810:8811,8813:8815,8821,8823,8834:8835)
soft.other.contingent.codes2 <- 9260
soft.other.contingent.codes2.sites <- c(0:399,470:759)
soft.other.contingent.codes3.sites <- c(0:399,470:639,659:768,809)
soft.other.contingent.codes3 <- 9365
soft.other.contingent.codes4 <- 9364
soft.other.contingent.codes4.sites <- c(0:399,470:639,659:699,739:768,809)
soft.other.contingent.codes5 <- 8963
soft.other.contingent.codes5.sites <- c(0:399,470:639,659:699,739:768,809)
soft.other.contingent.codes6 <- c(8800:8805,8830)
soft.other.contingent.codes6.sites <- c(0:399,440:768,809)
soft.other.contingent.codes7 <- c(9180,9210,9220,9240)
soft.other.contingent.codes7.sites <- 490:499
soft.other.contingent.codes8 <- 9240
soft.other.contingent.codes8.sites <- 720
intra.gct.codes <- c(9060:9065,9080:9084,9070,9072,9071,9100,9085,9101)
intra.gct.codes.sites <- c(700:729,751:753)
extra.gct.codes <- c(9060:9065,9080:9084,9070,9072,9071,9100,9103:9104,9085,9101,9102,9105)
extra.gct.codes.sites <- c(0:559,570:619,630:699,739:750,754:768,809)
gonad.gct.codes <- c(9060:9065,9080:9084,9090:9091,9070,9072,9071,9100,9085,9101,9073)
gonad.gct.codes.sites <- c(569,620:629)
other.unspec.non.contingent.codes <- c(8936,8971:8973,8930:8935,8950,8951,8974,8981,9050:9055,9110,8441:8444,8450,8451,8460:8473,8590:8671)
other.unspec.contingent.codes1 <- 9363
other.unspec.contingent.codes1.sites <- c(0:399,470:759)
other.unspec.contingent.codes2 <- 8000:8005
other.unspec.contingent.codes2.sites <- c(0:218,239:399,420:559,570:619,630:639,659:699,739:750,754:809)
other.unspec.contingent.codes3 <- c(8010:8041,8050:8075,8082,8120:8122,8130:8141,8143,8190:8201,8210,8211,8221:8241,8244:8246,8260:8263,8290,8310,8313,8320,8323,8380:8384,
                                    8430,8440,8480:8490,8504,8510,8550,8560:8573,9000,9014,9015,8000:8005) 
other.unspec.contingent.codes3.sites <- c(569,620:629)
other.unspec.contingent.codes4 <- 9080
other.unspec.contingent.codes4.sites <- 770
epithe.non.contingent.codes <- c(8370:8375,8330:8337,8340:8347,8350,8720:8780,8790)
epithe.contingent.codes1 <- c(8010:8084,8120:8157,8190:8264,8290,8310,8313:8315,8320:8325,8360,8380:8384,8430:8440,8452:8454,8480:8586,8588:8589,8940:8941,8983,9000,9010:9016,9020,9030)
epithe.contingent.codes1.sites <- c(79:89,180,182:189,199,209,210,218,181,340:349,379,500:509,530:539,670:679,690:699,760:768,809,0:69,90:109,129:179,239:339,480:488,510:529,540:549,559,570:619,630:639,
                                    659:669,680:689,700:729,750:759)
epithe.contingent.codes2 <- c(8010:8041,8050:8075,8310,8320,8323,8430)
epithe.contingent.codes2.sites <- c(739,110:119,440:449)
epithe.contingent.codes3 <- c(8082,8120:8122,8130:8141,8190,8200:8201,8211,8230,8231,8244:8246,8260:8263,8290,8440,8480:8481,8510,8560:8573)
epithe.contingent.codes3.sites <- 739
epithe.contingent.codes4 <- c(8082:8083,8120:8122,8130:8141,8190,8200:8201,8211,8230:8231,8244:8246,8260:8263,8290,8440,8480:8481,8500:8576)
epithe.contingent.codes4.sites <- 110:119
epithe.contingent.codes5 <- c(8078,8082,8090:8110,8140,8143,8147,8190,8200,8240,8246:8247,8260,8390:8420,8480,8542,8560,8570:8573,8940,8941)
epithe.contingent.codes5.sites <- 440:449
epithe.contingent.codes6 <- 8071
epithe.contingent.codes6.sites <- 411

cancer.codes <- list(all.non.contingent.codes = all.non.contingent.codes, 
                  all.contingent.codes = all.contingent.codes,
                  all.contingent.codes.sites = all.contingent.codes.sites,
                  aml.codes = aml.codes,
                  other.hem.codes = other.hem.codes,
                  hl.codes = hl.codes,
                  nhl.non.contingent.codes = nhl.non.contingent.codes,
                  nhl.contingent.codes = nhl.contingent.codes,
                  nhl.contingent.codes.sites = nhl.contingent.codes.sites,
                  lym.other.codes = lym.other.codes,
                  ependymoma.codes = ependymoma.codes,
                  astro.non.contingent.codes = astro.non.contingent.codes,
                  astro.contingent.codes = astro.contingent.codes,
                  astro.contingent.codes.sites = astro.contingent.codes.sites,
                  medullo.codes = medullo.codes,
                  pnet.codes = pnet.codes,
                  cns.other.non.contingent.codes = cns.other.non.contingent.codes,
                  cns.other.contingent.codes1 = cns.other.contingent.codes1,
                  cns.other.contingent.codes1.sites = cns.other.contingent.codes1.sites,
                  cns.other.contingent.codes2 = cns.other.contingent.codes2,
                  cns.other.contingent.codes2.sites = cns.other.contingent.codes2.sites,
                  cns.other.contingent.codes3 = cns.other.contingent.codes3,
                  cns.other.contingent.codes3.sites = cns.other.contingent.codes3.sites,
                  neuroblast.codes = neuroblast.codes,
                  pns.other.non.contingent.codes = pns.other.non.contingent.codes,
                  pns.other.contingent.codes = pns.other.contingent.codes,
                  pns.other.contingent.codes.sites = pns.other.contingent.codes.sites,
                  retino.codes = retino.codes,
                  nephro.codes = nephro.codes,
                  renal.other.non.contingent.codes = renal.other.non.contingent.codes,
                  renal.other.contingent.codes = renal.other.contingent.codes,
                  renal.other.contingent.codes.sites = renal.other.contingent.codes.sites,
                  hepatoblast.codes = hepatoblast.codes,
                  hepatic.other.non.contingent.codes = hepatic.other.non.contingent.codes,
                  hepatic.other.contingent.codes = hepatic.other.contingent.codes,
                  hepatic.other.contingent.codes.sites = hepatic.other.contingent.codes.sites,
                  osteo.codes = osteo.codes, 
                  osteo.codes.sites = osteo.codes.sites,
                  ewing.contingent.codes1 = ewing.contingent.codes1, 
                  ewing.contingent.codes1.sites = ewing.contingent.codes1.sites, 
                  ewing.contingent.codes2 = ewing.contingent.codes2, 
                  ewing.contingent.codes2.sites = ewing.contingent.codes2.sites,
                  bone.other.non.contingent.codes = bone.other.non.contingent.codes, 
                  bone.other.contingent.codes1 = bone.other.contingent.codes1, 
                  bone.other.contingent.codes1.sites =  bone.other.contingent.codes1.sites, 
                  bone.other.contingent.codes2 = bone.other.contingent.codes2,
                  bone.other.contingent.codes2.sites = bone.other.contingent.codes2.sites,
                  erms.codes = erms.codes,
                  arms.codes = arms.codes,
                  rms.other.codes = rms.other.codes,
                  soft.other.non.contingent.codes = soft.other.non.contingent.codes, 
                  soft.other.contingent.codes1 = soft.other.contingent.codes1, soft.other.contingent.codes1.sites = soft.other.contingent.codes1.sites, 
                  soft.other.contingent.codes2 = soft.other.contingent.codes2, soft.other.contingent.codes2.sites = soft.other.contingent.codes2.sites,
                  soft.other.contingent.codes3 = soft.other.contingent.codes3, soft.other.contingent.codes3.sites = soft.other.contingent.codes3.sites, soft.other.contingent.codes4 = soft.other.contingent.codes4, soft.other.contingent.codes4.sites = soft.other.contingent.codes4.sites,
                  soft.other.contingent.codes5 = soft.other.contingent.codes5, soft.other.contingent.codes5.sites = soft.other.contingent.codes5.sites, soft.other.contingent.codes6 = soft.other.contingent.codes6, soft.other.contingent.codes6.sites = soft.other.contingent.codes6.sites,
                  soft.other.contingent.codes7 = soft.other.contingent.codes7, soft.other.contingent.codes7.sites = soft.other.contingent.codes7.sites,
                  soft.other.contingent.codes8 = soft.other.contingent.codes8, soft.other.contingent.codes8.sites = soft.other.contingent.codes8.sites,
                  intra.gct.codes = intra.gct.codes, intra.gct.codes.sites = intra.gct.codes.sites,
                  extra.gct.codes = extra.gct.codes, extra.gct.codes.sites = extra.gct.codes.sites,
                  gonad.gct.codes = gonad.gct.codes, gonad.gct.codes.sites = gonad.gct.codes.sites,
                  other.unspec.non.contingent.codes = other.unspec.non.contingent.codes, other.unspec.contingent.codes1 = other.unspec.contingent.codes1, other.unspec.contingent.codes1.sites = other.unspec.contingent.codes1.sites, 
                  other.unspec.contingent.codes2 = other.unspec.contingent.codes2, other.unspec.contingent.codes2.sites = other.unspec.contingent.codes2.sites,
                  other.unspec.contingent.codes3 = other.unspec.contingent.codes3, other.unspec.contingent.codes3.sites = other.unspec.contingent.codes3.sites,
                  other.unspec.contingent.codes4 = other.unspec.contingent.codes4,
                  other.unspec.contingent.codes4.sites = other.unspec.contingent.codes4.sites,
                  epithe.contingent.codes1 = epithe.contingent.codes1, epithe.contingent.codes1.sites = epithe.contingent.codes1.sites, epithe.contingent.codes2 = epithe.contingent.codes2, epithe.contingent.codes2.sites = epithe.contingent.codes2.sites,
                  epithe.contingent.codes3 = epithe.contingent.codes3, epithe.contingent.codes3.sites = epithe.contingent.codes3.sites, epithe.contingent.codes4 = epithe.contingent.codes4, epithe.contingent.codes4.sites = epithe.contingent.codes4.sites,
                  epithe.contingent.codes5 = epithe.contingent.codes5, epithe.contingent.codes5.sites = epithe.contingent.codes5.sites, 
                  epithe.non.contingent.codes = epithe.non.contingent.codes,
                  epithe.contingent.codes6 = epithe.contingent.codes6, epithe.contingent.codes6.sites = epithe.contingent.codes6.sites)

save(cancer.codes, file = 'Z:/Jeremy/GOBACK/Datasets/iccc.codes.to.dx.mappings.v20171018.rdata')


