#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#'                                       
#'                                    ARKANSAS
#'                                       
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

# Prep environment --------------------------------------------------------
load("C:/Users/schraw/Downloads/Updated GOBACK files/arkansas.birthdefects.v20170826.1.rdata")
load("C:/Users/schraw/Downloads/Updated GOBACK files/arkansas.nobirthdefects.v20170826.1.rdata")

require(dplyr)
require(stringi)
require(stringr)
require(readstata13)



# User-defined functions --------------------------------------------------

#' Takes argument X: name of a data frame, and renames the columns the way I like.
restring.columns <- function(x){
  str_replace_all(tolower(colnames(x)),'_','.')
}

#' Takes arguments X: name of a data frame and Y: a vector of columns, and prints 
#' their structures and their unique values.
check.var.levels <- function(x, y){
  tmp <- x[, colnames(x) %in% y]
  apply(tmp[,-1], 2, function(x){
    str(x)
    unique(x)
  })
}

#' Takes argument X: name of a data frame, and selects the 78 variables in ar.aug.
get.ar.cols <- function(x){
  select(x,  infant.birth.weight ,                  
         infant.gestational.age ,                     
         infant.sex ,         
         maternal.race.ethnicity ,                      
         maternal.age.at.delivery ,                       
         maternal.education ,                 
         plurality ,        
         cancer ,     
         dx.date ,      
         age.in.months.at.cancer.diagnosis ,                                
         age.in.months.at.death ,      
         age.in.months.at.end.of.follow.up.if.alive,
         histology3 ,         
         primarysit ,         
         siterecode ,         
         behavior3 ,        
         laterality ,         
         grade ,    
         dxconfirm ,        
         anencephalus ,           
         spina.bifida.without.anencephalus ,                                
         hydrocephalus.without.spina.bifida ,                                 
         encephalocele ,            
         microcephalus ,            
         holoproencephaly ,               
         anophthalmia.microphthalmia ,                          
         congenital.cataract ,                  
         aniridia ,       
         anotia.microtia ,              
         common.truncus ,             
         transposition.of.great.arteries.any.type ,                                       
         ventricular.septal.defect ,                        
         atrial.septal.defect ,                   
         atrioventricular.septal.defect ,                             
         pulmonary.valve.atresia.and.stenosis ,                                   
         ebstein.s.anomaly ,                
         aortic.valve.stenosis ,                    
         hypoplastic.left.heart.syndrome ,                              
         patent.ductus.arteriosus ,                       
         coarctation.of.aorta ,                   
         total.anomalous.pulmonary.venous.return.tapvr ,                                            
         single.ventricle ,               
         interrupted.aortic.arch ,                      
         double.outlet.right.ventricle.dorv ,                                 
         cleft.palate.without.cleft.lip ,                             
         cleft.lip.only.without.cleft.palate ,                                  
         cleft.lip.with.cleft.palate ,                          
         choanal.atresia ,              
         esophageal.atresia.tracheoesophageal.fistula ,                                           
         rectal.large.intestinal.or.anal.atresia.stenosis ,                                               
         pyloric.stenosis ,               
         hirshsprung.s.disease.congenital.megacolon ,                                         
         biliary.atresia ,              
         small.intestinal.atresia.stenosis ,                                
         renal.agenesis.hypoplasia ,                        
         bladder.exstrophy ,                
         obstructive.genitourinary.defect ,                               
         hypospadias ,          
         epispadias ,         
         congenital.posterior.urethral.valves ,                                   
         cloacal.exstrophy ,                
         limb.deficiencies ,                
         gastroschisis ,            
         omphalocele ,          
         congenital.hip.dislocation ,                         
         diaphragmatic.hernia ,                   
         craniosynostosis ,               
         clubfoot ,       
         trisomy.13.patau.syndrome ,                        
         down.syndrome.trisomy.21 ,                       
         trisomy.18.edwards.syndrome ,                          
         turner.syndrome.trisomy.45 ,                         
         deletion.22.q11.2.di.george ,                          
         fetus.newborn.affected.by.maternal.alcohol.use ,                                             
         amniotic.bands ,             
         any.defect ,         
         bdate ,    
         origin.set)}

#' Takes arguments X: data frame and column to operate on and Y: a missing value code.
compute.m.edu2 <- function(x, y){
  factor(
          ifelse(x < 12, 1,
               ifelse(x == 12, 2,
                    ifelse(x %in% 12:98, 3,
                         ifelse(x == y, 7, NA)))),
          levels = c(1:3,7),
          labels = c('< HS', 'HS', '> HS', 'Uknown'))
}

#' Takes argument X: a variable in a data frame which should have certain values 
#' recoded as NA, and Y: a vector of numeric values which should be recoded as NA.
set.to.na <- function(x, y){
  ifelse(x %in% y, NA, x)
}

#' Computes the birth defects groupings variables.  Takes arguemnt X: a data frame
#' and Y: a vector of column indices to sum across.
compute.defect <- function(x, y){
  ifelse(rowSums(x[,y]) >= 1, 1, 0)
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

#' As above, but takes up to 4 arguments and requires 3.
populate.chd.var <- function(w,x,y,z){
  if(missing(y) & missing(z)){
    ifelse((w | x) == 1, 1, 0)
  }
  if(missing(z) & !missing(y)){
    ifelse((w | x | y) == 1, 1, 0)
  }
  else{
    ifelse((x | y | z | w) == 1, 1, 0)
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



# Data cleaning -----------------------------------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Clean, compute and rename columns to prepare for bind.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

#' Clean up the column names.
colnames(ar.jul) <- restring.columns(ar.jul)
colnames(ar.jul) <- str_replace_all(colnames(ar.jul), '[.]$','')

colnames(ar.aug) <- restring.columns(ar.aug)
colnames(ar.aug) <- str_replace_all(colnames(ar.aug), '[.]$','')
colnames(ar.aug) <- str_replace_all(colnames(ar.aug), '[.]{2}','.')

#' Which are conserved, which need to be added to July data?
col.same <- intersect(colnames(ar.aug), colnames(ar.jul))
col.diff <- setdiff(colnames(ar.aug), colnames(ar.jul))
col.diff2 <- setdiff(colnames(ar.aug), colnames(ar.jul))

#' Add birth defects columns to non-birth defects children and set to 0.
ar.jul <- rename(ar.jul, maternal.race.ethnicity = maternal.race, dx.date = dxdate)
ar.jul[ , 20:75] <- 0
names <- c(colnames(ar.jul[ , 1:19]), colnames(ar.aug[,20:75]))
colnames(ar.jul) <- names
ar.jul$any.defect <- 0

ar.jul$cancer <- ifelse(ar.jul$primarysit == "" & ar.jul$grade == "", 0, 1)

#' The new AR data has calculated age at DX and does not provide birth date.  
#' Will need to calculate age at DX in months for AR July data.
#' Will also need an empty variable for age in months at death, which is 
#' given in the August but not July data.
ar.jul$dx.date <- ifelse(ar.jul$dx.date == "", NA, ar.jul$dx.date)
ar.jul$dx.date <- as.Date(ar.jul$dx.date)
ar.jul$age.in.months.at.cancer.diagnosis <- ifelse(is.na(ar.jul$dx.date), NA, (ar.jul$dx.date - ar.jul$bdate)/30)
ar.jul$age.in.months.at.death <- as.numeric(NA)

#' Need to calculate age in months at end of follow up for children without cancer.
ar.jul$age.in.months.at.end.of.follow.up.if.alive <- ifelse(ar.jul$cancer == 0, (as.Date('2011-12-31') - ar.jul$bdate)/30, NA)

#' Some cleaning and restructuring in the August data.
ar.aug <- select(ar.aug, -any.defect.computed)
ar.aug$bdate <- as.Date(NA)
ar.aug <- rename(ar.aug, dx.date = diagnosis.date)
ar.aug$dx.date <- as.Date(ar.aug$dx.date, '%m/%d/%Y')

#' A variable for which dataset a row originated from.
ar.jul$origin.set <- 'July'
ar.aug$origin.set <- 'August'



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Bind birth defects and non-birth defects columns back together.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
tmp <- get.ar.cols(ar.jul)

#' Columns are in identical order.
identical(colnames(ar.aug), colnames(tmp))
col.same <- intersect(colnames(ar.aug), colnames(ar.jul))
col.diff <- setdiff(colnames(ar.aug), colnames(ar.jul))
col.diff2 <- setdiff(colnames(ar.aug), colnames(ar.jul))

ar <- rbind(ar.aug, tmp)

rm(tmp, names, col.diff, col.diff2, col.same)

save(ar, file = 'C:/Users/schraw/Downloads/Updated GOBACK files/arkansas.v20170828.1.rdata')





# Demographic variables ---------------------------------------------------
load("C:/Users/schraw/Downloads/Updated GOBACK files/arkansas.v20170828.1.rdata")

#' Generate a unique study ID, verify that none are duplicated, reorder columns.
ar$studyid <- paste0('AR',runif(629120,1,1000000),stri_rand_strings(629120,2,"[A-Z]"))
ar$studyid <- str_replace(ar$studyid, '[.]','')

ar <- ar[!duplicated(ar$studyid), ]
ar <- ar[, c(79,1:78)]

ar$state <- 'AR'

#' Paternal variables are all missing for AR.
ar$f.edu2 <- as.numeric(NA)
ar$f.race <- as.numeric(NA)
ar$f.age <- as.numeric(NA)

#' Rename and recode what's already been computed.
ar <- rename(ar, sex = infant.sex, m.race = maternal.race.ethnicity, m.age = maternal.age.at.delivery,
             gest.age = infant.gestational.age)

#' Because nothing is ever simple, old data has sex given as 1, 2 or 9 and new data is M or F.
#' Recode to numeric and check distributions against March data.
ar$sex <- as.numeric(ifelse(ar$sex %in% c('M','1'), 1,
                            ifelse(ar$sex %in% c('F','2'), 2,
                                   ifelse(ar$sex == '9', 9, ar$sex))))

ar.mar <- read.csv('C:/Users/schraw/Downloads/GOBACK_Arkansas_20170329.csv', header = TRUE, stringsAsFactors = FALSE)
ar.mar.stata <- read.dta13('C:/Users/schraw/Downloads/GOBACK_Arkansas_20170329.csv.dta')

table(ar$sex)
table(ar.mar$Infant.sex)

save(ar, file = 'C:/Users/schraw/Downloads/Updated GOBACK files/arkansas.v20170828.2.rdata')

#' Maternal race-ethnicity. In July data, 0 = Asian, 1 = white, 2 = black, 3 = AIAN, 4 = Chinese, 5 = Japanese, 6 = Hawaiian,
#' 7 = Other, 8 = Filipino, 9 = Unknown.
ar$m.race <- ifelse(ar$m.race == 'Hispanic', 'hisp(1)', 
                    ifelse(ar$m.race %in% c('1','White'), 'white(2)', ar$m.race))
ar$m.race <- ifelse(ar$m.race %in% c('2', 'Black'), 'black(3)',
                    ifelse(ar$m.race %in% c('0','Asian or Pacific Islander','4','5','6','8'), 'asian(4)', ar$m.race))
ar$m.race <- ifelse(ar$m.race == '7', 'other(5)',
                    ifelse(ar$m.race %in% c('3','American Indian or Alaska Native'), 'aian(6)', ar$m.race))
ar$m.race <- ifelse(ar$m.race %in% c('9','Other/Unknown', NA, ""), 'unknown(7)', ar$m.race)
ar$m.race <- as.numeric(
                ifelse(ar$m.race == 'hisp(1)', 1,
                    ifelse(ar$m.race == 'white(2)', 2,
                           ifelse(ar$m.race == 'black(3)', 3, 
                                  ifelse(ar$m.race == 'asian(4)', 4,
                                         ifelse(ar$m.race == 'other(5)', 5, 
                                                ifelse(ar$m.race == 'aian(6)', 6,
                                                       ifelse(ar$m.race == 'unknown(7)', 7, ar$m.race))))))))

#' Accounting for the fact that 23 NA values are not displayed in the ar.mar table, these totals add up.
#' table(ar.mar$Maternal.race.ethnicity)
#' table(ar$m.race)

#' Birth year (for children without birth defects) and plurality.
ar$birth.yr <- as.numeric(ifelse(is.na(ar$bdate), NA, substr(ar$bdate, 1, 4)))

ar$plu <- ifelse(ar$plurality %in% c('Singleton','1'),1,
                 ifelse(ar$plurality %in% c('Twin','2'),2,
                        ifelse(ar$plurality %in% c('Triplet','3'),3,
                               ifelse(ar$plurality == '4', 4, 
                                      ifelse(ar$plurality == '5', 5,
                                             ifelse(ar$plurality %in% c('Other Multiple Triplet','9','0'),9,9))))))
ar$plu <- ifelse(is.na(ar$plu), 9, ar$plu)

#' Need to figure out what to do about birth weight, which is provided categorically
#' for children in the birth defects set.
ar$birth.wt <- ifelse(ar$origin.set == 'July' & (!is.na(ar$infant.birth.weight)), as.numeric(ar$infant.birth.weight), NA)
ar$birth.wt <- set.to.na(ar$birth.wt, 0)

ar$birth.wt.cat <- ifelse(ar$origin.set == 'August' & ar$infant.birth.weight == '2500g+*', 'Normal or High BW',
                               ifelse(ar$origin.set == "July" & ar$birth.wt > 2499, 'Normal or High BW',
                                      ifelse(ar$origin.set == 'August' & ar$infant.birth.weight == '1500-2499g', "Low Birthweight", 
                                             ifelse(ar$origin.set == 'July' & ar$birth.wt %in% 1500:2499, 'Low Birthweight',
                                                    ifelse(ar$origin.set == 'July' & ar$birth.wt %in% 400:1499, '400-1499g',
                                                           ifelse(ar$origin.set == 'July' & ar$birth.wt < 400, '<400g', 99999))))))

ar$gest.age <- set.to.na(ar$gest.age, c(0,99))

save(ar, file = 'C:/Users/schraw/Downloads/Updated GOBACK files/arkansas.v20170828.3.rdata')

#' Maternal education is an integer where 99 codes for missing.
ar$maternal.education <- ifelse(is.na(ar$maternal.education), 99, ar$maternal.education)
ar$m.edu2 <- compute.m.edu2(ar$maternal.education, 99)

save(ar, file = 'C:/Users/schraw/Downloads/Updated GOBACK files/arkansas.v20170829.1.rdata')




# Person-years ------------------------------------------------------------
load("C:/Users/schraw/Downloads/Updated GOBACK files/arkansas.v20170829.1.rdata")

#' 'August' is equivalent to having a birth defect.
#' 'July' is equivalent to not having having a birth defect.
ar$person.yrs <- ifelse(ar$origin.set == "August" & ar$cancer == 0 & is.na(ar$age.in.months.at.death), ar$age.in.months.at.end.of.follow.up.if.alive/12,
                        ifelse(ar$origin.set == 'August' & ar$cancer == 0 & !is.na(ar$age.in.months.at.death), ar$age.in.months.at.death/12,
                                ifelse(ar$origin.set == 'August' & ar$cancer == 1, ar$age.in.months.at.cancer.diagnosis/12,
                                       ifelse(ar$origin.set == 'July' & ar$cancer == 0, (as.Date('2011-12-31') - ar$bdate)/365,
                                              ifelse(ar$origin.set == 'July' & ar$cancer == 1, (ar$dx.date - ar$bdate)/365, NA)))))

#' Looks good.
print(ar[1:20,c(89,79,9,10,11,12,13,78)])
print(ar[400000:400020,c(89,79,9,10,11,12,13,78)])

tmp <- arrange(ar, desc(cancer))
print(tmp[1:20,c(89,79,9,10,11,12,13,78)])
print(tmp[1000:1020,c(89,79,9,10,11,12,13,78)])

tmp <- arrange(ar, desc(person.yrs))
print(tmp[1:20,c(89,79,9,10,11,12,13,78)])

tmp[1:100,c(32,22,10:13,21)]
tmp[101:200,c(32,22,10:13,21)]

#' We have a single person who shows up as being diagnosed before birth.
#' Study ID is AR535467515269719ZR.
tmp <- arrange(ar, person.yrs)
print(tmp[1:20,c(1,89,79,9,10,11,12,13,78)])
print(tmp[tmp$studyid == 'AR535467515269719ZR', ])

summary(ar$person.yrs)
table(is.na(ar$person.yrs))

#' Although its not impossible this would happen, its more likely an error.
#' Remove this individual.
ar <- ar[ar$studyid != 'AR535467515269719ZR', ]

save(ar, file = 'C:/Users/schraw/Downloads/Updated GOBACK files/arkansas.v20170829.2.rdata')




# Checking for anomalies in birth cert vars -------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Check the birth certificate variables for anomalies.
#' 
#' I actually don't want to spend much time cleaning these right now.  It would be more
#' efficient to identify them, establish some rules and deal with them all at once 
#' when we bind the datasets together.
#' 
#' Make a note of these issues and set them aside fow now.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
load("C:/Users/schraw/Downloads/Updated GOBACK files/arkansas.v20170829.2.rdata")

birth.cert.vars <- c('studyid', 'state','sex','m.race','m.edu','m.edu2','m.age','birth.wt','gest.age','birth.yr','plu','f.race',
                     'f.edu','f.edu2','f.age', 'person.yrs')
check.var.levels(ar, birth.cert.vars)

#' Seeing 0's, 99's and 999's for some of the numeric variables: maternal age, birthweight, gestational age.
#' Gestational age has some values of 2, 8, 12, 16 and 91.
#' Maternal age has some unlikely values 391, 372, 394, 10, 11, 6 and in the range 54:67.
tmp <- arrange(ar, gest.age)
tmp <- arrange(ar, desc(gest.age))
tmp <- filter(ar, gest.age < 22)
tmp <- arrange(ar, birth.wt)
tmp <- arrange(ar, desc(birth.wt))

print(head(select(tmp, studyid, gest.age, birth.wt, infant.birth.weight, birth.wt.cat, origin.set, any.defect, age.in.months.at.death), 100))




# Birth defects variables -------------------------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' AR provided data on birth defects by name rather than BPA code.  
#' We need to make sure that these match the standard defects 
#' Tiffany computed.  Philip should have final sayon any recoding here.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
load("Z:/Jeremy/GOBACK/Datasets/Arkansas/arkansas.v20170829.2.rdata")
load("C:/Users/schraw/Downloads/Jeremy/Jeremy/R data files/Datasets/Michigan/michigan.rdata")

tiff.defects <- colnames(mi[177:277])
rm(mi)

tiff.defects <- tiff.defects[!grepl('_Other_', tiff.defects)]
tiff.defects <- tolower(tiff.defects)
tiff.defects <- str_replace_all(tiff.defects, '_','.')

#' Rename those which are the same minus capitalization and periods.
tmp <- colnames(ar[21:76])
tmp <- str_replace_all(tmp,'[.]','')
names(ar) <- c(colnames(ar[1:20]),tmp,colnames(ar[77:89]))

#' Rename the ones that have already been computed under different titles.
ar <- rename(ar, 
             any.birthdefect = any.defect,
             spinabifida.wo.anencephaly = spinabifidawithoutanencephalus,
             hydrocephalus.wo.spinabifida = hydrocephaluswithoutspinabifida,
             anotia.microtia = anotiamicrotia,
             anopthalmos.micropthalmos = anophthalmiamicrophthalmia,
             transpositionofgreatvessels = transpositionofgreatarteriesanytype,
             endocardialcushiondefect = atrioventricularseptaldefect,
             pulmvalveatresiaandstenosis = pulmonaryvalveatresiaandstenosis,
             ebsteinanomaly = ebsteinsanomaly,
             patentductusarteriosis = patentductusarteriosus,
             esophagealatre.tracheofist = esophagealatresiatracheoesophagealfistula,
             rectal.largeintestatresia.sten = rectallargeintestinaloranalatresiastenosis,
             hirshsprungdisease = hirshsprungsdiseasecongenitalmegacolon,
             renalagenesis.hypoplasia = renalagenesishypoplasia,
             obstructivegenitourinarydefects= obstructivegenitourinarydefect,
             trisomy13 = trisomy13patausyndrome,
             downsyndrome = downsyndrometrisomy21,
             trisomy18 = trisomy18edwardssyndrome,
             fetalalcoholsyndrome = fetusnewbornaffectedbymaternalalcoholuse,
             cleftpalate.wo.cleft.lip = cleftpalatewithoutcleftlip,
             cleftpalateandcleftlip = cleftlipwithcleftpalate,
             limb.deficiencies.unspecified = limbdeficiencies,
             di.george.syndrome = deletion22q112digeorge,
             turner.syndrome = turnersyndrometrisomy45)

#' Interrupted aortic arch is a special, complete case of coarctation (narrowing) of the aorta.
ar$coarctationofaorta <- compute.defect(ar, c(41,44))

#' Compute merged gastroschis.omphalocele, hypospadias.epispadias and cleft.lip variables.
ar$gastroschisis.omphalocele <- compute.defect(ar, 64:65)
ar$hypospadias.epispadias <- compute.defect(ar, 59:60)
ar$cleftlip.w.and.wo.cleftpalate <- compute.defect(ar, 47:48)
ar$matexposuresaffectingfetus <- ifelse(ar$fetalalcoholsyndrome == 1, 1, 0)

ar$defect.total <- rowSums(ar[,c(21:76)])

ar <- ar[ , -c(59, 60, 64, 65)]

#' Initalize empty variables for birth defects that are not provided.
ar$tetralogyoffallot <- as.numeric(NA)
ar$trivalveatresiaandstenosis <- as.numeric(NA)
ar$pulmonaryarteryanomalies <- as.numeric(NA)
ar$lungagenesis.hypoplasia <- as.numeric(NA)
ar$syphilis <- as.numeric(NA)
ar$otherinfectionsperinatalperiod <- as.numeric(NA)
ar$other.unspeccongenitalanomalies <- as.numeric(NA)
ar$infectiouscondperinatalperiod <- as.numeric(NA)
ar$familial.congenitalneoplasms <- as.numeric(NA)
ar$endocrine.metabolicdisorders <- as.numeric(NA)
ar$diseasesoftheblood <- as.numeric(NA)
ar$otherdiseaseseye <- as.numeric(NA)
ar$hearingdeficiency <- as.numeric(NA)
ar$conganomaliesintegument <- as.numeric(NA)
ar$upperlimbreductiondeformities <- as.numeric(NA)
ar$lowerlimbreductiondeformities <- as.numeric(NA)
ar$otherdiseasescns.pns <- as.numeric(NA)
ar$otherdiseasesheart.circsystem <- as.numeric(NA)
ar$otherdiseasesgastrosystem <- as.numeric(NA)
ar$otherdiseasesgeni.urinsystem <- as.numeric(NA)
ar$othermuscskelsystemdiseases <- as.numeric(NA)

#' Compute composite birth defects variables.
ar$conganomalies.cns <- compute.defect(ar, 21:26)
ar$conganomalies.eye <- compute.defect(ar, 27:29)
ar$oral.clefts <- compute.defect(ar, 45:47)
ar$conganomaliesdigestivesystem <- compute.defect(ar, c(50:55,60))
ar$conganomaliesgenitalandurinary <- compute.defect(ar, c(56:59,87))
ar$conganomaliesmusculoskelsys <- compute.defect(ar, c(61:65,86))
ar$chromosomalanomalies <- compute.defect(ar, 66:70)
ar$overall.chrom1 <- compute.defect(ar, 66:70)

ar$conganomalies.respsys <- ifelse(ar$choanalatresia == 1, 1, 0)
ar$conganomalies.ear.face.neck <- ifelse(ar$anotia.microtia == 1, 1, 0)

def.diff <- setdiff(tiff.defects, colnames(ar[-c(1:20,74:85)]))
print(def.diff)

#' This is the best we can do with the birth defects for now.
save(ar, file = 'C:/Users/schraw/Downloads/Updated GOBACK files/arkansas.v20170831.1.rdata')






# General cancer variables ------------------------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Clean up morphology and site codes; compute cancer1, laterality1, 
#' cancertime1 and number diagnoses.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
load("Z:/Jeremy/GOBACK/Datasets/Arkansas/arkansas.v20170831.1.rdata")
load("Z:/Jeremy/GOBACK/Datasets/iccc.codes.to.dx.mappings.v20171018.rdata")

ar <- rename(ar, morphology = histology3, site.code = primarysit)
ar$site.code <- str_replace(ar$site.code, '^C', "")
ar$site.code <- ifelse(ar$site.code == '', NA, ar$site.code)
ar$site.code <- as.numeric(ar$site.code)

ar$morphology <- ifelse(ar$morphology == '', NA, ar$morphology)
ar$morphology <- as.numeric(ar$morphology)

#' Compute cancer1.
ar$cancer1 <- compute.cancer1(ar$morphology, ar$site.code)
ar$cancer1 <- compute.cancer2(ar$morphology, ar$site.code, ar$cancer1)

#' Old method of computing cancer1.  Deprecated.  Kept in for reference.
#' ar <- left_join(ar, site.morph.codes, by = c('site.code','morphology'))

#' Laterality. In the July data non-cancer cases have laterality represented as an empty
#' string; in the August data they have it represented as NA. 
ar$laterality1 <- as.numeric(ifelse(ar$laterality == "", NA, ar$laterality))

#' Cancertime1.  Computed by subtracting DOB from date of DX for July data, and by using 
#' age.in.months.at.cancer.diagnosis for August data.
ar$cancertime1 <- ifelse(ar$cancer == 1 & ar$origin.set == 'July', ((ar$dx.date - ar$bdate)/365)*12,
                         ifelse(ar$cancer == 1 & ar$origin.set == 'August', ar$age.in.months.at.cancer.diagnosis,
                                ifelse(ar$cancer == 0, NA, 99999)))

#' Number of diagnoses.  A moot point for AR, which only reported the first DX.
ar$num.diagnoses <- ifelse(ar$cancer == 1, 1,
                           ifelse(ar$cancer == 0, 0, 99999))

save(ar, file = 'Z:/Jeremy/GOBACK/Datasets/Arkansas/arkansas.v20170908.1.rdata')




# Individual cancer variables ---------------------------------------------
load("Z:/Jeremy/GOBACK/Datasets/Arkansas/arkansas.v20170908.1.rdata")

#' A list of all the named cancer diagnoses.  Exclude NA.
l <- c(unique(ar$cancer1))
l <- l[-1]

#' Initialize an empty variable in the ar data frame for every item in the list.
#' All values for all cancers set to zero.
#' Requires the gen.cancer.var function.
for (i in l){
  x <- (paste0('ar$',i))
  ar[i] <- gen.cancer.var(x)
}

#' This is fine for the non-cancer kids.
#' For the cancer kids, for every column I created, go through every row in the 
#' dataset and if the value of the cancer1 variable is the same as the name of
#' the column, mark it 1.
tmp <- filter(ar, cancer == 0)
tmp2 <- filter(ar, cancer == 1)

for (i in 126:155){
  can.name <- (names(tmp2)[i])
  for (j in tmp2){
    tmp2[,i] <- ifelse(tmp2$cancer1 == can.name, 1, 0)
  }
}

#' Confirm that the number of 1's for each cancer correspond to the frequencies
#' of their diagnosis per the cancer1 variable.
for (i in 126:155){
  print(names(tmp2)[i])
  print(table(tmp2[i]))
}

#' Bind the data together.
ar <- arrange(rbind(tmp, tmp2), studyid)

save(ar, file = 'Z:/Jeremy/GOBACK/Datasets/Arkansas/arkansas.v20170908.2.rdata')




# Generate [cancer.category].any variables --------------------------------
load("Z:/Jeremy/GOBACK/Datasets/Arkansas/arkansas.v20170908.2.rdata")

ar$leu.any <- populate.cancer.var(ar$all, ar$aml, ar$leu.other)
ar$lym.any <- populate.cancer.var(ar$hl, ar$nhl, ar$lym.other)
ar$cns.any <- populate.cancer.var(ar$astro, ar$medullo, ar$cns.other)
ar$pns.any <- populate.cancer.var(ar$neuro, ar$pns.other)
ar$renal.any <- populate.cancer.var(ar$nephro, ar$renal.other)
ar$hepatic.any <- populate.cancer.var(ar$hepato, ar$hepatic.other)
ar$bone.any <- populate.cancer.var(ar$osteo, ar$ewing, ar$bone.other)
ar$rms.any <- populate.cancer.var(ar$erms, ar$arms, ar$rms.other)
ar$soft.any <- populate.cancer.var(ar$rms.any, ar$soft.other)
ar$gct.any <- populate.cancer.var(ar$gct.extra, ar$gct.gonad, ar$gct.intra)

save(ar, file = 'Z:/Jeremy/GOBACK/Datasets/Arkansas/arkansas.v20170908.3.rdata')





# Re-order columns --------------------------------------------------------

#' Move all birth defects columns so they are contiguous.
#' Move all demographics columns so they are contiguous.
ar <- ar[, c(1:20, 74:85, 21:73, 86:121, 122:165)]

save(ar, file = 'Z:/Jeremy/GOBACK/Datasets/Arkansas/arkansas.v20170912.1.rdata')





# Updates to birth defects coding -----------------------------------------

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2017.09.12.
#' 
#' I am revisiting the remaining birth defects variables now that we have settled on our
#' final rules.  Many changes/decisions were made in our meeting today:
#' 
#' - Keep holoprosencephaly, TAPVR, single ventricle, DORV, clubfoot, craniosynostosis, 
#'   small intestinal atresia, Turner Syndrome and DiGeorge Syndrome as their own 
#'   variables.
#' - Make an APVR category for CHDs.  TAPVR falls under this category.
#' - Make a conotruncal category for CHDs.  
#' - Make a septal defect category.  All VSD and ASD fall under this category.
#' - Split hypospadias and epispadias into two separate variables.
#' - If a child has amniotic bands, mark any.defect as 1, but drop the amniotic bands 
#'   variable (data is already coded this way).
#' - Fold congenital posterior urethral valves into obstructive.genitourinary.defect.
#' - Fold cloacal exstrophy into digestive.system.other.major.
#' - For AR data, common truncus is synonymous with both truncus arteriosis and IAA type
#'   B, and should be grouped with conotruncal CHDs.  The IAA variable is taken to refer
#'   specifically to types A or C and should be grouped with LVOT defects.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

#' I removed the hypospadias and epispadias variables previously.  Load the last dataset
#' they can be found in.
load("Z:/Jeremy/GOBACK/Datasets/Arkansas/arkansas.v20170829.2.rdata")
tmp <- ar

#' Load the current data.
load("Z:/Jeremy/GOBACK/Datasets/Arkansas/arkansas.v20170912.1.rdata")

ar <- left_join(ar, select(tmp, studyid, hypospadias, epispadias), by = 'studyid')

#' Philip is confident that the 6,535 children marked as 0 for any birth defect in the 
#' AR August data do in fact have some defect.
ar$any.birthdefect <- with(ar, ifelse(origin.set == 'August' & any.birthdefect == 0, 1, any.birthdefect))

#' No data is provided or can be discerned for minor defects in these children.
ar$minor.status <- NA
ar$minordefect.total <- NA

#' If any named defect is considered major, then major defect total would be the same as 
#' defect total in this data set.
ar$majordefect.total <- ar$defect.total

#' Fold CPUV into obstructive genitourinary defects.
ar$obstructivegenitourinarydefects <- ifelse(ar$congenitalposteriorurethralvalves == 1 & ar$obstructivegenitourinarydefects == 0, 
                                             1, ar$obstructivegenitourinarydefects)






# Generate [bodypart].other.[major or minor] variables --------------------
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

for (i in newcol){
  ar[, i] <- NA
}

ar$digestivesystem.other.major <- ifelse(ar$cloacalexstrophy == 1, 1, ar$digestivesystem.other.major)





# Generate level 3 CHD groups ---------------------------------------------
ar$conotruncal.defects <- populate.chd.var(ar$commontruncus, ar$transpositionofgreatvessels, ar$doubleoutletrightventricledorv)
ar$avsd <- ifelse(ar$endocardialcushiondefect == 1, 1, 0)
ar$apvr <- ifelse(ar$totalanomalouspulmonaryvenousreturntapvr == 1, 1, 0)
ar$lvot.defects <- populate.chd.var(ar$hypoplasticleftheartsyndrome, ar$interruptedaorticarch, ar$coarctationofaorta, ar$aorticvalvestenosis)
#' An ugly hack: the function requires 3 arguments so supply the same one multiple times when there are 2 subcategories of CHD in a group.
ar$rvot.defects <- populate.chd.var(ar$pulmvalveatresiaandstenosis, ar$ebsteinanomaly, ar$ebsteinanomaly)
ar$septal.defects <- populate.chd.var(ar$ventricularseptaldefect, ar$atrialseptaldefect, ar$atrialseptaldefect)

save(ar, file = 'Z:/Jeremy/GOBACK/Datasets/Arkansas/arkansas.v20170913.1.rdata')





# Standardize variable names ----------------------------------------------
load("Z:/Jeremy/GOBACK/Datasets/Arkansas/arkansas.v20170913.1.rdata")

ar <- rename(ar, overall.chrom = overall.chrom1,
             holoprosencephaly = holoproencephaly,
#             anopthalmus.micropthalmus = anopthalmos.micropthalmos,
#             congenital.cataract = congenitalcataract,
             lung.agenesis.hypoplasia = lungagenesis.hypoplasia,
             choanal.atresia = choanalatresia,
             cleft.palate.wo.cleft.lip = cleftpalate.wo.cleft.lip,
             cleft.lip.w.and.wo.cleft.palate = cleftlip.w.and.wo.cleftpalate,
             cleft.palate.and.cleft.lip = cleftpalateandcleftlip,
             conganomalies.digestivesystem = conganomaliesdigestivesystem,
             esophageal.atre.tracheofist = esophagealatre.tracheofist,
             pyloric.stenosis = pyloricstenosis,
             hirshsprung.disease = hirshsprungdisease,
             biliary.atresia = biliaryatresia,
             small.intestinal.atresia = smallintestinalatresiastenosis,
             conganomalies.genitalandurinary = conganomaliesgenitalandurinary,
             renal.agenesis.hypoplasia = renalagenesis.hypoplasia,
             bladder.exstrophy = bladderexstrophy,
             obstructive.genitourinary.defects = obstructivegenitourinarydefects,
             conganomalies.musculoskelsys = conganomaliesmusculoskelsys,
             limb.deformities.unspecified = limb.deficiencies.unspecified,
             upper.limb.reduction.deformities = upperlimbreductiondeformities,
             lower.limb.reduction.deformities = lowerlimbreductiondeformities,
             congenital.hip.dislocation = congenitalhipdislocation,
             diagphragmatic.hernia = diaphragmatichernia,
             conganomalies.integument = conganomaliesintegument,
             down.syndrome = downsyndrome,
             double.outlet.right.ventricle = doubleoutletrightventricledorv,
             total.anomalous.pulmonary.venous.return = totalanomalouspulmonaryvenousreturntapvr,
             interrupted.aortic.arch.type.a.or.c = interruptedaorticarch,
             common.truncus = commontruncus,
             aortic.valve.stenosis = aorticvalvestenosis,
             tetralogy.of.fallot = tetralogyoffallot,
             transposition.of.greatvessels = transpositionofgreatvessels)



# Clean the oral clefts variables -----------------------------------------

#' If a child is flagged as having both cleft palate with cleft lip & cleft palate without cleft lip, set 
#' cleft.palate.wo.cleft.lip to 0.
ar$cleft.palate.wo.cleft.lip <- ifelse(ar$cleft.palate.and.cleft.lip == 1 & ar$cleft.palate.wo.cleft.lip == 1,
                                       0, ar$cleft.palate.wo.cleft.lip)

#' If a child is indicated as having both cleft lip w/ or w/o cleft palate and as having cleft palate but not cleft lip,
#' recode cleft palate without cleft lip 0 and set cleft.palate.and.cleft.lip to 1.
ar$cleft.palate.and.cleft.lip <- ifelse(ar$cleft.lip.w.and.wo.cleft.palate == 1 & ar$cleft.palate.wo.cleft.lip == 1,
                                       0, ar$cleft.palate.and.cleft.lip)
ar$cleft.palate.wo.cleft.lip <- ifelse(ar$cleft.lip.w.and.wo.cleft.palate == 1 & ar$cleft.palate.wo.cleft.lip == 1,
                                       0, ar$cleft.palate.wo.cleft.lip)
ar$oral.clefts <- compute.defect(ar, 47:49)



# Last minute computation of missing variables ----------------------------
ar$conganomalies.heart.circsys <- compute.defect(ar, 43:57)



# Extract the standardized variables --------------------------------------
ar.dem <- get.standard.dem.vars(ar)
ar.def <- get.standard.def.vars(ar)
ar.can <- get.standard.cancer.vars(ar)

ar <- left_join(ar.dem, ar.def, by = 'studyid')
ar <- left_join(ar, ar.can, by = 'studyid')

save(ar, file = 'Z:/Jeremy/GOBACK/Datasets/Arkansas/arkansas.v20170913.2.rdata')

# Recover missing birthweights in BD cohort -------------------------------

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2017.11.16.
#' 
#' For reasons that I don't claim to understand, I decided that I could probably recover
#' the original continuous birthweight data from the set of Arkansas kids with birth
#' defects by creating a sort of pseudo-id variable.  I did this by concatenating several
#' extant variables into something that was ALMOST unique.  This allowed me to recover
#' birth weights on ~21,500/23,000 of kids with birth defects.
#' 
#' It was kind of a clusterfuck.                                       
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
load("Z:/Jeremy/GOBACK/Datasets/Arkansas/arkansas.raw.data.march.rdata")
load("Z:/Jeremy/GOBACK/Datasets/Arkansas/arkansas.birthdefects.v20170826.1.rdata")
load("Z:/Jeremy/GOBACK/Datasets/Arkansas/arkansas.v20170828.3.rdata")

pattern <- '[[:digit:]]+[.]{1}[[:digit:]]{3}'

colnames(ar.mar) <- restring.columns(ar.mar)
colnames(ar.aug) <- restring.columns(ar.aug)

#' Generate pseudo-ID for the old data.
tmp <- ar.mar[!is.na(ar.mar$any.defect), ]
tmp$infant.gestational.age <- ifelse(tmp$infant.gestational.age == 0, NA, tmp$infant.gestational.age)
tmp <- arrange(tmp, infant.gestational.age)
for (i in 13:68){
  tmp[,i] <- ifelse(tmp[,i] == 'Yes', 1, 0)
}
rm(i)
tmp$new.age <- ifelse(!is.na(tmp$age.in.months.at.end.of.follow.up.if.alive), tmp$age.in.months.at.end.of.follow.up.if.alive + 7.528764, tmp$age.in.months.at.end.of.follow.up.if.alive)
tmp$new.age.short <- stringr::str_extract(tmp$new.age, pattern)
tmp$infant.birth.weight <- set.to.na(tmp$infant.birth.weight, c(0,99))
tmp$infant.sex <- as.numeric(ifelse(tmp$infant.sex %in% c('M','1'), 1,
                                    ifelse(tmp$infant.sex %in% c('F','2'), 2, tmp$infant.sex)))
tmp <- arrange(tmp, infant.gestational.age)
tmp$hash <- as.character(tmp$infant.gestational.age)
for (i in c(3,5,6,72)){
  tmp$hash <- paste0(tmp$hash, tmp[,i])
}

indices <- which(duplicated(tmp$hash))
indices <- c(indices, indices - 1)

tmp <- filter(tmp, !(row.names(tmp) %in% as.character(indices)))
tmp <- tmp[!duplicated(tmp$hash), ]
tmp <- select(tmp, hash, infant.birth.weight, age.in.months.at.cancer.diagnosis, age.in.months.at.death, age.in.months.at.end.of.follow.up.if.alive)
tmp <- arrange(tmp, hash)

#' Generate pseudo-ID for the new data.
ar$age.short <- stringr::str_extract(ar$age.in.months.at.end.of.follow.up.if.alive, pattern)
ar <- arrange(filter(ar, origin.set == 'August'), gest.age)
ar <- arrange(ar, gest.age)
ar$hash <- as.character(ar$gest.age)
for (i in c(4,6,7,88)){
  ar$hash <- paste0(ar$hash, ar[,i])
}

tmp2 <- arrange(ar, hash)
indices <- which(duplicated(tmp2$hash))
indices <- c(indices, indices - 1)

tmp2 <- filter(tmp2, !(row.names(tmp2) %in% as.character(indices)))
tmp3 <- tmp2[!duplicated(tmp2$hash), ]
tmp2 <- select(tmp2, hash, studyid)
tmp2 <- left_join(tmp2, tmp, by = 'hash')
tmp2 <- rename(tmp2, new.birth.wt = infant.birth.weight)

recovered.ar.birthweights <- left_join(tmp2, select(ar, studyid, birth.wt.cat), by = 'studyid')
save(recovered.ar.birthweights, file = 'Z:/Jeremy/GOBACK/Datasets/Arkansas/ar.recovered.birthweights.v20180302.rdata')

rm(ar, ar.aug, ar.mar, recovered.ar.birthweights, tmp, tmp2, i, indices, pattern)
