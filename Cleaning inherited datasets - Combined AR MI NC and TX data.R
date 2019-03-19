#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#'                                       
#'                              ALL STATES
#'                                         
#' Ready to merge all four states together.
#' 
#' As much as possible, I refrained from doing data cleaning on each state
#' separately.  Much of it is done here on the combined file.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

# User-defined functions --------------------------------------------------
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
         soft.other,    gct.intra,     gct.extra,     gct.gonad,     gct.any,       epithe,        other.any,     ependymoma,    pnet,
         runif,         birth.char.flag)  
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



# Load in and merge data --------------------------------------------------
load("W:/Old_genepi2/Jeremy/GOBACK/Datasets/North Carolina/nc.v20171107.3.rdata")
load("W:/Old_genepi2/Jeremy/GOBACK/Datasets/Combined Arkansas Michigan and Texas/ar.mi.tx.v20171106.1.rdata")

goback <- rbind(armitx, nc)

rm(armitx, nc)
gc()

save(goback, file = 'W:/Old_genepi2/Jeremy/GOBACK/Datasets/GOBACK.v20171107.1.rdata')



# Clean up AR birth defects variables -------------------------------------
load('./GOBACK.v20171107.1.rdata')

ar <- filter(goback, state == 'AR')

#' Issue is that individual birth defects variables are coded either 0 or 1.
#' Should be NA if the child has a birth defect, but not the index defect.
for (i in c(23:30,32:36,38:40,46,47,49,51:54,56,57,59:66,68:71,73,74,76:83,85:91,93:103,105:112)){
ar[,i] <- ifelse(ar[,i] == 1, 1,
                   ifelse(ar[,i] == 0 & ar$any.birthdefect == 0, 0,
                                     ifelse(ar[,i] == 0 & ar$any.birthdefect == 1, NA, 2)))
}

#' Birth defects categories variables must now be updated.
update.defect <- function(constituent.defect.indices){
                  ifelse(rowSums(constituent.defect.indices, na.rm = TRUE) >= 1, 1,
                         ifelse(ar$any.birthdefect == 1 & rowSums(constituent.defect.indices) == 0, NA,
                                ifelse(ar$any.birthdefect == 0, 0, 2)))
}

ar$conganomalies.cns <- update.defect(ar[,23:30])
ar$conganomalies.eye <- update.defect(ar[,32:36])
ar$conganomalies.ear.face.neck <- update.defect(ar[,38:40])
ar$conganomalies.heart.circsys <- update.defect(ar[,42:66])
ar$conganomalies.respsys <- update.defect(ar[,68:71])
ar$oral.clefts <- update.defect(ar[,73:74])
ar$conganomalies.digestivesystem <- update.defect(ar[,76:83])
ar$conganomalies.genitalandurinary <- update.defect(ar[,85:91])
ar$conganomalies.musculoskelsys <- update.defect(ar[,93:102])
ar$chromosomalanomalies <- update.defect(ar[,105:111])

goback <- filter(goback, state != 'AR')
goback <- rbind(goback, ar)

rm(ar, i)
gc()

save(goback, file = './goback.v20171107.2.rdata')



# Flag implausible birthweight-gestational age combinations ---------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' I previously noted that some kids had implausible birthweights for their
#' given gestational ages, implying that one or both is incorrect.
#' 
#' One previous study compared several methods of tagging implausible 
#' birthweight-gestational age combinations, but stopped short of 
#' recommending one in particular.
#' 
#' Using one of the methods in this reference, we will flag kis with a BW
#' more than 4 SDs from the gestational age- and sex-specific median.
#' 
#' Median birthweights in GOBACK are similar to recent US data overall
#' (see Talge NM et al, "United States Birth Weight Reference Corrected for
#' Implausible Gestational Age Estimates"), and all sex-gestational age 
#' combinations have at >1,000 observations so I will use the median BWs 
#' from the GOBACK dataset itself to make these flags.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(dplyr)

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/goback.v20171107.2.rdata')

goback$birth.wt <- ifelse(goback$birth.wt == 9999, NA, goback$birth.wt)

print(ggplot(data=goback, aes(x = birth.wt)) +
        geom_histogram(color = 'red', fill = 'white'))

tmp <- aggregate(birth.wt ~ gest.age + sex, data = goback, median)
tmp2 <- aggregate(birth.wt ~ gest.age + sex, data = goback, sd)

birthweight.norms <- left_join(tmp, tmp2, by = c('sex','gest.age'))
birthweight.norms <- filter(birthweight.norms, sex != 9)
birthweight.norms <- rename(birthweight.norms, median.birth.wt = birth.wt.x, sd.birth.wt = birth.wt.y)
birthweight.norms$upper.limit <- birthweight.norms$median.birth.wt + 4*(birthweight.norms$sd.birth.wt) 
birthweight.norms$lower.limit <- ifelse((birthweight.norms$median.birth.wt - 4*birthweight.norms$sd.birth.wt) < 0, 0, birthweight.norms$median.birth.wt - 4*birthweight.norms$sd.birth.wt)

save(birthweight.norms, file = './birthweight.norms.rdata')

goback <- left_join(goback, birthweight.norms, by = c('sex','gest.age'))
goback$birth.char.flag <- ifelse(goback$birth.wt < goback$lower.limit | goback$birth.wt > goback$upper.limit, 1, 0)
goback <- select(goback, -median.birth.wt, -sd.birth.wt, -upper.limit, -lower.limit)

save(goback, file = './goback.v20171108.1.rdata')

rm(tmp, tmp2, birthweight.norms)

#' How many birth defects and cancer cases are we dropping by excluding these kids?
table(goback$birth.char.flag, goback$any.birthdefect, useNA = 'always') #' Excluding 2,416 with any defect, ~ 0.4% of the exposed kids.
gmodels::CrossTable(goback$birth.char.flag, goback$cancer) #' Excluding 38 cancer cases, ~0.25% of of the total kids with cancer.



# Clean up MI maternal demographic vars -----------------------------------
load("./goback.v20171108.1.rdata")

mi <- select(filter(goback, state == 'MI'), -m.edu2, -m.race)
mi <- arrange(mi, studyid)
goback <- filter(goback, state != "MI")

#' Upated MI data.
mi.new <- read.csv(file = 'C:/Users/schraw/Desktop/mi.updated.raceetheduc.csv', header = TRUE, stringsAsFactors = FALSE)
mi.new <- rename(select(mi.new, studyid, m_race, m_edu2), m.race = m_race, m.edu2 = m_edu2)
mi.new$studyid <- paste0('mi',mi.new$studyid)

mi <- left_join(mi, mi.new, by = 'studyid')
rm(mi.new)
gc()

mi.dem <- get.standard.dem.vars(mi)
mi.def <- get.standard.def.vars(mi)
mi.can <- get.standard.cancer.vars(mi)

tmp <- left_join(mi.dem, mi.def, by = 'studyid')
mi <- left_join(tmp, mi.can, by = 'studyid')
rm(mi.dem, mi.def, mi.can, tmp)
gc()

goback <- rbind(goback, mi)
rm(mi)
gc()

save(goback, file = './goback.v20171115.1.rdata')



# Clean up AR maternal education variable ---------------------------------
load('./goback.v20171115.1.rdata')

ar <- filter(goback, state == 'AR')
ar$m.edu2 <- ifelse(ar$m.edu2  == '< HS', 1, 
                        ifelse(ar$m.edu2 == 'HS', 2, 
                               ifelse(ar$m.edu2 == '> HS', 3,
                                      ifelse(ar$m.edu2 == 'Uknown', 7, ar$m.edu2))))

goback <- filter(goback, state != 'AR')
goback$m.edu2 <- as.numeric(goback$m.edu2)
goback$m.edu2 <- ifelse(goback$m.edu2 == 5, 1, 
                        ifelse(goback$m.edu2 == 6, 3, 
                               ifelse(goback$m.edu2 == 7, 2, 
                                      ifelse(goback$m.edu2 == 8, 7, goback$m.edu2))))
goback <- rbind(goback, ar)

save(goback, file = './goback.v20171115.2.rdata')
rm(ar)
gc()



# Convert select demographic variables to factors -------------------------
load('./goback.v20171115.2.rdata')

goback$m.race <- factor(ifelse(goback$m.race == 1, 1,
                            ifelse(goback$m.race == 2, 2,
                                   ifelse(goback$m.race == 3, 3,
                                          ifelse(goback$m.race == 4, 4, 
                                                 ifelse(goback$m.race == 5, 5,
                                                        ifelse(goback$m.race == 6, 6, 
                                                               ifelse(goback$m.race == 7, 7, NA))))))),
                     levels = 1:7, 
                     labels = c('Hispanic','NHW','NHB','Asian','Other','AIAN','Unknown'))

save(goback, file = './goback.v20171115.3.rdata')



# Clean up birthweight category -------------------------------------------
load('./goback.v20171115.3.rdata')

goback$birth.wt.cat <- ifelse(goback$birth.wt.cat == '<400 grams' | goback$birth.wt.cat == '<400g', '<400',
                           ifelse(goback$birth.wt.cat == '400-1499 grams' | goback$birth.wt.cat == '400-1499g', '400-1499',
                                  ifelse(goback$birth.wt.cat == '99999' | goback$birth.wt.cat == '', 'Unknown', 
                                         ifelse(goback$birth.wt.cat == 'Low Birthweight', '1500-2499',
                                                ifelse(goback$birth.wt.cat == 'Normal or High Birthweight' | goback$birth.wt.cat == 'Normal or High BW', '>2499',goback$birth.wt.cat)))))
goback$birth.wt.cat <- ifelse(is.na(goback$birth.wt.cat), 'Unknown',goback$birth.wt.cat)

save(goback, file = './goback.v20171116.1.rdata')



# Replace birthweight with recovered continuous values --------------------
load('./goback.continuous.birthweights.rdata')
load("./goback.v20171116.1.rdata")

goback <- select(goback, -birth.wt)
goback <- left_join(goback, goback.bws, by = 'studyid')
goback <- goback[, c(1:6,159,7:158)]

rm(goback.bws)

#' Inspect the new birth weight variable.  Looks consistent with what was computed before.
print(goback[1:200, 7:8])
print(goback[100000:100200, 7:8])

goback$bw.flag <- ifelse(goback$birth.wt.cat == '>2499' & goback$birth.wt <= 2499, 1, 0)
table(goback$birth.wt.flag)
goback$bw.flag <- ifelse(goback$birth.wt.cat == '1500-2499' & (goback$birth.wt < 1500 | goback$birth.wt > 2499), 1, 0)
table(goback$birth.wt.flag)
goback$bw.flag <- ifelse(goback$birth.wt.cat == '400-1499' & (goback$birth.wt < 400 | goback$birth.wt > 1499), 1, 0)

save(goback, file = './goback.v20171117.1.rdata')



# Merge other.minor and other.major ---------------------------------------

load('./goback.v20171117.1.rdata')

#' Clean up other.major and other.minor vars in AR in preparation.
#' If state is Arkansas and child has no defect, set to 0.  Otherwise, NA.
for (i in c(29,30,35,36,39,40,65,66,70,71,82,83,90,91,101,102,110,111)){
  print(names(goback[i]))
  print(table(goback$state, goback[,i], useNA = 'always'))
}
  
for (i in c(29,30,35,36,39,40,65,66,70,71,82,83,90,91,101,102,110,111)){
  goback[,i] <- ifelse(goback$any.birthdefect == 0 & goback$state == 'AR', 0, goback[,i])
}

goback$conganomalies.cns.other <- compute.other.var(goback$conganomalies.cns.other.major, goback$conganomalies.cns.other.minor)
goback$conganomalies.eye.other <- compute.other.var(goback$conganomalies.eye.other.major, goback$conganomalies.eye.other.minor)
goback$conganomalies.ear.face.neck.other <- compute.other.var(goback$conganomalies.ear.face.neck.other.major, goback$conganomalies.ear.face.neck.other.minor)
goback$conganomalies.heart.circsys.other <- compute.other.var(goback$heart.circsys.other.major, goback$heart.circsys.other.minor)
goback$conganomalies.respsys.other <- compute.other.var(goback$respsys.other.major, goback$respsys.other.minor)
goback$conganomalies.digestivesystem.other <- compute.other.var(goback$digestivesystem.other.major, goback$digestivesystem.other.minor)
goback$conganomalies.genitalandurinary.other <- compute.other.var(goback$genitalandurinary.other.major, goback$genitalandurinary.other.minor)
goback$conganomalies.musculoskelsys.other <- compute.other.var(goback$musculoskelsys.other.major, goback$musculoskelsys.other.minor)
goback$conganomalies.chromosomalanomalies.other <- compute.other.var(goback$chromosomalanomalies.other.major, goback$chromosomalanomalies.other.minor)

goback <- select(goback, -conganomalies.cns.other.major, -conganomalies.cns.other.minor, -conganomalies.eye.other.major, -conganomalies.eye.other.minor,
                 -conganomalies.ear.face.neck.other.major, -conganomalies.ear.face.neck.other.minor, - respsys.other.major, -respsys.other.minor,
                 -digestivesystem.other.major, -digestivesystem.other.minor, -genitalandurinary.other.major, -genitalandurinary.other.minor,
                 -musculoskelsys.other.major, -musculoskelsys.other.minor, -chromosomalanomalies.other.major, -chromosomalanomalies.other.minor,
                 -heart.circsys.other.major, -heart.circsys.other.minor)

goback <- goback[, c(1:28,143,29:32,144,33,34,145,35:58,146,59:61,147,62:71,148,72:77,149,78:86,150,87:93,151,94:142)]

#' Forgot to update other and unspecified variable.
goback$other.unspeccongenitalanomalies <- ifelse(goback$any.birthdefect == 0 & goback$state == 'AR', 0, goback$other.unspeccongenitalanomalies)

save(goback, file = 'goback.v20171130.1.rdata')

#' Crosstabs for the new variables.
#' Philip wanted to review each body part by state to make sure none where too out of whack.
sink(file = 'C:/Users/schraw/Desktop/updated.bd.crosstabs.txt')
for (i in 22:103){
  print(names(goback[i]))
  print(CrossTable(goback$state, goback[,i], prop.c = FALSE, prop.t = FALSE, prop.chisq = FALSE))
}
sink()



# Drop redundant CHD level 3 variables ------------------------------------

load('goback.v20171130.1.rdata')

#' Each of these has only one defect under the category.  No real need to keep.
goback <- select(goback, -avsd, -apvr)

save(goback, file = 'goback.v20171201.1.rdata')



# Update conganomalies.[bodypart] variables -------------------------------

load('goback.v20171201.1.rdata')

#' Let's check that the general body part variables equal the sum of their parts.
goback$conganomalies.cns <- with(goback, ifelse(rowSums(goback[,23:29], na.rm = TRUE) > 0, 1, goback$conganomalies.cns))
goback$conganomalies.eye <- ifelse(rowSums(goback[,31:34], na.rm = TRUE) > 0, 1, goback$conganomalies.eye)
goback$conganomalies.ear.face.neck <- ifelse(rowSums(goback[,36:37], na.rm = TRUE) > 0, 1, goback$conganomalies.ear.face.neck)
goback$conganomalies.heart.circsys <- ifelse(rowSums(goback[,c(40:44,46:49,51,52,54:60)], na.rm = TRUE) > 0, 1, goback$conganomalies.heart.circsys)
goback$conganomalies.respsys <- ifelse(rowSums(goback[,62:64], na.rm = TRUE) >0, 1, goback$conganomalies.respsys)
goback$conganomalies.digestivesystem <- ifelse(rowSums(goback[,69:75], na.rm = TRUE) > 0, 1, goback$conganomalies.digestivesystem)
goback$conganomalies.genitalandurinary <- ifelse(rowSums(goback[,77:82], na.rm = TRUE) > 0, 1, goback$conganomalies.genitalandurinary)
goback$conganomalies.musculoskelsys <- ifelse(rowSums(goback[,84:92], na.rm = TRUE) > 0, 1, goback$conganomalies.musculoskelsys)
goback$chromosomalanomalies <- ifelse(rowSums(goback[,95:100], na.rm = TRUE) > 0, 1, goback$chromosomalanomalies)

save(goback, file = 'goback.v20171201.2.rdata')

#' Re-gen crosstabs.
sink(file = 'C:/Users/schraw/Desktop/updated.bd.crosstabs.txt')
for (i in 22:103){
  print(names(goback[i]))
  print(CrossTable(goback$state, goback[,i], prop.c = FALSE, prop.t = FALSE, prop.chisq = FALSE))
}
sink()



# Review MI any.birthdefect variable --------------------------------------

load('goback.v20171201.2.rdata')
load("W:/Old_genepi2/Jeremy/GOBACK/Datasets/birth.defects.codes.rdata")

#' There are too many children in MI who are tagged as having any birth defect.
#' More specifically, 8.7%.  This drops to 6.2% if you only count kids who have one or
#' more of the defects we have variables for. 
CrossTable(goback$state, goback$any.birthdefect)
CrossTable(goback$state, goback$minor.status)

goback$defect.sum <- rowSums(goback[,22:101], na.rm = TRUE)
goback$defect.flag <- ifelse(goback$defect.sum > 0, 1, 0)
CrossTable(goback$state, goback$defect.flag, prop.chisq = FALSE, prop.c = FALSE)

#' 54,500 kids from MI with defects don't have an ICD9 code in the 740:759 range. 
tmp <- defect.codes
tmp <- tmp[grepl('^mi', tmp$studyid), ]
tmp$flag <- 0

for (i in 2:67){
  tmp$flag <- ifelse(floor(tmp[,i]) %in% 740:759, 1, tmp$flag)
}
table(tmp$flag, useNA = 'always')

tmp <- filter(tmp, flag == 0)
print(tmp[1:100,2:25])

#' The 6.2% of MI kids from GOBACK with a defect variable equal to 1 have at least one code
#' in the 740:759 range.  Conclude that this is not a coding issue.  It is something in the 
#' state data: more codes captured, more aggressive follow up, not all births from state 
#' present in the data etc.
tmp <- filter(filter(goback, any.birthdefect == 1), state == 'MI')
tmp2 <- filter(filter(goback, defect.flag == 1), state == 'MI')
tmp2 <- c(tmp2$studyid)
tmp2 <- defect.codes[defect.codes$studyid %in% tmp2, ]
tmp2$flag <- 0
for (i in 2:67){
  tmp2$flag <- ifelse(floor(tmp2[,i]) %in% 740:759, 1, tmp2$flag)
}
table(tmp2$flag, useNA = 'always')



# Clean blank AR birth defects variables ----------------------------------

load('goback.v20171201.2.rdata')
load("birth.defects.codes.rdata")

#' First off, clean some blank BD variables in the AR data.  If no defect, set to 0.
for (i in c(41,57,59,61,62,85,86,90,101)){
  print(names(goback[i]))
  print(table(goback$state, goback[,i], useNA = 'always'))
}

set.to.zero <- function(bdvar){
  ifelse(goback$any.birthdefect == 0 & goback$state == 'AR', 0, bdvar)
}

for (i in c(41,57,59,61,62,85,86,90,101)){
  goback[,i] <- set.to.zero(goback[,i])
}

save(goback, file = 'goback.v20171201.3.rdata')



# Clean gastroschisis/omphalocele -----------------------------------------

#' Philip wants gastroschisis and omphalocele separated into two different
#' variables.  In addition, it was coded incorrectly in TX and NC.  I will
#' update these two.  Tiffany will update MI.  I will go back to the AR raw
#' data from August and pull the original variables.  

load('goback.v20171201.3.rdata')
load("birth.defects.codes.rdata")
load("W:/Old_genepi2/Jeremy/GOBACK/Datasets/Arkansas/arkansas.v20170828.2.rdata")

#' North Carolina.
tmp <- defect.codes[grepl('^nc', defect.codes$studyid), ]
tmp$gast.flag <- 0
for (i in 2:67){
  tmp$gast.flag <- ifelse(tmp[,i] %in% c(756710:756719), 1, tmp$gast.flag)
}
table(tmp$gast.flag, useNA = 'always')
tmp <- filter(tmp, gast.flag == 1)
gast.ids.nc <- c(tmp$studyid)

tmp <- defect.codes[grepl('^nc', defect.codes$studyid), ]
tmp$om.flag <- 0
for (i in 2:67){
  tmp$om.flag <- ifelse(tmp[,i] %in% c(756700:756709), 1, tmp$om.flag)
}
table(tmp$om.flag, useNA = 'always')
tmp <- filter(tmp, om.flag == 1)
om.ids.nc <- c(tmp$studyid)

#' Texas.
#' The only two codes in the range for gastroschisis and omphalocele are 756.70 & 756.71.
tmp <- defect.codes[grepl('^tx', defect.codes$studyid), ]
tmp <- tmp[!duplicated(tmp$studyid), ]
tmp$gast.flag <- 0
for (i in 2:67){
  tmp$gast.flag <- ifelse(tmp[,i] %in% 756.71, 1, tmp$gast.flag)
}
table(tmp$gast.flag, useNA = 'always')
tmp <- filter(tmp, gast.flag == 1)
gast.ids.tx <- c(tmp$studyid)

tmp <- defect.codes[grepl('^tx', defect.codes$studyid), ]
tmp <- tmp[!duplicated(tmp$studyid), ]
tmp$om.flag <- 0
for (i in 2:67){
  tmp$om.flag <- ifelse(tmp[,i] %in% 756.70, 1, tmp$om.flag)
}
table(tmp$om.flag, useNA = 'always')
tmp <- filter(tmp, om.flag == 1)
om.ids.tx <- c(tmp$studyid)

#' Michigan uses ICD9 codes instead of BPA codes.
#' 756.73 = gastroschisis.
#' 756.72 = omphalocele.
#' These numbers are much lower than before.  But it doesn't appear that there
#' are any extra codes in this range that I'm not capturing, and Tiffany got the
#' same numbers.
tmp <- defect.codes[grepl('^mi', defect.codes$studyid), ]
tmp <- tmp[!duplicated(tmp$studyid), ]
tmp$gast.flag <- 0
for (i in 2:67){
  tmp$gast.flag <- ifelse(tmp[,i] %in% 756.73, 1, tmp$gast.flag)
}
table(tmp$gast.flag, useNA = 'always')
tmp <- filter(tmp, gast.flag == 1)
gast.ids.mi <- c(tmp$studyid)

tmp <- defect.codes[grepl('^mi', defect.codes$studyid), ]
tmp <- tmp[!duplicated(tmp$studyid), ]
tmp$om.flag <- 0
for (i in 2:67){
  tmp$om.flag <- ifelse(tmp[,i] %in% 756.72, 1, tmp$om.flag)
}
table(tmp$om.flag, useNA = 'always')
tmp <- filter(tmp, om.flag == 1)
om.ids.mi <- c(tmp$studyid)

#' Arkansas.
tmp <- select(ar, studyid, origin.set, gastroschisis, omphalocele)
tmp$gastroschisis <- ifelse(tmp$origin.set == 'August' & tmp$gastroschisis == 0, NA, tmp$gastroschisis)
tmp$omphalocele <- ifelse(tmp$origin.set == 'August' & tmp$omphalocele == 0, NA, tmp$omphalocele)

tmp2 <- filter(tmp, gastroschisis == 1)
gast.ids.ar <- c(tmp2$studyid)

tmp2 <- filter(tmp, omphalocele == 1)
om.ids.ar <- c(tmp2$studyid)

gast.ids <- c(gast.ids.ar, gast.ids.mi, gast.ids.nc, gast.ids.tx)
om.ids <- c(om.ids.ar, om.ids.mi, om.ids.nc, om.ids.tx)

rm(gast.ids.ar, gast.ids.mi, gast.ids.nc, gast.ids.tx,
   om.ids.ar, om.ids.mi, om.ids.nc, om.ids.tx,
   tmp, tmp2, ar)

goback$gastroschisis <- as.numeric(NA)
goback$gastroschisis <- ifelse(goback$any.birthdefect == 0, 0, 
                               ifelse(goback$studyid %in% gast.ids, 1, goback$gastroschisis))
table(goback$state, goback$gastroschisis, useNA = 'always')

goback$omphalocele <- as.numeric(NA)
goback$omphalocele <- ifelse(goback$any.birthdefect == 0, 0, 
                               ifelse(goback$studyid %in% om.ids, 1, goback$omphalocele))
table(goback$state, goback$omphalocele, useNA = 'always')

#' Kids who were coded 1 under the old variable, but NA under the new one, should have 
#' other musculoskeletal defects.  All kids in a vector resulting from filtering on these
#' conditions should have conganomalies.musculoskelsys.other updated to 1.
tmp <- filter(goback, gastroschisis.omphalocele == 1)
tmp <- tmp[is.na(tmp$gastroschisis), ]
tmp <- tmp[is.na(tmp$omphalocele), ]
tmp <- c(tmp$studyid)

table(goback$state, goback$conganomalies.musculoskelsys.other, useNA = 'always')
goback$conganomalies.musculoskelsys.other <- ifelse(goback$studyid %in% tmp, 1, goback$conganomalies.musculoskelsys.other)
table(goback$state, goback$conganomalies.musculoskelsys.other, useNA = 'always')

goback <- select(goback, -gastroschisis.omphalocele)
goback <- goback[,c(1:86,149,150,87:148)]

save(goback, file = 'goback.v20171206.1.rdata')



# Verify overall.chrom ----------------------------------------------------

#' For subsequent modeling steps the exposed children will be split into sets
#' with and without chromosomal defects.  It would be good to verify that this
#' variable holds up in light of all the changes we've made.
load('goback.v20171206.1.rdata')
load("birth.defects.codes.rdata")

table(goback$state, goback$overall.chrom, useNA = 'always')
unique(goback$overall.chrom)

#' First problem: In NC, its counting the number of chromosomal defects?  
#' Elsewhere, dichotomous.
goback$chrom.flag <- ifelse(goback$overall.chrom > 0, 1, 0)
table(goback$state, goback$chrom.flag, useNA = 'always')

#' Do those numbers jive with a search through the defects codes?
#' In the BPA system, chromosomal anomalies fall under 758.000:758.999.
tmp <- defect.codes[grepl('^nc', defect.codes$studyid), ]
tmp$flag <- 0
for (i in 2:67){
  tmp$flag <- ifelse(tmp[,i] %in% 758000:758999, 1, tmp$flag)
}
table(tmp$flag, useNA = 'always')
tmp <- filter(tmp, flag == 1)

ids1 <- c(tmp$studyid)

ids2 <- filter(filter(goback, overall.chrom > 0), state == 'NC')
ids2 <- c(ids2$studyid)

#' They're close, but I find a few extra kids in the defects codes data frame with chromosomal anomalies.
ids.diff <- setdiff(ids1, ids2)




# Verify overall.chrom: NC ------------------------------------------------

#' Down.
tmp <- defect.codes[grepl('^nc', defect.codes$studyid), ]
tmp$down.flag <- 0
for (i in 2:67){
  tmp$down.flag <- ifelse(tmp[,i] %in% 758000:758099, 1, tmp$down.flag)
}
table(tmp$down.flag, useNA = 'always')
tmp <- filter(tmp, down.flag == 1)
tmp <- c(tmp$studyid)
goback$down.syndrome <- ifelse(goback$studyid %in% tmp, 1, goback$down.syndrome)

#' Patau.
tmp <- defect.codes[grepl('^nc', defect.codes$studyid), ]
tmp$patau.flag <- 0
for (i in 2:67){
  tmp$patau.flag <- ifelse(tmp[,i] %in% 758100:758199, 1, tmp$patau.flag)
}
table(tmp$patau.flag, useNA = 'always')
tmp <- filter(tmp, patau.flag == 1)
tmp <- c(tmp$studyid)
goback$trisomy13 <- ifelse(goback$studyid %in% tmp, 1, goback$trisomy13)
goback$overall.chrom <- ifelse(goback$studyid %in% tmp, 1, goback$overall.chrom)

#' Edwards.
tmp <- defect.codes[grepl('^nc', defect.codes$studyid), ]
tmp$ed.flag <- 0
for (i in 2:67){
  tmp$ed.flag <- ifelse(tmp[,i] %in% 758200:758299, 1, tmp$ed.flag)
}
table(tmp$ed.flag, useNA = 'always')
tmp <- filter(tmp, ed.flag == 1)
tmp <- c(tmp$studyid)
goback$trisomy18 <- ifelse(goback$studyid %in% tmp, 1, goback$trisomy18)
goback$overall.chrom <- ifelse(goback$studyid %in% tmp, 1, goback$overall.chrom)

#' DiGeorge.
tmp <- defect.codes[grepl('^nc', defect.codes$studyid), ]
tmp$dig.flag <- 0
for (i in 2:67){
  tmp$dig.flag <- ifelse(tmp[,i] %in% 758370, 1, tmp$dig.flag)
}
table(tmp$dig.flag, useNA = 'always')
tmp <- filter(tmp, dig.flag == 1)
tmp <- c(tmp$studyid)
goback$di.george.syndrome <- ifelse(goback$studyid %in% tmp, 1, goback$di.george.syndrome)
goback$overall.chrom <- ifelse(goback$studyid %in% tmp, 1, goback$overall.chrom)

#' Turner.
tmp <- defect.codes[grepl('^nc', defect.codes$studyid), ]
tmp$turner.flag <- 0
for (i in 2:67){
  tmp$turner.flag <- ifelse(tmp[,i] %in% 758600:758690, 1, tmp$turner.flag)
}
table(tmp$turner.flag, useNA = 'always')
tmp <- filter(tmp, turner.flag == 1)
tmp <- c(tmp$studyid)
goback$turner.syndrome <- ifelse(goback$studyid %in% tmp, 1, goback$turner.syndrome)
goback$overall.chrom <- ifelse(goback$studyid %in% tmp, 1, goback$overall.chrom)



# Verify overall.chrom: TX ------------------------------------------------

#' Turner.  Looks good.
tmp <- defect.codes[grepl('^tx', defect.codes$studyid), ]
tmp$turner.flag <- 0
for (i in 2:67){
  tmp$turner.flag <- ifelse(tmp[,i] %in% c(758.600, 758.610, 758.611, 758.612, 758.613, 758.614, 758.615, 758.690), 1, tmp$turner.flag)
}
table(tmp$turner.flag, useNA = 'always')

#' DiGeorge.  No additions.
tmp <- defect.codes[grepl('^tx', defect.codes$studyid), ]
tmp$dig.flag <- 0
for (i in 2:67){
  tmp$dig.flag <- ifelse(tmp[,i] %in% 758.370, 1, tmp$dig.flag)
}
table(tmp$dig.flag, useNA = 'always')
tmp <- filter(tmp, dig.flag == 1)
tmp <- c(tmp$studyid)
goback$di.george.syndrome <- ifelse(goback$studyid %in% tmp, 1, goback$di.george.syndrome)
goback$overall.chrom <- ifelse(goback$studyid %in% tmp, 1, goback$overall.chrom)

#' Down.
tmp <- defect.codes[grepl('^tx', defect.codes$studyid), ]
tmp$dig.flag <- 0
for (i in 2:67){
  tmp$dig.flag <- ifelse(tmp[,i] %in% c(758.000, 758.01, 758.02, 758.03, 758.04, 758.09), 1, tmp$dig.flag)
}
table(tmp$dig.flag, useNA = 'always')



# Verify overall.chrom: MI ------------------------------------------------

#' Down.  Have about 80 extra cases in MI.
tmp <- defect.codes[grepl('^mi', defect.codes$studyid), ]
tmp$down.flag <- 0
for (i in 2:67){
  tmp$down.flag <- ifelse(tmp[,i] %in% c(758.0,758.1,758.2,758.3,758.31, 758.32, 758.33, 758.39,758.40,
                                         758.5, 758.6, 758.7, 758.8, 758.81, 758.89, 758.9), 1, tmp$down.flag)
}
table(tmp$down.flag, useNA = 'always')
tmp <- filter(tmp, down.flag == 1)
tmp <- c(tmp$studyid)
goback$down.syndrome <- ifelse(goback$studyid %in% tmp, 1, goback$down.syndrome)

#' There are no extra codes I'm missing.
load("W:/Old_genepi2/Jeremy/GOBACK/Datasets/list.of.unique.mi.bd.codes.rdata")
tmp <- filter(codes, codes >= 758 & codes < 759)
print(tmp)

tmp <- filter(filter(goback, down.syndrome == 1), state == 'MI')
tmp <- c(tmp$studyid)
tmp <- defect.codes[defect.codes$studyid %in% tmp, ]

#' This file has a record for EVERY child in Michigan.
#' Q: Do these mismatched kids even have a defect code?
#' A: Hard to tell, still can't reconcile the IDs.
load("W:/Old_genepi2/Jeremy/GOBACK/Datasets/Michigan/mi.birthdefects.codes.all.children.rdata")
mi.filt <- arrange(mi.filt, studyid)

tmp <- mi.filt[grepl('^mi', mi.filt$studyid), ]
tmp$down.flag <- 0
for (i in 2:67){
  tmp$down.flag <- ifelse(tmp[,i] %in% 758.0, 1, tmp$down.flag)
}
table(tmp$down.flag, useNA = 'always')

#' Let's just make sure that all kids tagged as having a given chromosomal defect do have one.
#' They all do.
tmp <- filter(filter(goback, overall.chrom == 1), state == 'MI')
tmp <- c(tmp$studyid)
tmp <- defect.codes[defect.codes$studyid %in% tmp, ]
tmp$flag <- 0
for (i in 2:67){
  tmp$flag <- ifelse(tmp[,i] %in% c(758.0,758.1,758.2,758.3,758.31, 758.32, 758.33, 758.39,758.40,
                                         758.5, 758.6, 758.7, 758.8, 758.81, 758.89, 758.9), 1, tmp$flag)
}
table(tmp$flag, useNA = 'always')
# Split into chromosomal and non-chromosomal sets -------------------------
load('goback.v20171206.1.rdata')

#' Re-order cancer variables alphabetically and move the '.any' variables
#' to the end.  
#' Required for some of code in the regression models to run correctly, due
#' to the order in which the levels of the cancer1 variable are printed.
goback <- goback[,c(1:107, 
                    108, 109, 135, 116, 133, 119, 146, 144, 134, 131, 141, 142, 140,
                    129, 127, 112, 111, 115, 117, 124, 120, 113, 130, 145, 147, 122, 
                    126, 123, 137, 139, 110, 114, 118, 121, 125, 128, 132, 136, 138, 
                    143, 148:150)]

save(goback, file = 'goback.v20171211.1.rdata')

chrom <- filter(goback, any.birthdefect == 1 & overall.chrom > 0)
no.chrom <- filter(goback, any.birthdefect == 1 & overall.chrom == 0)
control <- filter(goback, any.birthdefect == 0 & overall.chrom == 0)

rm(goback)
gc()

goback.nochrom <- rbind(no.chrom, control)
save(goback.nochrom, file = 'goback.no.chrom.v20171211.1.rdata')
rm(goback.nochrom, no.chrom)
gc()

goback.chrom <- rbind(chrom, control)
save(goback.chrom, file = 'goback.chrom.v20171211.1.rdata')
rm(goback.chrom, chrom, control)
gc()



# Right-censor data at 18: non-chromosomal --------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Philip and I met with Susan Hilsenbeck on 1-10-2018 to discuss moving
#' forward with Cox PH as the model of choice, given that the PH 
#' assumption appears to be largely upheld in children under 18 years of
#' age.
#' 
#' She agreed this was reasonable.  Before we can do the modeling, however,
#' we need to recode the cancer variables such that they are 1 if the child
#' was diagnosed before 18 and 0 otherwise.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

load('goback.no.chrom.v20171211.1.rdata')

goback.nochrom$dxby18 <- ifelse(goback.nochrom$cancer == 1 & goback.nochrom$person.yrs <= 18, 1, 0)

#' Crosstabs indicate there are very few children diagnosed with cancer at > 18.00 years in
#' the non-chromosomal set.
sink(file = 'C:/Users/schraw/Desktop/cancer.dx.by.18rs.txt')
for (i in 108:137){
  print(names(goback.nochrom[i]))
  print(gmodels::CrossTable(goback.nochrom$dxby18, goback.nochrom[,i], prop.chisq = FALSE, prop.t = FALSE))
}
sink() 

#' Remove cancer cases DX'd over 18.
tmp <- filter(goback.nochrom, cancer == 1 & dxby18 == 0)
tmp <- c(tmp$studyid)
goback.nochrom <- goback.nochrom[!(goback.nochrom$studyid %in% tmp), ]

#' Set max follow-up to 18 in non-cancer kids.
goback.nochrom$person.yrs <- ifelse(goback.nochrom$cancer == 0 & goback.nochrom$person.yrs > 18, 18, goback.nochrom$person.yrs)

aggregate(person.yrs ~ cancer + state, data = goback.nochrom, max)

save(goback.nochrom, file = 'goback.no.chrom.v20180113.1.rdata')

rm(i, tmp, goback.nochrom)



# Right-censor data at 18: chromosomal ------------------------------------

load('goback.chrom.v20171211.1.rdata')

goback.chrom$dxby18 <- ifelse(goback.chrom$cancer == 1 & goback.chrom$person.yrs <= 18, 1, 0)

sink(file = 'C:/Users/schraw/Desktop/cancer.dx.by.18rs.chrom.txt')
for (i in 108:137){
  print(names(goback.chrom[i]))
  print(gmodels::CrossTable(goback.chrom$dxby18, goback.chrom[,i], prop.chisq = FALSE, prop.t = FALSE))
}
sink() 

tmp <- filter(goback.chrom, cancer == 1 & dxby18 == 0)
tmp <- c(tmp$studyid)
goback.chrom <- goback.chrom[!(goback.chrom$studyid %in% tmp), ]

#' Set max follow-up to 18 in non-cancer kids.
goback.chrom$person.yrs <- ifelse(goback.chrom$cancer == 0 & goback.chrom$person.yrs > 18, 18, goback.chrom$person.yrs)

aggregate(person.yrs ~ cancer + state, data = goback.chrom, max)

save(goback.chrom, file = 'goback.chrom.v20180113.1.rdata')

rm(i, tmp, goback.chrom)
# Fix an issue with cancer variables in NC --------------------------------

load('goback.no.chrom.v20180113.1.rdata')

tmp <- filter(goback.nochrom, state == 'NC')
goback.nochrom <- filter(goback.nochrom, state != 'NC')

for (i in 108:147){
  print(table(tmp[,i], useNA = 'always'))
}

#' Individual cancer variables are either 1 or NA.  Never 0.
for (i in 108:147){
  tmp[ ,i] <- ifelse(is.na(tmp[ ,i]) & tmp$cancer == 0, 0, tmp[, i])
}

goback.nochrom <- rbind(goback.nochrom, tmp)
rm(i, tmp)

save(goback.nochrom, file = 'goback.no.chrom.v20180117.1.rdata')
rm(goback.nochrom); gc()

load('goback.chrom.v20180113.1.rdata')

tmp <- filter(goback.chrom, state == 'NC')
goback.chrom <- filter(goback.chrom, state != 'NC')

for (i in 108:147){
  print(table(tmp[,i], useNA = 'always'))
}

#' Individual cancer varaibles are either 1 or NA.  Never 0.
for (i in 108:147){
  tmp[ ,i] <- ifelse(is.na(tmp[ ,i]) & tmp$cancer == 0, 0, tmp[, i])
}

goback.chrom <- rbind(goback.chrom, tmp)
rm(i, tmp)

save(goback.chrom, file = 'goback.chrom.v20180117.1.rdata')
rm(goback.chrom); gc()



# Update birthweight category ---------------------------------------------

#' Non-chromosomal set.
load('goback.no.chrom.v20180117.1.rdata')

codes <- data.frame(state = c('NC','TX','MI','AR'),
                    state.num = 1:4)

goback.nochrom <- filter(goback.nochrom, sex != 9)
goback.nochrom <- left_join(goback.nochrom, codes, by = 'state')

#' I realize I never updated the birthweight category variable after I 
#' recovered the original continuous AR birthweight data.
goback.nochrom$birth.wt.cat <- factor(
  ifelse(goback.nochrom$birth.wt >= 2500 & goback.nochrom$birth.wt < 4000, 0,
         ifelse(goback.nochrom$birth.wt >= 4000, 1, 2)),
  levels = c(0:2),
  labels = c('NBW','HBW','LBW'))

save(goback.nochrom, file = 'goback.no.chrom.v20180122.1.rdata')

#' Chromosomal set.
load('goback.chrom.v20180117.1.rdata')

codes <- data.frame(state = c('NC','MI','TX','AR'),
                    state.num = 1:4)

goback.chrom <- filter(goback.chrom, sex != 9)
goback.chrom <- left_join(goback.chrom, codes, by = 'state')

goback.chrom$birth.wt.cat <- factor(
  ifelse(goback.chrom$birth.wt >= 2500 & goback.chrom$birth.wt < 4000, 0,
         ifelse(goback.chrom$birth.wt >= 4000, 1, 2)),
  levels = c(0:2),
  labels = c('NBW','HBW','LBW'))

save(goback.chrom, file='goback.chrom.v20180122.1.rdata')

# Drop kids with unknown sex from combined dataset ------------------------

require(dplyr)

#' For consistency with the chromosomal and non-chromosomal sets,
#' drop children with sex == 9 from the combined dataset.
load("W:/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/goback.v20171211.1.rdata")

goback <- filter(goback, sex != 9)

save(goback, file = 'goback.v20180125.1.rdata')



# Merge updated chromosomal and non-chromosomal sets ----------------------

load("W:/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/goback.no.chrom.v20180122.1.rdata")

controls <- filter(goback.nochrom, any.birthdefect == 0)
cases.nonchrom <- goback.nochrom[goback.nochrom$any.birthdefect == 1, ]
rm(goback.nochrom); gc()

load("goback.chrom.v20180122.1.rdata")

cases.chrom <- goback.chrom[goback.chrom$any.birthdefect == 1, ]
rm(goback.chrom); gc()
cases <- rbind(cases.chrom, cases.nonchrom); rm(cases.chrom, cases.nonchrom); gc()

goback <- rbind(cases, controls)
rm(cases, controls); gc()

goback <- goback[!duplicated(goback$studyid), ]
save(goback, file = 'goback.v20180216.1.rdata')



# Remove NA values in TX non-BD kids --------------------------------------

#' Texas kids with no birth defects have NA values for defect.total.
#' Convert to zero.

#' Non-chromosomal set.
setwd('W:/Old_genepi2/Jeremy/GOBACK/Datasets/')
load('./goback.no.chrom.v20180122.1.rdata') 

table(goback.nochrom$state, goback.nochrom$defect.total, useNA = 'always')
goback.nochrom$defect.total <- ifelse(goback.nochrom$state == 'TX' & is.na(goback.nochrom$defect.total), 0, goback.nochrom$defect.total)
table(goback.nochrom$state, goback.nochrom$defect.total, useNA = 'always')

save(goback.nochrom, file = './goback.nochrom.v20180419.rdata')

rm(goback.nochrom); gc()

#' Chromosomal set.
load('./goback.chrom.v20180122.1.rdata')

table(goback.chrom$state, goback.chrom$defect.total, useNA = 'always')
goback.chrom$defect.total <- ifelse(goback.chrom$state == 'TX' & is.na(goback.chrom$defect.total), 0, goback.chrom$defect.total)
table(goback.chrom$state, goback.chrom$defect.total, useNA = 'always')

save(goback.chrom, file = './goback.chrom.v20180419.rdata')

rm(goback.chrom); gc()

#' Overall set.
load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/goback.v20180216.1.rdata')

table(goback$state, goback$defect.total, useNA = 'always')
goback$defect.total <- ifelse(goback$state == 'TX' & is.na(goback$defect.total), 0, goback$defect.total)
table(goback$state, goback$defect.total, useNA = 'always')

save(goback, file = './goback.v20180419.rdata')

rm(goback); gc()




# Remove NA values from majordefect.total ---------------------------------

load('./goback.nochrom.v20180419.rdata')

table(goback.nochrom$majordefect.total, goback.nochrom$any.birthdefect, useNA = 'ifany')
goback.nochrom$majordefect.total <- ifelse(is.na(goback.nochrom$majordefect.total) & goback.nochrom$any.birthdefect == 0, 0, goback.nochrom$majordefect.total)
table(goback.nochrom$majordefect.total, goback.nochrom$any.birthdefect, useNA = 'ifany')

save(goback.nochrom, file = './goback.nochrom.v20180507.rdata')

rm(list = ls()); gc()

load('./goback.chrom.v20180419.rdata')

table(goback.chrom$majordefect.total, goback.chrom$any.birthdefect, useNA = 'ifany')
goback.chrom$majordefect.total <- ifelse(is.na(goback.chrom$majordefect.total) & goback.chrom$any.birthdefect == 0, 0, goback.chrom$majordefect.total)
table(goback.chrom$majordefect.total, goback.chrom$any.birthdefect, useNA = 'ifany')

save(goback.chrom, file = './goback.chrom.v20180507.rdata')

rm(list = ls()); gc()

load('./goback.v20180419.rdata')

table(goback$majordefect.total, goback$any.birthdefect, useNA = 'ifany')
goback$majordefect.total <- ifelse(is.na(goback$majordefect.total) & goback$any.birthdefect == 0, 0, goback$majordefect.total)
table(goback$majordefect.total, goback$any.birthdefect, useNA = 'ifany')

save(goback, file = './goback.v20180507.rdata')

rm(list = ls()); gc()



# Replace 9999 with NA for plurality --------------------------------------

#' Actually, I suspect 9 is also a missing value code in one or more states.
load('./goback.nochrom.v20180507.rdata')

table(goback.nochrom$plu, useNA = 'ifany')
goback.nochrom$plu <- ifelse(goback.nochrom$plu == 9999 | goback.nochrom$plu == 9, NA, goback.nochrom$plu)
table(goback.nochrom$plu, useNA = 'ifany')

save(goback.nochrom, file = './goback.nochrom.v20180517.rdata')

rm(list = ls()); gc()

load('./goback.chrom.v20180507.rdata')

table(goback.chrom$plu, useNA = 'ifany')
goback.chrom$plu <- ifelse(goback.chrom$plu == 9999 | goback.chrom$plu == 9, NA, goback.chrom$plu)
table(goback.chrom$plu, useNA = 'ifany')

save(goback.chrom, file = './goback.chrom.v20180517.rdata')

rm(list = ls()); gc()

load('./goback.v20180507.rdata')

table(goback$plu, useNA = 'ifany')
goback$plu <- ifelse(goback$plu == 9999 | goback$plu == 9, NA, goback$plu)
table(goback$plu, useNA = 'ifany')

save(goback, file = './goback.v20180517.rdata')

rm(list = ls()); gc()



# Replace NA values for NC cancer variables -------------------------------

#' NC Children diagnosed with cancer have an NA value for any specific cancer
#' variable if they were not diagnosed with THAT cancer.  Replace with zeroes.

load('./goback.v20180517.rdata')

for (i in 108:147){
  goback[,i] <- ifelse(goback$state == 'NC' & is.na(goback[,i]), 0, goback[,i])
}

save(goback, file = './goback.v20180530.1.rdata')

rm(list = ls()); gc()

load('./goback.nochrom.v20180517.rdata')

for (i in 108:147){
  goback.nochrom[,i] <- ifelse(goback.nochrom$state == 'NC' & is.na(goback.nochrom[,i]), 0, goback.nochrom[,i])
}

save(goback.nochrom, file = './goback.nochrom.v20180530.1.rdata')

rm(list = ls()); gc()

load('./goback.chrom.v20180517.rdata')

for (i in 108:147){
  goback.chrom[,i] <- ifelse(goback.chrom$state == 'NC' & is.na(goback.chrom[,i]), 0, goback.chrom[,i])
}

save(goback.chrom, file = './goback.chrom.v20180530.1.rdata')

rm(list = ls()); gc()



# Remove laterality values in children with no cancer ---------------------

#' A number of TX children without any cancer have values for laterality.
#' I can't find any record of these kids in the cancer registry data.
#' I am prepared to say we should just remove the values and replace with
#' NA.

load('./goback.v20180530.1.rdata')

goback$laterality1 <- ifelse(goback$state == 'TX' & goback$cancer == 0, NA, goback$laterality1)

save(goback, file = './goback.v20180530.2.rdata')

rm(list = ls()); gc()

load('./goback.nochrom.v20180530.1.rdata')

goback.nochrom$laterality1 <- ifelse(goback.nochrom$state == 'TX' & goback.nochrom$cancer == 0, NA, goback.nochrom$laterality1)

save(goback.nochrom, file = './goback.nochrom.v20180530.2.rdata')

rm(list = ls()); gc()

load('./goback.chrom.v20180530.1.rdata')

goback.chrom$laterality1 <- ifelse(goback.chrom$state == 'TX' & goback.chrom$cancer == 0, NA, goback.chrom$laterality1)

save(goback.chrom, file = './goback.chrom.v20180530.2.rdata')

rm(list = ls()); gc()



# Flag syndromic defects in the non-chromosomal set: generate list --------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.06.06.
#' 
#' I realized where performing analyses for BTEC that children with 
#' certain syndromes are still included in the non-chromosomal set, 
#' because the filtering that was performed was done to remove codes in the 
#' 758 range.
#' 
#' We want to identify kids with certain codes:
#' BPA 237.700, ICD9 237.7: neurofibromatosis
#' BPA 759.500, ICD9 759.5: tuberous sclerosis
#' BPA 759.800: possibly Costello or Noonan
#' BPA 279.110: phenotypic DX of DiGeorge (though these should already be 
#'              excluded)
#' BPA 759.840: possibly Rubinstein-Taybi
#' BPA 759.870: possibly Beckwith-Wiedemann 
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(dplyr)

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.v20180606.rdata')
load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.v20180606.rdata')

tsc.txnc <- filter(bd.codes.txnc, bpa1 == '759.500')
for (i in 3:67){
  tmp <- filter(bd.codes.txnc, bd.codes.txnc[, i] == '759.900')
  tsc.txnc <- rbind(tsc.txnc, tmp)
}
tsc.txnc <- tsc.txnc[!duplicated(tsc.txnc$studyid), ]
tsc.mi <- filter(bd.codes.mi, icd9cod1 == '759.5')
for (i in 3:25){
  tmp <- filter(bd.codes.mi, bd.codes.mi[, i] == '759.9')
  tsc.mi <- rbind(tsc.mi, tmp)
}
tsc.mi <- tsc.mi[!duplicated(tsc.mi$studyid), ]
tsc.ids <- c(tsc.mi$studyid, tsc.txnc$studyid)



neurofib.txnc <- filter(bd.codes.txnc, bpa1 == '237.700')
for (i in 3:67){
  tmp <- filter(bd.codes.txnc, bd.codes.txnc[, i] == '237.700')
  neurofib.txnc <- rbind(neurofib.txnc, tmp)
}
neurofib.txnc <- neurofib.txnc[!duplicated(neurofib.txnc$studyid), ]
neurofib.mi <- filter(bd.codes.mi, icd9cod1 == '237.7')
for (i in 3:25){
  tmp <- filter(bd.codes.mi, bd.codes.mi[, i] == '237.7')
  neurofib.mi <- rbind(neurofib.mi, tmp)
}
neurofib.mi <- neurofib.mi[!duplicated(neurofib.mi$studyid), ]
neurofib.ids <- c(neurofib.mi$studyid, neurofib.txnc$studyid)



cost.noon.txnc <- filter(bd.codes.txnc, bpa1 == '759.800')
for (i in 3:67){
  tmp <- filter(bd.codes.txnc, bd.codes.txnc[, i] == '759.800')
  cost.noon.txnc <- rbind(cost.noon.txnc, tmp)
}
cost.noon.txnc <- cost.noon.txnc[!duplicated(cost.noon.txnc$studyid), ]
cost.noon.ids.txnc <- cost.noon.txnc$studyid



digeorge.txnc <- filter(bd.codes.txnc, bpa1 == '279.110')
for (i in 3:67){
  tmp <- filter(bd.codes.txnc, bd.codes.txnc[, i] == '279.110')
  digeorge.txnc <- rbind(digeorge.txnc, tmp)
}
digeorge.txnc <- digeorge.txnc[!duplicated(digeorge.txnc$studyid), ]
digeorge.mi <- filter(bd.codes.mi, icd9cod1 == '279.11')
for (i in 3:25){
  tmp <- filter(bd.codes.mi, bd.codes.mi[, i] == '279.11')
  digeorge.mi <- rbind(digeorge.mi, tmp)
}
digeorge.mi <- digeorge.mi[!duplicated(digeorge.mi$studyid), ]
digeorge.ids <- c(digeorge.mi$studyid, digeorge.txnc$studyid)



rub.tay.txnc <- filter(bd.codes.txnc, bpa1 == '759.840')
for (i in 3:67){
  tmp <- filter(bd.codes.txnc, bd.codes.txnc[, i] == '759.840')
  rub.tay.txnc <- rbind(rub.tay.txnc, tmp)
}
rub.tay.txnc <- rub.tay.txnc[!duplicated(rub.tay.txnc$studyid), ]
rub.tay.ids.txnc <- rub.tay.txnc$studyid



bws.txnc <- filter(bd.codes.txnc, bpa1 == '759.870')
for (i in 3:67){
  tmp <- filter(bd.codes.txnc, bd.codes.txnc[, i] == '759.870')
  bws.txnc <- rbind(bws.txnc, tmp)
}
bws.txnc <- bws.txnc[!duplicated(bws.txnc$studyid), ]
bws.ids.txnc <- bws.txnc$studyid



syndromes.mi <- filter(bd.codes.mi, icd9cod1 == '759.89')
for (i in 3:25){
  tmp <- filter(bd.codes.mi, bd.codes.mi[, i] == '759.89')
  syndromes.mi <- rbind(syndromes.mi, tmp)
}
syndromes.mi <- syndromes.mi[!duplicated(syndromes.mi$studyid), ]
syndrome.ids.mi <- syndromes.mi$studyid



syndrome.ids <- list(tuberous.sclerosis = tsc.ids, 
                     neurofibromatosis = neurofib.ids, 
                     digeorge = digeorge.ids,
                     costello.noonan.txnc = cost.noon.ids.txnc, 
                     rubinstein.taybi.txnc = rub.tay.ids.txnc,
                     beckwith.wiedemann.tx.nc = bws.ids.txnc, 
                     cot.noon.rub.tay.bws.mi = syndrome.ids.mi)

save(syndrome.ids, file = 'W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/list.of.syndromic.kids.in.tx.mi.nc.rdata')



# Flag syndromic defects in the non-chromosomal set: count occurre --------

load("W:/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/goback.nochrom.v20180530.2.rdata")
load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/list.of.syndromic.kids.in.tx.mi.nc.rdata')

syndromic.birth.defects.cases <- as.data.frame(matrix(nrow = 1, ncol = 36))
df.names <- c('syndrome', 'num.cases', 'mi.cases', 'nc.cases', 'tx.cases', 'num.comorbid.cases', names(goback.nochrom[108:137]))
df.names <- df.names[c(1:29,31:36,30)]
names(syndromic.birth.defects.cases) <- df.names
syndromic.birth.defects.cases$syndrome <- as.character(syndromic.birth.defects.cases$syndrome)

for (i in 1:7){

tmp <- syndrome.ids[[i]]

index.syndrome <- names(syndrome.ids[i])

cases <- filter(goback.nochrom, studyid %in% tmp)
cases.mi <- filter(cases, state == 'MI')
cases.nc <- filter(cases, state == 'NC')
cases.tx <- filter(cases, state == 'TX')

comorbid.cases <- filter(cases, cancer == 1)

new.syndrome <- data.frame (syndrome = index.syndrome,
                            num.cases = nrow(cases),
                            mi.cases = nrow(cases.mi),
                            nc.cases = nrow(cases.nc),
                            tx.cases = nrow(cases.tx),
                            num.comorbid.cases = nrow(comorbid.cases),
                            all = nrow(filter(comorbid.cases, all == 1)),
                            aml = nrow(filter(comorbid.cases, aml == 1)),
                            arms = nrow(filter(comorbid.cases, arms == 1)),
                            astro = nrow(filter(comorbid.cases, astro == 1)),
                            bone.other = nrow(filter(comorbid.cases, bone.other == 1)),
                            cns.other = nrow(filter(comorbid.cases, cns.other == 1)),
                            ependymoma = nrow(filter(comorbid.cases, ependymoma == 1)),
                            epithe = nrow(filter(comorbid.cases, epithe == 1)),
                            erms = nrow(filter(comorbid.cases, erms == 1)),
                            ewing = nrow(filter(comorbid.cases, ewing == 1)),
                            gct.extra = nrow(filter(comorbid.cases, gct.extra == 1)),
                            gct.gonad = nrow(filter(comorbid.cases, gct.gonad == 1)),
                            gct.intra = nrow(filter(comorbid.cases, gct.intra == 1)),
                            hepatic.other = nrow(filter(comorbid.cases, hepatic.other == 1)),
                            hepato = nrow(filter(comorbid.cases, hepato == 1)),
                            hl = nrow(filter(comorbid.cases, hl == 1)),
                            leu.other = nrow(filter(comorbid.cases, leu.other == 1)),
                            lym.other = nrow(filter(comorbid.cases, lym.other == 1)),
                            medullo = nrow(filter(comorbid.cases, medullo == 1)),
                            nephro = nrow(filter(comorbid.cases, nephro == 1)),
                            neuro = nrow(filter(comorbid.cases, neuro == 1)),
                            nhl = nrow(filter(comorbid.cases, nhl == 1)),
                            osteo = nrow(filter(comorbid.cases, osteo == 1)),
                            pnet = nrow(filter(comorbid.cases, pnet == 1)),
                            pns.other = nrow(filter(comorbid.cases, pns.other == 1)),
                            renal.other = nrow(filter(comorbid.cases, renal.other == 1)),
                            retino = nrow(filter(comorbid.cases, retino == 1)),
                            rms.other = nrow(filter(comorbid.cases, rms.other == 1)),
                            soft.other = nrow(filter(comorbid.cases, soft.other == 1)),
                            other.any = nrow(filter(comorbid.cases, other.any == 1)))
                            

syndromic.birth.defects.cases <- rbind(syndromic.birth.defects.cases, new.syndrome)

}

syndromic.birth.defects.cases <- syndromic.birth.defects.cases[2:8, ]

write.csv(syndromic.birth.defects.cases, file = 'W:/Old_genepi2/Jeremy/GOBACK/R outputs/Syndromic cases/syndromic.case.counts.csv', row.names = FALSE)



# Compute new variables for TSC and NF cases ------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.06.11.
#' 
#' I discovered that children with tuberous sclerosis (TSC) and 
#' neurofibromatosis (NF) were still in the non-chromosomal dataset.
#' 
#' In conversation with Philip we decided that these children should not be
#' analyzed with this group, as they have known cancer predisposition 
#' syndromes.  Thus, I will move them to the chromosomal set after 
#' creating some new variables to reflect the new data partition:
#' - any.genetic.anomaly: an indicator variable taking value 1 if the child
#'   has either a chromosomal anomaly or a single-gene syndrome.
#' - single.gene.anomaly: an indicator variable taking value 1 if the 
#'   child has a DX of either TSC or NF, and 0 otherwise.
#' - tsc: an indicator variable taking value 1 if the child has a BPA or 
#'   ICD9 code for TSC (759.900 and 759.9, respectively).
#' - nf: an indicator variable taking value 1 if the child has a BPA or 
#'   ICD9 code for neurofibromatosis (237.700 and 237.7, respectively).
#'
#' Per the usual organization, each of these variables will take value NA
#' if the child does NOT have that defect, but does have one or more other
#' defects.
#'
#' In addition to this, we will completely exclude children who have the
#' code for the clinical diagnosis of DiGeorge syndrome, but do not have
#' the code for the accompanying genetic diagnosis.
#' 
#' Lastly, I will update each of the old birth defects variables such that 
#' if they currently have value 0 and the child has TSC or NF, they are NA.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(dplyr)

load("W:/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/goback.v20180530.2.rdata")
load("W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/list.of.syndromic.kids.in.tx.mi.nc.rdata")
load("W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.v20180606.rdata")
load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.transpose.v20180614.rdata')

#' Initialize and compute new BD variables.
goback[ ,153:156] <- 0
names(goback) <- c(names(goback)[1:152],'any.genetic.anomaly','single.gene.anomaly','tsc','nf')

goback$single.gene.anomaly <- ifelse(goback$studyid %in% syndrome.ids[[1]] | goback$studyid %in% syndrome.ids[[2]], 1, goback$single.gene.anomaly)
goback$single.gene.anomaly <- ifelse(goback$any.birthdefect == 1 & goback$single.gene.anomaly == 0, NA, goback$single.gene.anomaly)

goback$any.genetic.anomaly <- ifelse(goback$single.gene.anomaly == 1 | goback$chromosomalanomalies == 1, 1, goback$any.genetic.anomaly)
goback$any.genetic.anomaly <- ifelse(goback$any.birthdefect == 1 & goback$any.genetic.anomaly == 0, NA, goback$any.genetic.anomaly)

goback$tsc <- ifelse(goback$studyid %in% syndrome.ids[[1]], 1, 0)
goback$tsc <- ifelse(goback$tsc == 0 & goback$any.birthdefect == 1, NA, goback$tsc)

goback$nf <- ifelse(goback$studyid %in% syndrome.ids[[2]], 1, 0)
goback$nf <- ifelse(goback$nf == 0 & goback$any.birthdefect == 1, NA, goback$nf)

#' Definitively identify DiGeorge syndrome cases without genetic DX.
#' This distinction can only be established in the TX and NC data.
#' Exclude TX and NC children who have only a clinical diagnosis of DiGeorge.
tmp <- filter(bd.codes.txnc.transpose, (`279.110` == 1 | `279.118` == 1) & is.na(`758.370`))

goback <- filter(goback, !(studyid %in% tmp$studyid))

goback <- goback[, c(1:94,153:156,95:152)]
names(goback)[99] <- 'any.chromosomal.anomaly'

save(goback, file = 'W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20180611.rdata')

rm(list = ls()); gc()



# Replace 99 with NA for missing paternal age in TX -----------------------

#' During the initial phase of the DS-ALL BD project, I noticed that missing paternal ages in TX may be coded 99.
load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20180611.rdata')

goback$f.age <- ifelse(goback$f.age == 99, NA, goback$f.age)

save(goback, file = 'W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20180711.rdata')

rm(list = ls()); gc()



# Repair DiGeorge and compute 13q deletion variables ----------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.08.24.
#' 
#' We have shared a draft of the manuscript with external co-authors and 
#' Sonja Rasmussen questioned the DiGeorge-retinoblastoma association,
#' raising the possibility that these were misclassified cases of 13q 
#' deletion.
#' 
#' In fact, review of the birth defects codes for comorbid DiGeorge-cancer
#' cases revealed significant problems with the coding for this variable,
#' which likely stem from 1) too many BPA codes being called DiGeorge in the
#' .do file that someone originally wrote to compute the dummy variables, and 
#' 2) these BPA codes all being mapped to the uninformative ICD9 code for
#' other autosomal deletions, 758.39.
#' 
#' Shared findings with Philip. We will repair the DiGeorge variable as 
#' follows:
#' 1) compute a variable for 22q11 deletion in TX and NC kids.  Only 
#' children with the code 758.370 should be considered affected. We will
#' call this del.22q.
#' 2) Compute a variable for DiGeorge syndrome in all states. Any TX and 
#' NC kids with either or both of the codes 758.370/279.110 are considered
#' affected, as are any MI kids with the ICD9 code 279.11. We will just 
#' have to take AR at their word that any kid they said has DiGeorge does.
#' we will call this variable digeorge.syndrome.
#' 
#' We will also compute a variable for 13q deletion in TX and NC.
#' The code for this is 758.330.  We will call the variable del.13q.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(dplyr)

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20180711.rdata')
load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.transpose.v20180614.rdata')

#' Compute del.22q in TX and NC kids. Happily, there are no synonymous NC codes to contend with.
del.22q.ids <- select(filter(bd.codes.txnc.transpose, `758.370` == 1), studyid)

goback$del.22q <- ifelse(goback$studyid %in% del.22q.ids$studyid, 1, 
                         ifelse(goback$studyid %in% del.22q.ids$studyid & goback$any.birthdefect == 0, 2, # Checking for errors.  Should be empty.
                                ifelse(!(goback$studyid %in% del.22q.ids$studyid) & goback$any.birthdefect == 0, 0, NA)))

#' Compute del.13q in TX and NC kids. Again, no NC synonyms to contend with.
del.13q.ids <- select(filter(bd.codes.txnc.transpose, `758.330` == 1), studyid)

goback$del.13q <- ifelse(goback$studyid %in% del.13q.ids$studyid, 1, 
                         ifelse(goback$studyid %in% del.13q.ids$studyid & goback$any.birthdefect == 0, 2, # Checking for errors.  Should be empty.
                                ifelse(!(goback$studyid %in% del.13q.ids$studyid) & goback$any.birthdefect == 0, 0, NA)))

#' Remove DiGeorge and other chromosomal anomalies variables, re-order remaining coumns.
goback <- goback[, c(1:103,157,158,106:156)]

save(goback, file = 'W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20180829.rdata')

rm(list = ls()); gc()



# Correct a misclassified child from MI -----------------------------------

require(dplyr)

#' A microscale fix. One MI child with an unspecified chromosomal anomaly was considered 
#' to have zero major birth defects. Update to 1. 
load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20180829.rdata')

table(goback$majordefect.total, goback$any.genetic.anomaly, useNA = 'ifany')

goback$majordefect.total <- ifelse((!is.na(goback$any.genetic.anomaly) & goback$any.genetic.anomaly == 1 & goback$majordefect.total < 1), 1, goback$majordefect.total)

save(goback, file = paste0('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v',Sys.Date(),'.rdata'))

rm(list = ls()); gc()



# Remove obsolete variables -----------------------------------------------

require(dplyr)

load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20180829.rdata")

goback <- select(goback, 
                 -cancertime1, #' Contains errors, replaced by person.yrs.
                 -dxby18) #' Redundant, what with person.yrs being around.

save(goback, file = 'W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20190318.rdata')

# Split dataset into new chromosomal and non-chromosomal sets -------------

require(dplyr)

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v')

chrom <- filter(goback, any.birthdefect == 1 & any.genetic.anomaly == 1)
no.chrom <- filter(goback, any.birthdefect == 1 & is.na(any.genetic.anomaly))
control <- filter(goback, any.birthdefect == 0)

rm(goback); gc()

goback.nochrom <- rbind(no.chrom, control)
save(goback.nochrom, file = 'W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.nochrom.v20180829.rdata')
rm(goback.nochrom, no.chrom); gc()

goback.chrom <- rbind(chrom, control)
save(goback.chrom, file = 'W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.chrom.v20180829.rdata')

rm(list = ls()); gc()





