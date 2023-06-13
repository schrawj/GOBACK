#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Authored: 2020.01.17.
#'
#' Last updated: 2020.01.23.
#' 
#' Amanda's Oklahoma data have been cleaned and are ready to be harmonized. 
#' 
#' This new data cleaning script will pick up where the previous 4 state
#' data cleaning script left off, first by appending the OK data and then
#' as a record of edits to the resulting files.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

#' TODO: Determine whether pulmonary artery anomalies should be counted as RVOT defects.
#' TODO: Figure out why there are some children in the birth.defects.results data frame who aren't in the GOBACK data frame.

# Bind in OK data ---------------------------------------------------------

require(dplyr)

load('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20191114.rdata')
load('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/Oklahoma/ok.v20191206.rdata')

goback <- bind_rows(goback, ok.clean); rm(ok.clean); gc()

save(goback,
     file = '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20200117.rdata')

rm(list = ls()); gc()



# Insert the new birth defects variables ----------------------------------

require(dplyr); require(stringr)

load('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20200117.rdata')
load('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/computed.bd.vars.ncoktx.v20200117.rdata')

#' Birth defects variables in the GOBACK file need to be edited to match the new ones.
goback$cloacal.exstrophy <- ifelse(goback$any.birthdefect == 0, 0, NA)

#' Search variables up until about whehere cancer variables start in order to find "other unspecified" defects variables 
#' (and in order to NOT find "other cancer" variables). 
bd.names <- subset(names(goback[1:110]), str_detect(names(goback[1:110]), 'other'))

goback <- goback[, !(names(goback) %in% bd.names)]
goback <- rename(goback,
                 any.birth.defect = any.birthdefect,
                 number.defects = defect.total,
                 number.major.defects = majordefect.total,
                 number.minor.defects = minordefect.total,
                 spina.bifida.without.anencephalus = spinabifida.wo.anencephaly,
                 hydrocephalus.without.spina.bifida = hydrocephalus.wo.spinabifida,
                 reduction.deformities.of.brain = holoprosencephaly,
                 anophthalmia.microphthalmia =  anopthalmos.micropthalmos,
                 congenital.cataract = congenitalcataract,
                 anotia.microtia.or.anom.causing.hearing.impairment = anotia.microtia,
                 conganomalies.heart.circ.sys = conganomalies.heart.circsys,
                 transposition.of.great.vessels = transposition.of.greatvessels,
                 atrioventricular.septal.defect = endocardialcushiondefect,
                 hypoplastic.left.heart.syndrome = hypoplasticleftheartsyndrome,
                 interrupted.aortic.arch = interrupted.aortic.arch.type.a.or.c,
                 coarctation.of.aorta = coarctationofaorta,
                 pulm.valve.atresia.stenosis = pulmvalveatresiaandstenosis,
                 ebstein.anomaly = ebsteinanomaly,
                 pulm.artery.anomalies = pulmonaryarteryanomalies,
                 ventricular.septal.defect = ventricularseptaldefect,
                 asd.or.pfo = atrialseptaldefect,
                 single.ventricle = singleventricle,
                 tricuspid.valve.atresia.stenosis = trivalveatresiaandstenosis,
                 patent.ductus.arteriosus = patentductusarteriosis,
                 conganomalies.resp.sys = conganomalies.respsys,
                 cleft.palate.alone = cleft.palate.wo.cleft.lip,
                 cleft.lip.w.or.wo.cleft.palate = cleft.lip.w.and.wo.cleft.palate,
                 conganomalies.digestive.sys = conganomalies.digestivesystem,
                 esophageal.atresia.te.fistula = esophageal.atre.tracheofist,
                 rectal.large.instestinal.atresia.stenosis = rectal.largeintestatresia.sten,
                 hirschsprung.disease = hirshsprung.disease,
                 small.intestinal.atresia.stenosis = small.intestinal.atresia,
                 conganomalies.genitourinary = conganomalies.genitalandurinary,
                 renal.agenesis.or.hypoplasia =  renal.agenesis.hypoplasia,
                 conganomalies.musculoskeletal.sys = conganomalies.musculoskelsys,
                 limb.reduction.deformities = limb.deformities.unspecified,
                 diaphragmatic.hernia = diagphragmatic.hernia,
                 any.single.gene.anomaly = single.gene.anomaly,
                 neurofibromatosis = nf,                               
                 tuberous.sclerosis = tsc,
                 trisomy.13 = trisomy13,
                 trisomy.18  = trisomy18,
                 trisomy.21 = down.syndrome,                          
                 deletion.13q = del.13q,
                 deletion.22q  = del.22q)

#' For AR/MI, not much data manipulation is needed at the moment. Drop obsolete variables.
armi <- select(filter(goback, state %in% c('AR', 'MI')), -minor.status, -overall.chrom)

#' For NC/OK/TX, current birth defects variables will be overwritten. Note the number of rows for unaffected children.
goback <- select(goback, -(minor.status:deletion.13q), -cloacal.exstrophy)
goback <- filter(goback, state %in% c('OK', 'NC', 'TX'))

rows <- nrow(filter(goback, any.birth.defect == 0))

bd.names <- subset(names(birth.defects.results), str_detect(names(birth.defects.results), 'canonical', negate = T))

birth.defects.results <- birth.defects.results[ , names(birth.defects.results) %in% bd.names]
birth.defects.results <- birth.defects.results[!duplicated(birth.defects.results$studyid), ]
birth.defects.results <- filter(birth.defects.results, studyid != 'txNA')

#' Generate a data frame of zeroes for the unaffected children. Bind in rows for affected ones.
new.bd.data <- data.frame(matrix(nrow = rows, ncol = ncol(birth.defects.results)))
new.bd.data[,1] <- c(filter(goback, any.birth.defect == 0))$studyid

for (i in 2:ncol(new.bd.data)) { new.bd.data[, i] <- 0 }

names(new.bd.data) <- names(birth.defects.results)

new.bd.data <- bind_rows(birth.defects.results, new.bd.data)

#' Add birth defects data; combine rows with AR and MI; rearrange columns as previously.
goback <- left_join(select(goback, -any.birth.defect), new.bd.data, by = 'studyid')
goback <- bind_rows(goback, armi)
goback <- select(goback, studyid:person.yrs, runif:state.num, 
                 any.birth.defect:conganomalies.heart.circ.sys, conotruncal.defects, common.truncus:total.anomalous.pulmonary.venous.return,
                 lvot.defects, hypoplastic.left.heart.syndrome:aortic.valve.stenosis, rvot.defects, 
                 pulm.valve.atresia.stenosis:pulm.artery.anomalies, septal.defects, ventricular.septal.defect:deletion.22q,
                 cancer:gct.any)
goback <- rename(goback, m.edu = m.edu2)

save(goback, 
     file = '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20200123.rdata')

# Fix NC ethnicity variable -----------------------------------------------

require(dplyr)

load('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20200123.rdata')

goback$m.race <- as.numeric(goback$m.race)
goback$m.race <- ifelse(goback$m.race == 6 & goback$state == 'NC', 1, goback$m.race)
goback$m.race <- factor(goback$m.race, labels = c('Hispanic','NHW','NHB','Asian','Other','AIAN','Unknown'))

save(goback, file = '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20201209.rdata')

rm(list = ls()); gc()

# Correct OK cancer1 variable ---------------------------------------------

require(tidyverse)

setwd('//smb-main.ad.bcm.edu//genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/')

load('goback.v20201209.rdata')

goback %<>% mutate(cancer1 = ifelse(cancer1 == '', NA, cancer1))

saveRDS(goback, 'goback.v20210616.rds')

# Split dataset into new chromosomal and non-chromosomal sets -------------

require(tidyverse)

setwd('//smb-main.ad.bcm.edu//genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/')

goback <- readRDS('goback.v20210616.rds')

chrom <- filter(goback, any.birth.defect == 1 & any.genetic.anomaly == 1)
no.chrom <- filter(goback, any.birth.defect == 1 & is.na(any.genetic.anomaly))
control <- filter(goback, any.birth.defect == 0)

rm(goback); gc()

goback.nochrom <- bind_rows(no.chrom, control)
saveRDS(goback.nochrom, file = 'goback.nochrom.v20210616.rds')
rm(goback.nochrom, no.chrom); gc()

goback.chrom <- bind_rows(chrom, control)
saveRDS(goback.chrom, file = 'goback.chrom.v20210616.rds')

rm(list = ls()); gc()
