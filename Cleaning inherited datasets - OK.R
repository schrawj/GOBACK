#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Authored: 2019.12.05.
#' 
#' Last updated: 2019.12.05.
#' 
#' JP actually cleaned the Oklahoma GOBACK data to match the main file's 
#' specifications. This script is just to run some checks on those data
#' to look for any problems.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/Oklahoma/OK.clean.jp.11132019.2.rdata")

ok.clean <- select(ok_clean, -bw.flag); rm(ok_clean)

# Combine with GOBACK data and compare variables state-by-state -----------

require(dplyr); require(ggplot2); require(descr)

load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/Oklahoma/OK.clean.jp.11132019.2.rdata")
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20191114.rdata")

goback <- bind_rows(goback, ok.clean)


goback <- bind_rows(goback, ok.clean)

rm(ok.clean); gc()

con.vars <- c('birth.wt','gest.age','birth.yr','f.age','runif','person.yrs','defect.total','majordefect.total', 'minordefect.total')

defects <- c(16, 22:106)



for (i in con.vars){
  
  print('`````````````````````')
  print(i)
  print('`````````````````````')
  
  plot <- ggplot(data = goback) +
    geom_histogram(aes(x = goback[, i]), color = 'red', fill = 'grey') + 
    labs(x = i) + 
    facet_wrap(~state, scales = 'free')
  
  print(plot)
  
  print(aggregate(goback[, i] ~ state, data = goback, summary))
  
}

out <- as.data.frame(matrix(nrow = 0, ncol = 16))

for (i in defects){
  
  tab <- CrossTable(goback$state, goback[,i], dnn = c('state',names(goback)[i]))
  
  new.out <- data.frame(defect = names(goback[i]),
                        
                        ar.count = tab$tab['AR', 2],
                        ar.percent = round(tab$prop.row['AR', 2]*100, 2),
                        ar.prev.per.ten.k = tab$tab['AR', 2]/( (tab$tab['AR', 1] + tab$tab['AR', 2])/10000 ),
                        
                        mi.count = tab$tab['MI', 2],
                        mi.percent = round(tab$prop.row['MI', 2]*100, 2),
                        mi.prev.per.ten.k = tab$tab['MI', 2]/( (tab$tab['MI', 1] + tab$tab['MI', 2])/10000 ),
                        
                        nc.count = tab$tab['NC', 2],
                        nc.percent = round(tab$prop.row['NC', 2]*100, 2),
                        nc.prev.per.ten.k = tab$tab['NC', 2]/( (tab$tab['NC', 1] + tab$tab['NC', 2])/10000 ),
                        
                        ok.count = tab$tab['OK', 2],
                        ok.percent = round(tab$prop.row['OK', 2]*100, 2),
                        ok.prev.per.ten.k = tab$tab['OK', 2]/( (tab$tab['OK', 1] + tab$tab['OK', 2])/10000 ),
                        
                        tx.count = tab$tab['TX', 2],
                        tx.percent = round(tab$prop.row['TX', 2]*100, 2),
                        tx.prev.per.ten.k = tab$tab['TX', 2]/( (tab$tab['TX', 1] + tab$tab['TX', 2])/10000 )
  )
  
  out <- rbind(out, new.out)
  
}

for (i in defects){
  
  print(table(goback$any.birthdefect, goback[, i], useNA = 'ifany'))
  
}

write.csv(out,
          file = '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/R outputs/bd.prevalence.state.by.state.comparison.v20191204.csv',
          row.names = F)

#' I previously reviewed six-number summaries and histograms for demographic variables and did not find any errors or concerns.
#' I generated a .csv file with counts, percentages, and birth prevalence rates per 10,000 for each birth defect by state,
#' which is here: \\smb-main.ad.bcm.edu\genepi2\Old_genepi2\Jeremy\GOBACK\R outputs\bd.prevalence.state.by.state.comparison.v20191204,
#' and manually flagged any that varied in OK or in other states for further review.

#' Confirm that BD variables are structured properly with respect to coding of 0, 1, and NA values.
bd.vars <- c(16, 22:107)

for (i in bd.vars){
  
  print(table(goback$state, goback[,i], useNA = 'always', dnn = c('State',names(goback)[i])))
  
}



# Check on coding of some specific birth defects --------------------------

require(stringr); require(dplyr)

load('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txncok.transpose.v20191205.rdata')

names <- names(bd.codes.ok.nc.tx.transpose)

#' Code ranges for CNS defects.
pattern <- '^74[012]'

#' Code range for anecephalus.
pattern <- '^740.0'

#' Code range for spina bifida.
pattern <- '^741'

#' Code range for hydrocephalus.
pattern <- '^742.3'

#' Code range for ToF.
pattern <- '^745.2'

#' Code range for HLHS.
pattern <- '^746.7'

columns <- which(str_detect(names, pattern))

tmp <- filter(mutate(bd.codes.ok.nc.tx.transpose, 
                     has.code = as.logical(rowSums(bd.codes.ok.nc.tx.transpose[columns] > 0, na.rm = T))),
              has.code == T)

table(substr(tmp$studyid,1,2))

#' I count ~2,800 CNS defects in the OK data whereas JP has ~1,600.
#' Anencephalus is only off by 1.



# Check cancer variables --------------------------------------------------

require(dplyr); require(stringr)

load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/iccc.codes.to.dx.mappings.v20171018.rdata")
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/Oklahoma/oklahoma.raw.data.rdata")
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/Oklahoma/OK.clean.jp.11132019.2.rdata")

ok.clean <- select(ok_clean, -bw.flag); rm(ok_clean)

#' Three legacy user-defined functions from the NC data cleaning. Use these to determine diagnosis based on ICDO3 codes.
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
                                                                                                                                                                                                             ifelse(hist %in% cancer.codes$ewing.contingent.codes2 & site %in% cancer.codes$ewing.contingent.codes2.sites, 'ewing', 'other.any'))))))))))))))))))))))))))))))
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
populate.cancer.var <- function(x,y,z){
  if(missing(z)){
    ifelse((x | y) == 1, 1, 0)
  }
  else{
    ifelse((x | y | z) == 1, 1, 0)
  }
}

ok.cancer <- mutate(filter(ok, cancer == '1'),
                    cancer = as.numeric(cancer))

ok.cancer$studyid <- paste0('ok',ok.cancer$randID)

#' For children with non-malignant brain tumors, these are recorded in the 'site/hist/beh4' columns. Exclude.
ok.cancer <- filter(ok.cancer, site1 != '')

#' Format site and histology variables as required by the functions.
ok.cancer <- mutate(ok.cancer,
                       site1 = as.numeric(str_remove(ok.cancer$site1, 'C')),
                       hist1 = as.numeric(substr(ok.cancer$hist1,1,4)))

ok.cancer$cancer1 <- compute.cancer1(ok.cancer$hist1, ok.cancer$site1)
ok.cancer$cancer1 <- compute.cancer2(ok.cancer$hist1, ok.cancer$site1, ok.cancer$cancer1)

cancers <- unique(ok.cancer$cancer1)

#' Join the new cancer variable with the OK data and replace the old cancer1 variable.
ok.clean <- left_join(ok.clean, select(ok.cancer, studyid, cancer1), by = 'studyid')
ok.clean <- select(mutate(ok.clean, cancer1 = ifelse(is.na(cancer1.y), '', cancer1.y)),
                   -cancer1.x, -cancer1.y)

#' Populate the dummy variables for each specific diagnosis based on the values of cancer1.
for (i in cancers){
  
  ok.clean[,i] <- ifelse(ok.clean$cancer1 == i, 1, 0)
  
}

#' Populate the cancer category variables based on the values of the specific cancer dummy variables.
ok.clean$leu.any <- populate.cancer.var(ok.clean$all, ok.clean$aml, ok.clean$leu.other)
ok.clean$lym.any <- populate.cancer.var(ok.clean$hl, ok.clean$nhl, ok.clean$lym.other)
ok.clean$cns.any <- populate.cancer.var(ok.clean$astro, ok.clean$medullo, ok.clean$cns.other)
ok.clean$pns.any <- populate.cancer.var(ok.clean$neuro, ok.clean$pns.other)
ok.clean$renal.any <- populate.cancer.var(ok.clean$nephro, ok.clean$renal.other)
ok.clean$hepatic.any <- populate.cancer.var(ok.clean$hepato, ok.clean$hepatic.other)
ok.clean$bone.any <- populate.cancer.var(ok.clean$osteo, ok.clean$ewing, ok.clean$bone.other)
ok.clean$rms.any <- populate.cancer.var(ok.clean$erms, ok.clean$arms, ok.clean$rms.other)
ok.clean$soft.any <- populate.cancer.var(ok.clean$rms.any, ok.clean$soft.other)
ok.clean$gct.any <- populate.cancer.var(ok.clean$gct.extra, ok.clean$gct.gonad, ok.clean$gct.intra)

save(ok.clean, file = '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/Oklahoma/ok.v20191206.rdata')

rm(list = ls()); gc()