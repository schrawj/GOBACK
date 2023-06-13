#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' A template for coding cancer diagnoses to the GOBACK specifications 
#' using typical cancer registry data.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/iccc.codes.to.dx.mappings.v20171018.rdata")
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/Oklahoma/oklahoma.raw.data.rdata")
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/Oklahoma/OK.clean.jp.11132019.2.rdata")

ok.clean <- select(ok_clean, -bw.flag); rm(ok_clean)

# Define 3 convenience functions for computing cancer diagnosis -----------

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

populate.cancer.var <- function(x,y,z){
  if(missing(z)){
    ifelse((x | y) == 1, 1, 0)
  }
  else{
    ifelse((x | y | z) == 1, 1, 0)
  }
}



# Compute cancer diagnoses in the OK data ---------------------------------

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