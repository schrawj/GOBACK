# Filter down to specific BD-CC pairs (revised 20180223) ------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.06.14.
#' 
#' Philip's algorithm for selecting top results.
#' 1. Filter out non-specific models.
#' 2. Perform B-H adjustment only on the set of specific models.
#' 
#' DO exclude septal defects.  
#' DO NOT exclude RVOT or LVOT defects.
#' DO NOT exclude non-RMS soft tissues sarcomas (i.e., 'soft.other').
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

load('Z:/Jeremy/GOBACK/Datasets/goback.cox.ph.results.v20180124.1.rdata')

goback.coxmodels <- read.csv(file = 'Z:/Jeremy/GOBACK/R outputs/goback.coxph.models.v20180612.csv', header = FALSE, stringsAsFactors = FALSE)
names(goback.coxmodels) <- c('defect','cancer','hr','ci.lower','ci.upper','p.val.coef','n.comorbid')

#' P-values < 1.11 * 10^-16 are listed as zero.  Impute this value.
goback.coxmodels$p.val.coef <- ifelse(goback.coxmodels$p.val.coef < 1.11e-16, 1.11e-16, goback.coxmodels$p.val.coef)

pat1 <- 'conganomalies.'; pat2 <- 'any.chromosomal.anomaly' 
pat3 <- '.any'; pat4 <- '.other' 
pat5 <- 'septal.defects'; pat6 <- 'soft.other'

goback.coxmodels <- goback.coxmodels[!grepl(pat1, goback.coxmodels$defect), ]
goback.coxmodels <- goback.coxmodels[!grepl(pat2, goback.coxmodels$defect), ]
goback.coxmodels <- goback.coxmodels[!grepl(pat3, goback.coxmodels$cancer), ]

#' Grab soft.other models, set them aside, them bind them back in.
tmp <- goback.coxmodels[grepl(pat6, goback.coxmodels$cancer), ]
goback.coxmodels <- goback.coxmodels[!grepl(pat4, goback.coxmodels$cancer), ]
goback.coxmodels <- rbind(goback.coxmodels, tmp)
goback.coxmodels <- goback.coxmodels[!grepl(pat5, goback.coxmodels$defect), ]

#' Compute Bonferroni-Holm adjusted p-values.
goback.coxmodels <- arrange(goback.coxmodels, p.val.coef)
goback.coxmodels$rank <- 1:nrow(goback.coxmodels)
goback.coxmodels$crit.val.bonferroni <- 0.05/(nrow(goback.coxmodels)-goback.coxmodels$rank+1)
goback.coxmodels$delta <- goback.coxmodels$p.val.coef - goback.coxmodels$crit.val.bonferroni
goback.coxmodels <- filter(goback.coxmodels, delta < 0)

top.hits <- arrange(goback.coxmodels, defect, cancer)
save(top.hits, file = 'Z:/Jeremy/GOBACK/Datasets/Expanded datasets/goback.coxph.top.hits.v20180612.rdata')
write.csv(top.hits, file = 'Z:/Jeremy/GOBACK/R outputs/goback.coxph.top.hits.v20180612.csv', row.names = FALSE)

rm(list = ls()); gc()

