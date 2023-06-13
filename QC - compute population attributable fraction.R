#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2019.02.28.
#' 
#' Compute population attributable fraction where: 
#' 1. The outcome is any cancer before age 18 and the exposure is any 
#' chromosomal anomaly, single gene disorder, or structural birth defect.
#' 2. The outcome is any cancer before age 18 and the exposure is any 
#' chromosomal anomaly or single gene defect.
#' 3. The outcome is any cancer before age 18 and the exposure is any 
#' non-chromosomal structural birth defect.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20190529.2.rdata')

# Subsetting operations for parts 2 and 3 ---------------------------------

#' The code in this section is not required if all types of birth defects (chromosomal, single gene, non-chromosomal) are considered.
#' Otherwise, run either sections 1 or 2, and section 3.

#' 1: Remove children with single gene disorders.
ids <- c(subset(goback, goback$single.gene.anomaly == 1)$studyid)

#' 2: Remove children with single gene disorders or chromsomal anomalies.
ids <- c(subset(goback, goback$any.genetic.anomaly == 1)$studyid)

#' 3: In either case, run...
ids <- setdiff(goback$studyid, ids)
goback <- subset(goback, goback$studyid %in% ids)



# Compute PAF -------------------------------------------------------------

require(gmodels)

tab <- CrossTable(goback$any.birthdefect, goback$cancer)

#' Calculate cumulative incidence of cancer in exposed and unexposed children.
ci.bd <- tab$prop.row[2,2]
ci.no.bd <- tab$prop.row[1,2]

#' Attributable proportion in the exposed group.
attrib.exposed <- ((ci.bd - ci.no.bd)/ci.bd)*100

#' Compute the proportion of cases exposed.
prop.cases.exposed <- tab$t[2,2]/(tab$t[1,2]+tab$t[2,2])

#' Compute population attributable fraction.
paf <- prop.cases.exposed * attrib.exposed

rm(list = ls()); gc()



# Compute PAF in Paul Fisher's study --------------------------------------

#' From Paul Fisher's paper, including chromosomal and non-chromosomal:
paf.data <- matrix(data = c(222, 4647, 65363, 3151617), nrow=2,ncol=2)

ci.bd <- paf.data[1,1]/sum(paf.data[1,])
ci.no.bd <- paf.data[2,1]/sum(paf.data[2,])
attrib.exposed <- ((ci.bd - ci.no.bd)/ci.bd)
prop.cases.exposed <- paf.data[1,1]/sum(paf.data[,1])
paf <- prop.cases.exposed * attrib.exposed * 100
