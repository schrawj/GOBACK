
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.06.14.
#' 
#' Compute Benjamini-Hochberg FDR for Cox model results.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

# Load in model outputs, compute FDR --------------------------------------

require(dplyr)

goback.coxmodels <- read.csv(file = 'Z:/Jeremy/GOBACK/R outputs/goback.coxph.models.v20180612.csv', header = FALSE, stringsAsFactors = FALSE)
names(goback.coxmodels) <- c('defect','cancer','hr','ci.lower','ci.upper','p.val.coef','n.comorbid')

#' P-values < 1.11 * 10^-16 are listed as zero.  Impute this value.
goback.coxmodels$p.val.coef <- ifelse(goback.coxmodels$p.val.coef < 1.11e-16, 1.11e-16, goback.coxmodels$p.val.coef)

#' Compute Benjamini-Hochberg FDR.  
goback.coxmodels <- arrange(goback.coxmodels, p.val.coef) 
goback.coxmodels$j <- 1:nrow(goback.coxmodels) #' Ranks for p-values
m <- nrow(goback.coxmodels) #' Total number of tests
delta <- 0.05 #' Desired family-wise error rate
goback.coxmodels$bh.fdr <- (goback.coxmodels$j/m)*delta #' Compute value for FDR statistic at delta = 0.05
#' Find greatest rank for which p.val.coef is less than (j/m)*delta.
#' This test and all tests ranking below are considered significant.
goback.coxmodels$fdr.flag <- (goback.coxmodels$bh.fdr - goback.coxmodels$p.val.coef) > 0 

write.csv(goback.coxmodels, file = 'goback.coxph.models.with.fdr.v20180612.csv', row.names = FALSE)

save(goback.coxmodels, file = 'Z:/Jeremy/GOBACK/Datasets/Expanded datasets/goback.coxph.results.v20180612.rdata')