#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Authored: 2019.06.18.
#' 
#' Last updated: 2019.06.18.
#' 
#' Retrieve SEER summary stage at diagnosis from the GOBACK cancer registry
#' datasets. We will compare staging at diagnosis according to birth 
#' defects status in the descriptive epi ("Diagnostic Profiles") paper.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------


# Load in and save raw data -----------------------------------------------

require(haven); require(dplyr)

setwd('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Birth Defects-Childhood Cancer Projects/GOBACK Project (Tiffany)/For Jeremy/Datasets with SEER Summary Stages/')

mi.can <- read_dta(file = 'MI_cancer_062017_rg.dta')

nc.can <- filter(read_dta(file = 'NC_cancer_rg.dta'), cancer == '1')

tx.can <- read_dta(file = 'tcr_17_007_july20 (STATA).dta')

cancer.data <- list(mi.can = mi.can, nc.can = nc.can, tx.can = tx.can)

save(cancer.data, file = '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/raw.cancer.registry.data.v20190618.rdata')



# Recover SEER summary stage variables ------------------------------------

load('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/raw.cancer.registry.data.v20190618.rdata')
