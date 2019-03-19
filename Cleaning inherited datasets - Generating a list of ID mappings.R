#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2019.02.20.
#' 
#' It would be really nice if I had a data frame that mapped IDs between
#' all the various datasets I'm working with for GOBACK, such as the 
#' analytic file, the raw data, the registry-specifc identifiers, and the 
#' recruitment spreadsheets from Dani and the coordinators.
#' 
#' At some point, I started building this. I have no idea where that code 
#' is, but I'm trying to move all this into this script so it's easier to
#' find and maintain. 
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

# Recruitment IDs, NC comorbid cases --------------------------------------

require(xlsx); require(dplyr)

load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/linked.registry.ids.v20180908.rdata")

goback.ids$recruitment.id <- as.numeric(NA)

#' This has paired recruitment-analytic file IDs for comorbid NC children EXCEPT for those born in 2011.
nc.recruitment.ids <- read.xlsx(file = 'W:/Old_genepi2/Jeremy/GOBACK/Family-based cohort/nc.comorbid.cases.demographics.DLM.xlsx',
                                sheetIndex = 2, header = TRUE, stringsAsFactors = FALSE)
nc.recruitment.ids <- rename(select(nc.recruitment.ids, bcertno, ncid.AD1.AE77),
                             recruitment.id = bcertno, studyid = ncid.AD1.AE77)

#' These files have those two IDs for children born in 2011.
nc.2011.ncid <- read.csv(file = 'W:/Old_genepi2/Jeremy/GOBACK/Family-based cohort/nc.2011.children.ncid.and.codes.csv',
                         stringsAsFactors = FALSE, header = TRUE)

nc.2011.bcertno <- read.csv(file = 'W:/Old_genepi2/Jeremy/GOBACK/Family-based cohort/nc.2011.children.bcert.and.codes.csv',
                            stringsAsFactors = FALSE, header = TRUE)

names(nc.2011.ncid)[2:15] <- names(nc.2011.bcertno)[2:15]

#' Small file. Quickly review and remove incorrect possible matches.
nc.2011 <- left_join(nc.2011.bcertno, nc.2011.ncid, by = c('sex','behavior_code_icdo3','histology_type_icdo3','primary_site'))
nc.2011 <- nc.2011[c(1:5,8:15), ]
nc.2011 <- rename(select(nc.2011, bcertno, ncid),
                  recruitment.id = bcertno, studyid = ncid)

nc.recruitment.ids <- rbind(nc.recruitment.ids, nc.2011)
nc.recruitment.ids$studyid <- paste0('nc',nc.recruitment.ids$studyid)

goback.ids$recruitment.id <- ifelse(goback.ids$studyid %in% nc.recruitment.ids$studyid, nc.recruitment.ids$recruitment.id, goback.ids$recruitment.id)

save(goback.ids, file = paste0('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/linked.registry.ids.v',Sys.Date(),'.rdata'))

rm(list = ls()); gc()



# Recruitment IDs, TX comorbid cases --------------------------------------

require(dplyr)

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/linked.registry.ids.v2019-02-20.rdata')

tx.recruit.ids <- read.csv(file = 'W:/Old_genepi2/Jeremy/GOBACK/Family-based cohort/tx.contact.and.dx.info.w.maternal.education.v20190130.csv',
                           header = TRUE, stringsAsFactors = F)
tx.recruit.ids <- tx.recruit.ids[!is.na(tx.recruit.ids$patientid), ]

#' Currently works for 281 of 306 kids with birthdates post-1998.
goback.ids$recruitment.id <- ifelse(goback.ids$cancer.registry.id %in% tx.recruit.ids$patientid, tx.recruit.ids$bc_link, goback.ids$recruitment.id)

save(goback.ids, file = paste0('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/linked.registry.ids.v',Sys.Date(),'.2.rdata'))

rm(list = ls()); gc()



# Add birth defects and cancer status -------------------------------------

require(dplyr)

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/linked.registry.ids.v2019-02-20.2.rdata')
load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20180829.rdata')

goback.ids <- left_join(goback.ids, 
                        select(goback, studyid, any.birthdefect, cancer),
                        by = 'studyid')
goback.ids <- rename(goback.ids, birth.defect = any.birthdefect)

save(goback.ids, file = paste0('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/linked.registry.ids.v',Sys.Date(),'.3.rdata'))

rm(list = ls()); gc()