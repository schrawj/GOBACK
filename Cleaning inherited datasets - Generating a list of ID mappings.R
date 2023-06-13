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



# Recruiment IDs, subset of MI cases who matched --------------------------

require(dplyr); require(xlsx)

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/linked.registry.ids.v2019-02-20.2.rdata')

mi.recruit.ids <- read.xlsx(file = 'W:/Old_genepi2/Jeremy/GOBACK/Family-based cohort/mi.recruiting.matches.DLM.xlsx',
                            stringsAsFactors = F, sheetIndex = 1)

#' Checked and confirmed that there are no duplicated study ID's among the 549 who we matched.
#' Only need to filter on whether missing studyid.
mi.recruit.ids <- subset(mi.recruit.ids, !is.na(mi.recruit.ids$studyid))
mi.recruit.ids <- rename(mi.recruit.ids, recruitment.id = recruiting.id)

#' Left join MI IDs to GOBACK, copy recruitment IDs to existing variable, drop duplicated one.
goback.ids <- left_join(goback.ids, 
                        select(mi.recruit.ids, studyid, recruitment.id),
                        by = 'studyid')
goback.ids$recruitment.id.x <- ifelse(goback.ids$state == 'MI' & is.na(goback.ids$recruitment.id.x), goback.ids$recruitment.id.y, goback.ids$recruitment.id.x)
goback.ids <- rename(select(goback.ids, -recruitment.id.y), recruitment.id = recruitment.id.x)

save(goback.ids, file = 'W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/linked.registry.ids.v20190319.rdata')



# Add birth defects and cancer status -------------------------------------

require(dplyr)

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/linked.registry.ids.v2019-02-20.2.rdata')
load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20180829.rdata')

goback.ids <- left_join(goback.ids, 
                        select(goback, studyid, any.birthdefect, cancer),
                        by = 'studyid')
goback.ids <- rename(goback.ids, birth.defect = any.birthdefect)

save(goback.ids, file = paste0('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/linked.registry.ids.v',Sys.Date(),'.3.rdata'))

recruiting <- subset(goback.ids, !is.na(goback.ids$recruitment.id))

write.csv(recruiting, file = 'W:/Old_genepi2/Jeremy/GOBACK/Family-based cohort/goback.studyid.to.recruiting.id.mapping.csv', row.names = F)

rm(list = ls()); gc()



# Generate list of matched comorbid cases ---------------------------------

require(dplyr)

#' Went back to the raw data for Texas. Think there may have been mismatches.
load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/linked.registry.ids.v20190319.rdata')

family.candidates <- read.csv(file = 'W:/Old_genepi2/Jeremy/GOBACK/Family-based cohort/tx.contact.and.dx.info.w.maternal.education.v20190130.csv',
                              stringsAsFactors = F)
demo <- as.data.frame(haven::read_dta(file = 'C:/Users/schraw/Downloads/demo.dta'))
demo$birthID <- paste0('tx',demo$birthID)

tcr <- as.data.frame(haven::read_dta(file = 'C:/Users/schraw/Downloads/tcr_17_007_july20 (STATA).dta'))
tcr$birthid <- paste0('tx',tcr$birthid)

goback.to.bdr <- left_join(select(demo, CASE_ID, birthID),
                           goback.ids,
                           by = c('birthID' = 'studyid'))

goback.to.tcr <- left_join(select(tcr, patientid, birthid),
                           goback.ids,
                           by = c('birthid' = 'studyid'))

id.mapping <- left_join(select(goback.to.tcr, birthid, patientid),
                        select(goback.to.bdr, birthID, CASE_ID),
                        by = c('birthid' = 'birthID'))
id.mapping <- subset(id.mapping, !is.na(id.mapping$patientid) & !is.na(id.mapping$CASE_ID))
id.mapping <- subset(id.mapping, !duplicated(id.mapping$birthid))
id.mapping <- subset(id.mapping, id.mapping$patientid %in% family.candidates$patientid)
id.mapping$state <- 'TX'
id.mapping <- left_join(id.mapping,
                        select(family.candidates, patientid, bc_link),
                        by = 'patientid')
id.mapping <- rename(id.mapping, studyid = birthid, bd.registry.id = CASE_ID, cancer.registry.id = patientid, recruitment.id = bc_link)
id.mapping <- id.mapping[, c(1,4,3,2,5)]

#' Generate a list of IDs for comorbid cases. 
comorbid.ids <- subset(goback.ids, !is.na(goback.ids$recruitment.id) & goback.ids$state != 'TX')
comorbid.ids <- rbind(id.mapping, comorbid.ids)

save(comorbid.ids, file = 'W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/linked.ids.for.family.cohort.v20190319.rdata')

#' Update TX comorbid cases in goback.ids file.
keep.ids <- setdiff(goback.ids$studyid, comorbid.ids$studyid)
goback.ids <- subset(goback.ids, goback.ids$studyid %in% keep.ids)
goback.ids <- rbind(goback.ids, comorbid.ids)
goback.ids <- subset(goback.ids, !duplicated(goback.ids$studyid))
goback.ids <- arrange(goback.ids, studyid)

save(goback.ids, file = 'W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded Datasets/linked.registry.ids.v20190319.2.rdata')

#' Pull maternal education, age, race-ethnicity from GOBACK and send to Dani.
rm(demo, family.candidates, goback.to.bdr, goback.to.tcr, id.mapping, tcr); gc()

load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20190318.rdata")

dani <- left_join(comorbid.ids, select(goback, studyid, m.age, m.race, m.edu2, f.age), by = 'studyid')

write.csv(dani, file = 'W:/Old_genepi2/Jeremy/GOBACK/Family-based cohort/goback.studyid.to.recruitingid.mapping.csv', row.names = F)

rm(list = ls()); gc()
