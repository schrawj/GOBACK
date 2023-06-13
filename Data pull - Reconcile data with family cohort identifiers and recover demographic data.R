#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2019.01.30.
#' 
#' For TX, recover maternal education from 1996-1999 birth certificate files.
#' 
#' For NC, recover year of DX and age at DX, race-ethnicity, sex, and
#' maternal education.   
#'    
#' share with Dani by appending this to the info she sent for those cases.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

# Get Texas info ----------------------------------------------------------

require(foreign); require(xlsx); require(dplyr)

old.bcs <- read.dbf(file = 'W:/Old_genepi2/Birth Defects-Childhood Cancer Projects/GOBACK Study Documents and MISC/TCR data/birth_1996_1999_complete/birth_1996_1999_complete.dbf')
names(old.bcs) <- tolower(names(old.bcs))
old.bcs$bc_link <- as.character(old.bcs$bc_link)

dani <- read.xlsx(file = 'W:/Old_genepi2/Jeremy/GOBACK/Family-based cohort/TX contact and dx info MA.xlsx',
                  sheetName = 'Contact & Dx Info',
                  rowIndex = 1:360,
                  header = TRUE,
                  stringsAsFactors = FALSE)

names(dani) <- tolower(names(dani))

dani <- left_join(dani, select(old.bcs, bc_link, bs_bday, b_month, bs_byear, b_m_educ), by = 'bc_link')

write.csv(dani, file = 'W:/Old_genepi2/Jeremy/GOBACK/Family-based cohort/tx.contact.and.dx.info.w.maternal.education.v20190130.csv',
          row.names = FALSE)

rm(list = ls()); gc()




# Get NC info -------------------------------------------------------------

require(haven); require(xlsx); require(dplyr)

dani <- read.xlsx(file = 'W:/Old_genepi2/Jeremy/GOBACK/Family-based cohort/RAW_NC contact and dx info.xlsx',
                  header = TRUE, sheetName = 'Dx', stringsAsFactors = FALSE) #' Downloaded from the 'GOBACK Paper' Box folder on 1/30/2019.

nc <- read_sas(data_file = "W:/Old_genepi2/Birth Defects-Childhood Cancer Projects/GOBACK/North Carolina Data/ALSF files from NC to BCM/Comorbid patients eligible for follow-up/nc_comordid_forbcm.sas7bdat")
nc.linked <- read_sas(data_file = 'W:/Old_genepi2/Birth Defects-Childhood Cancer Projects/GOBACK/North Carolina Data/ALSF files from NC to BCM/Association analyses/nc_linked_forbcm.sas7bdat')

names(nc) <- tolower(names(nc))
names(nc.linked) <- tolower(names(nc.linked))

request <- select(nc, bcertno, date_of_diagnosis, age_at_diagnosis, racegroup, sex_ccr, meduc)
request$year.of.diagnosis <- substr(request$date_of_diagnosis, 1, 4)

lev <- c(1:9,NA)
lab <- c('8th grade or less','9th-12th grade no diploma','HS diploma or GED','Some college',"Associate's degree",
         "Bachelor's degree","Master's","Doctorate or professional",'Unknown')

request$meduc <- factor(request$meduc,
                           levels = lev,
                           labels = lab)

#' It's never easy. The meduc variable in this dataset is only populated for kids born in 2011. 
#' In some other datasets from NC there were two maternal education variables, presumably because they changed the way it was recorded.
#' Let's find Dani's comorbid cases for whom meduc was not recovered, and try to match them to other datasets where possible.

search.space <- subset(nc.linked, nc.linked$comorbid_flag == 1) #' The total population of comorbid cases in NC.
search.space$yearbth <- as.numeric(search.space$yearbth)
search.space <- subset(search.space, search.space$yearbth != 2011) #' Kids born this year ought not to be candidates.
search.space <- rename(search.space, birth.year = yearbth)

unmet.request <- data.frame(subset(request, is.na(request$meduc))) #' All of Dani's kids who are still missing maternal education.
unmet.request <- rename(unmet.request, sex = sex_ccr)
unmet.request <- left_join(unmet.request, select(nc, bcertno, mage, plural, dx1:dx18)) #' Some additional covariates that will hopefully help us locate matches.
unmet.request$sex <- as.numeric(unmet.request$sex)
unmet.request$mage <- as.numeric(unmet.request$mage)
unmet.request$birth.year <- as.numeric(paste0('20',substr(unmet.request$bcertno,1,2))) 

#' Join by race, birth year, sex, and maternal age. Almost sufficient for 1:1.
#' 11 duplicates can hopefully be resolved by review of BD codes.
unmet.request <- as.data.frame(
                                left_join(unmet.request,
                                          select(search.space, ncid, birth.year, meduc, sex, mage, racegroup, dx1:dx50),
                                          by = c('racegroup','birth.year','sex','mage')))

names(unmet.request)[11:28] <- paste0(rep('dx',18),1:18,'.dani')
names(unmet.request)[31:80] <- paste0(rep('dx',50),1:50,'.goback')

met.request <- subset(request, !is.na(request$meduc))
met.request <- rename(select(met.request, -date_of_diagnosis),
                      race.ethnicity = racegroup, 
                      sex = sex_ccr, 
                      maternal.education = meduc)

unmet.request <- rename(select(unmet.request, -date_of_diagnosis, -meduc.x),
                        maternal.education = meduc.y,
                        race.ethnicity = racegroup)

#' The children born in 2011 for whom I could recover meduc by matching exactly on the bcertno ID variable.
write.xlsx(met.request, 
           file = 'W:/Old_genepi2/Jeremy/GOBACK/Family-based cohort/nc.comorbid.cases.demographics.v20190130.xlsx',
           sheetName = 'Exact matches')

#' The children for whom I had to match by demographic factors and birth defects diagnoses.
write.xlsx(unmet.request, file = 'W:/Old_genepi2/Jeremy/GOBACK/Family-based cohort/nc.comorbid.cases.demographics.v20190130.xlsx',
           sheetName = 'Probabilistic matches', append = TRUE)

rm(list = ls()); gc()



# Get MI info -------------------------------------------------------------

require(xlsx); require(dplyr)

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20180829.rdata')
load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.v20180712.rdata')
load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/cancer.codes.v20180227.1.rdata')

goback <- filter(goback, cancer == 1 & any.birthdefect == 1 & state == 'MI'); gc()

cancer.codes <- subset(cancer.codes, grepl('mi', substr(cancer.codes$studyid,1,2)) == TRUE)

#' The cancer and birth defects code for comorbid cases, hopefully.
#' With covariates from GOBACK.
codes.mi <- left_join(filter(bd.codes.mi, studyid %in% cancer.codes$studyid),
                         cancer.codes,
                         by = 'studyid')

codes.mi <- left_join(codes.mi,
                      select(goback, studyid, sex, birth.yr, m.edu2, m.race, m.age, person.yrs),
                      by = 'studyid')

#' Both files downloaded from Dani's 'GOBACK Paper' Box folder on 2/1/2019.
dani.mi.diagnostics <- read.xlsx('W:/Old_genepi2/Jeremy/GOBACK/Family-based cohort/MI Recruiting with diagnositcs.xlsx',
                                 sheetIndex = 1, header = TRUE, stringsAsFactors = FALSE)
dani.mi.recruitment <- read.xlsx('W:/Old_genepi2/Jeremy/GOBACK/Family-based cohort/MI Recruitment summary file 4-19-2018.xlsx',
                                 sheetIndex = 1, colIndex = 1:33, header = TRUE, stringsAsFactors = FALSE)

names(dani.mi.diagnostics) <- tolower(names(dani.mi.diagnostics))
names(dani.mi.recruitment) <- tolower(names(dani.mi.recruitment))
dani.mi.recruitment$recruiting.id <- as.character(dani.mi.recruitment$recruiting.id)

#' Number of rows do not match, no duplicates in either file. Join and hope for the best?
dani <- left_join(dani.mi.recruitment, dani.mi.diagnostics, by = c("recruiting.id" = 'recruiting'))

#' Try to sync up variable types and formats.
codes.mi$site_code1 <- paste0('C',codes.mi$site_code1)
codes.mi$birth.yr <- as.character(codes.mi$birth.yr)
codes.mi$sex <- as.character(codes.mi$sex)
codes.mi$morph31 <- as.character(codes.mi$morph31)
codes.mi$behavior1 <- as.character(codes.mi$behavior1)
codes.mi <- rename(codes.mi, cell.behavior = behavior1, icd.o.iii.cell.type = morph31, icd.o.iii.site = site_code1)
codes.mi$hispanic <- ifelse(codes.mi$m.race == 'Hispanic', 'Y', "N")
codes.mi$m.age.group <- with(codes.mi, ifelse(m.age < 20, '1',
                                       ifelse(m.age %in% 20:24, '2', 
                                       ifelse(m.age %in% 25:29, '3', 
                                       ifelse(m.age %in% 30:34, '4',
                                       ifelse(m.age >= 35, '5', 'NA'))))))
codes.mi$m.race.dani <- with(codes.mi, ifelse(m.race %in% c('Hispanic','NHW'), '1',
                                       ifelse(m.race == 'NHB', '2',
                                       ifelse(m.race == 'Unknown', '9', '3'))))
names(codes.mi)[2:25] <- paste0(names(codes.mi)[2:25], '.goback')

dani <- rename(dani, hispanic = hispanic.y.yes.n.no.u.unknown, 
               m.race = momther.s.race.1.white.2.black.3.other.9.unknown,
               m.age.group = mother.s.age.group.1..20.2.20.24.3.25.29.4.30.34.5.35.,
               birth.yr = bxyear)
names(dani)[38:61] <- paste0(names(dani)[38:61],'.dani')

#' Recover year of cancer DX if possible.
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/Michigan/mi.raw.data.v20170818.rdata")

mi <- select(mi, studyid, cancer_yr1)
mi$studyid <- paste0('mi',mi$studyid)

codes.mi <- left_join(codes.mi,mi, by = 'studyid')
codes.mi <- rename(codes.mi, cancer.diagnosis.year = cancer_yr1)
codes.mi$cancer.diagnosis.year <- as.character(codes.mi$cancer.diagnosis.year)

#' Careful with maternal race. It may not match exactly.
dani.matches <- left_join(dani, 
                          codes.mi, 
                          by = c('sex',
                                 'icd.o.iii.site',
                                 'cell.behavior',
                                 'icd.o.iii.cell.type',
                                 'm.age.group',
                                 'birth.yr',
                                 'cancer.diagnosis.year'))
dani.matches <- rename(select(dani.matches, - m.race.dani), 
                       m.race.dani = m.race.x,
                       hispanic.dani = hispanic.x,
                       hispanic.goback = hispanic.y,
                       m.race.goback = m.race.y,
                       maternal.education = m.edu2)
dani.matches$maternal.education <- with(dani.matches, ifelse(maternal.education == 1, 'Less than HS diploma',
                                                      ifelse(maternal.education == 2, 'HS diploma',
                                                      ifelse(maternal.education == 3, 'Greater than HS diploma', NA))))

write.xlsx(dani.matches, file = 'W:/Old_genepi2/Jeremy/GOBACK/Family-based cohort/mi.recruiting.matches.v20190201.xlsx', 
           row.names = FALSE, showNA = FALSE)



# Revisit NC: Recover NCID for 2011-born children -------------------------

#' Dani wants NCID for the 13 exact matches. Don't have it. Can get it. 
require(haven); require(xlsx); require(dplyr)

dani <- read.xlsx(file = 'W:/Old_genepi2/Jeremy/GOBACK/Family-based cohort/RAW_NC contact and dx info.xlsx',
                  header = TRUE, sheetName = 'Dx', stringsAsFactors = FALSE) #' Downloaded from the 'GOBACK Paper' Box folder on 1/30/2019.

nc <- read_sas(data_file = "W:/Old_genepi2/Birth Defects-Childhood Cancer Projects/GOBACK/North Carolina Data/ALSF files from NC to BCM/Comorbid patients eligible for follow-up/nc_comordid_forbcm.sas7bdat")
nc.linked <- read_sas(data_file = 'W:/Old_genepi2/Birth Defects-Childhood Cancer Projects/GOBACK/North Carolina Data/ALSF files from NC to BCM/Association analyses/nc_linked_forbcm.sas7bdat')

names(nc) <- tolower(names(nc))
names(nc.linked) <- tolower(names(nc.linked))

search.space.2011 <- subset(nc.linked, nc.linked$comorbid_flag == 1 & search.space.2011$yearbth == 2011) #' The total population of comorbid cases in NC.
search.space.2011$yearbth <- as.numeric(search.space.2011$yearbth)
search.space.2011$sex <- as.character(search.space.2011$sex)

request <- nc
request$year.of.diagnosis <- substr(request$date_of_diagnosis, 1, 4)
request <- subset(request, request$yob == '2011')

lev <- c(1:9,NA)
lab <- c('8th grade or less','9th-12th grade no diploma','HS diploma or GED','Some college',"Associate's degree",
         "Bachelor's degree","Master's","Doctorate or professional",'Unknown')

request$meduc <- factor(request$meduc,
                        levels = lev,
                        labels = lab)
request$monbth <- substr(request$date_of_birth, 5, 6)
request <- dplyr::rename(request, primary_site_1 = primary_site, histologic_type_icdo3_1 = histologic_type_icdo3)


request <- left_join(select(request, bcertno, sex, primary_site_1, histologic_type_icdo3_1),
                     select(search.space.2011, ncid, sex, primary_site_1, histologic_type_icdo3_1),
                     by = c('sex','primary_site_1','histologic_type_icdo3_1'))

t(head(request))
t(head(search.space.2011))



# Revisit MI --------------------------------------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2019.02.19.
#' 
#' I was not able to match all Dani's MI kids during the first pass.
#' 
#' Revisit the MI data, taking a closer look at the comorbid cases who did 
#' not match the first time. 
#' 
#' The previous match required exact match on sex, ICDO3 cell type, 
#' behavior, and site, birth year, cancer diagnosis year, and maternal age
#' as a categorical variable.
#' 
#' I manually reviewed codes for a few instances where my info didn't 
#' match hers. Of the first two I found, one was mismatched on ICDO3 
#' histology and the other on site.
#' 
#' This leaves me with a couple options. We could try matching based only on sex, birth 
#' and diagnosis years, and maternal age group. This would be easiest.
#' Second, we could try ranking potential matches based on how many of their ICDO3 codes match.
#' Third, we could try cleaning up the birth defects codes, or at least matching based on number 
#' of birth defects. The order they're listed in probably the order I'd be inclined to try them
#' in.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(xlsx); require(dplyr); require(stringr)

#' The most up-to-date version of the MI family cohort info. Has studyid for the kids we could match the first time.
dani.round.two <- read.xlsx(file = 'W:/Old_genepi2/Jeremy/GOBACK/Family-based cohort/mi.recruiting.matches.DLM.xlsx',
                            sheetIndex = 1, stringsAsFactors = FALSE)

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20180829.rdata')
load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.v20180712.rdata')
load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/cancer.codes.v20180227.1.rdata')
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/Michigan/mi.raw.data.v20170818.rdata")

matched.ids <- c(subset(dani.round.two$studyid, !is.na(dani.round.two$studyid)))

#' Filter down to comorbid MI cases, then to filter both datasets to the unmatched kids.                 
goback <- filter(goback, cancer == 1 & any.birthdefect == 1 & state == 'MI'); gc()
goback <- goback[!goback$studyid %in% matched.ids, ]
goback$m.age.group <- with(goback, ifelse(m.age < 20, '1',
                                   ifelse(m.age %in% 20:24, '2', 
                                   ifelse(m.age %in% 25:29, '3', 
                                   ifelse(m.age %in% 30:34, '4',
                                   ifelse(m.age >= 35, '5', 'NA'))))))
goback$sex <- as.character(goback$sex)
goback$birth.yr <- as.character(goback$birth.yr)

unmatched.ids <- c(goback$studyid)

dani.round.two <- dani.round.two[is.na(dani.round.two$studyid), ]
dani.round.two <- dani.round.two[ , 1:61]

bd.codes.mi <- bd.codes.mi[bd.codes.mi$studyid %in% unmatched.ids, ]

cancer.codes <- subset(cancer.codes, cancer.codes$studyid %in% unmatched.ids)

mi$studyid <- paste0('mi',mi$studyid)
mi <- filter(select(mi, studyid, cancer_yr1), studyid %in% unmatched.ids)
mi <- rename(mi, yeardx = cancer_yr1)

codes.mi <- left_join(bd.codes.mi, select(goback, studyid, sex, birth.yr, m.age.group), by = 'studyid')
codes.mi <- left_join(codes.mi, mi, by = 'studyid')
codes.mi <- left_join(codes.mi, select(cancer.codes, studyid, morph31, site_code1, behavior1), by = 'studyid')
codes.mi$yeardx <- as.character(codes.mi$yeardx)


#' Choose one of the following two approaches.

#' One: match by sex, birth year, and maternal age group.
matches <- left_join(dani.round.two, codes.mi, by = c('sex','birth.yr', 'm.age.group'))

#' Two: match by sex, birth year, and year of cancer diagnosis.
#' dani.round.two <- rename(select(dani.round.two, -yeardx), yeardx = cancer.diagnosis.year)
#' matches <- left_join(dani.round.two, codes.mi, by = c('sex','birth.yr','yeardx'))


#rm(bd.codes.mi, cancer.codes, codes.mi, dani.round.two, goback, mi, matched.ids, unmatched.ids); gc()

matches$site_code1 <- paste0('C',matches$site_code1)
matches$behavior1 <- as.character(matches$behavior1)

#' Find first and last column in the icd9cod[X].dani sequence. Column indices vary based on matching method.
first <- grep('icd9cod.+dani', names(matches))[1]
last <- rev(grep('icd9cod.+dani', names(matches)))[1]

for (i in first:last){
  matches[,i] <- ifelse(matches[, i] == "", NA, matches[,i])
}

matches$diff.num.defects <- abs(rowSums(!is.na(matches[,38:61])) - rowSums(!is.na(matches[,63:86])))

matches$morph.match <- as.numeric(NA)
matches$behavior.match <- as.numeric(NA)
matches$site.match <- as.numeric(NA)

for (i in 1:nrow(matches)){
  
  matches$morph.match[i] <- as.numeric(identical(matches$morph31[i], matches$icd.o.iii.cell.type[i]))
  matches$behavior.match[i] <- as.numeric(identical(matches$cell.behavior[i], matches$behavior1[i]))
  matches$site.match[i] <- as.numeric(identical(matches$icd.o.iii.site[i], matches$site_code1[i]))
  
}

#' Find first and last in the ICDO3 code match true/false columns. Column indices vary based on matching method.
first <- grep('morph.match', names(matches))
last <- grep('site.match', names(matches))

matches$icd3o.match.count <- rowSums(matches[first:last])

matches <- arrange(matches, recruiting.id, diff.num.defects, desc(icd3o.match.count))

#' Path to the file for potential matches based on sex, year of birth, and maternal age.
#' A low success rate for matching,
write.xlsx(matches, file = 'W:/Old_genepi2/Jeremy/GOBACK/Family-based cohort/mi.recruiting.matches.round.two.xlsx', row.names = FALSE, showNA = FALSE)

#' Path to the file for potential matches based on sex, year of birth, and year of diagnosis.
#' Unfortunately, this approach is not yielding matches in instances where the first approached also failed.
write.xlsx(matches, file = 'W:/Old_genepi2/Jeremy/GOBACK/Family-based cohort/mi.recruiting.matches.round.two.matched.by.yeardx.xlsx', row.names = F, showNA = F)

rm(list = ls()); gc()





# Get Texas info, updated recruitment sheet -------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2019.03.27.
#' 
#' Dani found and sent some newer recruitment sheet for TX cases.
#' 
#' It looks better organized, has both patientid and case_ID, and has about
#' 100 more cases. 
#' 
#' After joining these files based on cancer registry ID, i did some
#' QC:
#' - 765 of the 847 rows were matched to a record in my GOBACK file.
#' - Birth year from the recruiting sheet matched well with birth year and
#'   person years from my data.
#' - ICDO3 codes matched well.
#' - Child's race-eth from the recruiting sheet and maternal race-eth from 
#'   my file matched pretty well.
#'   
#' These are good matches and it looks like a better success rate than
#' previously. 
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(dplyr); require(xlsx); require(stringr)

load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/linked.registry.ids.v20190319.2.rdata")

recruit <- read.xlsx(file = 'W:/Old_genepi2/Jeremy/GOBACK/Family-based cohort/TX REGISTRY DATA UPDATED 15-008_forMA.xlsx',
                     sheetName = 'Sheet1', stringsAsFactors = F)

tx.ids <- select(filter(goback.ids, state == 'TX'), -recruitment.id)
goback.ids <- filter(goback.ids, state != 'TX')

#' recruit$case_id should correspond to BD registry id, 
#' recruit$patientid should correspond to cancer registry id.
tx.ids <- left_join(tx.ids, 
                    select(recruit, patientid, case_id),
                    by = c('cancer.registry.id' = 'patientid'))

tx.ids$recruitment.id <- ifelse(!is.na(tx.ids$case_id), tx.ids$cancer.registry.id, as.character(NA))
tx.ids <- select(tx.ids, -case_id)

goback.ids <- rbind(tx.ids, goback.ids)

save(goback.ids, file = 'W:/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/linked.registry.ids.v20190327.rdata')

rm(list = ls()); gc()
