#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Dani wants maternal education, race, and age at DX for MI children and 
#' maternal education for TX children in the family-based cohort.
#' 
#' TX is easy but in MI the identifiers don't overlap between her set and
#' mine. Code below is pulls this info for TX, and for the MI kids that 
#' were relatively easy to map based on their birth defects and cancer 
#' codes.
#' 
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------


require(xlsx); require(dplyr); require(stringr)


# Read in and clean files Dani sent ---------------------------------------

mi.demo <- read.xlsx(file = 'W:/Old_genepi2/Jeremy/GOBACK/Manuscripts/Recruitment Paper/MI Recruitment summary file 4-19-2018.xlsx',
                     sheetIndex = 1, colIndex = 1:29, stringsAsFactors = FALSE, header = TRUE)
names(mi.demo) <- tolower(names(mi.demo))
mi.demo$recruiting.id <- as.character(mi.demo$recruiting.id)

mi.codes <- read.xlsx(file = 'W:/Old_genepi2/Jeremy/GOBACK/Manuscripts/Recruitment Paper/MI Recruiting with diagnositcs to philip.xlsx',
                      sheetIndex = 1, stringsAsFactors = FALSE, header = TRUE)

dani <- left_join(mi.codes, mi.demo, by = c('recruiting' = 'recruiting.id'))

#' Remove leading and trailing zeroes.
for (i in 6:29){
  dani[,i] <- str_remove(dani[,i], '[0]+$')
  dani[,i] <- str_remove(dani[,i], '^[0]+')
}

#' Sort codes in decreasing order. Remove any duplicated codes and pad with empty character strings.
for (i in 1:nrow(dani)){
  sorted.codes <- sort(as.character(dani[i,6:29]), decreasing = TRUE)
  sorted.codes <- subset(sorted.codes, !duplicated(sorted.codes) | str_length(sorted.codes) == 0)
  
  len.diff <- 24 - length(sorted.codes)
  sorted.codes <- c(sorted.codes, rep("", len.diff))
  
  dani[i,6:29] <- sorted.codes
}


# Read in and clean my MI data --------------------------------------------

load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.v20180712.rdata")
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/cancer.codes.v20180227.1.rdata")
load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20180829.rdata')
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/Michigan/mi.raw.data.v20170818.rdata")

mi$studyid <- paste0('mi',mi$studyid)

cancer.codes$state <- substr(cancer.codes$studyid, 1, 2)
cancer.codes <- filter(cancer.codes, state == 'mi')

mi.comorbid <- filter(bd.codes.mi, studyid %in% cancer.codes$studyid)
mi.comorbid <- left_join(mi.comorbid, select(cancer.codes, -state), by = 'studyid')
mi.comorbid <- left_join(mi.comorbid, select(goback, studyid, sex, birth.yr, m.race, person.yrs))
#' In one instance, cancer_yr2 is less than cancer_yr1.
mi.comorbid <- left_join(mi.comorbid, select(mi, studyid, cancer_yr1, cancer_yr2, cancer_yr3, cancer_yr4), by = 'studyid')

mi.comorbid <- rename(mi.comorbid, bxyear = birth.yr, icd.o.iii.site = site_code1, icd.o.iii.cell.type = morph31, cell.behavior = behavior1)
mi.comorbid$icd.o.iii.site <- paste0('C',mi.comorbid$icd.o.iii.site)
mi.comorbid$sex <- as.character(mi.comorbid$sex)
columns <- c('sex','icd.o.iii.cell.type','icd.o.iii.site','cell.behavior')

#' Reformat character classes for important columns to match Dani's.
for (i in columns){
  mi.comorbid[,i] <- as.character(mi.comorbid[,i])
}

#' Reformat ICD9 codes to match Dani's.
for (i in 2:25){
  mi.comorbid[,i] <- as.character(mi.comorbid[,i])
  mi.comorbid[,i] <- str_remove(mi.comorbid[,i], '[[:punct:]]')
  mi.comorbid[,i] <- ifelse(is.na(mi.comorbid[,i]), '', mi.comorbid[,i])
}

#' Sort codes to match Dani's. Remove duplicate codes and pad with empty strings.
for (i in 1:nrow(mi.comorbid)){
  
  sorted.codes <- sort(as.character(mi.comorbid[i,2:25]), decreasing = TRUE)
  sorted.codes <- subset(sorted.codes, !duplicated(sorted.codes) | str_length(sorted.codes) == 0)
  
  len.diff <- 24 - length(sorted.codes)
  sorted.codes <- c(sorted.codes, rep("", len.diff))
  
  mi.comorbid[i,2:25] <- sorted.codes
  
}



# Attempt to join the two -------------------------------------------------

dani <- arrange(dani, desc(icd9cod1))
mi.comorbid <- arrange(mi.comorbid, desc(icd9cod1))

rm(goback, bd.codes.mi, cancer.codes, mi.codes, mi.demo, mi, columns, i, sorted.codes)

years <- sort(unique(dani$bxyear))

tmp <- filter(dani, bxyear == years[1])
tmp2 <- filter(mi.comorbid, bxyear == years[1])

tmp <- left_join(tmp, tmp2, by = c('icd9cod1','icd9cod2','icd9cod3','icd9cod4', 'icd9cod5', 
                                   'icd.o.iii.site','icd.o.iii.cell.type','cell.behavior'))
matched <- tmp[!is.na(tmp$studyid),]

for (i in years[2:22]){
  
  tmp <- dani[dani$bxyear == i, ]
  tmp2 <- mi.comorbid[mi.comorbid$bxyear == i, ]
  
  tmp <- left_join(tmp, tmp2, by = c('icd9cod1','icd9cod2','icd9cod3','icd9cod4', 'icd9cod5', 
                                     #'icd9cod6', 'icd9cod7', 'icd9cod8', 'icd9cod9', 'icd9cod10',
                                     'icd.o.iii.site','icd.o.iii.cell.type','cell.behavior'))
  
  matched <- rbind(matched, tmp[!is.na(tmp$studyid),])
  
}

map <- select(matched, recruiting, studyid)

map <- left_join(map, select(goback, studyid, m.edu2, person.yrs, m.race), by = 'studyid')
map <- rename(map, age.at.cancer.dx = person.yrs, maternal.education = m.edu2)

write.csv(map, file = 'W:/Old_genepi2/Jeremy/GOBACK/Manuscripts/Recruitment Paper/mi.covariates.for.dani.csv', row.names = FALSE)

rm(list = ls()); gc()



# Read in TX data and add maternal education ------------------------------

load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20180829.rdata")
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/linked.registry.ids.v20180908.rdata")

goback <- filter(goback, state == 'TX')
goback.ids <- filter(goback.ids, state == 'TX')

#' Note: I manually moved patientid up to the second column in this Excel file.
tx <- read.xlsx(file = 'W:/Old_genepi2/Jeremy/GOBACK/Manuscripts/Recruitment Paper/TX contact and dx info MA.xlsx',
                sheetIndex = 1, rowIndex = 1:360, header = TRUE, stringsAsFactors = FALSE)
tx <- select(tx, bc_link, patientid, histbehav_desc1)
tx <- rename(tx, studyid = bc_link, cancer.registry.id = patientid)
tx$studyid <- paste0('tx', tx$studyid)

tmp2 <- filter(goback.ids, state == 'TX')
tmp2 <- tmp2[!is.na(tmp2$cancer.registry.id), ]

tmp <- left_join(tx, tmp2, by = 'cancer.registry.id')
tmp <- left_join(tmp, select(goback, studyid, m.edu2), by = c('studyid.y' = 'studyid'))
tmp <- select(rename(tmp, bc_link = studyid.x, patientid = cancer.registry.id, jeremy.studyid = studyid.y),
              bc_link, patientid, jeremy.studyid, bd.registry.id, m.edu2)

write.csv(tmp, file = 'W:/Old_genepi2/Jeremy/GOBACK/Manuscripts/Recruitment Paper/tx.maternal.education.for.dani.csv',
          row.names = FALSE)







tmp <- tmp[!is.na(tmp$state), ]
tmp <- left_join(tx, tmp2, by = 'studyid')
tmp <- left_join(tx, goback, by = 'studyid')




# Scratch paper -----------------------------------------------------------

tmp <- tmp[!is.na(tmp$studyid),]



#' For now, remove less useful columns. Will make comparison easier.
dani <- dani[,c(1,6:29,36,37,40,52:57)]
mi.comorbid <- mi.comorbid[,c(1:26,31,36,46:50)]

not.matched <- dani[!(dani$recruiting %in% matched$recruiting), ]

candidates <- filter(mi.comorbid, !(studyid %in% matched$studyid))

not.matched[35:85,c(1:9,27,31:33)]
candidates[1:50,c(1:9,29,26:28)]
candidates[51:100,c(1:9,29,26:28)]

tmp <- left_join(not.matched, candidates, by = c('icd9cod1','icd9cod2','icd9cod3','icd9cod4', 'icd9cod5', 
                                                'icd9cod6', 'icd9cod7', 'icd9cod8', 
                                                'icd.o.iii.site','icd.o.iii.cell.type','cell.behavior'))



#' Generates 4 new rows. 2 matches per row in a few cases?
dani <- left_join(dani, mi.comorbid, by = c('icd9cod1','icd9cod2','icd9cod3','icd9cod4','icd.o.iii.site','icd.o.iii.cell.type','cell.behavior'))


codes <- unique(dani$icd9cod1)
print(sort(codes))


mi.comorbid[1:25,c(1:10,26:28)]
 
code.sort <- mi.comorbid

for (i in 1:nrow(code.sort)){
  sorted.codes <- sort(as.character(code.sort[i,2:25]))
  code.sort[i,2:25] <- sorted.codes
}



tmp <- filter(dani, cancer.diagnosis.year == 2002)
tmp <- arrange(tmp, icd9cod1)

tmp2 <- filter(mi.comorbid, cancer_yr1 == 2002)
tmp2 <- arrange(tmp2, icd9cod1)

tmp <- left_join(tmp, tmp2, by = c('icd9cod1','icd9cod2','icd9cod3','icd9cod4','icd.o.iii.site','icd.o.iii.cell.type','cell.behavior'))

tmp[,1:12]
tmp2[,1:10]

search <- subset(goback.ids, !is.na(goback.ids$cancer.registry.id))
search <- subset(search, grepl('tx2005', search$studyid))
search <- arrange(search, studyid)

tx.2005 <- subset(tx, grepl('tx2005', tx$studyid))
tx.2005 <- arrange(tx.2005, studyid)

head(tx.2005)

pattern <- '32'
find <- subset(search, grepl(pattern, search$studyid))

load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/Texas/tx.raw.data.date.cancer.dx.v20171018.rdata")

tmp <- left_join(tx, tx.datedx, by = c('cancer.registry.id' = 'patientid'))
tmp$studyid <- paste0('tx', tmp$birthid)
tmp <- select(tmp, patientid, studyid)
tmp <- left_join(tmp, select(goback, studyid, m.edu2), by = 'studyid')
