#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Generates a data frame of all the BPA codes for kids in TX and NC, and 
#' another with the ICD codes for kids in MI.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(tidyverse); require(magrittr); require(haven)

#' TODO: Append MI children born after 9/2015 using ICD10 codes.
#' TODO: Append corrected three-digit ICD9 codes for MI kids.

# Load TX data, create BD codes data frame --------------------------------

load('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/Texas/tx.raw.data.rdata')
load("Z:/Jeremy/GOBACK/Datasets/Texas/tx.raw.data.bpa.codes.v20171106.rdata")
load('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/Texas/tx.raw.data.bpa.codes.v20171106.rdata')

tx.raw <- select(tx.raw, birthID, CASE_ID)
tx.sel <- left_join(tx.bd, tx.raw, by = 'CASE_ID')
tx.sel$studyid <- paste0('tx',tx.sel$birthID)
tx.sel <- tx.sel[,c(70, 3:68)]

for (i in 2:67){
  tx.sel[, i] <- ifelse(str_length(tx.sel[, i]) == 2, paste0('0',tx.sel[, i],'.000'),
                        ifelse(str_length(tx.sel[, i]) == 3,paste0(tx.sel[, i],'.000'),
                               ifelse(str_length(tx.sel[, i]) == 5, paste0(tx.sel[, i], '00'), 
                                      ifelse(str_length(tx.sel[, i]) == 6, paste0(tx.sel[, i], '0'), 
                                             ifelse(str_length(tx.sel[, i]) == 7, tx.sel[, i], tx.sel[, i])))))
}

rm(tx.raw, tx.bd)

load("Z:/Jeremy/GOBACK/Datasets/Old Datasets/North Carolina/nc.birth.defect.data.rdata")

#' Select only populated defects columns.
nc.sel <- select(nc.bd, NCID, DX1:DX47)
nc.sel$NCID <- paste0('nc',nc.sel$NCID)

for (i in 2:48){
  nc.sel[,i] <- ifelse(nc.sel[,i] == "", NA, nc.sel[,i]) 
}

for (i in 49:67){
  nc.sel[,i] <- as.character(NA)
}

for (i in 2:67){
  nc.sel[ ,i] <- ifelse(str_length(nc.sel[,i]) == 6, paste0(substr(nc.sel[,i], 1, 3), '.', substr(nc.sel[,i], 4, 6)), nc.sel[,i])
}

names(nc.sel) <- names(tx.sel)

nc.sel <- nc.sel[!is.na(nc.sel$bpa1), ]

bd.codes.txnc <- rbind(tx.sel, nc.sel)

rm(nc.bd, tx.sel, nc.sel)

load("Z:/Jeremy/GOBACK/Datasets/Michigan/mi.birthdefects.codes.rdata")

bd.codes.mi <- mi.bd[, c(3, 110:133)]
bd.codes.mi <- bd.codes.mi[!is.na(bd.codes.mi$ICD9COD1), ]
names(bd.codes.mi) <- tolower(names(bd.codes.mi))
bd.codes.mi$studyid <- paste0('mi',bd.codes.mi$studyid)

rm(mi.bd)

save(bd.codes.txnc, file = 'Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.v20180606.rdata')
save(bd.codes.mi, file = 'Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.v20180606.rdata')

rm(list = ls());gc()

#' Addendum: this is code taken from the data cleaning script in the DS-ALL BD project that adds back in
#' MI kids with cancer, who were not previously included in this data frame.
require(readstata13); require(dplyr)

load("Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.transpose.v20180614.rdata")
load("Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.v20180606.rdata")

#' A file from Tiffany that has the BD codes for kids with cancer.
mi.bd.cancercases <- read.dta13('Z:/Jeremy/DS-ALL BD project/Datasets/Raw datasets/mi_cancer_long_jeremy_071118.dta', convert.underscore = TRUE)
mi.bd.cancercases <- subset(mi.bd.cancercases, !duplicated(mi.bd.cancercases$studyid))
mi.bd.cancercases$studyid <- paste0('mi',mi.bd.cancercases$studyid)
mi.bd.cancercases <- mi.bd.cancercases[ , c(2,109:132)]
names(mi.bd.cancercases) <- tolower(names(mi.bd.cancercases))

bd.codes.mi <- rbind(bd.codes.mi, mi.bd.cancercases)
bd.codes.mi <- subset(bd.codes.mi, !duplicated(bd.codes.mi$studyid))

#' Remove new rows for cancer cases with no birth defects.
bd.codes.mi$rowsum <- rowSums(bd.codes.mi[2:25], na.rm = TRUE)
bd.codes.mi <- subset(bd.codes.mi, bd.codes.mi$rowsum > 0)
bd.codes.mi <- bd.codes.mi[, 1:25]

save(bd.codes.mi, file = 'Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.v20180712.rdata')

# Scratch paper -----------------------------------------------------------


#' Actually a terrible name for this function: it computes whether there is ANY defect.
compute.num.defect <- function(data.frame){ifelse(rowSums(data.frame[,2:67], na.rm = TRUE) >= 1, 1, 0)}

#' Sorts the vectors of unique codes in ascending order.
sort.codes <- function(dataframe){
  dataframe <- data.frame(codes = dataframe)
  dataframe <- arrange(dataframe, codes)
  dataframe <- c(dataframe$codes)
}

# Merge state-level data --------------------------------------------------

load("Z:/Jeremy/GOBACK/Datasets/Texas/tx.raw.data.bpa.codes.v20171106.rdata")
load("Z:/Jeremy/GOBACK/Datasets/Michigan/mi.birthdefects.codes.rdata")
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/North Carolina/nc.birth.defect.data.rdata")
load("Z:/Jeremy/GOBACK/Datasets/Texas/tx.raw.data.rdata")

mi.filt <- select(mi.bd, studyid, ICD9COD1:ICD9COD24)
nc.filt <- select(nc.bd, NCID, DX1:DX47)

#' Pad the length of short NC BD codes.
nc.filt %<>% 
  mutate(across(starts_with('DX'), ~ ifelse(nchar(.x) == 4, paste0(.x, '00'),
                                     ifelse(str_length(.x) == 5, paste0(.x, '0'), .x))))

tx.ids <- select(tx.raw, birthID, CASE_ID)
tx.filt <- left_join(tx.bd, tx.ids, by = 'CASE_ID')
tx.filt <- tx.filt[,c(69, 3:68)]

rm(tx.bd, tx.ids, tx.raw, mi.bd, bc.bd)

for (i in 2:48){
  nc.filt[,i] <- ifelse(nc.filt[,i] == "", NA, nc.filt[,i]) 
  nc.filt[,i] <- as.numeric(nc.filt[,i])
}

for (i in 49:67){
  nc.filt[,i] <- as.numeric(NA)
}

for (i in 26:67){
  mi.filt[,i] <- as.numeric(NA)
}

newcols <- c('studyid',rep(paste0('defectcode',1:66)))

colnames(nc.filt) <- newcols
colnames(mi.filt) <- newcols
colnames(tx.filt) <- newcols

tx.filt$studyid <- paste0('tx',tx.filt$studyid)
nc.filt$studyid <- paste0('nc',nc.filt$studyid)
mi.filt$studyid <- paste0('mi',mi.filt$studyid)

mi.filt$num.defect <- compute.num.defect(mi.filt)
nc.filt %<>% 
  mutate(num.defect = if_all(starts_with('defectcode'), ~ is.na(.x)))
tx.filt$num.defect <- compute.num.defect(tx.filt)

defect.codes <- rbind(tx.filt, mi.filt, nc.filt)
defect.codes <- filter(defect.codes, num.defect == 1)
defect.codes <- defect.codes[!duplicated(defect.codes$studyid), ]

rm(mi.filt, nc.filt, tx.filt, i, newcols)

save(defect.codes, file = 'Z:/Jeremy/GOBACK/Datasets/birth.defects.codes.rdata')

# Inspect MI defects: Do all MI kids have a defect? -----------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2017.11.08.
#' 
#' Should be able to use this data frame to investigate whether kids who
#' are in that discordant group for the new and old any.defect variables
#' actually do or do not have a defect.
#' 
#' Michigan has ~10% of their dataset flagged as having a defect, which 
#' seems tremendously unlikely.  Look into the status of these codes in 
#' the MI data.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

load("Z:/Jeremy/GOBACK/Datasets/goback.v20171108.1.rdata")

goback$computed.any.defect <- ifelse(rowSums(goback[,22:112], na.rm = TRUE) >= 1, 1, 0)  

tmp <- filter(goback, computed.any.defect == 0 & any.birthdefect == 1)
tmp <- select(tmp, studyid)
tmp <- left_join(tmp, defect.codes, by = 'studyid')
tmp$computed.any.defect <- compute.num.defect(tmp)
tmp <- left_join(tmp, select(goback, studyid, state), by = 'studyid')

mi.ids.old <- filter(tmp, state == 'MI')
mi.ids.old <- c(mi.ids.old$studyid)

tmp <- filter(tmp, computed.any.defect == 1)

mi.ids.new <- filter(tmp, state == 'MI')
mi.ids.new <- c(mi.ids.new$studyid)

no.defects <- setdiff(mi.ids.old, mi.ids.new)

#' No record of these IDs in the defect codes data frame.
#' But they ARE valid studyids, as evidenced by the 2nd statement.
tmp <- defect.codes[defect.codes$studyid %in% no.defects, ]
tmp <- goback[goback$studyid %in% no.defects, ]

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Looks like 69,265 of the 76,236 kids tagged with the old any defect 
#' variable but not the new one DO have a birth defect.  
#' 
#' Most of the rest are from AR, and we already addressed this issue. 
#' Philip is confident these kids have SOME defect, and we are going to 
#' leave them classified this way.
#' 
#' However, there are 1,018 kids in MI where the situation is the same: 
#' they're tagged with the old any defect variable, but don't have any 
#' actual defect codes. Will point this out to Philip.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

mi <- filter(goback, state == 'MI')
mi <- filter(mi, any.birthdefect == 1)
mi <- left_join(select(mi, studyid), defect.codes, by = 'studyid')
mi$num.defect <- compute.num.defect(mi)
tmp <- c(mi$studyid)
mi <- filter(mi, num.defect == 1)
tmp2 <- c(mi$studyid)

#' A vector of MI kids marked as having a birth defect in GOBACK but for which no evidence of a defect exists.
no.defect <- setdiff(tmp, tmp2)
save(no.defect, file = 'Z:/Jeremy/GOBACK/Datasets/Michigan/potential.bd.misclassifications.rdata')

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' This is still insufficient to explain why 10% of MI has at least one 
#' birth defect code.  Two things to consider:
#' 
#' - Is MI classifying certain things as birth defects that TX and NC are 
#' not?  Could evaluate this by comparing the differences of all the unique
#' codes in MI as compared to the other two states.
#' 
#' - Why are the prevalences of the 'other.major' defects categories so 
#' high in MI?  Is there misclassification here?
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

# Inspect MI defects: Is MI including different codes? --------------------

load("Z:/Jeremy/GOBACK/Datasets/birth.defects.codes.rdata")
load("Z:/Jeremy/GOBACK/Datasets/goback.v20171108.1.rdata")

#' Generate a list of all the unique codes from MI.
pattern <- '^mi'

tmp <- defect.codes[grep(pattern, defect.codes$studyid), ]

mi.codes <- unique(tmp$defectcode1)

for (i in 3:67){
  new.codes <- unique(tmp[,i])
  diff <- setdiff(new.codes, mi.codes)
  codes <- c(mi.codes, diff)
}

#' Generate a list of all the unique codes from TX.
pattern <- '^tx'
tmp <- defect.codes[grep(pattern, defect.codes$studyid), ]

tx.codes <- unique(tmp$defectcode1)

for (i in 3:67){
  new.codes <- unique(tmp[,i])
  diff <- setdiff(new.codes, tx.codes)
  tx.codes <- c(tx.codes, diff)
}

tx.codes <- as.character(tx.codes)

#' Generate a list of all the unique codes from NC.
pattern <- '^nc'
tmp <- defect.codes[grep(pattern, defect.codes$studyid), ]

nc.codes <- unique(tmp$defectcode1)

for (i in 3:67){
  new.codes <- unique(tmp[,i])
  diff <- setdiff(new.codes, nc.codes)
  nc.codes <- c(nc.codes, diff)
}

rm(new.codes, diff, tmp, pattern, i)

tx.codes <- sort.codes(tx.codes)
mi.codes <- sort.codes(mi.codes)
nc.codes <- sort.codes(nc.codes)

mi.v.nc <- setdiff(mi.codes, nc.codes)
mi.v.tx <- setdiff(mi.codes, tx.codes)

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' There is definitely variation in how many codes are recorded in each
#' state.  TX has 1,005 codes, NC has 2,241 and MI has 1,974.
#' 
#' MI is also definitely including codes that are outside the ranges NC and
#' TX are.  Likely this explains the increased prevalence of defects in MI?
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

#' TX appears to have every code we are interested in and ONLY codes we are interested in.
#' What would be the frequency of birth defects in MI if we picked out only the kids whose codes were in TX?
pattern <- '^mi'

#' Cut out NA.
codes <- tx.codes[-1005]

tmp <- defect.codes[grep(pattern, defect.codes$studyid), ]
cases <- tmp[tmp[,2] %in% codes, ]

for (i in 3:67){
  new.cases <- tmp[tmp[,i] %in% codes, ]
  cases <- rbind(cases, new.cases)
  cases <- cases[!duplicated(cases$studyid), ]
}

#' If we do this, we cut the number of kids with birth defects down from 234k to 177k.
#' That's a more reasonable number, but still high.



# Inspect MI defects: other.major variables -------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' The other mystery here is why so many kids in MI fall into these
#' 'other.major' categories.
#' 
#' As an example, there are 2x as many kids in MI as compared to TX with 
#' 'other major congenital anomalies of the eye' despite the fact that the
#' TX dataset is more than twice as large.
#' 
#' What sorts of codes are these MI children flagged with?
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
load("Z:/Jeremy/GOBACK/Datasets/birth.defects.codes.rdata")
load("Z:/Jeremy/GOBACK/Datasets/goback.v20171108.1.rdata")

tmp <- filter(goback, state == 'MI' & conganomalies.eye.other.major == 1)
eye.other.ids <- c(tmp$studyid)
tmp <- defect.codes[defect.codes$studyid %in% eye.other.ids, ]

codes <- unique(tmp$defectcode1)

for (i in 3:67){
  new.codes <- unique(tmp[,i])
  diff <- setdiff(new.codes, codes)
  codes <- c(codes, diff)
}

codes <- sort.codes(codes)

print(tmp[1:100,2:26])

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' At a glance, all these kids do appear to have codes for the eye 
#' (range 743).  Let's confirm this.
#' 
#' The code below confirms that all kids who are flagged as having an 
#' other.major eye defect do in fact have it.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

eye.other.codes <- c(743.20,743.21,743.22,
                     743.35,743.36,743.37,743.39, #' Certain subtypes of congenital cataract should be included with other major.
                     743.40,743.41,743.42,743.43,743.44,743.45,743.46,743.47,743.48,743.49,
                     743.50,743.51,743.52,743.53,743.54,743.55,743.56,743.57,743.58,743.59,
                     743.60,743.61,743.62,743.63,743.64,743.65,743.66,743.69,
                     743.80,743.90)

tmp <- filter(goback, state == 'MI' & conganomalies.eye.other.major == 1)
tmp2 <- c(tmp$studyid)
tmp <- defect.codes[defect.codes$studyid %in% tmp2, ]

confirmed <- tmp[tmp[,2] %in% eye.other.codes, ]

for (i in 3:25){
  tmp2 <- tmp[tmp[,i] %in% eye.other.codes, ]
  confirmed <- rbind(confirmed, tmp2)
  confirmed <- confirmed[!duplicated(confirmed$studyid), ]
}

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' If they also have other eye defects, are they flagged for those?
#' 
#' Using anopthalmos and congenital cataracts as examples, yes.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
anopthalmos.codes <- c(743.00, 743.03, 743.06)

tmp <- filter(goback, state == 'MI' & conganomalies.eye.other.major == 1)
tmp2 <- c(tmp$studyid)
tmp <- defect.codes[defect.codes$studyid %in% tmp2, ]

anop <- tmp[tmp[,2] %in% anopthalmos.codes, ]
for (i in 3:25){
  tmp2 <- tmp[tmp[,i] %in% anopthalmos.codes, ]
  anop <- rbind(anop, tmp2)
  anop <- anop[!duplicated(anop$studyid), ]
}

tmp <- c(anop$studyid)
tmp <- goback[goback$studyid %in% tmp, ]
table(tmp$anopthalmos.micropthalmos)

cat.codes <- c(743.30,743.31,743.32,743.33,743.34)

tmp <- filter(goback, state == 'MI' & conganomalies.eye.other.major == 1)
tmp2 <- c(tmp$studyid)
tmp <- defect.codes[defect.codes$studyid %in% tmp2, ]

cat <- tmp[tmp[,2] %in% cat.codes, ]
for (i in 3:25){
  tmp2 <- tmp[tmp[,i] %in% cat.codes, ]
  cat <- rbind(cat, tmp2)
  cat <- cat[!duplicated(cat$studyid), ]
}

tmp <- c(cat$studyid)
tmp <- goback[goback$studyid %in% tmp, ]
table(tmp$congenitalcataract)

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Let's try another other.major variable just to cut down on the chances
#' this one is an aberration.
#' 
#' Heart and circulatory system other major is also very enriched in the
#' MI data.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
tmp <- filter(goback, state == 'MI' & heart.circsys.other.major == 1)
tmp2 <- c(tmp$studyid)
tmp <- defect.codes[defect.codes$studyid %in% tmp2, ]

codes <- unique(tmp$defectcode1)

for (i in 3:25){
  new.codes <- unique(tmp[,i])
  diff <- setdiff(new.codes, codes)
  codes <- c(codes, diff)
}

codes <- sort.codes(codes)

print(tmp[1:100,2:26])
print(tmp[101:200,1:26])

heart.other.codes <- c(745.70, 745.80, 745.90, 746.00, 746.09, 746.40, 746.50, 746.60,
                       746.80, 746.81, 746.82, 746.83, 746.84, 746.85, 746.86, 746.87, 746.89, 746.90,
                       747.20, 747.21, 747.22, 747.29,
                       747.40, 747.41, 747.42, 747.49, 
                       747.50,
                       747.60, 747.69, 747.80, 747.83, 747.89, 747.90)

tmp <- filter(goback, state == 'MI' & heart.circsys.other.major == 1)
heart.other.ids <- c(tmp$studyid)
tmp <- defect.codes[defect.codes$studyid %in% heart.other.ids, ]

confirmed <- tmp[tmp[,2] %in% heart.other.codes, ]
for (i in 3:25){
  tmp2 <- tmp[tmp[,i] %in% heart.other.codes, ]
  confirmed <- rbind(confirmed, tmp2)
  confirmed <- confirmed[!duplicated(confirmed$studyid), ]
}

confirmed.ids <- c(confirmed$studyid)
unconfirmed.ids <- setdiff(heart.other.ids, confirmed.ids)

tmp <- defect.codes[defect.codes$studyid %in% unconfirmed.ids, ]

print(tmp[1:100,1:25])
print(tmp[101:200,1:25])

print(goback[goback$studyid == 'mi2001048286', ])
print(goback[goback$studyid == 'mi2001282210', ])
print(goback[goback$studyid == 'mi2000795719', ])
print(goback[goback$studyid == 'mi2002417880', ])
print(goback[goback$studyid == 'mi2002412748', ])

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' In the case of other major CHDs, I can't square up the current 
#' frequencies of that variable with the ICD codes.
#' 
#' There are 13,075 MI kids in GOBACK with chd.other.major == 1. 
#' 
#' 13,018 of them appear in the MI birth defects file.
#' 
#' As far as I can tell, 12,458 of them DO have other major defects.
#' 
#' However, there are also clearly errors: for example, 747.11 is IAA, but
#' it isn't being used to assign kids to that category and instead is 
#' showing up as other major.  Only a fraction of the IAA kids are flagged
#' as having it.
#' 
#' VSD and endocardial cushion defect are both correct, though.  I think 
#' IAA is problematic because it has that 747.11 cateogry that didn't get
#' accounted for in the code.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
iaa.codes <- c(747.10, 747.11)

tmp <- filter(goback, state == 'MI' & heart.circsys.other.major == 1)
tmp2 <- c(tmp$studyid)
tmp <- defect.codes[defect.codes$studyid %in% tmp2, ]

iaa <- tmp[tmp[,2] %in% iaa.codes, ]
for (i in 3:25){
  tmp2 <- tmp[tmp[,i] %in% iaa.codes, ]
  iaa <- rbind(iaa, tmp2)
  iaa <- iaa[!duplicated(iaa$studyid), ]
}

tmp <- c(iaa$studyid)
tmp <- goback[goback$studyid %in% tmp, ]
table(tmp$interrupted.aortic.arch.type.a.or.c, useNA = 'always')

#' Try another.
cushion.codes <- c(745.60,745.61,745.63,745.68,745.69)

tmp <- filter(goback, state == 'MI' & heart.circsys.other.major == 1)
tmp2 <- c(tmp$studyid)
tmp <- defect.codes[defect.codes$studyid %in% tmp2, ]

cushion <- tmp[tmp[,2] %in% cushion.codes, ]
for (i in 3:25){
  tmp2 <- tmp[tmp[,i] %in% cushion.codes, ]
  cushion <- rbind(cushion, tmp2)
  cushion <- cushion[!duplicated(cushion$studyid), ]
}

tmp <- c(cushion$studyid)
tmp <- goback[goback$studyid %in% tmp, ]
table(tmp$endocardialcushiondefect, useNA = 'always')

vsd.codes <- c(745.40)

tmp <- filter(goback, state == 'MI' & heart.circsys.other.major == 1)
tmp2 <- c(tmp$studyid)
tmp <- defect.codes[defect.codes$studyid %in% tmp2, ]

vsd <- tmp[tmp[,2] %in% vsd.codes, ]
for (i in 3:25){
  tmp2 <- tmp[tmp[,i] %in% vsd.codes, ]
  vsd <- rbind(vsd, tmp2)
  vsd <- vsd[!duplicated(vsd$studyid), ]
}

tmp <- c(vsd$studyid)
tmp <- goback[goback$studyid %in% tmp, ]
table(tmp$ventricularseptaldefect, useNA = 'always')

rm(vsd,tmp, tmp2, vsd.codes, i)

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' So what have we learned.
#' 
#' There are 1,018 kids from Michigan marked as having a birth defect
#' who have no diagnostic codes.
#' 
#' Many MI codes fall outside the ranges in the 'six digit codes for 
#' reportable birth defects' document.  If we restrict to the codes found
#' in the TX dataset, the number of MI kids with defects drops from 
#' 234,000 to 177,000.
#' 
#' The other major categories are enriched in MI data, and other minor 
#' categories are generally depleted.  Some of the kids currently in these
#' categories may not belong there, but this is insufficient to explain 
#' the preponderance of birth defects in the MI dataset.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------




# Could other major defects be UNDERdiagnosed in other states? ------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' If it's the case that all or most MI kids flagged as other major do 
#' belong in that category, but the numbers are higher than in other 
#' states, could these defects not be captured correctly in other states?
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
load("Z:/Jeremy/GOBACK/Datasets/birth.defects.codes.rdata")
load("Z:/Jeremy/GOBACK/Datasets/goback.v20171108.1.rdata")

eye.other.codes <- c(743.20,743.21,743.22,
                     743.35,743.36,743.37,743.39, #' Certain subtypes of congenital cataract should be included with other major.
                     743.40,743.41,743.42,743.43,743.44,743.45,743.46,743.47,743.48,743.49,
                     743.50,743.51,743.52,743.53,743.54,743.55,743.56,743.57,743.58,743.59,
                     743.60,743.61,743.62,743.63,743.64,743.65,743.66,743.69,
                     743.80,743.90)

tmp <- filter(goback, state == 'TX')

table(tmp$conganomalies.eye.other.major, useNA = 'always')

tmp2 <- c(tmp$studyid)
tmp <- defect.codes[defect.codes$studyid %in% tmp2, ]

confirmed <- tmp[tmp[,2] %in% eye.other.codes, ]

for (i in 3:67){
  tmp2 <- tmp[tmp[,i] %in% eye.other.codes, ]
  confirmed <- rbind(confirmed, tmp2)
  confirmed <- confirmed[!duplicated(confirmed$studyid), ]
}

print(confirmed[1:100,2:18])
print(confirmed[101:200,2:18])
print(confirmed[15908:16108,2:18])

tmp <- c(confirmed$studyid)
tmp <- goback[goback$studyid %in% tmp, ]
table(tmp$conganomalies.eye.other.major, useNA = 'always')

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Well well...going by the codes, it appears there would be 16,000 kids in 
#' TX diagnosed with these anomalies.  
#' 
#' This is more in line with what I would expect (slightly more than 2x 
#' the number in MI).  Let's see if this is also the situation for 
#' CHD other major.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
heart.other.codes <- c(745.70, 745.80, 745.90, 746.00, 746.09, 746.40, 746.50, 746.60,
                       746.80, 746.81, 746.82, 746.83, 746.84, 746.85, 746.86, 746.87, 746.89, 746.90,
                       747.20, 747.21, 747.22, 747.29,
                       747.40, 747.41, 747.42, 747.49, 
                       747.50,
                       747.60, 747.69, 747.80, 747.83, 747.89, 747.90)

tmp <- filter(goback, state == 'TX')
table(tmp$heart.circsys.other.major, useNA = 'always')

tmp2 <- c(tmp$studyid)
tmp <- defect.codes[defect.codes$studyid %in% tmp2, ]

confirmed <- tmp[tmp[,2] %in% heart.other.codes, ]

for (i in 3:67){
  tmp2 <- tmp[tmp[,i] %in% heart.other.codes, ]
  confirmed <- rbind(confirmed, tmp2)
  confirmed <- confirmed[!duplicated(confirmed$studyid), ]
}

print(confirmed[1:100,2:18])
print(confirmed[101:200,2:18])
print(confirmed[28660:28760,2:18])

tmp <- c(confirmed$studyid)
tmp <- goback[goback$studyid %in% tmp, ]
table(tmp$conganomalies.eye.other.major, useNA = 'always')

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' For CHD other major, I come away with a list of ~29,000 TX kids who
#' should be flagged 1, compared to the 7,253 who are flagged that way now.
#' 
#' So contrary to my initial hypothesis, these are not OVERrepresented in 
#' Michigan, they're UNDERrepresented everywhere else.  Mind = blown.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------







# Append Oklahoma data ----------------------------------------------------

require(dplyr)

load('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.v20180606.rdata')
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/Oklahoma/oklahoma.raw.data.rdata")

ok$studyid <- paste0('ok',ok$randID)

ok <- select(filter(ok, bd == 1), studyid, CODE30:CODE37)

#' Standardize the names and values in each BD column.
for (i in 2:ncol(ok)){
  
  ok[,i] <- ifelse(ok[,i] == '.', NA, ok[,i])
  
  ok[,i] <- ifelse(nchar(ok[,i]) == 6, paste0(substr(ok[,i],1,3),'.',substr(ok[,i],4,6)),
            ifelse(nchar(ok[,i]) == 5, paste0(substr(ok[,i],1,3),'.',substr(ok[,i],4,5),'0'), NA))

}

#' For each row, sort the codes in ascending order.
for (i in 1:nrow(ok)){
  
 sorted.codes <- c(ok[i,'studyid'],
                   sort(paste(ok[i,2:ncol(ok)])))
 
 ok[i,] <- sorted.codes
 
}

#' Append empty columns to match the existing BD data frame and name them.
ok[,(ncol(ok)+1):ncol(bd.codes.txnc)] <- as.character(NA)

names(ok)[2:ncol(ok)] <- paste0('bpa',1:(ncol(ok)-1))

bd.codes.ok.nc.tx <- bind_rows(bd.codes.txnc, ok)

save(bd.codes.ok.nc.tx, 
     file = '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txncok.v20191205.rdata')


# Append MA data ----------------------------------------------------------

require(tidyverse)

load('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txncok.v20191205.rdata')

ma.bd <- read_dta('All BD.dta') %>% 
  mutate(studyid = paste0('ma', KIDUID)) %>% 
  as.data.frame()

#' Select the study ID column and the columns containing BPA codes.
columns <- names(ma.bd) %>% subset(grepl('DxCode',.)) %>% c('studyid',.)

ma.bd  %<>% select(all_of(columns))

#' For each row, sort the codes in ascending order.
#' Pad codes with trailing zeroes based on string length.
for (i in 1:nrow(ma.bd)){
  
  sorted.codes <- unlist(c(ma.bd[i,'studyid'],
                         sort(paste(ma.bd[i,2:ncol(ma.bd)]))))
  
  sorted.codes <- ifelse(nchar(sorted.codes) == 3, paste0(sorted.codes,'.000'), 
                  ifelse(nchar(sorted.codes) == 5, paste0(sorted.codes, '00'),
                  ifelse(nchar(sorted.codes) == 6, paste0(sorted.codes, '0'), sorted.codes)))
  
  ma.bd[i,columns] <- sorted.codes
  
}

#' Append empty columns to match the existing BD data frame and name them.
ma.bd[,(ncol(ma.bd)+1):ncol(bd.codes.txnc)] <- as.character(NA)

names(ma.bd)[2:ncol(ma.bd)] <- paste0('bpa',1:(ncol(ma.bd)-1))

bd.codes.merge <- bind_rows(bd.codes.ok.nc.tx, ma.bd)

saveRDS(bd.codes.merge, 
        '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.merge.v20210806.rds')

# Append updated MI data --------------------------------------------------

#' New Michigan data.
mi <- readRDS('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/Michigan/mi.raw.data.v20220511.rds')

#' ICD9 codes were used to record birth defect diagnoses until 10/2015.
#' For now, I will append older data with ICD9 codes.
mi %<>% 
  select(STUDYID, DOB, BDefect_ICD9COD1:BDefect_ICD9COD30) %>% 
  filter(DOB < '2015-10-01') %>% 
  mutate(STUDYID = paste0('mi', STUDYID))

names(mi) %<>% tolower() %<>% str_replace_all('_', '.')

#' Kids w/o ICD9 codes are represented by either "" or ".". Set all to "".
mi %<>% mutate(across(bdefect.icd9cod1:bdefect.icd9cod30, ~ str_replace(.x, '\\Q.\\E', '')))

#' One or two digit codes would not be valid. Set to "".
mi %<>% mutate(across(bdefect.icd9cod1:bdefect.icd9cod30, ~ ifelse(str_length(.x) < 3, '', .x)))

#' Drop kids without any birth defects (i.e., kids for whom every icd9 code == "").
mi %<>% 
  mutate(flag = as.logical(if_all(bdefect.icd9cod1:bdefect.icd9cod30, ~ifelse(.x == '', 1, 0)))) %>% 
  filter(flag == F) %>% 
  select(-flag)

#' Change " to NA to match existing data.
mi %<>% mutate(across(bdefect.icd9cod1:bdefect.icd9cod30, ~ ifelse(.x == "", NA, .x)))

#' Standardize the format of BD codes: all five digit character vectors with a period after the third.
mi %<>% mutate(across(bdefect.icd9cod1:bdefect.icd9cod30, ~ ifelse(nchar(.x) == 3, paste0(.x, '.00'),
                                                            ifelse(nchar(.x) == 4, paste0(substr(.x,1,3), '.', substr(.x,4,4)),
                                                            ifelse(nchar(.x) == 5, paste0(substr(.x,1,3), '.', substr(.x,4,5)),
                                                            .x)))))

#' Drop kids without at least one defect in the ranges that contribute to our computed BD variables.
mi %<>% 
  mutate(flag = ifelse(if_any(bdefect.icd9cod1:bdefect.icd9cod30, 
                                  ~str_detect(.x, str_c('^216', '^237.7', '^7[45]', sep = '|'))), 1, 0)) %>% 
  filter(flag == 1)

#' Sort BD columns such that each child's codes run from 217 through 759, with NAs at the end.
#' This is slow, since I'm updating every cell of a large matrix individually.

#' Names for the thirty BD columns.
my.cols <- paste0('bdefect.icd9cod',1:30)

for (i in 1:nrow(mi)){
  
  # A sorted vector of the ith child's BD codes, with NAs at the end.
  index.codes <- mi %>% select(bdefect.icd9cod1:bdefect.icd9cod30) %>% slice(i) %>% c() %>% unlist() %>% sort(na.last = T)
  
  #' Replace BD variables for the ith row with the sorted vector, so that they will be sorted.
  for (j in 1:length(my.cols)){
    
    mi[i, my.cols[j]] <- index.codes[j]
    
  }

}

mi %<>% rename_with(~paste0('icd9code',1:30), .cols = starts_with('bdefect')) %>% 
  mutate(state = 'MI',
         coding.system = 'ICD-9') %>% 
  select(studyid, state, coding.system, dob, icd9code1:icd9code30)

saveRDS(mi, 
        '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/birth_defects_codes_icd9_states.v20220609.rds')

rm(list = ls()); gc()

# Append updated TX data --------------------------------------------------

merged.codes <- readRDS('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.merge.v20210806.rds')

merged.codes %<>% 
  filter(!str_detect(studyid, '^tx')) %>% 
  select(-c(bpa64:bpa66))

#' A flag indicating whether a child has at least one phenotype we're counting towards our list of birth defects.
#' Remove kids without any birth defects that count towards our computed variables.
merged.codes %<>%  
  mutate(flag = ifelse(if_any(bpa1:bpa63, ~str_detect(.x, str_c('^216', '^228.0', '^237.7', '^7[45]', sep = '|'))), 1, 0)) %>% 
  filter(!is.na(flag)) %>% 
  select(-flag)

#' Read birth defects file. Looks like possible/probable cases have already been excluded. There are no duplicate IDs.
bds <- read_dta('//smb-main.ad.bcm.edu/genepi3/TiffanyChambers/GOBACK Data/Texas 2021 Data/Stata/Birth Defects Data_.dta')

bds %<>% 
  select(STUDYID, starts_with('bpa')) %>% 
  mutate(across(bpa1:bpa63, ~ as.character(.x)))

#' Coerce codes to six digits with a decimal after three.
recode <- function(x){
  
  ifelse(nchar(x) == 2, paste0('0', x, '.000'), 
  ifelse(nchar(x) == 3, paste0(x, '.000'),
  ifelse(nchar(x) == 5, paste0(x, '00'),
  ifelse(nchar(x) == 6, paste0(x, '0'), x))))
  
}

bds %<>% mutate(across(bpa1:bpa63, ~recode(.x)))

#' A flag indicating whether a child has at least one phenotype we're counting towards our list of birth defects.
#' Remove kids without any birth defects that count towards our computed variables.
bds %<>% 
  mutate(flag = ifelse(if_any(bpa1:bpa63, ~str_detect(.x, str_c('^216', '^237.7', '^7[45]', sep = '|'))), 1, 0)) %>% 
  filter(!is.na(flag)) %>% 
  select(-flag)

#' Normally I would sort BD columns from numerically lowest to numerically highest at this point, but it has already been done.
bds %<>% 
  filter(!is.na(bpa1)) %>% 
  mutate(state = 'TX',
         coding.system = 'CDC-BPA',
         STUDYID = paste0('tx', STUDYID)) %>% 
  rename(studyid = STUDYID) %>% 
  select(studyid, state, coding.system, bpa1:bpa63)

#' Same formatting for existing BD codes data frame.
merged.codes %<>% 
  mutate(state = toupper(substr(studyid,1,2)),
         coding.system = 'CDC-BPA') %>% 
  select(studyid, state, coding.system, bpa1:bpa63)

#' Join the two together.
merged.codes %<>% bind_rows(bds)

saveRDS(merged.codes,
        '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/birth_defects_codes_bpa_states.v20220602.rds')

# Corrections to NC codes -------------------------------------------------

codes <- readRDS('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/birth_defects_codes_bpa_states.v20220602.rds')

#' NC codes were not coerced to the uniform length. Some codes are 4 or 5 digits. None are 6 digits or <4 digits.
nc <- codes %>% 
  filter(substr(studyid,1,2) == 'nc') %>% 
  mutate(across(starts_with('bpa'), ~ ifelse(nchar(.x) == 4, paste0(substr(.x,1,3),'.',substr(.x,4,4),'00'),
                                      ifelse(nchar(.x) == 5, paste0(substr(.x,1,3),'.',substr(.x,4,5),'0'), .x))))

codes %<>% filter(state != 'NC')

codes %<>% bind_rows(nc)  

saveRDS(codes, 
        '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/birth_defects_codes_bpa_states.v20220603.rds')

# Discount hemangiomas ----------------------------------------------------

bds <- readRDS('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/birth_defects_codes_bpa_states.v20220603.rds')

#' Previously, we counted hemangiomas as birth defects. 
#' In MA, NC, and OK, children with hemangiomas but no other birth defects were included in the birth defects codes data frame.
#' Remove them to harmonize data structure with newew TX and MI data.
bds %<>% 
  mutate(flag = ifelse(if_any(bpa1:bpa63, ~str_detect(.x, str_c('^216', '^237.7', '^7[45]', sep = '|'))), 1, 0)) %>% 
  filter(!is.na(flag)) %>% 
  select(-flag)

saveRDS(bds,
        '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/birth_defects_codes_bpa_states.v20220604.rds')

# Sratch paper ------------------------------------------------------------
