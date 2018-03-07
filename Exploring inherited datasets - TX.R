#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2017.10.03.
#' 
#' After our last group meeting for GOBACK, Tiffany cleaned up the cancer codes and 
#' computed the new and revised birth defects variables we discussed.
#' 
#' She uploaded two final datasets for review.  Load them in and inspect.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
require(dplyr)

restring.columns <- function(x){
  stringr::str_replace_all(tolower(colnames(x)),'_','.')
}

tx.raw <- read.csv('Z:/Jeremy/GOBACK/Datasets/Texas/tx.raw.data.csv', header = TRUE, stringsAsFactors = FALSE)
save(tx.raw, file = 'Z:/Jeremy/GOBACK/Datasets/Texas/tx.raw.data.rdata')



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Earlier today (2017.10.04) I identified several issues with the MI data that likely
#' apply to TX as well.  
#' 
#' Go through the TX data.  Determine whether which of the above issues applies here, and
#' whether any issues exist in TX that do not in MI.  Develop a plan for fixing them and
#' joining the data in as few steps as possible.
#' 
#' m.age should have values < 13 or > 50 set to missing.
#' 
#' There is no birth.wt.cat variable.
#' 
#' gest.age needs values < 22 or > 44 set to missing.
#' 
#' plu numeric instead of a factor.  Has values 6, 7, 8 that are not standard.
#' 
#' f.race does not include any values = 1 (Hispanic).  No ethnicity for father?
#' 
#' There is no f.edu2 variable.  This should have been in the raw data except for 2004.
#' 
#' f.age has values as low as 10 and as high as 90.
#' 
#' personyears is useless: it is the number of years between the birth.yr and 2014.
#' 
#' cancertime1 is computed, presumably in months.  I don't know what it was computed from, 
#' however, so I don't trust it.
#' 
#' Any.birthdefect only takes values 1 and NA.  Should be 0 if there is no defect.
#' 
#' Individual birth defects variables are NA for children with no birth defects.
#' Should be zero.
#' 
#' minor.status levels given as character instead of numeric.
#' 
#' overall.chrom levels given as character instead of numeric.
#' 
#' majordefect.total and minordefect.total are both NA for children without defects.
#' 
#' num.diagnoses is NA for any child without cancer.
#' 
#' cancer1 has only a subset of diagnoses.  Likely not computed correctly.  Issue 
#' appears to be precisely the same as in MI.
#' 
#' cancer1 should be 'no cancer' for children never DX'd.
#' 
#' The inidvidual cancer variables are only coded as 1 or NA.  Set to 0 if no cancer.
#' 
#' Each of the individual birth defects variable should be set to 0 if the child has no 
#' defects.
#' 
#' Three birth defects and six cancer variables have non-standard names.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
tx <- tx.raw
colnames(tx) <- restring.columns(tx)
rm(tx.raw)

#' Create studyid variable.
tx$studyid <- paste0(tx$state, tx$X)

save(tx, file = './Texas/tx.v20171004.1.rdata')



# Demographics variables --------------------------------------------------
str(tx$state)

tmp <- tx[!duplicated(tx$X), ]

str(tx$sex)
table(tx$sex)

str(tx$m.race)
table(tx$m.race)


str(tx$m.edu)
table(tx$m.edu)

str(tx$m.age)
summary(tx$m.age)
unique(tx$m.age)

str(tx$birth.wt)
summary(tx$birth.wt)

str(tx$gest.age)
summary(tx$gest.age)
print(unique(tx$gest.age))

str(tx$birth.year)
table(tx$birth.year)

str(tx$plu)
table(tx$plu)

str(tx$f.race)
table(tx$f.race)

str(tx$f.age)
summary(tx$f.age)
unique(tx$f.age)
# Person-years and cancertime ------------------------------------------------------------
str(tx$person.years)
table(tx$person.years)

tx[1:100,c(197,11,13)]

tmp <- filter(tx, cancer == 1)
tmp[1:100,c(197,11,13,153)]

table(is.na(tx$cancertime1))

#' We have the same problem as in MI.
#' Person years only gives us time from birth to end of study (2014 in this case)?
table(tx$birth.year, tx$person.years)





# Birth defects variables -------------------------------------------------
str(tx$any.birthdefect)
unique(tx$any.birthdefect)

str(tx$minor.status)
unique(tx$minor.status)

str(tx$overall.chrom)
unique(tx$overall.chrom)

str(tx$majordefect.total)
table(tx$majordefect.total)
table(is.na(tx$majordefect.total))

str(tx$minordefect.total)
table(tx$minordefect.total)
table(is.na(tx$minordefect.total))

table(tx$majordefect.total, tx$any.birthdefect)

#' Verified that every child flagged by any.birthdefect has at least one major or minor defect.  
#' They do.
tmp <- tx
tmp$any.defect <- ifelse(tx$majordefect.total >= 1 | tx$minordefect.total >= 1, 1, 0)
table(tmp$any.defect, tmp$any.birthdefect)

tmp <- filter(tx, any.birthdefect == 1)
print(tmp[1:100, 27:152])

tmp <- tx[is.na(tx$any.birthdefect), ]
print(tmp[1:100, 26:50])

#' Is cleftpalateandcleftlip a catch-all?  Or is the co-occurrence of those two.
tmp <- filter(tx, cleftpalateandcleftlip == 1)
print(tmp[1:100,83:88])
tmp$comorbid.cleft <- rowSums(tmp[,84:85], na.rm = TRUE)
table(tmp$comorbid.cleft)

#' We have 3,043 cases of cleft.palate.wo.cleft.lip who also have cleft.palate.and.cleft.lip.
#' The only way this makes sense is if cleft.palate.and.cleft.lip is a catch-all category.
table(tmp$cleft.palate.wo.cleft.lip)
tmp2 <- filter(tmp, comorbid.cleft == 2)
print(tmp2[,83:88])



# Cancer variables (other than cancertime) --------------------------------
str(tx$cancer)
table(tx$cancer)
table(is.na(tx$cancer))

str(tx$num.diagnoses)
table(tx$num.diagnoses)

str(tx$cancer1)
unique(tx$cancer1)

str(tx$laterality1)
table(tx$laterality1)




# Trying to locate a file with f.edu and f.race ---------------------------
tx.bc <- read.csv(file = 'C:/Users/schraw/Desktop/tx.bc.csv', header = TRUE, stringsAsFactors = FALSE)
tx.bc <- readstata13::read.dta13('Z:/Birth Defects-Childhood Cancer Projects/GOBACK/Texas Birth Record Data/birthrecords.dta')
tx.bc.five <- readstata13::read.dta13('Z:/Birth Defects-Childhood Cancer Projects/Updated Texas Data/Tiffany/TX05.dta')
tx.bc <- readstata13::read.dta13('Z:/Birth Defects-Childhood Cancer Projects/Updated Texas Data/Tiffany/TXall.dta')
tx.bc <- haven::read_sas(data_file = 'Z:/Birth Defects-Childhood Cancer Projects/Updated Texas Data/Tiffany/txall_addedvar.sas7bdat')


# Trying to locate dates for cancer dates, codes --------------------------
tx.can <- readstata13::read.dta13('Z:/Birth Defects-Childhood Cancer Projects/GOBACK Project (Tiffany)/Texas/TX_Cancer_added can variables.dta')
colnames(tx.can) <- restring.columns(tx.can)
tx.can <- select(tx.can, birthid, morph31, site.code1, morph32, site.code2, morph33, site.code3, morph34, site.code4)
save(tx.can, file = 'Z:/Jeremy/GOBACK/Datasets/Texas/tx.cancer1.codes.rdata')
