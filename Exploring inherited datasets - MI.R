#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Michigan.
#' 
#' Load a .csv file I generated using the pandas package for Python.
#' This is a copy of the 'Final MI Merge Dataset.dta' file Tiffany put together.
#' Stata -> Python -> CSV -> R is actually more efficient than trying to read in the 
#' original stata file directly using the foreign, haven or readstata13 packages.
#' 
#' This file has all the birth defects variables and (presumably) all the cancer 
#' variables.  I'm still investigating exactly what cancer variables we want.
#' 
#' Prep environment.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
require(dplyr)

restring.columns <- function(x){
  stringr::str_replace_all(tolower(colnames(x)),'_','.')
}

setwd('Z:/Jeremy/GOBACK/Datasets/')
mi.raw <- read.csv('./Michigan/mi.raw.data.v20171003.csv', header = TRUE, stringsAsFactors = FALSE)
save(mi.raw, file = './Michigan/mi.raw.data.v20171003.rdata')



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Go through the standard variables.  Find any that are missing or not coded in the 
#' standard way.
#' 
#' It's likely that some or all of these will also apply to TX, which I have not looked
#' through yet.  For now, just ID any potential problems.  After reviewing TX, decide
#' which need to be solved now and which can be solved after MI is joined to TX.
#' 
#' m.edu2 is continuous.
#' 
#' m.age should have values < 13 or > 50 set to missing.
#' 
#' birth.wt is missing from all non-cancer cases.
#' 
#' birth.wt.cat has all missing birth.wt values coded as >2500 g.  In addition, the 
#' variable has a non-standard name and the levels are not coded according to the 
#' data dictionary.
#' 
#' gest.age needs values < 22 or > 44 set to missing.
#' 
#' plu likely just numeric instead of a factor.  Has value 6 that's not standard.
#' 
#' personyears is useless: it is the number of years between the birth.yr and 2011.
#' 
#' num.diagnoses not computed.  Only first DX recorded in MI data.  Will need to create.
#' 
#' cancer1 has only a subset of diagnoses.  Likely not computed correctly.  
#' 
#' cancer1 should be 'no cancer' for children never DX'd.
#' 
#' The individual cancer variables are only coded as 1 or NA.  Set to 0 if no cancer.
#' 
#' Each of the individual birth defects variable should be set to 0 if the child has no 
#' defects.
#' 
#' minor.status levels given as character instead of numeric.
#' 
#' overall.chrom has character levels instead of numeric.  It is also set to 'no' for 
#' everyone in the dataset.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
mi <- mi.raw
rm(mi.raw)

colnames(mi) <- restring.columns(mi)

str(mi$m.edu2)
table(mi$m.edu2)

str(mi$m.age)
table(mi$m.age)

str(mi$birth.wt)
summary(mi$birth.wt)

str(mi$gest.age)
table(mi$gest.age)

str(mi$birth.yr)
table(mi$birth.yr)

str(mi$plu)
table(mi$plu)

str(mi$f.race)
table(is.na(mi$f.race))
table(is.na(mi$f.edu2))
table(is.na(mi$f.age))

#' Remove a duplicated column (presumably a holdover from some prior table join) that is making dplyr angry.
mi <- mi[,-232]



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' The Michigan time variables are confusing me.
#' 
#' age.diag1: child's age, in years, at cancer diagnosis.  Computed from dxtime1?
#' 
#' cancertime: Approx. age in yrs at DX.  Computed by subtracting year of dx from 
#' birth year.  Used to replace age.diag1 in 4 instances where age.diag1 had illogical
#' values.
#' 
#' cancertime1: Age at diagnosis measured in months?  Values/12 align with age.diag1.
#' 
#' casurvival1: time, in months, from diagnosis to end of study or death?
#' 
#' dxtime1: presumed to be age in months at diagnosis.
#' 
#' personyears: nothing really.  Years between year of birth and 2011.
#' 
#' yeardiag1: year the child was diagnosed with cancer.
#' 
#' yeardx: unknown.
#' 
#' yeardx1: unknown.
#' 
#' The operations that led me to make these determinations are detailed below.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
str(mi$yeardx1)
table(mi$yeardx1) #' Not clear what yeardx and yeardx1 refer to.  Neither appears to have anything to do with cancer.
table(mi$yeardx)
table(mi$yeardiag1) # This is the year of cancer diagnosis.
table(is.na(mi$yeardx1))
table(is.na(mi$yeardx))
table(is.na(mi$yeardiag1))

tmp <- mi[mi$any.cancer == 1, ]
tmp[1:100, c(2,12,65,68,224)]

#' person.years is not calculated correctly.  It's 2011 - birth.yr.  
#' age.diag1 gives a reliable age at cancer diagnosis.
table(mi$person.years, mi$birth.yr)

tmp <- arrange(tmp, desc(person.years))
tmp[1:100, c(2,12,65,68,224)]
tmp[1000:1100, c(2,12,65,68,224)]

one <- c(tmp$age.diag1)
two <- c(tmp$yeardiag1 - tmp$birth.yr)
three <- two - one
table(three)

mi$calculated.person.yrs <- mi$yeardiag1 - mi$birth.yr
mi$person.yrs.diff <- abs(mi$calculated.person.yrs - mi$age.diag1)
table(mi$person.yrs.diff)

#' 4 observations that appear to have unreliable age.diag1.
tmp <- filter(mi, person.yrs.diff > 1)
tmp[ ,c(2,12,65,68,70:72,96,97,224,225)]

list.of.bad.agediag1 <- c(tmp$studyid)

mi$age.diag1 <- ifelse(mi$person.yrs.diff > 1, mi$cancertime,  mi$age.diag1)

tmp <- mi[mi$studyid %in% list.of.bad.agediag1, ]
tmp[ ,c(2,12,65,68,96,97,224,225,239,240)]

#' Problem solved.
mi$person.yrs.diff <- abs(mi$calculated.person.yrs - mi$age.diag1)
table(mi$person.yrs.diff)

tmp <- c(mi$calculated.person.yrs)
tmp2 <- c(mi$cancertime)
tmp3 <- tmp2 - tmp
table(tmp3)



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Inspect cancer variables.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
table(is.na(mi$cancer1))
table(mi$cancer1)

str(mi$cancertime)
str(mi$cancertime1)

mi[1:100, c(2,3,96,225,12,65,68,97,224)]



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Done for the day.  Save the data.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
save(mi, file = './Michigan/mi.v20171003.rdata')



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2017.10.04.
#' 
#' Continue through cancer variables.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
str(mi$laterality1)
unique(mi$laterality1)

unique(mi$cancer1)

#' Look through the individual cancer variables.
tmp[1:100,20:60]

#' Some have non-standard names.
mi <- rename(mi, leu.any = leu, lym.any = lym, pns.any = pns, rms.any = rms, 
             soft.any = soft, gct.any = gct)



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Birth defects variables.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

#' Some have non-standard names.
mi <- rename(mi, double.outlet.right.ventricle = dououtletrightvent,
                 rvot.defects = rvot,
                 lvot.defects = lvot)

#' Look through a set.
mi[1000:1100,98:203]

tmp <- filter(mi, any.birthdefect == 0)
tmp[1000:1100,98:203]

tmp <- filter(mi, any.birthdefect == 1)
tmp[1000:1100,98:203]

#' Review the level 3 grouping variables.
#' Conotruncal defects.
tmp <- filter(mi, conotruncal.defects == 1)
tmp$conotruncal.defects.sum <- rowSums(tmp[,126:129], na.rm = TRUE)
print(tmp[1:100,c(130,126:129,242)])
table(tmp$conotruncal.defects.sum)

#' AVSD.
tmp <- filter(mi, avsd == 1)
table(tmp$endocardialcushiondefect)

#' APVR.
tmp <- filter(mi, apvr == 1)
table(tmp$totalanompulmvenousreturn)

#' LVOT defects.
tmp <- filter(mi, lvot.defects == 1)
tmp$lvot.defects.sum <- rowSums(tmp[,135:138], na.rm = TRUE)
print(tmp[1:100,c(139,135:138,242)])
table(tmp$lvot.defects.sum)

#' Negative contol: All children w/o LVOT defects have NA values for these 4 defects?
tmp <- mi[is.na(mi$lvot.defects), ]
#' Yes.
for (i in 135:138){
  print(table(is.na(tmp[,i])))
}



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Save the data.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
save(mi, file = './Michigan/mi.v20171004.1.rdata')



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Need to update some of the cancer coding.  Some codes for gonadal carcinomas
#' and other and unspecified gonadal tumors were incorrectly mapped to gonadal germ cell
#' tumors.  
#' 
#' The latest round of files Tiffany passed to me don't have the diagnostic codes.
#' 
#' See if I can find them.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
require(readstata13)
mi.can <- read.dta13('Z:/Birth Defects-Childhood Cancer Projects/GOBACK Project (Tiffany)/Michigan/MI_cancer_062017_rg.dta', 
                     convert.factors = TRUE, convert.underscore = TRUE) 
mi.can <- mi.can[,c(1,52:55,68:71,124)]

#' Confirm studyids map to the same individuals by comparing birth weights.
tmp <- left_join(select(mi.can, studyid, birth.wt),
                 select(mi, studyid, birth.wt), by = 'studyid')
tmp$birth.wt.x <- as.numeric(tmp$birth.wt.x)
tmp$birth.wt.diff <- tmp$birth.wt.x - tmp$birth.wt.y
summary(tmp$birth.wt.diff)



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Need to recover birth defects diagnostic codes in order to compute a variable for
#' cleft.palate.and.cleft.lip.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
mi.bd <- read.csv('C:/Users/schraw/Desktop/mi.old.csv', header = TRUE, stringsAsFactors = FALSE)

str(mi.bd$F01)
table(mi.bd$F01)
table(mi.bd$F02)

table(mi.bd$F01, mi.bd$F02)

unique(mi.bd$ICD9COD1)

save(mi.bd, file = 'Z:/Jeremy/GOBACK/Datasets/Michigan/mi.birthdefects.codes.rdata')
