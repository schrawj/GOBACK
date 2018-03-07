
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#'                                       
#'                                       
#'                                       ARKANSAS
#'                                       
#'                                       
#' 2017.08.15.
#' 
#' In advance of my meeting with Heather on 8/17/17 I want to identify all the raw
#' files and finish my 'final data schema' document.
#' 
#' Through conversations with Ruthie, Heather and Tiffany I understand that AR, MI and NC
#' data are all in single files.  TX data is distributed across a number of files.
#' 
#' Old AR data has 629,120 rows.  Birth defects coded by name.
#' 
#' New AR data has 605,778 rows.  Children with birth defects are not included.
#' 
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Prep environment.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
require(dplyr)
require(readstata13)
require(haven)

setwd('C:/Users/schraw/Downloads/')

#' Load in old AR data.  These were received March 2017.
ar.old <- read.csv(file = './GOBACK_Arkansas_20170329.csv', header = TRUE, stringsAsFactors = FALSE)

#' Load in new AR data from the provided STATA file.  These were received July 2017.
ar.new <- read.dta13('./bth9511vs_update.dta', convert.underscore = TRUE)



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Overview of demographic variables.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
names(ar.old)
names(ar.new)

#' Check levels of the demographic variables.
unique(ar.old$Maternal.race.ethnicity)
unique(ar.new$Maternal.race)
unique(ar.new$Infant.sex)

#' There are no labels for the numeric values of maternal race.
#' Its unclear whether there never were any, or if there was a 
#' loss of fidelity with the import.  But the readstata13 command
#' to recover these does not work if they aren't factors.
get.origin.codes(ar.new$Maternal.race, labtab)

#' Maternal education takes some interesting values.  There are 0s, NAs, 46 and 99.
unique(ar.new$Maternal.education)
table(ar.new$Maternal.education)

#' Gestational age has both NA's and 99's as missing value codes.
summary(ar.new$Infant.gestational.age)

#' Birthweight has some outliers.  There are children with a birthweight of 0 and the max 
#' is 12360 (equivalent to 27 lbs).
#' Fortunately it looks like no one tried to use 99 as a missing value code.
summary(ar.new$Infant.birth.weight)
which(ar.new$Infant.birth.weight == 99)

#' Cancer diagnoses are coded according to ICD-9 codes.
sites <- data.frame(site = unique(ar.new$primarysit))
sites <- arrange(sites, site)
print(sites)

print(ar.new[1:25,])
print(ar.old[1:25,])

#' Number of children w/o birth defect in old dataset = number in new dataset.
table(ar.old$Group)



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Overview of birth defects variables.
#' 
#' The birth defects in ar.old were provided by name.  These names correspond to the 
#' diagnostic groupings in the 'bx defects code groups' Tiffany used to harmonize the 
#' other files.
#' 
#' I will need to harmonize the names and order of these variables with her datasets.  
#' For each diagnostic grouping in that Excel sheet, I will also need to generate: 
#' 
#'      1) a variable 'CongAnomalies_[body part]' == 1 if any of the defects in that 
#'      category are present, 0 if not, NA if any other defect present?
#'      
#'      2) a variable 'CongAnomalies_[body part]_Other_Minor' == 1 if another unspecified 
#'      defect present, 0 if not, NA if any other defect present?
#'      
#'      3) a variable 'CongAnomalies_[body part]_Other_Major' == 1 if another unspecified 
#'      defect present, 0 if not, NA if any other defect present?
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

#' Birth defect names from ar.old.
bd.names.ar <- colnames(ar.old[13:69])

#' How do these compare to the computed birth defects?
load("Y:/Jeremy Schraw/GOBACK project/Datasets/michigan.rdata")
bd.names.mi <- colnames(mi.bd.ca[177:277])

unique(mi.bd.ca$CongAnomalies_CNS)
unique(mi.bd.ca$CongAnomalies_CNS_Other_Major)
unique(mi.bd.ca$CongAnomalies_CNS_Other_Minor)

unique(mi.bd.ca$CongAnomalies_Eye)
unique(mi.bd.ca$CongAnomalies_Eye_Other_Major)
unique(mi.bd.ca$CongAnomalies_Eye_Other_Minor)

unique(mi.bd.ca$any_birthdefect)



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Overview of cancer variables.
#' 
#' The variable 'primarysit' is the ICD-O-3.1 topographical code.
#' 
#' The variable 'HISTOLOGY3' is the ICD-O-3.1 morphological code.
#' 
#' Unclear what 'siterecode' represents. 
#' 5 digit numbers taking values between 20000-40000.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

#' Alive/dead?
unique(ar.new$vstatus)

#' Grade and laterality are numbers in character format.  These were probably factors for 
#' which the formatting was lost.
str(ar.new)

#' They have a wider range of values than expected...must represent bilaterality or 
#' uncertainty in laterality/grade.
unique(ar.new$laterality)
unique(ar.new$grade)

#' Looks like these two are the primary cancer dx variables.
unique(ar.new$primarysit)
unique(ar.new$HISTOLOGY3)

#' Levels for other cancer variables.
unique(ar.new$dxconfirm)
unique(ar.new$reportsrce)
unique(ar.new$BEHAVIOR3)
unique(ar.new$siterecode)



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Trying to determine what dates are available, if person-time is calculated, or if it
#' is calculable.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

#' For some reason this is missing in all but 1000 subjects.
ar.new$datelast <- ifelse(ar.new$datelast == "", NA, ar.new$datelast)
table(is.na(ar.new$datelast))

table(is.na(ar.old))



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2017.09.08.
#' 
#' Review birthweight, gestational age and maternal age variables in advance of meeting
#' next Monday.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
require(dplyr)
require(ggplot2)

load("Z:/GOBACK/Jeremy/Datasets/Arkansas/arkansas.v20170831.1.rdata")

tmp <- arrange(ar, birth.wt)
tmp <- arrange(ar, gest.age)
tmp <- arrange(ar, m.age)

tmp <- arrange(ar, desc(birth.wt))
tmp <- arrange(ar, desc(gest.age))
tmp <- arrange(ar, desc(m.age))

print(tmp[1:100, c(1,75,6,3,82,83,12)])
print(tmp[101:200, c(1,75,6,3,82,83,12)])
print(tmp[201:300, c(1,75,3,6,82,83,12)])

print(ggplot(data = ar, aes(x = birth.wt)) + geom_histogram(binwidth = 200, color = 'red', fill = 'white') + xlim(0,6000))

table(ar$gest.age)



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2017.09.11.
#' 
#' Trying to figure out what's up with these 6535 kids in the new data who don't have a 
#' defect recorded.
#' 
#' Going back to the original data.  Did they have these 'other' conditions in the 
#' 200.xxx range of the six-digit birth defect codes?
#' 
#' I guess not.  I don't know what's going on here.  It's all very suspicious.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

ar.mar <- read.csv(file = 'C:/Users/schraw/Downloads/GOBACK_Arkansas_20170329.csv', header = TRUE, stringsAsFactors = FALSE)
save(ar.mar, file = 'Z:/GOBACK/Jeremy/Datasets/Arkansas/arkansas.raw.data.march.rdata')

tmp <- filter(ar.mar, Any.defect == 0)
head(tmp, 5)

table(tmp$Atrial.septal.defect)
table(tmp$Atrioventricular.septal.defect)
table(tmp$Clubfoot)
table(tmp$Spina.bifida.without.anencephalus)
table(tmp$Hydrocephalus.without.Spina.Bifida)



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2017.09.11.
#' 
#' Now investigating whether the etiologic groups described for CHD in Botto et al 2007
#' can be applied to our data.  A possible limitation in our data relates to the 
#' fact that ICD9-based birth defects codes are insufficient to distinguish different
#' forms of interrupted aortic arch and ventricualr septal defect, which we would need
#' to be able to do.  Botto et al modified the codes to allow for this.  Which were we
#' provided from our states?
#' 
#' It won't be possible for Arkansas, which didn't provide any kind of code.
#' It won't be possible for Michigan, which didn't provide any kind of code.
#' It may be possible for North Carolina.
#' As of today I don't know whether its possible for TX, b/c I haven't looked over their
#' new data dump.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
tmp <- filter(nc.bdef, ventricularseptaldefect == 1)
studyid <- c(tmp$studyid)

tmp <- nc.raw[nc.raw$NCID %in% studyid, ]
tmp <- arrange(tmp[!duplicated(tmp$DX1), ], as.numeric(DX1))
tmp <- c(as.numeric(tmp$DX1))

#' Are the DX codes detailed enough?  Refer to Botto et al 2007 Appendix 3.
print(tmp)

#' Turns out they are detailed enough.  How about for IAA?
tmp <- filter(nc.bdef, coarctationofaorta == 1)

#' Unclear whether it will for IAA.  Have codes 747215 and 747216 but no code 747217 
#' (IAA Type B, should represent ~50% of IAA cases).
studyid <- c(tmp$studyid)