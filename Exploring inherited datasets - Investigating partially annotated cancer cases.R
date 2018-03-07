#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2017.08.22.
#' 
#' While reviewing the NC data I noticed an issue whereby for ~200 children, cancers will
#' take value 1 for a '.any' variable but not for any specific subtype for that category 
#' of tumors.
#' 
#' Presumably, a child with a value 1 for a '.any' variable should also either have a 
#' value 1 for one of the major subtypes, or if not, for the '.other' variable.
#' 
#' E.g., a child with 'leu.any' == 1 should also be flagged as having either all, aml
#' or other leukemia.
#' 
#' It would be useful to figure out what kinds of tumors these children have, and if we 
#' need to make any changes to how their cancer variables are coded.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Prep environment.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
require(dplyr)
require(stringr)
require(xlsx)
require(readstata13)

restring.columns <- function(x){
  str_replace_all(tolower(colnames(x)),'_','.')
}

#' For children with no cancers, set number of diagnoses equal to zero.
nc.cancer$num.diagnoses <- ifelse(is.na(nc.cancer$num.diagnoses), 0, nc.cancer$num.diagnoses)

#' Subset of NC children with exactly 1 cancer DX.  N = 1,305.
nc.cancer.sub2 <- filter(nc.cancer, num.diagnoses == 1)

#' Initialize an empty variable to hold the name of their DX.
nc.cancer.sub2$cancer1 <- NA

#' For each row, search through all specific cancer diagnoses (excluding '.any' variables).  
#' Write the name of their DX to the cancer1 variable.
nc.cancer.sub2$cancer1 <- apply(nc.cancer.sub2[, c(60,61,63:65,67:69,71:72,74:76,78,79,81,82,83,85,86,87,89,91:94,96:99)], 1, function(x){
  names(which(x == 1))
})  

#' Select pertinent columns.
tmp <- nc.cancer.sub2[,c(2,31,16,36,59:99,103)]

#' Change variable type to character; set to NA if missing.
tmp$cancer1 <- as.character(tmp$cancer1)
tmp$cancer1 <- ifelse(tmp$cancer1 == 'character(0)', NA, tmp$cancer1)

#' Arrange the data by cancer1; print last 50 rows.
tmp <- arrange(tmp, cancer1)
tail(tmp, 50)

#' The CDC publishes a list of ICD-O-3 histology and behavior codes and their accompanying descriptions.
#' https://seer.cancer.gov/icd-o-3/
seer.codes <- read.xlsx(file = 'C:/Users/schraw/Downloads/sitetype.icdo3.d20150918.xls', sheetIndex = 1, header = TRUE, 
                        colIndex = 3:6, stringsAsFactors = FALSE)

#' Filter down to one instance of every unique ICD-O-3 histology and behavior code combination.
seer.codes <- seer.codes[!duplicated(seer.codes$Histology.Behavior), ]
#' column names to lowercase, replace underscores with periods.
colnames(seer.codes) <- restring.columns(seer.codes)

save(seer.codes, file = 'Y:/Jeremy Schraw/GOBACK project/Datasets/icdo3.codes.descriptions.from.seer.rdata')

#' Create a histology.behavior column in tmp for a join.
tmp$histology.behavior <- as.character(paste0(tmp$morph31,'/',tmp$behavior.code.icdo3.1))

#' For every row in the NC cancer subset, add the corresponding diagnostic descriptions from the SEER ICD-O-3 sheet, 
#' based on a matching morphology/behavior code.
tmp <- left_join(tmp, seer.codes, by = 'histology.behavior')

#' Select pertinent columns.
tmp <- tmp[, c(1:4,47:50,46,5:45)]

#' Save as a .csv.
write.csv(tmp, file = 'C:/Users/schraw/Desktop/nc.cancers.csv', row.names = FALSE)

cancer1.codes <- tmp[, c(5, 8, 9)]
cancer1.codes <- cancer1.codes[!duplicated(cancer1.codes$histology.behavior), ]

save(cancer1.codes, file = 'Y:/Jeremy Schraw/GOBACK project/Datasets/icod3.codes.nc.rdata')



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2017.08.22.
#' 
#' Separately, we are trying to decide if we might use a laterality variable and are  
#' interested in which cancers laterality is relevant to, beyond retinoblastoma.
#' 
#' SEER requires laterality data for a variety of cancers, listed here by topography:
#' https://seer.cancer.gov/manuals/primsite.laterality.pdf
#' 
#' Using the AR data as a guide, which cancers show up with laterality values of 1 and 2,
#' which indicate right and left origin respectively?
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

#' Load in new AR data from the provided STATA file.  These were received July 2017.
ar.new <- read.dta13('./bth9511vs_update.dta', convert.underscore = TRUE)

#' Filter AR down to rows where tumor is marked as right or left origin.
tmp <- filter(ar.new, laterality == 1 | laterality == 2)

tmp$histology.behavior <- as.character(paste0(tmp$HISTOLOGY3,'/',tmp$BEHAVIOR3))

tmp <- left_join(tmp, seer.codes, by = 'histology.behavior')

tab <- table(tmp$histology.behavior.description, tmp$laterality)

write.csv(tab, file = 'C:/Users/schraw/Desktop/tumors.with.laterality.csv')



head(cancer1.codes)
