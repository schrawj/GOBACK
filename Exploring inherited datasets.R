#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Prep environment.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
setwd('C:/Users/schraw/Downloads/')

require(haven)

require(foreign)
require(readstata13)

require(dplyr)
require(stringr)



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Ruthie said MI and NC data housed in files labelled "xx_append.sas7bdat" are ready to 
#' go.  Read in some files and examine.  Rewrite column names in my preferred style.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

#' Using haven.
nc.append <- read_sas('./nc_append.sas7bdat')
mi.append <- read_sas('./mi_append.sas7bdat')

nc.append <- data.frame(nc.append)
mi.append <- data.frame(mi.append)

colnames(mi.append) <- tolower(colnames(mi.append))
str_locate(colnames(mi.append), '_')
colnames(mi.append) <- str_replace(colnames(mi.append), '_','.')

save(nc.append, file = 'Y:/Jeremy Schraw/BD-CC project/Datasets/nc.append.v20170808.rdata')
save(mi.append, file = 'Y:/Jeremy Schraw/BD-CC project/Datasets/mi.append.v20170808.rdata')

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' After playing around with various files for a while, it looks like the desired final 
#' names and codes for demographic variables are listed in an Excel workbook called 
#' 'TX_NC_MI'.  
#' 
#' The sheet 'Variables for stdrizing' gives the desired final name of a 
#' variable.  The columns to the right give the original names in the various datasets.
#'      green names are from the NC dataset.
#'  
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

print(nc.append[1:100, c(1,123:129)])



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2017.08.15.
#' 
#' In advance of my meeting with Heather this Thursday I want to identify all the raw
#' files and finish my 'final data schema' document.
#' 
#' Through conversations with Ruthie, Heather and Tiffany I understand that AR, MI and NC
#' data are all in single files.  TX data is distributed across a number of files.
#' 
#' Ruthie's handwritten notes suggest 5-6 birth record files for TX alone, and that 
#' cancer registry and birth registry data are in separate files.
#' 
#' I think I've already found the MI and NC data.  Hunt down the AR data and then start 
#' playing with the TX files.  
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Texas.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

#' This file was provided by Tiffany.  It contains birth defects for non-cancer patients.
tx.nc <- read.dta13('./TX_NC_W_Birth Defect Variables_merge.dta', convert.factors = TRUE, 
                    convert.underscore = TRUE)

#' The birth defects variables are coded 1 if yes, NA otherwise.  May need to revise this.
summary(tx.nc$CongAnomalies.CNS)
table(tx.nc$CongAnomalies.CNS)
summary(tx.nc$Encephalocele)

save(tx.nc, file = 'Y:/Jeremy Schraw/GOBACK project/Datasets/tx.nc.birthdefects.noncancer.subjects.rdata')

#' This file as the combined TX and NC cancer patient data.
#' Does not have birth defect data.
tx.nc.ca <- read.csv('./TX_NC_w_Cancer Variables_merge.csv', header = TRUE, stringsAsFactors = FALSE)
colnames(tx.nc.ca) <- str_replace_all(colnames(tx.nc.ca),'_','.')
save(tx.nc.ca, file='Y:/Jeremy Schraw/GOBACK project/Datasets/tx.nc.cancer.rdata')

#' Looks like the maximum number of cancer diagnoses in a single individual is 5 (ouch).
table(tx.nc.ca$ccr.numdiagnoses)

#' The 'cancer' variable has level 1 = 6896.  This must be the number of cancers in the dataset.
table(tx.nc.ca$cancer)

tx.ca <- filter(tx.nc.ca, state == 'TX')
nc.ca <- filter(tx.nc.ca, state == 'NC')

#' I think this is a file Ruthie generated with the standardized demographic variables
#' for TX children.
tx.demo <- read_sas('C:/Users/schraw/Downloads/tx_tiffany.sas7bdat')
colnames(tx.demo) <- tolower(colnames(tx.demo))
colnames(tx.demo) <- str_replace_all(colnames(tx.demo), "_", ".")

save(tx.demo,
     file = 'Y:/Jeremy Schraw/GOBACK project/Datasets/tx.demographics.rdata')

#' Investigate some of the variables

#' Looks like a unique ID.
str(tx.demo$bc.link)
table(duplicated(tx.demo$bc.link))

#' Looks like 5571 cancer cases in TX.
str(tx.demo$cancer)
table(tx.demo$cancer)

#' 5571 TX cancer cases + 1325 cancer cases in NC = 6896 cancers.  
#' Guess I solved that puzzle.
table(nc.ca$cancer)

#' Also appears to have all Ruthie's demographic variables EXCEPT m.edu and f.edu.
#' bc.link stands in for studyid.  birth.yr is named birth.year in this file.
tmp <- select(tx.demo, 
              bc.link, sex, m.race, m.edu2, 
              m.age, birth.wt, gest.age, birth.year,
              plu, f.race, f.edu2, f.age)



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' I cannot find a list of the specific cancer variables we want.  
#' 
#' Instead I have created two forensically by 
#'      1) comparing the common cancer variables from the MI and TX/NC datasets.
#'      2) Pulling all variable names from the MI dataset that are formatted as 
#'      'cancerxyz[1-4].
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
tx.cancers <- colnames(tx.nc.ca)
mi.cancers <- colnames(mi)
common.cancers <- intersect(tx.cancers, mi.cancers)

#' Drop non-cancer variables.
common.cancers <- common.cancers[-c(1,2,11,12,13,14:22)]

#' Was unclear what status was.  Think its alive/dead.  Definitely dichotomous.
commmon.cancers <- common.cancers[-9]

#' Looking at the 'morph3' variables.  These must be ICD-0-3 morphology codes.
str(mi$morph31)
summary(mi$morph31)

#' Looking at the 'site.code' variables.  These must be ICD-0-3 topography codes.
str(mi$site.code1)
summary(mi$site.code1)

#' Remove them for the moment.
common.cancers <- common.cancers[-(1:8)]

#' Pull the variable names for the cancers where serial diagnosis information was
#' computed.
pattern <- '[[:lower:]]+[1]$'
serial.cancers <- grep(pattern, colnames(mi.bd.ca), value = TRUE)
#' Filter non-cancer variables.
serial.cancers <- serial.cancers[-c(2:6,8:10,13:15,26,27,29,30)]

write.table(serial.cancers, file='C:/Users/schraw/Desktop/serial.cancers.txt',quote=FALSE)

