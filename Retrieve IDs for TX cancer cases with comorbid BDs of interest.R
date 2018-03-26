
# Recovered IDs from Excel sheets for comorbid cases ----------------------

require(xlsx)

setwd('Z:/Jeremy/GOBACK/R outputs/Cancer codes in comorbid cases/')

files <- Sys.glob('*.xlsx')

for (i in files){
  
  wb <- loadWorkbook(paste0('Z:/Jeremy/GOBACK/R outputs/Cancer codes in comorbid cases/',i))
  sheets <- getSheets(wb)
  sheets <- names(sheets)  
  
  for (j in sheets){
    
  tmp <- read.xlsx(file = paste0('Z:/Jeremy/GOBACK/R outputs/Cancer codes in comorbid cases/',i), sheetIndex = j, colIndex = 1, stringsAsFactors = FALSE)
  tmp <- c(tmp[,1])
  ids <- c(ids, tmp)
  
  }

}

rm(i, j, files, sheets, tmp, wb)

ids <- ids[grepl('^tx', ids, ignore.case = TRUE)]



# Remove leading text and convert to data frame ---------------------------

require(stringr); require(readstata13); require(dplyr)

ids <- str_replace(ids, 'tx','')

id.mappings <- data.frame(birthID = as.numeric(ids))

#' Read in demographic variables for TX BD cases.
demo <- read.dta13('C:/Users/schraw/Downloads/demo.dta')

id.mappings <- left_join(id.mappings, select(demo, birthID, CASE_ID), by = 'birthID')

#' Read in identifiers from the TX cancer registry.
tcr <- read.dta13('C:/Users/schraw/Downloads/tcr.dta')

id.mappings <- arrange(
                        rename(
                              left_join(id.mappings, select(tcr, patientid, sex, birthdate, birthid), by = c('birthID' = 'birthid')),
                        birthid = birthID, caseid = CASE_ID),
                      birthid)

#' No IDs missing anywhere.  
for (i in 1:4){
  print(table(is.na(id.mappings[,i])))
}

#' IDs are duplicated if they had more than cancer and more than one of the selected defects.
#' Remove duplicated rows.
id.mappings <- id.mappings[!duplicated(id.mappings$birthid), ]

write.csv(id.mappings, file = 'C:/Users/schraw/Desktop/tx.id.mappings.csv', row.names = FALSE)

rm(demo, tcr, i, ids); gc()
