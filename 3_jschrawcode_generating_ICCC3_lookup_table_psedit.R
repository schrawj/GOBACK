#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Pagna sent this on 2019.12.11, saying it was an edited version of the 
#' code I wrote to generate an ICCC3 lookup table for childhood cancer
#' diagnoses.
#' 
#' I haven't reviewed what's different about it in any level of detail.
#' 
#' He said he got the same results as when using mine. So although they're 
#' not all that orthogonal, I suppose that's a good sign.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

library(tidyr); library(dplyr); library(data.table)
workingdir <- 'Z:/PagnaSok/12_COG_CCRN/morphology_topography_reference/'

#' For web scraping.
require(xml2); require(rvest)
#' For data manipulation.
require(stringr); require(dplyr)

page <- read_html('https://seer.cancer.gov/iccc/iccc3_ext.html')
codes <- html_table(page, fill = T)[[1]]

#' Extract all Roman numerals from the site group column.
codes$Group <- str_extract(codes$`Site Group`, '^[IVX]+')
icccgroup <- codes %>% filter(is.na(Group)==F) %>% select(`Site Group`, Group)

codes <- codes %>% mutate(Group = replace(Group,
                                          Group %in% icccgroup$Group,
                                          icccgroup$`Site Group`))

#' If there's a non-missing value in the new "group" column, record that value.
#' Fill in all missing instances with that value until a new non-missing value is encountered.
for (i in 1:nrow(codes)) {

  if (is.na(codes$Group[i]) == F){

    group <- as.character(codes$Group[i])

  }

  else {

    codes$Group[i] <- group

  }
}

codes <- filter(codes, grepl('^\\d', codes$`ICD-O-3 Histology (Type)`))

codes$`ICD-O-3 Histology (Type)` <- str_replace_all(codes$`ICD-O-3 Histology (Type)`, '-', ':')
codes$`ICD-O-2/3 Site` <- str_replace_all(codes$`ICD-O-2/3 Site`, '-', ':')
codes$`ICD-O-2/3 Site` <- str_remove_all(codes$`ICD-O-2/3 Site`, 'C')

for (i in 1:nrow(codes)){

  #' Generate a vector of histology codes.
  index.hist.codes <- codes$`ICD-O-3 Histology (Type)`[i]
  #index.hist.codes <- unlist(regmatches(index.hist.codes, gregexpr('(?<!:)[[:digit:]]{4}(?!:)', perl =T, index.hist.codes)))

  index.hist.codes.character <- c(unlist(strsplit(index.hist.codes, ",")))
  index.hist.codes <- as.numeric(unlist(strsplit(index.hist.codes, ",")))
  index.hist.codes <- index.hist.codes[complete.cases(index.hist.codes)]

  n = sum(str_detect(index.hist.codes.character, ':'))

  if (n>0){

    for (j in 1:n){
      list <- index.hist.codes.character[str_detect(index.hist.codes.character, ':')][j]
      num1 <- as.numeric(str_split(list,':')[[1]][1])
      num2 <- as.numeric(str_split(list,':')[[1]][2])
      num <- num1:num2

      index.hist.codes <- c(index.hist.codes, num)

    }

  }

  index.hist.codes <- paste(index.hist.codes, collapse=',')

  codes$`ICD-O-3 Histology (Type)`[i] <- index.hist.codes



  index.site.codes <- codes$`ICD-O-2/3 Site`[i]
  index.site.codes.character <- c(unlist(strsplit(index.site.codes, ',')))
  index.site.codes <- as.numeric(unlist(strsplit(index.site.codes, ",")))
  index.site.codes <- index.site.codes[complete.cases(index.site.codes)]

  n = sum(str_detect(index.site.codes.character, ':'))

  if (n>0){

    for (j in 1:n){
      list <- index.site.codes.character[str_detect(index.site.codes.character, ':')][j]
      num1 <- as.numeric(str_split(list,':')[[1]][1])
      num2 <- as.numeric(str_split(list,':')[[1]][2])
      num <- num1:num2

      index.site.codes <- c(index.site.codes, num)

    }

  }

  index.site.codes <- paste(index.site.codes, collapse=',')

  codes$`ICD-O-2/3 Site`[i] <- index.site.codes

}

codes2 <- codes %>% mutate(`ICD-O-3 Histology (Type)` = strsplit(as.character(`ICD-O-3 Histology (Type)`), ",")) %>% unnest(`ICD-O-3 Histology (Type)`)
codes2 <- codes2 %>% mutate(`ICD-O-2/3 Site` = strsplit(as.character(`ICD-O-2/3 Site`), ',')) %>% unnest(`ICD-O-2/3 Site`)
codes2$`ICD-O-2/3 Site` <- ifelse(nchar(codes2$`ICD-O-2/3 Site`) == 1, paste0('C00',codes2$`ICD-O-2/3 Site`),
                                      ifelse(nchar(codes2$`ICD-O-2/3 Site`) == 2, paste0('C0',codes2$`ICD-O-2/3 Site`),
                                             paste0('C',codes2$`ICD-O-2/3 Site`)))

fwrite(codes2, paste0(workingdir,'4_jschrawcode_iccc_table_psedit.csv'))
