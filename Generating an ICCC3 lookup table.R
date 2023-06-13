
# For Pagna ---------------------------------------------------------------

#' For web scraping.
require(xml2); require(rvest); require(magrittr)

#' For data manipulation.
require(tidyverse)

#page <- read_html('https://seer.cancer.gov/iccc/iccc3_ext.html')
page <- read_html('https://seer.cancer.gov/iccc/iccc-iarc-2017.html')

codes <- html_table(page)[[1]]

names(codes) %<>% tolower() %>% str_replace_all(' ', '.') %>% str_replace_all('-','.')

#' Fails to record site group for one row in particular.
new.codes <- data.frame(site.group = rep('(a.4) Lymphoid leukemia, NOS', 2),
                        icd.o.3.histology = c('9591','9820'),
                        icd.o.3.primary.site = c('420, 421, 423, 424','000-809'),
                        icd.o.3.behavior = rep('3', 2),
                        extended.recode = rep('004',2),
                        regular.recode = rep('011',2))

codes <- bind_rows(slice(codes,1:8), new.codes, slice(codes, 10:nrow(codes)))

#' Extract all Roman numerals from the site group column.
codes$group <- str_extract(codes$site.group, '^[IVX]+')

#' Extract subgroup information in parentheses.
pattern <- '^[(].{1,4}[)]'

codes$subgroup <- str_remove_all(str_extract(codes$site.group, pattern), '[()]')

#' If there's a non-missing value in the new "group" column, record that value. 
#' Fill in all missing instances with that value until a new non-missing value is encountered.
for (i in 1:nrow(codes)) {
  
  if (is.na(codes$group[i]) == F){
    
    group <- as.character(codes$group[i])
    
  }
  
  else {
    
    codes$group[i] <- group
    
  }
}

#' This line will remove the header rows for groups (I-XII) and subgroups (e.g., a-f) 
#' that are interspersed throughout the table.
codes <- filter(codes, grepl('^\\d', codes$icd.o.3.histology))

codes$icd.o.3.histology <- str_replace_all(codes$icd.o.3.histology, '-', ':')

codes$icd.o.3.primary.site <- str_replace_all(codes$icd.o.3.primary.site, '-', ':')

#' Initialize object to hold long-format codes.
iccc.lookup <- matrix(nrow = 0, ncol = 6)

for (i in 1:nrow(codes)){
  
  index.iccc.code <- stringr::str_trim(stringr::str_remove(codes$site.group[i], pattern), 'both')
  
  index.site.group <- codes$group[i]
  
  index.subgroup <- codes$subgroup[i]
  
  index.extended.class <- codes$extended.recode[i]

#' Generate a vector of histology codes.
#' Negative lookaheads and lookbehinds to separately find 4 digit strings which are and are not adjacent to a colon. 
#' Those which are will have to be "flattened out" and appended to those which are not, so we create a second object to hold them.
  index.hist.codes <- codes$icd.o.3.histology[i]
  
  code.sequences <- c(unlist(regmatches(index.hist.codes, gregexpr('[[:digit:]]{4}:[[:digit:]]{4}', perl =T, index.hist.codes))))
  
  index.hist.codes <- as.numeric(c(unlist(regmatches(index.hist.codes, gregexpr('(?<!:)[[:digit:]]{4}(?!:)', perl =T, index.hist.codes)))))
  
  for (j in seq_along(code.sequences)) {

      start <- as.numeric(substr(code.sequences[j], 1, 4))
      
      stop <- as.numeric(substr(code.sequences[j], 6, 9))
      
      sequence <- start:stop

      index.hist.codes <- c(index.hist.codes, sequence)
      
    }
    
#' Generate a vector of site codes.
  index.site.codes <- codes$icd.o.3.primary.site[i]

  index.site.codes <- c(unlist(regmatches(index.site.codes, gregexpr('(?<!:)[[:digit:]]{3}(?!:)', perl =T, index.site.codes))),
                        unlist(str_match_all(index.site.codes, '[[:digit:]]{3}:[[:digit:]]{3}')))
  
  index.site.codes.numeric <- as.numeric()
  
  for (j in seq_along(index.site.codes)){
    
    if (str_length(index.site.codes)[j] == 7){
      
      start <- substr(index.site.codes[j], 1, 3)
      
      end <- substr(index.site.codes[j], 5, 7)
      
      sequence <- start:end
      
      index.site.codes.numeric <- c(index.site.codes.numeric, sequence)
      
    }

    else if (str_length(index.site.codes[j]) == 3){
      
      index.site.codes.numeric <- c(index.site.codes.numeric, as.numeric(index.site.codes[j]))
      
    }
  
  }
  
  new.codes <- data.frame(iccc.group = rep(index.site.group, times = length(index.hist.codes)*length(index.site.codes.numeric)),
                          iccc.subgroup = rep(index.subgroup, times = length(index.hist.codes)*length(index.site.codes.numeric)),
                          iccc.name = rep(index.iccc.code, times = length(index.hist.codes)*length(index.site.codes.numeric)),
                          iccc.extended.classification = rep(index.extended.class, each = length(index.site.codes.numeric)),
                          icdo3.hist = rep(index.hist.codes, each = length(index.site.codes.numeric)),
                          icdo3.site.code = rep(index.site.codes.numeric, times = length(index.hist.codes)))
  
  iccc.lookup <- rbind(iccc.lookup, new.codes)
  
}

iccc.lookup$icdo3.site.code <- as.character(iccc.lookup$icdo3.site.code)
iccc.lookup$icdo3.site.code <- ifelse(nchar(iccc.lookup$icdo3.site.code) == 1, paste0('C00',iccc.lookup$icdo3.site.code), 
                               ifelse(nchar(iccc.lookup$icdo3.site.code) == 2, paste0('C0',iccc.lookup$icdo3.site.code), 
                                                                               paste0('C',iccc.lookup$icdo3.site.code)))

#' Minor fix. Some Roman numerals retaind in name column.
iccc.lookup$iccc.name <- str_trim(str_remove(iccc.lookup$iccc.name, '^[IVX]{1,4}\\Q.\\E'), 'both')

save(iccc.lookup, file = '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/iccc.lookup.table.v20210826.rdata')

rm(list = ls()); gc()

# For Casey ---------------------------------------------------------------

#' Saved code sheet with .xlsx extension.
codes <- read.xlsx(file = 'C:/Users/schraw/Downloads/Copy of sitetype.icdo3.20180323.xlsx', sheetIndex = 1, stringsAsFactors = F)

codes$Site.recode.numeric <- str_remove_all(codes$Site.recode, 'C')
codes$Site.recode.numeric <- str_replace_all(codes$Site.recode.numeric, '-', ':')

#' For testing different unique cases. 
codes.subset <- codes[!duplicated(str_length(codes$Site.recode.numeric)), ]

#' Initalize an empty object to hold the long format codes data.
long.codes <- matrix(nrow = 0, ncol = 8)

for (i in 1:nrow(codes.subset)){
  
  index.codes <- codes.subset$Site.recode.numeric[i]
  
  new.values <- c(unlist(regmatches(index.codes, gregexpr('(?<!:)[[:digit:]]{3}(?!:)', perl =T, index.codes))),
                  unlist(str_match_all(index.codes, '[[:digit:]]{3}:[[:digit:]]{3}')))

  for (j in seq_along(new.values)){
    
    if (str_length(new.values)[j] == 7){
      
      start <- substr(new.values[j], 1, 3)
      
      end <- substr(new.values[j], 5, 7)
      
      sequence <- start:end
      
      new.data <- cbind(dplyr::slice(codes.subset, rep(i, times = length(sequence))), sequence)
      new.data <- rename(new.data, icdo3.site.code = sequence)
      
      long.codes <- rbind(long.codes, new.data)
      
    }
    
    else if (str_length(new.values[j]) == 3){
      
      new.data <- cbind(dplyr::slice(codes.subset, i), new.values[j])
      new.data <- dplyr::rename(new.data, icdo3.site.code = `new.values[j]`)
      
      long.codes <- rbind(long.codes, new.data)
      
    } 

  }
  
}

print(select(long.codes, -Site.recode, -Histology.Behavior.Description))

