

# For Pagna ---------------------------------------------------------------

#' For web scraping.
require(xml2); require(rvest)

#' For data manipulation.
require(xlsx); require(stringr); require(dplyr)

page <- read_html('https://seer.cancer.gov/iccc/iccc3_ext.html')

codes <- html_table(page, fill = T)[[1]]

codes <- filter(codes, grepl('^\\d', codes$`ICD-O-3 Histology (Type)`))

codes$`ICD-O-3 Histology (Type)` <- str_replace_all(codes$`ICD-O-3 Histology (Type)`, '-', ':')

codes$`ICD-O-2/3 Site` <- str_replace_all(codes$`ICD-O-2/3 Site`, '-', ':')
codes$`ICD-O-2/3 Site` <- str_remove_all(codes$`ICD-O-2/3 Site`, 'C')

#' Initialize object to hold long-format codes.
iccc.lookup <- matrix(nrow = 0, ncol = 3)

for (i in 1:nrow(codes)){
  
  index.iccc.code <- codes$`Site Group`[i]

#' Generate a vector of histology codes.
  index.hist.codes <- codes$`ICD-O-3 Histology (Type)`[i]
  index.hist.codes <- unlist(regmatches(index.hist.codes, gregexpr('(?<!:)[[:digit:]]{4}(?!:)', perl =T, index.hist.codes)))
  
#' Generate a vector of site codes.
  index.site.codes <- codes$`ICD-O-2/3 Site`[i]
  
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
  
  new.codes <- data.frame(iccc.site.group = rep(index.iccc.code, times = length(index.hist.codes)*length(index.site.codes.numeric)),
                          icd03.hist = rep(index.hist.codes, each = length(index.site.codes.numeric)),
                          icdo3.site.code = rep(index.site.codes.numeric, times = length(index.hist.codes)))
  
  iccc.lookup <- rbind(iccc.lookup, new.codes)
  
}

iccc.lookup$icdo3.site.code <- as.character(iccc.lookup$icdo3.site.code)
iccc.lookup$icdo3.site.code <- ifelse(nchar(iccc.lookup$icdo3.site.code) == 1, paste0('C00',iccc.lookup$icdo3.site.code), 
                               ifelse(nchar(iccc.lookup$icdo3.site.code) == 2, paste0('C0',iccc.lookup$icdo3.site.code), 
                                                                               paste0('C',iccc.lookup$icdo3.site.code)))

save(iccc.lookup, file = '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/iccc.lookup.table.v20191007.rdata')



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

