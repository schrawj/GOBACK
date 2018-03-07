require(xlsx)
require(stringr)

def.vars <- read.xlsx(file = 'Z:/GOBAcK/Jeremy/Data dictionaries and data definitions/Variable list.xlsx', sheetName = 'Birth defects variables',
                      header = FALSE, colIndex = 1)
def.vars$X1 <- str_replace_all(def.vars$X1, '_','.')
def.vars$X1 <- str_replace(def.vars$X1, pattern ='\\Q(\\E+NEW+\\Q)\\E', '')
def.vars$X1 <- trimws(def.vars$X1, 'right')
def.vars$X1 <- tolower(def.vars$X1)
def.vars$X1 <- str_replace(def.vars$X1, pattern = '$', ',')

tmp <- c(def.vars$X1)
print(tmp, quote = FALSE, row.names = FALSE)

list.of.standard.birth.defects.variables <- tmp
save(list.of.standard.birth.defects.variables, file = 'Z:/GOBACK/Jeremy/Datasets/list.of.standard.birth.defects.variables.rdata')

can.vars <- read.xlsx(file = 'Z:/GOBAcK/Jeremy/Data dictionaries and data definitions/Variable list.xlsx', sheetName = 'Cancer variables',
                      header = FALSE, colIndex = 1)
can.vars$X1 <- str_replace_all(can.vars$X1, '_','.')
can.vars$X1 <- tolower(can.vars$X1)
can.vars$X1 <- str_replace(can.vars$X1, pattern = '$', ',')

tmp <- c(can.vars$X1)
print(tmp, quote = FALSE, row.names = FALSE)

list.of.standard.cancer.variables <- tmp
save(list.of.standard.cancer.variables, file = 'Z:/GOBACK/Jeremy/Datasets/list.of.standard.cancer.variables.rdata')
