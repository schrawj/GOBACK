


# Prep environment --------------------------------------------------------

require(dplyr)



# Arkansas ----------------------------------------------------------------

load("Z:/Jeremy/GOBACK/Datasets/Arkansas/arkansas.v20170913.1.rdata")

mort <- function(dataframe, filtervar){
  tmp <- filter(dataframe, filtervar == 1)
  print(table(is.na(tmp$age.in.months.at.death)))
  print(summary(tmp$age.in.months.at.death))
  rm(tmp)
}

unique(ar$age.in.months.at.death)
summary(ar$age.in.months.at.death)

#' A few situations where we would expect high mortality.
mort(ar, ar$anencephalus)
mort(ar, ar$trisomy13)
mort(ar, ar$cancer)



# Michigan ----------------------------------------------------------------

tmp <- as.data.frame(haven::read_dta(file = 'C:/Users/schraw/Downloads/noncancerbx922011w-wodefects_2x.dta'))
tmp <- select(tmp, -AGE_DIAG)
tmp <- select(tmp, -yeardiag)
tmp <- select(tmp, -tumorid)
tmp <- tmp[,c(1:135,158:175)]


# North Carolina ----------------------------------------------------------

load("Z:/Jeremy/GOBACK/Datasets/North Carolina/nc.v20171107.2.rdata")

mort <- function(dataframe, filtervar){
  tmp <- filter(dataframe, filtervar == 1)
  print(table(is.na(tmp$agedth)))
  print(summary(tmp$agedth))
  rm(tmp)
}

mort(nc, nc$anencephalus)
mort(nc, nc$trisomy13)
mort(nc, nc$cancer)



# Texas -------------------------------------------------------------------

load("Z:/Jeremy/GOBACK/Datasets/Texas/tx.v20171019.2.rdata")

str(tx$dob)
str(tx$dxdate)

#' Attempt to recover DOD, at least for BD children.
demo <- as.data.frame(haven::read_dta(file = 'C:/Users/schraw/Downloads/demo.dta'))
demo <- demo[!duplicated(demo$CASE_ID), ]
demo <- demo[!duplicated(demo$birthID), ]

#' Set this aside.  We may need it.
save(demo, file='Z:/Jeremy/GOBACK/Datasets/Texas/tx.raw.data.dod.in.birth.defects.cases.rdata')

#' Any DOD info in non-BD kids?
birth1113 <- as.data.frame(haven::read_dta(file = 'C:/Users/schraw/Downloads/deid_birth1113.dta'))



