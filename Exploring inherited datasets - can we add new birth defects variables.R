#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Philip is interested in adding variables for several additional birth defects.  These
#' are all available in the NC data for example, at the following incidence: 
#' 
#' craniosynostosis, 756.000, .005, .010, .020, .030, (~300 cases);
#' clubfoot, 754.730 (~259 cases); 
#' Turner Syndrome (trisomy 45), 758.690 (12 cases); 
#' DiGeorge Syndrome (deletion of 22q11.2), 279.110 and 758.380 (code for deletion 22q)
#'      (14 cases DiGeorge, 35 cases del22q11); 
#' small intestinal atresia, 751.100, .110, .120, .190, .195 (~265 cases).
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
require(dplyr)

#' Load raw NC data and processed birth defects data.
load("Z:/GOBACK/Jeremy/Datasets/North Carolina/nc.raw.data.rdata")
load("C:/Users/schraw/Downloads/Jeremy/Jeremy/R data files/Datasets/North Carolina/nc.birth.defect.data.v20170821.2.rdata")

names(nc.bdef)

#' They are not computed in the processed file.  Do we have codes for them in the raw data?
nc.raw$DX1 <- ifelse(nc.raw$DX1 == "", NA, nc.raw$DX1)
nc.raw$DX1 <- as.numeric(nc.raw$DX1)
nc.raw <- arrange(nc.raw, DX1)

tmp <- filter(nc.raw, DX1 %in% 756000:756040)
table(tmp$DX1)

tmp <- filter(nc.raw, DX1 %in% 751100:751200)
table(tmp$DX1)

tmp <- filter(nc.raw, DX1 == 754730)

tmp <- filter(nc.raw, DX1 %in% 758000:759000)
table(tmp$DX1)

#' Code for DiGeorge.
tmp <- filter(nc.raw, DX1 == 279110)
table(tmp$DX1)

#' Code for confirmed del22q11.
tmp <- filter(nc.raw, DX1 == 758380)
table(tmp$DX1)



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Philip is interested in categorizing individual CHDs into broader groups per methods
#' that have been published before.
#' 
#' For VSDs and IAA, there is some ambiguity as to what category these might fall under
#' depending on the specifics of the diagnosis.  This is not captured well in standard
#' six digit birth defects codes but is captured under some expanded coding systems used
#' by cardiologists.
#' 
#' Search for codes 747216 and 747517, which distinguish IAA types A or C from type B.
#' 
#' Search for codes in the range 74548x.  If only 745480 is present, we cannot 
#' distinguish subtypes of VSD.
#' 
#' In the NC raw data, we do have codes that specify membranous vs. other types of VSD.
#' 
#' In the NC raw data, we do have codes for IAA, NOS (19 cases) and IAA type A or C 
#' (6 cases).  There are no codes for IAA type B.  It's unclear whether this data is not 
#' available or whether IAA type B equals IAA, NOS - IAA type A/C (13 cases).  
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
tmp <- filter(nc.raw, DX1 %in% c(745480:745490, 746880, 745685))
tmp <- filter(nc.raw, DX1 %in% 747215:747217)
unique(tmp$DX1)

table(tmp$DX1)
