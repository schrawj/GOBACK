#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#'                                        NORTH CAROLINA
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
load("Y:/Jeremy Schraw/GOBACK project/Datasets/North Carolina/nc.cancer.data.v20170818.1.rdata")
require(dplyr)
require(haven)


#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Has NC provided the years that cancers were diagnosed? 
#' 
#' Variable cancer == 1 indicates cancer.
#' 
#' Columns 143:145 are dxdate[1-3].  All are missing.  
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
tmp <- filter(nc.cancer, cancer == 1)

print(tmp[1:50, c(2,59, 143:145),])

unique(tmp$dx.date1)
unique(tmp$dx.date2)
unique(tmp$dx.date3)



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2017.08.18.
#' 
#' I am in the process of harmonizing the NC data and its very unclear to me whether the 
#' age.at.diagnosis. variables are reliable.  
#' 
#' Age at dx n is frequently less than or equal to age at dx (n-1).
#' The only plausible explanations here are that they're either not in order, or each
#' variable is time from dx n-1 to dx n (the data dictionary does not describe them this
#' way)?
#' 
#' I don't think I can figure this out without loading in the raw data.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

#' Behold...the weirdness of the age.at.diagnosis.variables.
tmp <- nc.cancer[, c(2, 3, 10:15)]
tmp <- filter(tmp, ccr.numdiagnoses > 0)
tmp <- arrange(tmp, desc(ccr.numdiagnoses))
print(tmp)

#' Load in raw NC data.
nc.original <- read_sas('C:/Users/schraw/Downloads/nc_linked_forbcm.sas7bdat')
colnames(nc.original) <- restring.columns(nc.original)
names(nc.original)

nc.original <- select(nc.original, ncid, yearbth, monbth, ccr.numdiagnoses, age.at.diagnosis.1,
                      age.at.diagnosis.2, age.at.diagnosis.3, age.at.diagnosis.4, age.at.diagnosis.5)
tmp <- filter(nc.original, ccr.numdiagnoses > 0)
tmp <- arrange(tmp, desc(ccr.numdiagnoses))
tmp <- data.frame(tmp)
tmp[1:25,]

#' I'm no closer to understanding this.  I'm gonna save it and ask about it at the meeting on Monday.
save(tmp, file = 'C:/Users/schraw/Desktop/NC.cancers.rdata')
