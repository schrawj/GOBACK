#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2017.11.17.
#' 
#' Now that I've recovered continuous birthweights for AR birth defects 
#' kids, combine all GOBACK birthweights into a single data frame.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------


# Prep environment --------------------------------------------------------

require(dplyr)
require(ggplot2)
setwd('Z:/Jeremy/GOBACK/Datasets/')


# Load in and merge data --------------------------------------------------

load("Z:/Jeremy/GOBACK/Datasets/Texas/tx.v20171019.2.rdata")
tx <- select(tx, studyid, birth.wt)

load("Z:/Jeremy/GOBACK/Datasets/Michigan/mi.v20171006.3.rdata")
mi <- select(mi, studyid, birth.wt)
mi$studyid <- paste0('mi',mi$studyid)

load("Z:/Jeremy/GOBACK/Datasets/North Carolina/nc.v20171107.2.rdata")
nc <- select(nc, studyid, birth.wt)  

#' Recovered birthweights for birth defects children.
load("Z:/Jeremy/GOBACK/Datasets/Arkansas/ar.recovered.birthweights.rdata")
recovered.ar.birthweights <- rename(select(recovered.ar.birthweights, studyid, new.birth.wt), birth.wt = new.birth.wt)

#' Provided birthweights for non-birth defects children.
load("Z:/Jeremy/GOBACK/Datasets/Arkansas/arkansas.v20170913.1.rdata")
ar <- select(filter(ar, any.birthdefect == 0), studyid, birth.wt)

goback.bws <- rbind(tx, mi, nc, ar, recovered.ar.birthweights)

rm(ar, mi, nc, tx, recovered.ar.birthweights)

print(ggplot(data = goback.bws, aes(x=birth.wt)) + geom_histogram(color = 'red', fill = 'white'))

goback.bws$birth.wt <- ifelse(goback.bws$birth.wt == 9999, NA, goback.bws$birth.wt)

save(goback.bws, file = './goback.continuous.birthweights.rdata')
