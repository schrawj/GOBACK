
require(dplyr); require(descr); require(tictoc)
require(parallel); require(doParallel); require(foreach)

load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/computed.bd.vars.ncoktx.v20200109.rdata")
load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txncok.v20191205.rdata")

cl <- makePSOCKcluster(7)
registerDoParallel(cl)

vars <- c(2,6:ncol(birth.defects.results))

tic()
out <- foreach (i = vars, .combine = 'rbind') %dopar% {

  new.out <- data.frame(defect = names(birth.defects.results)[i], 
                        n = length(which(birth.defects.results[, i] == 1)))

}
toc()

write.csv(out,
          file = '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/R outputs/computed.bd.vars.ncoktx.counts.v20200109.csv',
          row.names = F)



# For checking new vs old variables ---------------------------------------

fetch <- function(defect) { c(filter(goback, goback[, defect] == 1)$studyid) }
show <- function(ids) { filter(bd.codes.ok.nc.tx, studyid %in% ids) }

holo <- fetch('holoprosencephaly')
holo <- show(holo)

