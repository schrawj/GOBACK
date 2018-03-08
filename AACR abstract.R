#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.03.08.
#' 
#' Code for AACR 2018 Annual Meeting abstract/poster.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------



# Cox models for abstract -------------------------------------------------

setwd('Z:/Jeremy/GOBACK/')
load('./Datasets/Combined Arkansas Michigan and Texas/ar.mi.tx.nochromosomaldefects.v20171106.1.rdata')

#' Chose to highlight astrocytoma and ependymoma as they were both associated with multiple birth defects.
#' Sensitivity analyses for hydrocephalus: exlude children < 1.
tmp <- glm(astro ~ hydrocephalus.wo.spinabifida, data = armitx.nochrom, family = binomial(link='logit'))
tmp <- glm(astro ~ hydrocephalus.wo.spinabifida, data = armitx.nochrom, subset = person.yrs > 1, family = binomial(link='logit'))

tmp <- glm(ependymoma ~ hydrocephalus.wo.spinabifida, data = armitx.nochrom, family = binomial(link='logit'))
tmp <- glm(ependymoma ~ hydrocephalus.wo.spinabifida, data = armitx.nochrom, subset = person.yrs > 1, family = binomial(link='logit'))

rm(tmp, armitx.nochrom); gc()



# Table 1 -----------------------------------------------------------------

require(gmodels); require(dplyr); require(tictoc)

load("Z:/Jeremy/GOBACK/Datasets/goback.no.chrom.v20180122.1.rdata")

#' Revise table 1 to refer only to children with non-chromosomal defects.
CrossTable(goback.nochrom$state, prop.t = FALSE, prop.chisq = FALSE)

vars <- c('any.birthdefect','cancer')

for (i in vars){
  CrossTable(goback.nochrom$state, goback.nochrom[, i], prop.t = FALSE, prop.chisq = FALSE)
}

states <- unique(goback.nochrom$state)

#' Experimenting with the tictoc package.
for (i in states){
  tic('loop')
  print(i)
  tic('filter')
  tmp <- dplyr::filter(goback.nochrom, state == i)
  toc()
  gmodels::CrossTable(tmp$any.birthdefect, tmp$cancer, prop.t = FALSE, prop.chisq = FALSE)
  toc()
}

table(goback.nochrom$any.birthdefect, goback.nochrom$cancer)

rm(goback.nochrom, tmp, i, states, vars); gc()
