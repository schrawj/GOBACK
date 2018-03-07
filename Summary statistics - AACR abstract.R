#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.01.09.
#' 
#' Summary statistics for 2018 AACR Annual Meeting abstract.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

#' Chose to highlight astrocytoma and ependymoma as they were both associated with multiple birth defects.
#' Sensitivity analyses for hydrocephalus: exlude children < 1.
tmp <- glm(astro ~ hydrocephalus.wo.spinabifida, data = armitx.nochrom, family = binomial(link='logit'))
tmp <- glm(astro ~ hydrocephalus.wo.spinabifida, data = armitx.nochrom, subset = person.yrs > 1, family = binomial(link='logit'))

tmp <- glm(ependymoma ~ hydrocephalus.wo.spinabifida, data = armitx.nochrom, family = binomial(link='logit'))
tmp <- glm(ependymoma ~ hydrocephalus.wo.spinabifida, data = armitx.nochrom, subset = person.yrs > 1, family = binomial(link='logit'))
