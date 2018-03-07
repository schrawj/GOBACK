setwd('Z:/Jeremy/GOBACK/Datasets/')
require(dplyr)


# ALL ---------------------------------------------------------------------
load('./goback.v20171115.3.rdata')
austin <- filter(goback, m.race %in% c('Hispanic','NHW'))
rm(goback)
gc()

tmp <- filter(austin, all == 1 | cancer == 0)
hisp <- filter(tmp, m.race == 'Hispanic')
nothisp <- filter(tmp, m.race == 'NHW')

rm(austin)
gc()

sink(file = 'C:/Users/schraw/Desktop/down-all crosstabs stratified by ethnicity.txt')

gmodels::CrossTable(hisp$down.syndrome, hisp$all, prop.t = FALSE, prop.chisq = FALSE, chisq = TRUE)
gmodels::CrossTable(nothisp$down.syndrome, nothisp$all, prop.t = FALSE, prop.chisq = FALSE, chisq = TRUE)
model <- summary(glm(all ~ down.syndrome + m.race + down.syndrome*m.race, data = tmp, family = binomial(link = 'logit')))
print(model)

sink()

rm(model, hisp, nothisp, tmp)
gc()



# AML ---------------------------------------------------------------------
load('./goback.v20171115.3.rdata')
austin <- filter(goback, m.race %in% c('Hispanic','NHW'))
rm(goback)
gc()

tmp <- filter(austin, aml == 1 | cancer == 0)
hisp <- filter(tmp, m.race == 'Hispanic')
nothisp <- filter(tmp, m.race == 'NHW')

rm(austin)
gc()

sink(file = 'C:/Users/schraw/Desktop/down-aml crosstabs stratified by ethnicity.txt')

gmodels::CrossTable(hisp$down.syndrome, hisp$aml, prop.t = FALSE, prop.chisq = FALSE, chisq = TRUE)
gmodels::CrossTable(nothisp$down.syndrome, nothisp$aml, prop.t = FALSE, prop.chisq = FALSE, chisq = TRUE)
model <- summary(glm(aml ~ down.syndrome + m.race + down.syndrome*m.race, data = tmp, family = binomial(link = 'logit')))
print(model)

sink()

rm(model, hisp, nothisp, tmp)
gc()