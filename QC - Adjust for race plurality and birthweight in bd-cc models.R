#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.10.22.
#' 
#' Reviewer comments on the NEJM draft were something of a mixed bag 
#' regarding how we should treat covariates. Reviewer 1 suggested we should
#' also adjust for race-ethnicity and plurality, and should include 
#' birthweight in more models. 
#' 
#' Reviewer 2 says adjusting for birthweight is unnecessary.
#' 
#' I will randomly select 10% of our  600 modeled associations, adjust them
#' for maternal race, birthweight, and singleton vs. twin or higher then 
#' compare the results with baseline. If they are not much different, we
#' move on. If they are different, Philip and I can discuss.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

# Import and choose models ------------------------------------------------

require(xlsx); require(dplyr); require(survival)

#' Choose 10% of models for sensitivity analysis.
models <- read.xlsx(file = 'W:/Old_genepi2/Jeremy/GOBACK/R outputs/goback.cox.models.v20180119.xlsx',
                    header = TRUE, stringsAsFactors = FALSE, sheetIndex = 1)
names(models) <- tolower(names(models))
#' Flag any defects that should be modeled in the syndromic dataset. For simplicity, will call these ineligible.
models$use.syndromic.cases <- ifelse(models$defect %in% c('down.syndrome','conganomalies.chromosomalanomalies.other','chromosomalanomalies',
                                                          'di.george.syndrome','trisomy18'), 1, 0)
models <- filter(models, use.syndromic.cases == 0)

set.seed(109)

models$random.uniform <- runif(nrow(models), 0, 1)
models <- arrange(models, random.uniform)
models <- select(models[1:60,], defect, cancer, hr, ci.lower, ci.upper, p.val.coef, random.uniform)

defects.non.syndromic <- c(models$defect)
cancers.non.syndromic <- c(models$cancer)



# Load in GOBACK data and compute new variables ---------------------------

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.nochrom.v20180829.rdata')

#' There are 192,370 births in NC to women of AIAN ancestry. WTF is that about?
table(goback.nochrom$m.race, goback.nochrom$state, useNA = 'ifany')

goback.nochrom$twin <- factor(ifelse(goback.nochrom$plu > 1, 1, 0),
                              levels = c(0,1),
                              labels = c('Singleton','Multiple'))

#' TODO: Why are there so many AIAN births in NC? Should this be Hispanic?
goback.nochrom$m.race.collapsed <- factor(ifelse(goback.nochrom$m.race %in% c('Hispanic','NHW'),0,
                                            ifelse(goback.nochrom$m.race == 'NHB', 1, 
                                               ifelse(goback.nochrom$m.race %in% c('Asian','Other','AIAN','Unknown'), 2, 3))),
                                          levels = c(0:2),
                                          labels = c('White','Black','Other'))

goback.nochrom <- goback.nochrom[,c(1,3,6,7,15:94,106,109,112:151,156:158)]



# Baseline models ---------------------------------------------------------

baseline.models.nonsyndromic <- data.frame(defect = as.character(),
                                           cancer = as.character(),
                                           defect.hazard.ratio = as.numeric(),
                                           defect.ci.lower = as.numeric(),
                                           defect.ci.upper = as.numeric(),
                                           n.comorbid = as.numeric())

for (i in 1:length(defects.non.syndromic)){
  
  print(paste('Computing model', i, ':' ,'state-, sex-, and maternal age-adjusted model for risk of', cancers.non.syndromic[i], 'in kids with', defects.non.syndromic[i]))
  
  goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                            cancer = goback.nochrom[,cancers.non.syndromic[i]],
                            defect = goback.nochrom[,defects.non.syndromic[i]],
                            sex = factor(goback.nochrom$sex,
                                         levels = c(1,2),
                                         labels = c('Male','Female')),
                            m.age = goback.nochrom$m.age,
                            state = goback.nochrom$state.num)
  
  cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
  cox.coef <- summary(cox)$coefficients

  estimates <- data.frame(defect = defects.non.syndromic[i], 
                          cancer = cancers.non.syndromic[i], 
                          defect.hazard.ratio = exp(cox.coef[1,1]), 
                          defect.ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                          defect.ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                          n.comorbid = as.numeric(table(goback.nochrom[,defects.non.syndromic[i]], goback.nochrom[,cancers.non.syndromic[i]])[2,2]))
  
  baseline.models.nonsyndromic <- rbind(baseline.models.nonsyndromic, estimates)
  
}

rm(cox.coef, estimates, goback.surv, models, cox, i)



# Also adjusted for race plurality and birthweight ------------------------

adjusted.models.nonsyndromic <- data.frame(defect = as.character(),
                                           cancer = as.character(),
                                           defect.hazard.ratio = as.numeric(),
                                           defect.ci.lower = as.numeric(),
                                           defect.ci.upper = as.numeric(),
                                           plu.hazard.ratio = as.numeric(),
                                           plu.ci.lower = as.numeric(),
                                           plu.ci.upper = as.numeric(),
                                           birthweight.p.value = as.numeric(),
                                           race.nhb.hazard.ratio = as.numeric(),
                                           race.nhb.ci.lower = as.numeric(),
                                           race.nhb.ci.upper = as.numeric(),
                                           race.other.hazard.ratio = as.numeric(),
                                           race.other.ci.lower = as.numeric(),
                                           race.other.ci.upper = as.numeric(),
                                           n.comorbid = as.numeric())

for (i in 1:length(defects.non.syndromic)){
  
  print(paste('Computing model', i, ':' ,'extensively adjusted model for risk of', cancers.non.syndromic[i], 'in kids with', defects.non.syndromic[i]))
  
  goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                            cancer = goback.nochrom[,cancers.non.syndromic[i]],
                            defect = goback.nochrom[,defects.non.syndromic[i]],
                            sex = factor(goback.nochrom$sex,
                                         levels = c(1,2),
                                         labels = c('Male','Female')),
                            m.age = goback.nochrom$m.age,
                            state = goback.nochrom$state.num,
                            plurality = goback.nochrom$twin,
                            race = goback.nochrom$m.race.collapsed,
                            birth.wt = goback.nochrom$birth.wt)
  
  cox <- coxph(Surv(time, cancer) ~ defect + sex + m.age + state + plurality + race + birth.wt, data = goback.surv)
  cox.coef <- summary(cox)$coefficients
  
  estimates <- data.frame(defect = defects.non.syndromic[i], 
                          cancer = cancers.non.syndromic[i], 
                          defect.hazard.ratio = exp(cox.coef[1,1]), 
                          defect.ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                          defect.ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                          plu.hazard.ratio = cox.coef[5,2],
                          plu.ci.lower = exp(cox.coef[5,1]-(1.96*cox.coef[5,3])),
                          plu.ci.upper = exp(cox.coef[5,1]+(1.96*cox.coef[5,3])),
                          birthweight.p.value = cox.coef[8,5],
                          race.nhb.hazard.ratio = cox.coef[6,2],
                          race.nhb.ci.lower = exp(cox.coef[6,1]-(1.96*cox.coef[6,3])),
                          race.nhb.ci.upper = exp(cox.coef[6,1]+(1.96*cox.coef[6,3])),
                          race.other.hazard.ratio = cox.coef[7,2],
                          race.other.ci.lower = exp(cox.coef[7,1]-(1.96*cox.coef[7,3])),
                          race.other.ci.upper = exp(cox.coef[7,1]+(1.96*cox.coef[7,3])),
                          n.comorbid = as.numeric(table(goback.nochrom[,defects.non.syndromic[i]], goback.nochrom[,cancers.non.syndromic[i]])[2,2]))
  
  adjusted.models.nonsyndromic <- rbind(adjusted.models.nonsyndromic, estimates)
  
}

rm(cox.coef, estimates, goback.surv, models, cox, i)



# Left join and write to file ---------------------------------------------

require(dplyr); require(xlsx)

baseline.models.nonsyndromic <- rename(baseline.models.nonsyndromic, 
                                       defect.hazard.ratio.baseline = defect.hazard.ratio, defect.ci.lower.baseline = defect.ci.lower, defect.ci.upper.baseline = defect.ci.upper)
adjusted.models.nonsyndromic <- rename(adjusted.models.nonsyndromic, 
                                       defect.hazard.ratio.adjusted = defect.hazard.ratio, defect.ci.lower.adjusted = defect.ci.lower, defect.ci.upper.adjusted = defect.ci.upper,
                                       n.comorbid.adjusted = n.comorbid)

models <- left_join(baseline.models.nonsyndromic, adjusted.models.nonsyndromic, by = c('defect','cancer'))
models$hazard.ratio.diff <- abs((models$defect.hazard.ratio.adjusted-models$defect.hazard.ratio.baseline)/models$defect.hazard.ratio.baseline)
models$hazard.ratio.changed.ten.percent <- ifelse(models$hazard.ratio.diff >= 0.1, 1, 0)

models <- models[, c(1:5,7:9,22,10:19,6,20)]

write.xlsx(models, file = 'W:/Old_genepi2/Jeremy/GOBACK/R outputs/goback.cox.ph.models.race.ethnicity.plurality.sensitivity.analysis.xlsx', row.names = FALSE)

