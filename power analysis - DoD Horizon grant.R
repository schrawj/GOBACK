#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.09.10.
#' 
#' Some code that will generate synthetic datasets that resemble GOBACK, 
#' and perform power analysis via simulation based on user-supplied 
#' parameters.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

# Determine what parameters to supply -------------------------------------

#' I need estimates of the birth prevalence of CHD variables among cases and controls 
#' and the distribution of maternal age to proceed.
require(gmodels); require(ggplot2)

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.v20180829.rdata')

CrossTable(goback$neuro)
CrossTable(goback$rvot.defects, goback$neuro, chisq = FALSE, prop.r = FALSE)
CrossTable(goback$lvot.defects, goback$neuro, chisq = FALSE, prop.r = FALSE)
CrossTable(goback$pulmvalveatresiaandstenosis, goback$neuro, chisq = FALSE, prop.r = FALSE)

print(ggplot(data = goback) + geom_density(aes(x = m.age)) + facet_wrap(~neuro))

aggregate(m.age ~ neuro, data = goback, mean, na.rm = TRUE)
aggregate(m.age ~ neuro, data = goback, sd, na.rm = TRUE)

rm(goback); gc()



# Power analysis ----------------------------------------------------------

require(tictoc); require(speedglm)

#' Supply parameters.
repetitions <- 10
significant <- matrix(nrow=repetitions, ncol=1)
n.cases <- 1000
n.controls <- 5000000
case.exposure.prob <- 0.007 # Consistently, exposure fraction among cases is 0.007.
case.maternal.age <- 27.5
control.exposure.prob <- 0.001 # Consistently, exposure fraction among controls is 0.001.
control.maternal.age <- 26.75
maternal.age.sd <- 6

#' For this simulation, I don't anticipate varying conditions for controls much.
#' It would be computationally wasteful to generate a new control dataset for 
#' every simulation.
tic()
controls <- data.frame(defect = rbinom(n.controls, 1, control.exposure.prob),
                       cancer = 0,
                       sex = rbinom(n.controls, 1, 0.5),
                       m.age = rnorm(n.controls, control.maternal.age, maternal.age.sd))

#' Run simulations.
for (i in 1:repetitions){
  
  cases <- data.frame(defect = rbinom(n.cases, 1, case.exposure.prob),
                      cancer = 1,
                      sex = rbinom(n.cases, 1, 0.5),
                      m.age = rnorm(n.cases, case.maternal.age, maternal.age.sd))
  
  sim.dat <- rbind(cases, controls)

  fit <- speedglm(cancer ~ defect + m.age, data = sim.dat, family = binomial(link = 'logit'))
  #fit <- summary(fit)$coefficients


  significant[i, ] <- fit
  
}

(table(significant[,1] <= 0.05)[2]/repetitions)*100
toc()

# Scratch paper -----------------------------------------------------------


repetitions <- 5000
significant <- matrix(nrow=repetitions, ncol=1)
fold.change <- log(1.25)
n <- 50
prop.mrd.neg.nci.hr <- 0.33
prop.mrd.pos.nci.hr <- 0.67
prop.mrd.neg.unfav.cyto <- 0.1
prop.mrd.pos.unfav.cyto <- 0.25

for (i in 1:repetitions){
  
  mrd.neg <- data.frame(sim.metab = rlnorm(n, meanlog = 0, sdlog= 0.5),
                        sim.nci = rbinom(n, 1, 0.33), 
                        sim.cyto = rbinom(n, 1, 0.1),
                        sim.immuno = rbinom(n, 1, 0.08),
                        mrd = 0)
  mrd.pos <- data.frame(sim.metab = rlnorm(n, meanlog = fold.change, sdlog = 0.75),
                        sim.nci = rbinom(n, 1, 0.66), 
                        sim.cyto = rbinom(n, 1, 0.25),
                        sim.immuno = rbinom(n, 1, 0.12),
                        mrd = 1)
  sim.dat <- rbind(mrd.neg, mrd.pos)
  
  fit <- glm(mrd ~ sim.metab + sim.nci + sim.cyto + sim.immuno, data = sim.dat, family = binomial (link = 'logit'))
  
  significant[i, ] <- summary(fit)$coefficients[1,4]
  
}

(table(significant[,1] <= 0.05)[2]/5000)*100


model <- function(x) {glm(cancer ~ defect + m.age, data = x, family = binomial(link = 'logit'))}
lapply(sim.dat, model)

predictors <- c(1,3,4)
out <- lapply(predictors, function(x) glm(sim.dat[, 2] ~ sim.dat[, x]))

detectCores()

lapply(sim.dat, function(x) glm(cancer ~ defect + m.age, data = x, family = binomial(link = 'logit')))
