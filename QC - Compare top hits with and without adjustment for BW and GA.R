#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.10.31.
#' 
#' Mixed messages from the reviewers about adjusting for Bw and GA.
#' 
#' One recommended, one did not. Co-authors largely support keeping these
#' in. Beth suggested doing both. Not sure I like that idea, but easy 
#' enough to write the code.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(xlsx); require(dplyr); require(survival)

#' Load in list of top hits.
top.hits <- read.csv(file = 'W:/Old_genepi2/Jeremy/GOBACK/R outputs/goback.cox.ph.top.hits.v20180612.csv',
                     stringsAsFactors = FALSE)

defects <- c(top.hits$defect)
cancers <- c(top.hits$cancer)

#' Some cleaning required. Exclude DiGeorge syndrome.
#' Separate remaining exposures into syndromic and non-syndromic.
index <- which(defects == 'di.george.syndrome')

defects <- defects[-index]
cancers <- cancers[-index]

syndromes <- c('down.syndrome','nf','trisomy18')

syndromic.defects <- defects[which(defects %in% syndromes)]
syndromic.cancers <- cancers[which(defects %in% syndromes)]

nonsyndromic.defects <- defects[which(!(defects %in% syndromes))]
nonsyndromic.cancers <- cancers[which(!(defects %in% syndromes))]

rm(defects, cancers, index, syndromes)



# Models for syndromic defects --------------------------------------------

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.chrom.v20180829.rdata')

adjusted.estimates <- data.frame(
                         defect = as.character(), 
                         cancer = as.character(), 
                         adjusted.hazard.ratio = as.numeric(), 
                         adjusted.ci.lower = as.numeric(), 
                         adjusted.ci.upper = as.numeric(),
                         adjusted.p.value.coef = as.numeric(),
                         num.comorbid = as.numeric())

for (i in 1:length(syndromic.defects)){
  
  if (syndromic.cancers[i] %in% c('all','nephro','hepato')){
  
    print(paste('Computing adjusted model', i, ':' ,'BW-adjusted model for risk of', syndromic.cancers[i], 'in kids with', syndromic.defects[i]))
    
    goback.surv <- data.frame(time = goback.chrom$person.yrs,
                              cancer = goback.chrom[,syndromic.cancers[i]],
                              defect = goback.chrom[,syndromic.defects[i]],
                              sex = factor(goback.chrom$sex,
                                           levels = c(1,2),
                                           labels = c('Male','Female')),
                              m.age = goback.chrom$m.age,
                              state = goback.chrom$state.num,
                              birth.wt = goback.chrom$birth.wt)
                              #,gest.age = goback.chrom$gest.age)
    
    cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state + birth.wt, data = goback.surv)
    cox.coef <- summary(cox)$coefficients
    
    new.estimates <- data.frame(
                                defect = syndromic.defects[i], 
                                cancer = syndromic.cancers[i], 
                                adjusted.hazard.ratio = exp(cox.coef[1,1]), 
                                adjusted.ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                                adjusted.ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                                adjusted.p.value.coef = cox.coef[1,5],
                                num.comorbid = as.numeric(table(goback.chrom[,syndromic.defects[i]], goback.chrom[,syndromic.cancers[i]])[2,2]))
    
    adjusted.estimates <- rbind(adjusted.estimates, new.estimates)
  
  }
  
  else{
    
    next
  }
  
}

unadjusted.estimates <- data.frame(
                                    defect = as.character(), 
                                    cancer = as.character(), 
                                    unadjusted.hazard.ratio = as.numeric(), 
                                    unadjusted.ci.lower = as.numeric(), 
                                    unadjusted.ci.upper = as.numeric(),
                                    unadjusted.p.value.coef = as.numeric())

for (i in 1:length(syndromic.defects)){
  
  if (syndromic.cancers[i] %in% c('all','nephro','hepato')){
  
    print(paste('Computing model', i, ':' ,'unadjusted model for risk of', syndromic.cancers[i], 'in kids with', syndromic.defects[i]))
    
    goback.surv <- data.frame(time = goback.chrom$person.yrs,
                              cancer = goback.chrom[,syndromic.cancers[i]],
                              defect = goback.chrom[,syndromic.defects[i]],
                              sex = factor(goback.chrom$sex,
                                           levels = c(1,2),
                                           labels = c('Male','Female')),
                              m.age = goback.chrom$m.age,
                              state = goback.chrom$state)
    
    cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
    cox.coef <- summary(cox)$coefficients
    
    new.estimates <- data.frame(
                                defect = syndromic.defects[i], 
                                cancer = syndromic.cancers[i], 
                                unadjusted.hazard.ratio = exp(cox.coef[1,1]), 
                                unadjusted.ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
                                unadjusted.ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
                                unadjusted.p.value.coef = cox.coef[1,5])
  }
  
  else{
    
    next
  }
                            
  unadjusted.estimates <- rbind(unadjusted.estimates, new.estimates)
  
}

chrom.estimates <- left_join(adjusted.estimates, unadjusted.estimates, by = c('defect','cancer'))

rm(adjusted.estimates, cox.coef, goback.chrom, goback.surv,new.estimates, unadjusted.estimates, cox, i, syndromic.cancers, syndromic.defects)



# Models for non-syndromic defects ----------------------------------------

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.nochrom.v20180829.rdata')

adjusted.estimates <- data.frame(
                                  defect = as.character(), 
                                  cancer = as.character(), 
                                  adjusted.hazard.ratio = as.numeric(), 
                                  adjusted.ci.lower = as.numeric(), 
                                  adjusted.ci.upper = as.numeric(),
                                  adjusted.p.value.coef = as.numeric(),
                                  num.comorbid = as.numeric())

#' Adjusted models.
#' Structure is a little more complicated.
#' If defect is PDA, ASD, or VSD, adjust for BW and GA.
#' Otherwise if cancer is ALL, hepato, or nephro adjust for BW.
#' Otherwise do not run model.
for (i in 1:length(nonsyndromic.defects)){
  
  if (nonsyndromic.defects[i] %in% c('atrialseptaldefect','patentductusarteriosis','ventricularseptaldefect')){
    
    print(paste('Computing adjusted model', i, ':' ,'BW- and GA-adjusted model for risk of', nonsyndromic.cancers[i], 'in kids with', nonsyndromic.defects[i]))
    
    goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                              cancer = goback.nochrom[,nonsyndromic.cancers[i]],
                              defect = goback.nochrom[,nonsyndromic.defects[i]],
                              sex = factor(goback.nochrom$sex,
                                           levels = c(1,2),
                                           labels = c('Male','Female')),
                              m.age = goback.nochrom$m.age,
                              state = goback.nochrom$state.num,
                              birth.wt = goback.nochrom$birth.wt,
                              gest.age = goback.nochrom$gest.age)
    
    cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state + birth.wt + gest.age, data = goback.surv)
    cox.coef <- summary(cox)$coefficients
    
    new.estimates <- data.frame(
      defect = nonsyndromic.defects[i], 
      cancer = nonsyndromic.cancers[i], 
      adjusted.hazard.ratio = exp(cox.coef[1,1]), 
      adjusted.ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
      adjusted.ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
      adjusted.p.value.coef = cox.coef[1,5],
      num.comorbid = as.numeric(table(goback.nochrom[,nonsyndromic.defects[i]], goback.nochrom[,nonsyndromic.cancers[i]])[2,2]))
    
    adjusted.estimates <- rbind(adjusted.estimates, new.estimates)
    
  }
  
  else if (nonsyndromic.cancers[i] %in% c('all','nephro','hepato')){
    
    print(paste('Computing adjusted model', i, ':' ,'BW-adjusted model for risk of', nonsyndromic.cancers[i], 'in kids with', nonsyndromic.defects[i]))
    
    goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                              cancer = goback.nochrom[,nonsyndromic.cancers[i]],
                              defect = goback.nochrom[,nonsyndromic.defects[i]],
                              sex = factor(goback.nochrom$sex,
                                           levels = c(1,2),
                                           labels = c('Male','Female')),
                              m.age = goback.nochrom$m.age,
                              state = goback.nochrom$state.num,
                              birth.wt = goback.nochrom$birth.wt)
    
    cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state + birth.wt, data = goback.surv)
    cox.coef <- summary(cox)$coefficients
    
    new.estimates <- data.frame(
      defect = nonsyndromic.defects[i], 
      cancer = nonsyndromic.cancers[i], 
      adjusted.hazard.ratio = exp(cox.coef[1,1]), 
      adjusted.ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
      adjusted.ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
      adjusted.p.value.coef = cox.coef[1,5],
      num.comorbid = as.numeric(table(goback.nochrom[,nonsyndromic.defects[i]], goback.nochrom[,nonsyndromic.cancers[i]])[2,2]))
    
    adjusted.estimates <- rbind(adjusted.estimates, new.estimates)
    
  }
  
  else{
    
    next
    
  }
  
}

unadjusted.estimates <- data.frame(
                                    defect = as.character(), 
                                    cancer = as.character(), 
                                    unadjusted.hazard.ratio = as.numeric(), 
                                    unadjusted.ci.lower = as.numeric(), 
                                    unadjusted.ci.upper = as.numeric(),
                                    unadjusted.p.value.coef = as.numeric())

for (i in 1:length(nonsyndromic.defects)){
  
  if (nonsyndromic.defects[i] %in% c('atrialseptaldefect','patentductusarteriosis','ventricularseptaldefect')){
    
    print(paste('Computing unadjusted model', i, ':' ,'for risk of', nonsyndromic.cancers[i], 'in kids with', nonsyndromic.defects[i]))
    
    goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                              cancer = goback.nochrom[,nonsyndromic.cancers[i]],
                              defect = goback.nochrom[,nonsyndromic.defects[i]],
                              sex = factor(goback.nochrom$sex,
                                           levels = c(1,2),
                                           labels = c('Male','Female')),
                              m.age = goback.nochrom$m.age,
                              state = goback.nochrom$state.num)
    
    cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
    cox.coef <- summary(cox)$coefficients
    
    new.estimates <- data.frame(
      defect = nonsyndromic.defects[i], 
      cancer = nonsyndromic.cancers[i], 
      unadjusted.hazard.ratio = exp(cox.coef[1,1]), 
      unadjusted.ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
      unadjusted.ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
      unadjusted.p.value.coef = cox.coef[1,5])
    
    unadjusted.estimates <- rbind(unadjusted.estimates, new.estimates)
    
  }
  
  else if (nonsyndromic.cancers[i] %in% c('all','nephro','hepato')){
    
    print(paste('Computing unadjusted model', i, ':' ,'for risk of', nonsyndromic.cancers[i], 'in kids with', nonsyndromic.defects[i]))
    
    goback.surv <- data.frame(time = goback.nochrom$person.yrs,
                              cancer = goback.nochrom[,nonsyndromic.cancers[i]],
                              defect = goback.nochrom[,nonsyndromic.defects[i]],
                              sex = factor(goback.nochrom$sex,
                                           levels = c(1,2),
                                           labels = c('Male','Female')),
                              m.age = goback.nochrom$m.age,
                              state = goback.nochrom$state.num)
    
    cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = goback.surv)
    cox.coef <- summary(cox)$coefficients
    
    new.estimates <- data.frame(
      defect = nonsyndromic.defects[i], 
      cancer = nonsyndromic.cancers[i], 
      unadjusted.hazard.ratio = exp(cox.coef[1,1]), 
      unadjusted.ci.lower = exp(cox.coef[1,1]-(1.96*cox.coef[1,3])), 
      unadjusted.ci.upper = exp(cox.coef[1,1]+(1.96*cox.coef[1,3])),
      unadjusted.p.value.coef = cox.coef[1,5])
    
    unadjusted.estimates <- rbind(unadjusted.estimates, new.estimates)
    
  }
  
  else{
    
    next
    
  }
  
}

no.chrom.estimates <- left_join(adjusted.estimates, unadjusted.estimates, by = c('defect','cancer'))

estimates <- rbind(chrom.estimates, no.chrom.estimates)
estimates <- estimates[,c(1:6,8:11,7)]

write.csv(estimates, file = 'W:/Old_genepi2/Jeremy/GOBACK/R outputs/goback.cox.top.hits.with.and.without.bw.ga.adjustment.v20181031.csv',
          row.names = FALSE)

rm(list = ls()); gc()
