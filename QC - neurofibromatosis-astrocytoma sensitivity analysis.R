#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.10.26.
#' 
#' This is me trying to figure out why the sensitivity analysis for the 
#' neurofibromatosis-astrocytoma association in NC and TX children aged 
#' >1 year is not performing the way it should.
#' 
#' Error messages are produced and the beta value is meaningless.
#' 
#' I have no idea why.
#' 
#' This code is not production quality and may contain errors or have 
#' elements out of sequence.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.chrom.v20180829.rdata')

#' An older version of the dataset, but it was the one the Cox models were generated from.
load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/goback.chrom.v20180611.rdata')

require(dplyr); require(survival)

#' Reproduce HR and CI without excluding cases < 1 year.
surv.data <- data.frame(studyid = goback.chrom$studyid, 
                        state = goback.chrom$state.num, 
                        sex = factor(goback.chrom$sex, 
                                     levels = c(1,2),
                                     labels = c('Male','Female')),
                        m.age = goback.chrom$m.age,
                        time = goback.chrom$person.yrs,
                        defect = goback.chrom$nf,
                        cancer = goback.chrom$astro)

cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = surv.data)
cox.coef <- summary(cox)$coefficients

#' Exclude infant cancers.
goback.chrom <- filter(goback.chrom, state %in% c('TX','NC'))
goback.chrom$state.num <- ifelse(goback.chrom$state.num == 3, 2, goback.chrom$state.num)
goback.chrom$exclude <- ifelse(goback.chrom$person.yrs <= 1, 1, 0)
goback.chrom <- filter(goback.chrom, exclude == 0)


surv.data <- data.frame(studyid = goback.chrom$studyid, 
                        time = goback.chrom$person.yrs,
                        defect = goback.chrom$nf,
                        cancer = goback.chrom$astro)
surv.data <- surv.data[!is.na(surv.data$defect), ]

tmp <- filter(goback.chrom, nf == 1 & astro == 1)
tmp <- filter(goback.chrom, nf == 1 & astro == 0)

cox <- coxph(Surv(time, cancer) ~ defect, data = surv.data)
cox <- coxph(Surv(person.yrs, astro) ~ nf, data = goback.chrom)
cox <- coxph(Surv(time, cancer) ~ defect + m.age + sex + state, data = surv.data)

model <- glm(cancer ~ )

cox.coef <- summary(cox)$coefficients
cox.coef

require(survminer)
fit <- survfit(Surv(time, cancer) ~ defect, data = surv.data)
ggsurvplot(fit, 
           risk.table = TRUE,
           ylim = c(0.9,1))

aggregate(time ~ defect + cancer, data = surv.data, summary)
