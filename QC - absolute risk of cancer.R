#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.10.23.
#' 
#' NEJM reviewers commented that they would like to know the absolute risk
#' of cancer in children with birth defects.
#' 
#' For the associations in table 3, we will provide:
#' - The proportion of children with that cancer who have that birth defect
#' - The proportion of children with that birth defect who go on to develop
#'   that cancer.
#'   
#' Also report absolute risk of the cancers in Figure 1 according to 
#' number of major birth defects.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(xlsx); require(tictoc)

#' Load in list of top hits.
top.hits <- read.csv(file = 'W:/Old_genepi2/Jeremy/GOBACK/R outputs/goback.cox.ph.top.hits.v20180612.csv',
                      stringsAsFactors = FALSE)

defects <- c(top.hits$defect)
cancers <- c(top.hits$cancer)

#' Some cleaning required. Exclude DiGeorge syndrome.
#' Separate remaining exposures into syndromic and non-syndromic.
index <- which(defects == 'di.george.syndrome')

defects <- defects[c(1:12,14:43)]
cancers <- cancers[c(1:12,14:43)]

syndromes <- c('down.syndrome','nf','trisomy18')

syndromic.defects <- defects[which(defects %in% syndromes)]
syndromic.cancers <- cancers[which(defects %in% syndromes)]

nonsyndromic.defects <- defects[which(!(defects %in% syndromes))]
nonsyndromic.cancers <- cancers[which(!(defects %in% syndromes))]

rm(defects, cancers, index, syndromes)

#' For syndromic birth defects.
load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.chrom.v20180829.rdata')

out <- data.frame(defect = as.character(),
                  cancer = as.numeric(),
                  n.with.defect = as.numeric(),
                  n.with.cancer = as.numeric(),
                  n.with.defect.and.cancer = as.numeric(),
                  pct.kids.w.defect.who.develop.cancer = as.numeric(),
                  pct.kids.w.cancer.who.have.defect = as.numeric())

for (i in 1:length(syndromic.defects)){
  
  tic()
  
  defect <- syndromic.defects[i]
  cancer <- syndromic.cancers[i]
  
  tab <- table(goback.chrom[,defect], goback.chrom[,cancer], useNA = 'ifany')
  
  n.with.defect <- sum(tab[2,])
  n.with.cancer <- sum(tab[,2])
  n.with.defect.and.cancer <- tab[2,2]
  pct.defect.with.cancer <- (n.with.defect.and.cancer/n.with.defect)*100
  pct.cancer.with.defect <- (n.with.defect.and.cancer/n.with.cancer)*100
  
  
  new.out <- data.frame(defect = defect,
                        cancer = cancer,
                        n.with.defect = n.with.defect,
                        n.with.cancer = n.with.cancer,
                        n.with.defect.and.cancer = n.with.defect.and.cancer,
                        pct.kids.w.defect.who.develop.cancer = pct.defect.with.cancer,
                        pct.kids.w.cancer.who.have.defect = pct.cancer.with.defect)
  
  out <- rbind(out, new.out)
  
  toc()
  
}

rm(goback.chrom); gc()

#' For non-syndromic birth defects.
load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.nochrom.v20180829.rdata')

for (i in 1:length(nonsyndromic.defects)){
  
  tic()
  
  defect <- nonsyndromic.defects[i]
  cancer <- nonsyndromic.cancers[i]
  
  tab <- table(goback.nochrom[,defect], goback.nochrom[,cancer], useNA = 'ifany')
  
  n.with.defect <- sum(tab[2,])
  n.with.cancer <- sum(tab[,2])
  n.with.defect.and.cancer <- tab[2,2]
  pct.defect.with.cancer <- (n.with.defect.and.cancer/n.with.defect)*100
  pct.cancer.with.defect <- (n.with.defect.and.cancer/n.with.cancer)*100
  
  
  new.out <- data.frame(defect = defect,
                        cancer = cancer,
                        n.with.defect = n.with.defect,
                        n.with.cancer = n.with.cancer,
                        n.with.defect.and.cancer = n.with.defect.and.cancer,
                        pct.kids.w.defect.who.develop.cancer = pct.defect.with.cancer,
                        pct.kids.w.cancer.who.have.defect = pct.cancer.with.defect)
  
  out <- rbind(out, new.out)
  
  toc()
  
}

write.xlsx(out, file = 'W:/Old_genepi2/Jeremy/GOBACK/R outputs/goback.cox.ph.top.hits.absolute.cancer.risk.xlsx', row.names = FALSE)



# According to number of birth defects ------------------------------------

require(xlsx)

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.nochrom.v20180829.rdata')

goback.nochrom$majordefect.cat <- factor(ifelse(goback.nochrom$majordefect.total == 0, 0,
                                                ifelse(goback.nochrom$majordefect.total == 1, 1,
                                                       ifelse(goback.nochrom$majordefect.total == 2, 2, 
                                                              ifelse(goback.nochrom$majordefect.total == 3, 3,
                                                                     ifelse(goback.nochrom$majordefect.total >= 4, 4, NA))))),
                                         levels = c(0:4),
                                         labels = c('0', '1', '2', '3', '4 or more'))

#' Same new variables as in the cumulative incidence plots.
goback.nochrom$any.heme.cancer <- ifelse(goback.nochrom$leu.any == 1 | goback.nochrom$lym.any == 1, 1, 0)

solid.tumors <- unique(goback.nochrom$cancer1)
solid.tumors <- subset(solid.tumors, !(solid.tumors %in% c(NA, 'all','leu.other','aml','hl','nhl','cns.other','medullo','pnet','lym.other',
                                                           'gct.intra','astro','ependymoma')))

goback.nochrom$any.non.cns.solid.tumor <- ifelse(goback.nochrom$cancer1 %in% solid.tumors, 1, 0)

#' First iteration.
i <- 'cancer'

out <- data.frame(number.of.defects = c('0','1','2','3','4 or more'))

tab <- table(goback.nochrom$majordefect.cat, goback.nochrom[,i], dnn = c('number of major defects','number of cancer cases'))

out$cancer <- rep(i,5)
out$n.in.category <- c(sum(tab[1,]), sum(tab[2,]), sum(tab[3,]), sum(tab[4,]), sum(tab[5, ]))
out$n.with.cancer <- c(tab[1,2], tab[2,2], tab[3,2], tab[4,2], tab[5,2])
out$percent.developing.cancer <- (out$n.with.cancer/out$n.in.category)*100  
  
#' Vector of remaining outcome variables to loop over.
other.outcomes <- c('any.heme.cancer','cns.any','any.non.cns.solid.tumor')  
  
for (i in other.outcomes){
  
  new.out <- data.frame(number.of.defects = c('0','1','2','3','4 or more'))
  
  tab <- table(goback.nochrom$majordefect.cat, goback.nochrom[,i], dnn = c('number of major defects','number of cancer cases'))
  
  new.out$cancer <- rep(i,5)
  new.out$n.in.category <- c(sum(tab[1,]), sum(tab[2,]), sum(tab[3,]), sum(tab[4,]), sum(tab[5, ]))
  new.out$n.with.cancer <- c(tab[1,2], tab[2,2], tab[3,2], tab[4,2], tab[5,2])
  new.out$percent.developing.cancer <- (new.out$n.with.cancer/new.out$n.in.category)*100  
  
  out <- rbind(out, new.out)
  
}

write.xlsx(out, file = 'W:/Old_genepi2/Jeremy/GOBACK/R outputs/goback.cox.ph.top.hits.absolute.cancer.risk.xlsx',sheetName = 'Number of Defects',
           row.names = FALSE, append = TRUE)
