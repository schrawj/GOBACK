#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' GOBACK heatmaps, production version.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------



# Prep environment --------------------------------------------------------
require(dplyr); require(stringr); require(ggplot2); require(ggthemes)

setwd('Z:/Jeremy/GOBACK/Datasets/') 
load('goback.cox.ph.results.v20180124.1.rdata')



# Generate heatmap --------------------------------------------------------

#' Need an updated defects matrix.
defects <- c("conganomalies.cns", "conganomalies.eye", "conganomalies.ear.face.neck", "conganomalies.heart.circsys",    
             "conganomalies.respsys", "oral.clefts", "conganomalies.digestivesystem", "conganomalies.genitalandurinary",
             "conganomalies.musculoskelsys", "conganomalies.integument", "chromosomalanomalies")
cancers <- c('all','aml','hl','nhl','cns.any','astro','medullo','neuro','retino','nephro','hepato','bone.any','osteo','rms.any','soft.other','gct.any')

bd.cc.associations <- data.frame(defect = rep(defects, each = 16), cancer = rep(cancers, 11))
bd.cc.associations <- left_join(bd.cc.associations, goback.coxmodels, by = c('defect','cancer'))
bd.cc.associations <- rename(bd.cc.associations, hr = HR, ci.lower = CI.lower, ci.upper = CI.upper)

bd.cc.associations$signif.cat <- 9999
bd.cc.associations$signif.cat <- factor(
  ifelse(bd.cc.associations$signif.cat == 9999 & is.na(bd.cc.associations$hr), 6,
         ifelse(bd.cc.associations$ci.lower < 1 & bd.cc.associations$ci.upper > 1, 2,
                ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$ci.upper > 1 & bd.cc.associations$hr > 1 & bd.cc.associations$hr <= 5, 3,
                       ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$ci.upper > 1 & bd.cc.associations$hr > 5 & bd.cc.associations$hr <= 25, 4,
                              ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$ci.upper > 1 & bd.cc.associations$hr > 25, 5,
                                     ifelse(bd.cc.associations$ci.lower < 1 & bd.cc.associations$ci.upper < 1 & bd.cc.associations$hr < 1, 1, NA)))))),
  levels = c(1:6),
  labels = c('< 1','Null','1-5', '5.01-25','>25','Not Tested'))

#' Convert defect and cancer variables to factor form for prettier display.
bd.cc.associations$cancer <- factor(bd.cc.associations$cancer, 
                                    levels = c('all','aml','hl','nhl','cns.any','astro','medullo', 'retino',
                                               'neuro','nephro','hepato','bone.any','osteo','rms.any','soft.other','gct.any'),
                                    labels = c('ALL', 'AML', "Hodgkin's Lymphoma", "Non-Hodgkin's Lymphoma", 'Any CNS Tumor', 'Astrocytoma','Medulloblastoma', 'Retinoblastoma',
                                               'Neuroblastoma',"Wilm's Tumor",'Hepatoblastoma','Any Bone Tumor','Osteosarcoma','Rhabdomyosarcoma','Non-RMS Soft Tissue Sarcoma',
                                               'Any Germ Cell Tumor'))

bd.cc.associations$defect <- factor(bd.cc.associations$defect, 
                                    levels = c( "conganomalies.cns",               "conganomalies.eye",               "conganomalies.ear.face.neck",     "conganomalies.heart.circsys",    
                                                "conganomalies.respsys",           "oral.clefts",                     "conganomalies.digestivesystem",   "conganomalies.genitalandurinary",
                                                "conganomalies.musculoskelsys",    "conganomalies.integument",        "chromosomalanomalies"),
                                    labels = c('Any CNS Anomaly', 'Any Eye Anomaly', 'Any Ear, Face or Neck Anomaly', 'Any Heart or Circulatory System Anomaly',
                                               'Any Respiratory System Anomaly', 'Cleft Lip and Cleft Palate', 'Any Digestive System Anomaly', 'Any Genitourinary Anomaly',
                                               'Any Musculoskeletal Anomaly', 'Any Integument Anomaly', 'Any Chromosomal Anomaly'))

plot<-  ggplot(data = bd.cc.associations, aes(x = cancer, y = defect)) +
  geom_tile(aes(fill = signif.cat), color = 'black') + 
  scale_fill_manual(values = c('green','grey50','indianred1','firebrick1','darkred','white'), drop = FALSE) +
  guides(fill = guide_legend(title='Hazard Ratio')) +
  ggtitle('Figure 1. Associations Between Different Categories of Birth Defects and Childhood Cancers') +
  labs(x = "Cancer", y = "Anomaly") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

print(plot)




# Omit osteosarcoma (no associations tested) ------------------------------

defects <- c("conganomalies.cns", "conganomalies.eye", "conganomalies.ear.face.neck", "conganomalies.heart.circsys",    
             "conganomalies.respsys", "oral.clefts", "conganomalies.digestivesystem", "conganomalies.genitalandurinary",
             "conganomalies.musculoskelsys", "conganomalies.integument", "chromosomalanomalies")
cancers <- c('all','aml','hl','nhl','cns.any','astro','medullo','retino','neuro','nephro','hepato','bone.any','rms.any','soft.other','gct.any')

bd.cc.associations <- data.frame(defect = rep(defects, each = 15), cancer = rep(cancers, 11))
bd.cc.associations <- left_join(bd.cc.associations, goback.coxmodels, by = c('defect','cancer'))
bd.cc.associations <- rename(bd.cc.associations, hr = HR, ci.lower = CI.lower, ci.upper = CI.upper)

bd.cc.associations$signif.cat <- 9999
bd.cc.associations$signif.cat <- factor(
  ifelse(bd.cc.associations$signif.cat == 9999 & is.na(bd.cc.associations$hr), 6,
         ifelse(bd.cc.associations$ci.lower < 1 & bd.cc.associations$ci.upper > 1, 2,
                ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$ci.upper > 1 & bd.cc.associations$hr > 1 & bd.cc.associations$hr <= 5, 3,
                       ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$ci.upper > 1 & bd.cc.associations$hr > 5 & bd.cc.associations$hr <= 25, 4,
                              ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$ci.upper > 1 & bd.cc.associations$hr > 25, 5,
                                     ifelse(bd.cc.associations$ci.lower < 1 & bd.cc.associations$ci.upper < 1 & bd.cc.associations$hr < 1, 1, NA)))))),
  levels = c(1:6),
  labels = c('< 1','Null','1-5', '5.01-25','>25','Not Tested'))

#' A significance category variable with more gradations.
bd.cc.associations$signif.cat.fine <- 9999
bd.cc.associations$signif.cat.fine <- factor(
  
                                               ifelse(bd.cc.associations$signif.cat.fine == 9999 & is.na(bd.cc.associations$hr), 9,
                                               ifelse(bd.cc.associations$ci.lower < 1 & bd.cc.associations$ci.upper > 1, 2,
                                               ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$ci.upper > 1 & bd.cc.associations$hr > 1 & bd.cc.associations$hr <= 5, 3,
                                               ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$ci.upper > 1 & bd.cc.associations$hr > 5 & bd.cc.associations$hr <= 10, 4,
                                               ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$ci.upper > 1 & bd.cc.associations$hr > 10 & bd.cc.associations$hr <= 15, 5,
                                               ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$ci.upper > 1 & bd.cc.associations$hr > 15 & bd.cc.associations$hr <= 20, 6,
                                               ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$ci.upper > 1 & bd.cc.associations$hr > 20 & bd.cc.associations$hr <= 25, 7,
                                               ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$ci.upper > 1 & bd.cc.associations$hr > 25, 8,
                                               ifelse(bd.cc.associations$ci.lower < 1 & bd.cc.associations$ci.upper < 1, 1, NA))))))))),
                                      
                                      levels = c(1:9),
                                      labels = c('< 1','Null','1.01-5', '5.01-10','10.01-15','15.01-20','20.01-25','>25','Not Tested'))

#' Intermediate number of gradations.
bd.cc.associations$signif.cat.int <- 9999
bd.cc.associations$signif.cat.int <- factor(
  
ifelse(bd.cc.associations$signif.cat.int == 9999 & is.na(bd.cc.associations$hr), 7,
ifelse(bd.cc.associations$ci.lower < 1 & bd.cc.associations$ci.upper > 1, 2,
ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$ci.upper > 1 & bd.cc.associations$hr > 1 & bd.cc.associations$hr <= 5, 3,
ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$ci.upper > 1 & bd.cc.associations$hr > 5 & bd.cc.associations$hr <= 10, 4,
ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$ci.upper > 1 & bd.cc.associations$hr > 10 & bd.cc.associations$hr <= 20, 5,
ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$ci.upper > 1 & bd.cc.associations$hr > 20, 6,
ifelse(bd.cc.associations$ci.lower < 1 & bd.cc.associations$ci.upper < 1, 1, NA))))))),
  
levels = c(1:7),
labels = c('< 1','Null','1.01-5', '5.01-10','10.01-20','>20','Not Tested'))

#' Convert axis variables to factor form for prettier display.
bd.cc.associations$cancer <- factor(bd.cc.associations$cancer, 
                                    levels = c('all','aml','hl','nhl','cns.any','astro','medullo', 'retino',
                                               'neuro','nephro','hepato','bone.any','rms.any','soft.other','gct.any'),
                                    labels = c('ALL', 'AML', "Hodgkin's Lymphoma", "Non-Hodgkin's Lymphoma", 'Any CNS Tumor', 'Astrocytoma','Medulloblastoma', 'Retinoblastoma',
                                               'Neuroblastoma',"Wilm's Tumor",'Hepatoblastoma','Any Bone Tumor','Rhabdomyosarcoma','Non-RMS Soft Tissue Sarcoma',
                                               'Any Germ Cell Tumor'))

bd.cc.associations$defect <- factor(bd.cc.associations$defect, 
                                    levels = c( "conganomalies.cns",               "conganomalies.eye",               "conganomalies.ear.face.neck",     "conganomalies.heart.circsys",    
                                                "conganomalies.respsys",           "oral.clefts",                     "conganomalies.digestivesystem",   "conganomalies.genitalandurinary",
                                                "conganomalies.musculoskelsys",    "conganomalies.integument",        "chromosomalanomalies"),
                                    labels = c('Any CNS Anomaly', 'Any Eye Anomaly', 'Any Ear, Face or Neck Anomaly', 'Any Heart or Circulatory System Anomaly',
                                               'Any Respiratory System Anomaly', 'Oral Clefts', 'Any Digestive System Anomaly', 'Any Genitourinary Anomaly',
                                               'Any Musculoskeletal Anomaly', 'Any Integument Anomaly', 'Any Chromosomal Anomaly'))

#' Coarse color gradient for significance.
plot <-  ggplot(data = bd.cc.associations, aes(x = cancer, y = defect)) +
  geom_tile(aes(fill = signif.cat), color = 'black') + 
  scale_fill_manual(values = c('green','grey50','indianred1','firebrick1','darkred','white'), drop = FALSE) +
  guides(fill = guide_legend(title='Hazard Ratio')) +
  ggtitle('Figure 1. Associations Between Different Categories of Birth Defects and Childhood Cancers') +
  labs(x = "Cancer", y = "Anomaly") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
print(plot)

#' Fine color gradient for significance.
plot2 <-  ggplot(data = bd.cc.associations, aes(x = cancer, y = defect)) +
  geom_tile(aes(fill = signif.cat.fine), color = 'black') + 
  scale_fill_manual(values = c('green','grey50','yellow','indianred1','tomato','firebrick1','red3','darkred','white'), drop = FALSE) +
  guides(fill = guide_legend(title='Hazard Ratio')) +
  ggtitle('Figure 1. Associations Between Different Categories of Birth Defects and Childhood Cancers') +
  labs(x = "Cancer", y = "Anomaly") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
print(plot2)

#' Intermediate number of gradations, all red.
plot3 <-  ggplot(data = bd.cc.associations, aes(x = cancer, y = defect)) +
  geom_tile(aes(fill = signif.cat.int), color = 'black') + 
  scale_fill_manual(values = c('green','grey75','indianred1','firebrick3','red','darkred','white'), drop = FALSE) +
  guides(fill = guide_legend(title='Hazard Ratio')) +
  ggtitle('Figure 1. Associations Between Different Categories of Birth Defects and Childhood Cancers') +
  labs(x = "Cancer", y = "Anomaly") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
print(plot3)
