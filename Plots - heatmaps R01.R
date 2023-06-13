#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' The minimal production version of the code for generating GOBACK 
#' heatmaps.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(dplyr); require(ggplot2)

load('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/goback.coxph.results.v20180612.rdata')

#' Need an updated defects matrix.
#' Omit osteosarcoma from cancers, because no associations were tested.
defects <- c("conganomalies.cns", "conganomalies.eye", "conganomalies.ear.face.neck", "conganomalies.heart.circsys",    
             "conganomalies.respsys", "oral.clefts", "conganomalies.digestivesystem", "conganomalies.genitalandurinary",
             "conganomalies.musculoskelsys", "conganomalies.integument", "any.chromosomal.anomaly")
cancers <- c('all','aml','hl','nhl','cns.any','astro','medullo','retino','neuro','nephro','hepato','bone.any','rms.any','soft.other','gct.any')

bd.cc.associations <- data.frame(defect = rep(defects, each = 15), cancer = rep(cancers, 11))
bd.cc.associations <- left_join(bd.cc.associations, goback.coxmodels, by = c('defect','cancer'))

#' Convert axis variables to factor form for prettier display.
bd.cc.associations$cancer <- factor(bd.cc.associations$cancer, 
                                    levels = c('all','aml','hl','nhl','cns.any','astro','medullo', 'retino',
                                               'neuro','nephro','hepato','bone.any','rms.any','soft.other','gct.any'),
                                    labels = c('ALL', 'AML', "HL", "NHL", 'CNS', 'ASTRO','MB', 'RB',
                                               'NB',"NB",'HB','Bone','RMS','STS',
                                               'GCT'))

bd.cc.associations$defect <- factor(bd.cc.associations$defect, 
                                    levels = c( "conganomalies.cns",               "conganomalies.eye",               "conganomalies.ear.face.neck",     "conganomalies.heart.circsys",    
                                                "conganomalies.respsys",           "oral.clefts",                     "conganomalies.digestivesystem",   "conganomalies.genitalandurinary",
                                                "conganomalies.musculoskelsys",    "conganomalies.integument",        "any.chromosomal.anomaly"),
                                    labels = c('CNS', 'Eye', 'CRF', 'CHD',
                                               'RESP', 'Clefts', 'GI', 'GU',
                                               'MSK', 'Skin', 'CHR'))

#' Bins associations into the currently agreed upon categories based on their HR.
bd.cc.associations$signif.cat <- 9999
bd.cc.associations$signif.cat <- factor(
  ifelse(bd.cc.associations$signif.cat == 9999 & is.na(bd.cc.associations$hr), 6,
         ifelse(bd.cc.associations$ci.lower < 1 & bd.cc.associations$ci.upper > 1, 2,
                ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$hr > 1 & bd.cc.associations$hr <= 5, 3,
                       ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$hr > 5 & bd.cc.associations$hr <= 10, 4,
                              ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$hr > 10, 5,
                                     ifelse(bd.cc.associations$ci.upper < 1 & bd.cc.associations$hr < 1, 1, NA)))))),
  levels = c(1:6),
  labels = c('< 1','Null','1.01-5.00', '5.01-10.00','>10.00','Not Tested'))

p <-  ggplot(data = bd.cc.associations, aes(y = cancer, x = defect)) +
  
  geom_tile(aes(fill = signif.cat), color = 'black', width = 0.5, height = 0.5) + 

  scale_fill_manual(values = c('grey75','indianred1','firebrick3','darkred','white')) +
  
  guides(fill = guide_legend(title='Adjusted Hazard Ratio')) +
  
  theme_bw() +
  
  labs(x = 'ANOMALY', y = 'CANCER') +
  
  theme(axis.text.y = element_text(hjust = 1, size = 11, face = 'bold'),
        axis.text.x = element_text(size = 11, face = 'bold', angle = 90, hjust = 1, vjust = 0.5),
        #axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 13, face = 'bold'),
        legend.title = element_blank(),
        legend.text = element_text(size = 11, face = 'bold'),
        legend.direction = 'horizontal',
        legend.spacing.x = unit(0.1, 'cm'),
        legend.position = 'top',
        panel.grid = element_blank(),
        panel.border = element_blank())

print(p)
