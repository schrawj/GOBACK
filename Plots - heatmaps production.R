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
                                    labels = c('ALL', 'AML', "Hodgkin Lymphoma", "Non-Hodgkin Lymphoma", 'Any CNS Tumor', 'Astrocytoma','Medulloblastoma', 'Retinoblastoma',
                                               'Neuroblastoma',"Wilms Tumor",'Hepatoblastoma','Any Bone Tumor','Rhabdomyosarcoma','Non-RMS Soft Tissue Sarcoma',
                                               'Any Germ Cell Tumor'))

bd.cc.associations$defect <- factor(bd.cc.associations$defect, 
                                    levels = c( "conganomalies.cns",               "conganomalies.eye",               "conganomalies.ear.face.neck",     "conganomalies.heart.circsys",    
                                                "conganomalies.respsys",           "oral.clefts",                     "conganomalies.digestivesystem",   "conganomalies.genitalandurinary",
                                                "conganomalies.musculoskelsys",    "conganomalies.integument",        "any.chromosomal.anomaly"),
                                    labels = c('Any CNS Defect', 'Any Eye Defect', 'Any Ear, Face or Neck Defect', 'Any Heart or Circulatory System Defect',
                                               'Any Respiratory System Defect', 'Oral Clefts', 'Any Digestive System Defect', 'Any Genitourinary Defect',
                                               'Any Musculoskeletal Defect', 'Any Integument Defect', 'Any Chromosomal Anomaly'))

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

#' Intermediate number of gradations, all red.
#' Do not plot title, or force <1 category to display in legend.
#' Font sizes have been puffed up for large-format PPT prentations.
p <-  ggplot(data = bd.cc.associations, aes(x = cancer, y = defect)) +
             geom_tile(aes(fill = signif.cat), color = 'black') + 
             scale_fill_manual(values = c('grey75','indianred1','firebrick3','darkred','white')) +
             guides(fill = guide_legend(title='Hazard Ratio')) +
             labs(x = "Cancer", y = "Anomaly") +
             theme_bw() +
             theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14)) +
             theme(axis.title.x = element_text(face = 'bold', size = 22)) +
             theme(axis.text.y = element_text(size = 14)) +
             theme(axis.title.y = element_text(face = "bold", size = 22)) +
             theme(legend.title = element_text(face = 'bold', size = 22)) +
             theme(legend.text = element_text(size = 14))
print(p)

p <-  ggplot(data = bd.cc.associations, aes(x = cancer, y = defect)) +
      
      geom_tile(aes(fill = signif.cat), color = 'black') + 
  
      scale_fill_manual(values = c('grey75','indianred1','firebrick3','darkred','white')) +
  
      guides(fill = guide_legend(title='Adjusted Hazard Ratio')) +
  
      theme_bw() +
  
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = 'bold'),
            axis.text.y = element_text(size = 14, face = 'bold'),
            axis.title = element_blank(),
            legend.title = element_text(face = 'bold', size = 22),
            legend.text = element_text(size = 14, face = 'bold'),
            legend.position = 'top')

print(p)