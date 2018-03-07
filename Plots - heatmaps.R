#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' GOBACK heatmaps, full iterative version.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------



# Prep environment --------------------------------------------------------
require(dplyr); require(stringr); require(ggplot2); require(ggthemes)



# User-defined functions --------------------------------------------------

#' Generates heatmaps with an eye towards comprehensively displaying BD-CC associations.
#' Saves each plot to the working directory as a PDF at 1.75x magnification.
generate.heatmap.verbose <- function(vector.of.defects, text.size, plot.title, name.of.file){
  
  data <- filter(bd.cc.associations, birth.defect %in% vector.of.defects)
  
  min.or <- min(data$or, na.rm = TRUE)
  mean.or <- mean(data$or, na.rm = TRUE)
  max.or <- max(data$or, na.rm = TRUE)
  
  labels = paste(data$or.new,"\n", data$ci, sep = "")
  labels = ifelse(labels == "NA\n", NA, labels)
  
  plot <- ggplot(data = data, aes(x = birth.defect, y = cancer)) +
          geom_tile(aes(fill = OR),color = "white") +
          geom_text(aes(fill = data$or.new, label = labels),size = text.size) +
          scale_fill_gradient2(low = "green", mid = "white",high = "indianred1",midpoint = mean.or, na.value = "grey90", limits = c(min.or,max.or)) + 
          ggtitle(plot.title) +
          labs(y = "Cancer", x = "Defect Type") +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45,hjust = 1))

  print(plot)
  
  ggsave(name.of.file, scale = 1.75)
  
}

#' Generates heatmaps with an eye towards minimizing white space.
#' Saves each plot to the working directory as a PDF at 1.75x magnification.
generate.heatmap.succinct <- function(vector.of.defects, text.size, plot.title, name.of.file){
  
  data <- arrange(filter(bd.cc.associations, birth.defect %in% vector.of.defects), desc(cancer))
  data <- data[!is.na(data$or), ]
  
  min.or <- min(data$or, na.rm = TRUE)
  mean.or <- mean(data$or, na.rm = TRUE)
  max.or <- max(data$or, na.rm = TRUE)
  
  labels = paste(data$or.new,"\n", data$ci, sep = "")
  labels = ifelse(labels == "NA\n", NA, labels)
  
  plot<- ggplot(data = data, aes(x = birth.defect, y = cancer)) +
    geom_tile(aes(fill = OR),color = "black") +
    geom_text(aes(fill = data$or.new, label = labels),size = text.size) +
    scale_fill_gradient2(low = "green", mid = "white",high = "indianred1",midpoint = mean.or, na.value = "grey90", limits = c(min.or,max.or)) + 
    ggtitle(plot.title) +
    labs(y = "Cancer", x = "Defect Type") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45,hjust = 1))
  
  print(plot)
  
  ggsave(filename = name.of.file, scale = 1.75)
  
}

#' Succinct heatmaps with new color scheme.
generate.heatmap.new <- function(vector.of.defects, plot.title, name.of.file){
  
  data <- arrange(filter(bd.cc.associations, birth.defect %in% vector.of.defects), desc(cancer))
  data <- data[!is.na(data$or), ]
  
  min.or <- min(data$or, na.rm = TRUE)
  mean.or <- mean(data$or, na.rm = TRUE)
  max.or <- max(data$or, na.rm = TRUE)
  
  labels = paste(data$or.new,"\n", data$ci, sep = "")
  labels = ifelse(labels == "NA\n", NA, labels)
  
  plot <- ggplot(data = data, aes(x = birth.defect, y = cancer, label = or.new)) +
          geom_tile(aes(fill = signif.cat), color = 'black') + 
          scale_fill_manual(values = c('grey50','yellow','indianred1','firebrick')) +
          guides(fill = guide_legend(title='Odds Ratio')) +
          ggtitle(plot.title) +
          labs(y = "Cancer", x = "Defect Type") +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45,hjust = 1))
  
  print(plot)
  
  ggsave(filename = name.of.file, scale = 1.75)
  
}

# Generate BD-CC matrix to house regression results -----------------------
load("//discovery2/bcm-pedi-genepi2/Jeremy/GOBACK/Datasets/Combined Arkansas Michigan and Texas/ar.mi.tx.nochromosomaldefects.v20171025.1.rdata")

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' I will want a data frame with a row for every possible BD-CC pair to
#' house the data we import from our regression models.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

defects <- c(colnames(armitx.nochrom[22:111]))
cancers <- c(colnames(armitx.nochrom[118:157]))

bd.cc.associations <- data.frame(birth.defect = rep(defects, each = 40), cancer = rep(cancers, 90))

print(bd.cc.associations[1:200,])

#' Check only unique combinations incorporated.
bd.cc.associations$pair <- paste0(bd.cc.associations$birth.defect,'-',bd.cc.associations$cancer)
bd.cc.associations <- bd.cc.associations[!duplicated(bd.cc.associations$pair), ]
bd.cc.associations <- bd.cc.associations[,-3]

save(bd.cc.associations, file = '//discovery2/bcm-pedi-genepi2/Jeremy/GOBACK/Datasets/list.of.crude.bd.cc.associations.rdata')

# Heatmaps: generate input data -------------------------------------------

#' A data frame with a row for every possible bd-cc pair.
load("Z:/Jeremy/GOBACK/Datasets/list.of.crude.bd.cc.associations.rdata")

model.outputs <- read.csv(file = 'Z:/Jeremy/GOBACK/R outputs/BD-CC associations.csv', header = TRUE, sep = ',', stringsAsFactors = FALSE)
model.outputs <- arrange(model.outputs, birth.defect)

bd.cc.associations$birth.defect <- as.character(bd.cc.associations$birth.defect)
bd.cc.associations$cancer <- as.character(bd.cc.associations$cancer)
bd.cc.associations <- arrange(bd.cc.associations, birth.defect)
bd.cc.associations <- left_join(bd.cc.associations, model.outputs, by = c('birth.defect','cancer'))

rm(model.outputs)

save(bd.cc.associations, file = 'Z:/Jeremy/GOBACK/Datasets/crude.bd.cc.associations.v20171031.1.rdata')

#' Can't figure out how to tell R to round a number to a specific number of decimal places.
#' A workaround.  Convert CI variables to string, extract all digits before and the first two digits after the decimal.
#' This is equivalent to always rounding down at the second digit, but I'm at peace with that.
bd.cc.associations$or.new <- as.numeric(ifelse(is.na(bd.cc.associations$or), NA, 
                                               ifelse(!is.na(bd.cc.associations$or), str_extract(as.character(bd.cc.associations$or), '^\\d+[.]\\d{2}'), 'this is an error')))
bd.cc.associations$ci.lower.new <- ifelse(is.na(bd.cc.associations$or), NA, 
                                          ifelse(!is.na(bd.cc.associations$or), str_extract(as.character(bd.cc.associations$ci.lower), '^\\d+[.]\\d{2}'), 'this is an error'))
bd.cc.associations$ci.upper.new <- ifelse(is.na(bd.cc.associations$or), NA, 
                                          ifelse(!is.na(bd.cc.associations$or), str_extract(as.character(bd.cc.associations$ci.upper), '^\\d+[.]\\d{2}'), 'this is an error'))

bd.cc.associations$ci <- ifelse(is.na(bd.cc.associations$or), '',
                                ifelse(!is.na(bd.cc.associations$or), paste0(bd.cc.associations$ci.lower.new,'-',bd.cc.associations$ci.upper.new), '9999-9999'))

bd.cc.associations$OR <- as.numeric(ifelse(bd.cc.associations$ci.lower < 1 & bd.cc.associations$ci.upper < 1, bd.cc.associations$or.new,
                                           ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$ci.upper > 1, bd.cc.associations$or.new, ''))) 
bd.cc.associations$ci.new.signif <- ifelse(bd.cc.associations$ci.lower < 1 & bd.cc.associations$ci.upper < 1, bd.cc.associations$ci,
                                           ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$ci.upper > 1, bd.cc.associations$ci, '')) 
bd.cc.associations$ci.new.signif <- ifelse(is.na(bd.cc.associations$ci.new.signif), '', bd.cc.associations$ci.new.signif)

save(bd.cc.associations, file = 'Z:/Jeremy/GOBACK/Datasets/crude.bd.cc.associations.v20171031.2.rdata')



# Heatmaps: generate succinct plots ---------------------------------------
load("Z:/Jeremy/GOBACK/Datasets/crude.bd.cc.associations.v20171031.2.rdata")
load("Z:/Jeremy/GOBACK/Datasets/Michigan/mi.v20171006.4.rdata")

setwd('Z:/Jeremy/GOBACK/R outputs/Heatmaps/Succinct/')

#' Create labels
labels = paste(bd.cc.associations$or.new,"\n", bd.cc.associations$ci, sep = "")
labels = ifelse(labels == "NA\n", NA, labels)

#' Generate plots
generate.heatmap.succinct(names(mi[104:111]), 2.5, 'Chromosomal Anomalies and Childhood Cancer','chromosomal.anomalies.pdf')
generate.heatmap.succinct(names(mi[22:30]), 2.5, 'Congenital Anomalies of the CNS and Childhood Cancer', 'cns.defects.pdf')
generate.heatmap.succinct(names(mi[31:36]), 2.5, 'Congenital Anomalies of the Eye and Childhood Cancer', 'eye.defects.pdf')
generate.heatmap.succinct(names(mi[37:40]), 2.5, 'Congenital Anomalies of the Ear, Face and Neck and Childhood Cancer', 'ear.face.neck.defects.pdf')
generate.heatmap.succinct(names(mi[41:66]), 2.5, 'Congenital Anomalies of the Heart and Childhood Cancer', 'heart.defects.pdf')
generate.heatmap.succinct(names(mi[67:71]), 2.5, 'Congenital Anomalies of the Respiratory System and Childhood Cancer', 'respsys.defects.pdf')
generate.heatmap.succinct(names(mi[72:74]), 2.5, 'Orofacial Clefts and Childhood Cancer', 'orofacial.clefts.pdf')
generate.heatmap.succinct(names(mi[73:83]), 2.5, 'Congenital Anomalies of the Digestive System and Childhood Cancer', 'digestive.defects.pdf')
generate.heatmap.succinct(names(mi[84:91]), 2.5, 'Congenital Anomalies of the Genitourinary System and Childhood Cancer', 'genitourinary.defects.pdf')
generate.heatmap.succinct(names(mi[92:102]), 2.5, 'Congenital Anomalies of the Musculoskeletal System and Childhood Cancer', 'muscskel.defects.pdf')
generate.heatmap.succinct(names(mi[103]), 2.5, 'Congenital Anomalies of the Integument and Childhood Cancer', 'integument.defects.pdf')




# Heatmaps: generate verbose plots ----------------------------------------
load("Z:/Jeremy/GOBACK/Datasets/crude.bd.cc.associations.v20171031.2.rdata")
load("Z:/Jeremy/GOBACK/Datasets/Michigan/mi.v20171006.4.rdata")

setwd('Z:/Jeremy/GOBACK/R outputs/Heatmaps/Verbose/')

#' Create labels
labels = paste(bd.cc.associations$or.new,"\n", bd.cc.associations$ci, sep = "")
labels = ifelse(labels == "NA\n", NA, labels)

#' Generate plots
generate.heatmap.verbose(names(mi[104:111]), 2, 'Chromosomal Anomalies and Childhood Cancer','chromosomal.anomalies.pdf')
generate.heatmap.verbose(names(mi[22:30]), 2, 'Congenital Anomalies of the CNS and Childhood Cancer', 'cns.defects.pdf')
generate.heatmap.verbose(names(mi[31:36]), 2, 'Congenital Anomalies of the Eye and Childhood Cancer', 'eye.defects.pdf')
generate.heatmap.verbose(names(mi[37:40]), 2, 'Congenital Anomalies of the Ear, Face and Neck and Childhood Cancer', 'ear.face.neck.defects.pdf')
generate.heatmap.verbose(names(mi[41:66]), 1.9, 'Congenital Anomalies of the Heart and Childhood Cancer', 'heart.defects.pdf')
generate.heatmap.verbose(names(mi[67:71]), 2, 'Congenital Anomalies of the Respiratory System and Childhood Cancer', 'respsys.defects.pdf')
generate.heatmap.verbose(names(mi[72:74]), 2, 'Orofacial Clefts and Childhood Cancer', 'orofacial.clefts.pdf')
generate.heatmap.verbose(names(mi[73:83]), 2, 'Congenital Anomalies of the Digestive System and Childhood Cancer', 'digestive.defects.pdf')
generate.heatmap.verbose(names(mi[84:91]), 2, 'Congenital Anomalies of the Genitourinary System and Childhood Cancer', 'genitourinary.defects.pdf')
generate.heatmap.verbose(names(mi[92:102]), 2, 'Congenital Anomalies of the Musculoskeletal System and Childhood Cancer', 'muscskel.defects.pdf')
generate.heatmap.verbose(names(mi[103]), 2, 'Congenital Anomalies of the Integument and Childhood Cancer', 'integument.defects.pdf')



# Heatmaps: generate plots for special instances --------------------------

#' Verbose plot of all BD-CC associations, without text.
tmp <- bd.cc.associations
tmp$OR <- ifelse(tmp$OR > 25, 25, tmp$OR)

print(ggplot(data = tmp, aes(x = cancer, y = birth.defect)) +
  geom_tile(aes(fill = OR),color = "white") +
#  geom_text(aes(fill = data$or.new, label = labels),size = text.size) +
  scale_fill_gradient2(low = "green", mid = "grey90",high = "indianred1",midpoint = 1, na.value = "grey90", limits = c(0,25)) + 
  ggtitle('All Congenital Anomaly-Childhood Cancer Associations') +
  labs(y = "Congenital Anomaly", x = "Cancer", caption = 'To improve color contrast, maximum color saturation is reached at an OR of 25') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)))

#' Succinct plot of all BD-CC associations, without text.
tmp <- bd.cc.associations
tmp <- tmp[!is.na(tmp$OR), ]
tmp$OR <- ifelse(tmp$OR > 25, 25, tmp$OR)

print(ggplot(data = tmp, aes(x = cancer, y = birth.defect)) +
        geom_tile(aes(fill = OR),color = "black") +
        #  geom_text(aes(fill = data$or.new, label = labels),size = text.size) +
        scale_fill_gradient2(low = "green", mid = "white",high = "indianred1",midpoint = 12.5, na.value = "grey90", limits = c(0,25)) + 
        ggtitle('All Congenital Anomaly-Childhood Cancer Associations') +
        labs(y = "Congenital Anomaly", x = "Cancer", caption = 'To improve color contrast, maximum color saturation is reached at an OR of 25') +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45,hjust = 1)))

#' A recreation of the WA State heatmap
setwd('Z:/Jeremy/GOBACK/R outputs/Heatmaps')
generate.heatmap.verbose(names(mi[c(22,41,69,75,84,104,92,103)]), 2, 'Congenital Anomalies and Risk of Childhood Cancer', 'wa.state.verbose.pdf')
generate.heatmap.succinct(names(mi[c(22,41,69,75,84,104,92,103)]), 2.5, 'Congenital Anomalies and Risk of Childhood Cancer', 'wa.state.succinct.pdf')




# New heatmaps: recoding --------------------------------------------------



#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' A few changes requested at today's meeting.
#' 
#'  - Change heatmap color scheme.  Label associations by a four level 
#'  categorical variable: deep red OR > 25; light red OR 5-24; 
#'  yellow OR 1-4; gray effect is NS.
#'  
#'  - Generate a table with all CNS tumors: gct.intra, ependymoma, pnet,
#'  medullo, astro, cns.other. 
#'  
#'  - Philip specifically requested certain heatmaps be re-sent when 
#'  complete.  Refer to Google Keep for list.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
load("Z:/Jeremy/GOBACK/Datasets/crude.bd.cc.associations.v20171031.2.rdata")

#' Generate new labeling variable.
bd.cc.associations$signif.cat <- 9999
bd.cc.associations$signif.cat <- factor(
                                        ifelse(bd.cc.associations$signif.cat == 9999 & is.na(bd.cc.associations$or), 6,
                                               ifelse(bd.cc.associations$ci.lower < 1 & bd.cc.associations$ci.upper > 1, 2,
                                                      ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$ci.upper > 1 & bd.cc.associations$or > 1 & bd.cc.associations$or <= 5, 3,
                                                             ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$ci.upper > 1 & bd.cc.associations$or > 5 & bd.cc.associations$or <= 25, 4,
                                                                    ifelse(bd.cc.associations$ci.lower > 1 & bd.cc.associations$ci.upper > 1 & bd.cc.associations$or > 25, 5,
                                                                           ifelse(bd.cc.associations$ci.lower < 1 & bd.cc.associations$ci.upper < 1 & bd.cc.associations$or < 1, 1, NA)))))),
                                        levels = c(1:6),
                                        labels = c('< 1','Null','1-5', '5.01-25','>25','Not Tested'))

save(bd.cc.associations, file = 'Z:/Jeremy/GOBACK/Datasets/crude.bd.cc.assocations.v20171101.1.rdata')



# New heatmaps: all BD-CC associations ------------------------------------
setwd('Z:/Jeremy/GOBACK/R outputs/Heatmaps/')

vector.of.defects <- names(mi[22:111])

data <- arrange(filter(bd.cc.associations, birth.defect %in% vector.of.defects), desc(cancer))

vector.of.cancers <- c()

for (i in unique(data$cancer)){
  
  tmp <- data[data$cancer == i, ]
  flag <- unique(tmp$or)
  
  if (length(flag) == 1){
    
    next
    
  }
  
  else {
    
    vector.of.cancers <- c(vector.of.cancers, i)
  }
}

data <- data[data$cancer %in% vector.of.cancers, ]

defects <- c()

for (i in unique(data$birth.defect)){
  
  tmp <- data[data$birth.defect == i, ]
  flag <- unique(tmp$or)    
  
  if (length(flag) == 1){
    
    next
    
  }
  
  else {
    
    defects <- c(defects, i)
  }
}

data <- data[data$birth.defect %in% defects, ]

plot<-  ggplot(data = data, aes(x = cancer, y = birth.defect)) +
  geom_tile(aes(fill = signif.cat), color = 'black') + 
  scale_fill_manual(values = c('green','grey50','indianred1','firebrick1','darkred','white'), drop = FALSE) +
  guides(fill = guide_legend(title='Odds Ratio')) +
  ggtitle('All Congenital Anomaly-Childhood Cancer Associations') +
  labs(x = "Cancer", y = "Anomaly") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

print(plot)

ggsave(filename = 'all.bd.cc.associations.pdf', scale = 1.75)




# New heatmaps: CNS defects -----------------------------------------------


vector.of.defects <- names(mi[22:30])

data <- arrange(filter(bd.cc.associations, birth.defect %in% vector.of.defects), desc(cancer))

vector.of.cancers <- c()

for (i in unique(data$cancer)){
  
  tmp <- data[data$cancer == i, ]
  flag <- unique(tmp$or)
  
  if (length(flag) == 1){
    next
  }
  else {
    vector.of.cancers <- c(vector.of.cancers, i)
    
  }
}

data <- data[data$cancer %in% vector.of.cancers, ]

defects <- c()

for (i in unique(data$birth.defect)){
  
  tmp <- data[data$birth.defect == i, ]
  flag <- unique(tmp$or)    
  
  if (length(flag) == 1){
    next
  }
  else {
    defects <- c(defects, i)
  }
}

data <- data[data$birth.defect %in% defects, ]

plot<-  ggplot(data = data, aes(x = cancer, y = birth.defect)) +
  geom_tile(aes(fill = signif.cat), color = 'black') + 
  scale_fill_manual(values = c('green','grey50','indianred1','firebrick1','darkred','white'), drop = FALSE) +
  guides(fill = guide_legend(title='Odds Ratio')) +
  ggtitle('Congenital CNS Anomalies and the Risk of Childhood Cancers') +
  labs(x = "Cancer", y = "Anomaly") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

print(plot)

ggsave(filename = 'cns.defects.pdf', scale = 1.75)

# New heatmaps: CNS tumros ------------------------------------------------
vector.of.cancers <- c('astro','pnet','ependymoma','gct.intra','medullo','cns.other','cns.any')
data <- arrange(filter(bd.cc.associations, cancer %in% vector.of.cancers), desc(cancer))

defects <- c()

for (i in unique(data$birth.defect)){
  
  tmp <- data[data$birth.defect == i, ]
  flag <- unique(tmp$or)    
  
  if (length(flag) == 1){
    next
  }
  else {
    defects <- c(defects, i)
  }
}

data <- data[data$birth.defect %in% defects, ]

plot<-  ggplot(data = data, aes(x = cancer, y = birth.defect)) +
  geom_tile(aes(fill = signif.cat), color = 'black') + 
  scale_fill_manual(values = c('green','grey50','indianred1','firebrick1','darkred','white'), drop = FALSE) +
  guides(fill = guide_legend(title='Odds Ratio')) +
  ggtitle('Congenital Anomalies and the Risk of Childhood CNS Cancers') +
  labs(x = "Cancer", y = "Anomaly") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

print(plot)

ggsave(filename = 'cns.tumors.and.all.birth.defects.pdf', scale = 1.75)

# New heatmaps: genitourinary defects -------------------------------------

vector.of.defects <- names(mi[84:91])

data <- arrange(filter(bd.cc.associations, birth.defect %in% vector.of.defects), desc(cancer))

vector.of.cancers <- c()

for (i in unique(data$cancer)){
  
  tmp <- data[data$cancer == i, ]
  flag <- unique(tmp$or)
  
  if (length(flag) == 1){
    next
  }
  else {
    vector.of.cancers <- c(vector.of.cancers, i)
  }
}

data <- data[data$cancer %in% vector.of.cancers, ]

defects <- c()

for (i in unique(data$birth.defect)){
  
  tmp <- data[data$birth.defect == i, ]
  flag <- unique(tmp$or)    
  
  if (length(flag) == 1){
    next
  }
  else {
    defects <- c(defects, i)
  }
}

data <- data[data$birth.defect %in% defects, ]

plot<-  ggplot(data = data, aes(x = cancer, y = birth.defect)) +
  geom_tile(aes(fill = signif.cat), color = 'black') + 
  scale_fill_manual(values = c('green','grey50','indianred1','firebrick1','darkred','white'), drop = FALSE) +
  guides(fill = guide_legend(title='Odds Ratio')) +
  ggtitle('Congenital Genitourinary Anomalies and the Risk of Childhood Cancers') +
  labs(x = "Cancer", y = "Anomaly") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

print(plot)

ggsave(filename = 'genitourinary.defects.pdf', scale = 1.75)



# Figure 1 ----------------------------------------------------------------

setwd('Z:/Jeremy/GOBACK/Datasets/') 
load('figure1.v20180126.1.rdata')
load('goback.no.chrom.v20180122.1.rdata')

#' Need an updated defects matrix.
defects <- names(goback.nochrom[c(22,30,35,38,61,65,68,76,83,94,95)])
cancers <- c(colnames(goback.nochrom[138:147]))

bd.cc.associations <- data.frame(defect = rep(defects, each = 10), cancer = rep(cancers, 11))
bd.cc.associations <- left_join(bd.cc.associations, heatmap, by = c('defect','cancer'))
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

#' Convert axis variables to factor form for prettier display.
bd.cc.associations$cancer <- factor(bd.cc.associations$cancer, 
                                    levels = c("leu.any",     "lym.any",     "cns.any",     "pns.any",     "renal.any",   
                                               "hepatic.any", "bone.any",    "rms.any",     "soft.any",    "gct.any"),
                                    labels = c('Any Leukemia', 'Any Lymphoma', 'Any CNS Tumor', 'Any PNS Tumor', 'Any Renal Tumor',
                                               'Any Hepatic Tumor', 'Any Bone Tumor', 'Any Rhabdomyosarcoma', 'Any Soft Tissue Tumor', 
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



# Figure 1 revised --------------------------------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Philip requested changes to the cancer variables in the heatmap.
#' 
#' Remove any leukemia.  Replace with columns for ALL and AML.
#' Remove and lymphoma.  Replace with columns for HL and NHL.
#' Add columns for medulloblastoma and astrocytoma.
#' Remove any PNS tumor; replace with neuroblastoma.
#' Remove any renal tumor; replace with Wilms.
#' Remove any hepatic tumor; replace with hepatoblastoma.
#' Remove any soft tissue; replace with soft.other labelled as "Any NRSTS"
#' (non-rhabdomyosarcoma soft tissue sarcoma).
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

setwd('Z:/Jeremy/GOBACK/Datasets/') 
load('goback.no.chrom.v20180122.1.rdata')
load('goback.cox.ph.results.v20180124.1.rdata')

#' Need an updated defects matrix.
defects <- names(goback.nochrom[c(22,30,35,38,61,65,68,76,83,94,95)])
cancers <- c('all','aml','hl','nhl','cns.any','astro','medullo','neuro','nephro','hepato','bone.any','osteo','rms.any','soft.other','gct.any')

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

#' Convert axis variables to factor form for prettier display.
bd.cc.associations$cancer <- factor(bd.cc.associations$cancer, 
                                    levels = c('all','aml','hl','nhl','cns.any','astro','medullo',
                                               'neuro','nephro','hepato','bone.any','osteo','rms.any','soft.other','gct.any'),
                                    labels = c('ALL', 'AML', "Hodgkin's Lymphoma", "Non-Hodgkin's Lymphoma", 'Any CNS Tumor', 'Astrocytoma','Medulloblastoma',
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

#' No hits for osteosarcoma.  Remove.
defects <- names(goback.nochrom[c(22,30,35,38,61,65,68,76,83,94,95)])
cancers <- c('all','aml','hl','nhl','cns.any','astro','medullo','neuro','nephro','hepato','bone.any','rms.any','soft.other','gct.any')

bd.cc.associations <- data.frame(defect = rep(defects, each = 14), cancer = rep(cancers, 11))
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

#' Convert axis variables to factor form for prettier display.
bd.cc.associations$cancer <- factor(bd.cc.associations$cancer, 
                                    levels = c('all','aml','hl','nhl','cns.any','astro','medullo',
                                               'neuro','nephro','hepato','bone.any','rms.any','soft.other','gct.any'),
                                    labels = c('ALL', 'AML', "Hodgkin's Lymphoma", "Non-Hodgkin's Lymphoma", 'Any CNS Tumor', 'Astrocytoma','Medulloblastoma',
                                               'Neuroblastoma',"Wilm's Tumor",'Hepatoblastoma','Any Bone Tumor','Rhabdomyosarcoma','Non-RMS Soft Tissue Sarcoma',
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
