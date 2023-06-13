require(dplyr); require(stringr); require(circlize)

load('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded_datasets/goback.coxph.results.v20180612.rdata')

load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old_Datasets/goback.v20180611.rdata")

# Specify plot data and parameters ----------------------------------------

#' Create a vector of birth defect column names from an old GOBACK dataset.
defects <- names(goback)[22:106]

rm(goback); gc()

#' Restrict to individual named phenotypes.
defects <- defects %>% 
  subset(str_detect(., 'any.|conganomalies.|other.|single.gene', negate = T))

#' Move choanal atresia to MSK group.
defects <- defects[c(1:32,34:56,33,57:63)]

#' Plot Wilms tumor.
plot.data <- 
  
  data.frame(cancer = rep('nephro', length(defects)),
             defect = defects) %>% 
  
  left_join(select(goback.coxmodels, defect, cancer, hr, p.val.coef),
            by = c('defect' = 'defect', 'cancer' = 'cancer'))

plot.data$hr <- ifelse(is.na(plot.data$hr), 1, plot.data$hr)
plot.data$p.val.coef <- ifelse(is.na(plot.data$p.val.coef), 0.1, plot.data$p.val.coef)
plot.data$organ.system <- c(rep(0, 6), # CNS
                            rep(1, 4), # Eye/Ear
                            rep(2, 22), # Heart
                            rep(3, 3), # Clefts
                            rep(4, 6), # GI
                            rep(5, 5), # GU
                            rep(6, 10), # MSK
                            rep(7, 7))
plot.data$organ.system <- factor(plot.data$organ.system,
                                 labels = c('CNS','EYE/EAR','CARDIOVASCULAR','CLEFTS','GI','GU','MUSCULOSKELETAL','SYNDROMES'))

#' Some plotting parameters.
plot.colors <- c('black', rep('dodgerblue', 6), rep('limegreen',4), rep('red',22), rep('gold',3),
                        rep('dodgerblue', 6), rep('limegreen', 5), rep('red', 10), rep('gold',7))

col.fun <- colorRamp2(sqrt(range(plot.data$hr, na.rm = T)), c("#FFEEEE", "#FF0000"), transparency = 0.5)

plot.order <- c('nephro', plot.data$defect[1:(nrow(plot.data))])

# Generate plot -----------------------------------------------------------

svg('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/R outputs/Chord diagrams/wilms_tumor_20230217.svg', height = 8, width = 8)

#' A vector of spacings between different defects. 
#' Increases visual separation between A/M and co-occurring defects, and to a lesser extent, between categories of co-occurring defects.
circos.clear()

circos.par(gap.after = c(16,
                         rep(1.75, 5),
                         3,
                         rep(1.75, 3),
                         3,
                         rep(1.75, 21),
                         3,
                         rep(1.75, 2),
                         3,
                         rep(1.75, 5),
                         3,
                         rep(1.75, 4),
                         3,
                         rep(1.75, 9),
                         3,
                         rep(1.75, 6),
                         16),
           start.degree = 329)

chordDiagram(
  x = plot.data[plot.data$cancer == 'nephro', 1:3], 
  preAllocateTracks = list(track.height = 0.1),
  order = plot.order, 
  grid.col = plot.colors, 
  col = col.fun,
  directional = -1, 
  diffHeight = uh(5, 'mm'),
  direction.type = c('diffHeight','arrows'), 
  link.arr.type = 'big.arrow',
  annotationTrack = 'grid')

title("\n Hazard Ratio of Wilms Tumor among Children with Birth Defects")

highlight.sector(sector.index = c('nephro'), track.index = 1, text = 'WILMS TUMOR', 
                 facing = 'bending.inside', niceFacing = T,
                 col = 'black', text.col = 'white', font = 2, cex = 1.2)

highlight.sector(sector.index = plot.order[2:7], track.index = 1, text = 'CNS',
                 col = 'dodgerblue', text.col = 'white', font = 2, cex = 1)

highlight.sector(sector.index = plot.order[8:11], track.index = 1, text = 'EYE/EAR',
                 cex = 0.7, col = 'limegreen', text.col = 'white', facing = 'bending.inside', font = 2)

highlight.sector(sector.index = plot.order[12:33], 
                 track.index = 1, text = 'CARDIOVASCULAR', facing = 'bending.inside',
                 col = 'red', text.col = 'white', 
                 cex = 1, niceFacing = T, font = 2)

highlight.sector(sector.index = plot.order[34:36], 
                 track.index = 1, text = 'CLEFTS', facing = 'bending.inside',
                 col = 'gold', text.col = 'white', cex = 0.6, font = 2)

highlight.sector(sector.index = plot.order[37:42], 
                 track.index = 1, text = 'GI', cex = 1,
                 col = 'dodgerblue', text.col = 'white', niceFacing = T, font = 2)

highlight.sector(sector.index = plot.order[43:47], 
                 track.index = 1, text = 'GENITOURINARY', cex = 1, facing = 'bending.inside',
                 col = 'limegreen', text.col = 'white', niceFacing = T, font = 2)

highlight.sector(sector.index = plot.order[48:57], 
                 track.index = 1, text = 'MUSCULOSKELETAL', facing = 'bending.inside',
                 col = 'red', text.col = 'white', niceFacing = T, font = 2, cex = 0.9)

highlight.sector(sector.index = plot.order[58:64], 
                 track.index = 1, text = 'SYNDROMES', facing = 'bending.inside',
                 col = 'gold', text.col = 'white', cex = 0.9, niceFacing = F, font = 2)

dev.off()
