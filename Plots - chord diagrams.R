#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Authored: 2020.04.18.
#' 
#' Last updated: 
#' 
#' Generate chord diagrams for a few cancers. Hoping to use these in the 
#' CAC2 webinar on 4/22.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(dplyr); require(stringr); require(circlize)

load('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/goback.coxph.results.v20180612.rdata')

load("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/goback.v20180611.rdata")

#' Create a vector of birth defect column names from an old GOBACK dataset.
defects <- names(goback)[22:106]

rm(goback); gc()

#' Str_detect not vectorized? Supply each pattern individually.
defects <- defects %>% 
  subset(!str_detect(., 'any.')) %>% 
  subset(!str_detect(., 'conganomalies.')) %>% 
  subset(!str_detect(., 'other.')) %>% 
  subset(!str_detect(., 'single.gene'))

#' Move choanal atresia to MSK group.
defects <- defects[c(1:32,34:56,33,57:63)]

#' Let's plot ALL and hepatoblastoma.
plot.data <- 
  
  data.frame(cancer = rep(c('all','hepato'), each = length(defects)),
                        defect = rep(defects, 2)) %>% 
  
  left_join(select(goback.coxmodels, defect, cancer, hr, p.val.coef),
            by = c('defect' = 'defect', 'cancer' = 'cancer'))

plot.data$hr <- ifelse(is.na(plot.data$hr), 1, plot.data$hr)
plot.data$p.val.coef <- ifelse(is.na(plot.data$p.val.coef), 0.1, plot.data$p.val.coef)
plot.data$organ.system <- rep(c(rep(0, 6), # CNS
                                rep(1, 4), # Eye/Ear
                                rep(2, 22), # Heart
                                rep(3, 3), # Clefts
                                rep(4, 6), # GI
                                rep(5, 5), # GU
                                rep(6, 10), # MSK
                                rep(7, 7)), # Syndromes
                              times  = 2)
plot.data$organ.system <- factor(plot.data$organ.system,
                                 labels = c('CNS','EYE/EAR','CARDIOVASCULAR','CLEFTS','GI','GU','MUSCULOSKELETAL','SYNDROMES'))

plot.colors <- c('black', rep('dodgerblue', 6), rep('limegreen',4), rep('red',22), rep('gold',3),
                          rep('dodgerblue', 6), rep('limegreen', 5), rep('red', 10), rep('gold',7))

col.fun = colorRamp2(range(sqrt(adj.list$o.e.adj)), c("#FFEEEE", "#FF0000"), transparency = 0.5)



# ALL plot ----------------------------------------------------------------

#' Specify plot: ALL or hepatoblastoma.
plot.order <- c('all', plot.data$defect[1:(nrow(plot.data)/2)])

#' A vector of spacings between different defects. 
#' Increases visual separation between A/M and co-occurring defects, and to a lesser extent,
#' between categories of co-occurring defects.
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
  x = plot.data[plot.data$cancer == 'all', 1:3], 
  preAllocateTracks = list(track.height = 0.1),
  order = plot.order, 
  grid.col = plot.colors, 
  col = col.fun,
  directional = -1, 
  diffHeight = uh(5, 'mm'),
  direction.type = c('diffHeight','arrows'), 
  link.arr.type = 'big.arrow',
  annotationTrack = 'grid')

title("\n Adjusted Hazard Ratios for Acute Lymphoblastic Leukemia \nAmong Children with Different Birth Defects")

highlight.sector(sector.index = c('all'), track.index = 1, text = 'ACUTE LYMPHOBLASTIC LEUKEMIA', 
                 facing = 'bending.inside', niceFacing = T,
                 col = 'black', text.col = 'white', font = 2, cex = 1.2)

highlight.sector(sector.index = plot.order[2:7], track.index = 1, text = 'CNS',
                 col = 'dodgerblue', text.col = 'white', font = 2, cex = 1)

highlight.sector(sector.index = plot.order[8:11], track.index = 1, text = 'EYE/EAR',
                 cex = 0.75, col = 'limegreen', text.col = 'white', facing = 'bending.inside', font = 2)

highlight.sector(sector.index = plot.order[12:33], 
                 track.index = 1, text = 'CARDIOVASCULAR', facing = 'bending.inside',
                 col = 'red', text.col = 'white', 
                 cex = 1, niceFacing = T, font = 2)

highlight.sector(sector.index = plot.order[34:36], 
                 track.index = 1, text = 'CLEFTS', facing = 'bending.inside',
                 col = 'gold', text.col = 'white', cex = 0.75, font = 2)

highlight.sector(sector.index = plot.order[37:42], 
                 track.index = 1, text = 'GI', cex = 1,
                 col = 'dodgerblue', text.col = 'white', niceFacing = T, font = 2)

highlight.sector(sector.index = plot.order[43:47], 
                 track.index = 1, text = 'GU', cex = 1,
                 col = 'limegreen', text.col = 'white', niceFacing = T, font = 2)

highlight.sector(sector.index = plot.order[48:57], 
                 track.index = 1, text = 'MUSCULOSKELETAL', facing = 'bending.inside',
                 col = 'red', text.col = 'white', niceFacing = T, font = 2, cex = 1)

highlight.sector(sector.index = plot.order[58:64], 
                 track.index = 1, text = 'GENETIC SYNDROMES', facing = 'bending.inside',
                 col = 'gold', text.col = 'white', cex = 1, niceFacing = F, font = 2)



# Hepatoblastoma plot -----------------------------------------------------

#' Specify the order in which phenotypes are plotted around the circle.
plot.order <- c('hepato', plot.data$defect[((nrow(plot.data)/2)+1):nrow(plot.data)])

#' A vector of spacings between different defects. 
#' Increases visual separation between A/M and co-occurring defects, and to a lesser extent,
#' between categories of co-occurring defects.
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
             x = plot.data[plot.data$cancer == 'hepato', 1:3],
             preAllocateTracks = list(track.height = 0.1),
             order = plot.order, 
             grid.col = plot.colors, 
             col = col.fun,
             directional = -1, 
             diffHeight = uh(5, 'mm'),
             direction.type = c('diffHeight','arrows'), 
             link.arr.type = 'big.arrow',
             annotationTrack = 'grid')

title("\n Adjusted Hazard Ratios for Hepatoblastoma \nAmong Children with Different Birth Defects")

highlight.sector(sector.index = c('hepato'), track.index = 1, text = 'HEPATOBLASTOMA', 
                 facing = 'bending.inside', niceFacing = T,
                 col = 'black', text.col = 'white', font = 2, cex = 1.2)

highlight.sector(sector.index = plot.order[2:7], track.index = 1, text = 'CNS',
                 col = 'dodgerblue', text.col = 'white', font = 2, cex = 1)

highlight.sector(sector.index = plot.order[8:11], track.index = 1, text = '',
                 cex = 0.75, col = 'limegreen', text.col = 'white', facing = 'bending.inside', font = 2)

highlight.sector(sector.index = plot.order[12:33], 
                 track.index = 1, text = 'CARDIOVASCULAR', facing = 'bending.inside',
                 col = 'red', text.col = 'white', 
                 cex = 1, niceFacing = T, font = 2)

highlight.sector(sector.index = plot.order[34:36], 
                 track.index = 1, text = '', facing = 'bending.inside',
                 col = 'gold', text.col = 'white', cex = 0.8, font = 2)

highlight.sector(sector.index = plot.order[37:42], 
                 track.index = 1, text = 'GI', cex = 1,
                 col = 'dodgerblue', text.col = 'white', niceFacing = T, font = 2)

highlight.sector(sector.index = plot.order[43:47], 
                 track.index = 1, text = 'GU', cex = 1,
                 col = 'limegreen', text.col = 'white', niceFacing = T, font = 2)

highlight.sector(sector.index = plot.order[48:57], 
                 track.index = 1, text = 'MUSCULOSKELETAL', facing = 'bending.inside',
                 col = 'red', text.col = 'white', niceFacing = T, font = 2, cex = 0.9)

highlight.sector(sector.index = plot.order[58:64], 
                 track.index = 1, text = 'GENETIC SYNDROMES', facing = 'bending.inside',
                 col = 'gold', text.col = 'white', cex = 1, niceFacing = F, font = 2)
