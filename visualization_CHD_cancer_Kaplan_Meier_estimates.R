require(tidyverse); require(survival); require(survminer); require(ggsci); require(viridis)

#' Load current GOBACK file.
goback <- readRDS("//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/goback_20230523.rds")

#' Organ system-level major defects variables, apart from the heart.
#' Will be used to identify children with CHD complicated by the co-occurrence of other major anomalies.
other.systems <- goback %>% 
  select(-c(major.heart.circulatory.anomaly, number.major.defects)) %>% 
  names() %>% 
  str_subset('major')

#' Remove children with genetic or chromosomal syndromes or unknown sex, and right censor at 18.
goback <- goback %>% 
  filter(genetic.anomaly == 0 | is.na(genetic.anomaly), sex != 'Unknown') %>% 
  select(birth.defect, cancer, all:gct.any, person.yrs, sex, birth.wt, gest.age, state, m.age, major.heart.circulatory.anomaly, all_of(other.systems), number.major.defects) %>% 
  mutate(cancer = ifelse(cancer == 1 & person.yrs > 18, 0, cancer),
         person.yrs = ifelse(person.yrs > 18, 18, person.yrs),
         non.chd.defects = rowSums(.[, other.systems], na.rm = T),
         chd.plus = factor(ifelse(major.heart.circulatory.anomaly == 1 & non.chd.defects >= 1, 2,
                                  ifelse(major.heart.circulatory.anomaly == 1 & non.chd.defects == 0, 1, 
                                         ifelse(birth.defect == 0, 0, NA))),
                           labels = c('No birth defect', 'CHD only', 'CHD plus')))

#' Categories for number of major defects.
goback <- goback %>% 
  mutate(major.defects.cat = factor(ifelse(number.major.defects >= 4, 4, number.major.defects),
                                    labels = c('0', '1', '2', '3', '4 or more')))

# K-M curves for number of major defects ----------------------------------

fit <- surv_fit(Surv(person.yrs, cancer) ~ major.defects.cat, data = goback)

plot <- ggsurvplot(fit, 
                   data = goback, 
                   fun = 'event',
                   xlim = c(0,18),
                   break.x.by = 2,
                   ylim = c(0, 0.03),
                   axes.offset = T,
                   conf.int = F,
                   risk.table = F,
                   pval = T,
                   pval.coord = c(2.7, .015),
                   pval.method = T,
                   pval.method.coord = c(1, .015),
                   surv.scale = 'percent',
                   censor = F,
                   legend.labs = unique(levels(goback$major.defects.cat)),
                   legend.title = 'Number of major birth defects',
                   xlab = 'Time (yrs)',
                   ylab = 'Cumulative incidence',
                   ggtheme = theme_bw())

plot <- plot$plot +
  
  scale_color_viridis(discrete = T,
                      limits = c('4 or more', '3', '2', '1', '0')) + 
  
  theme(axis.text = element_text(size = 15, face = 'bold'),
        axis.title = element_text(size = 17, face = 'bold'),
        
        panel.grid = element_blank(),
        
        legend.position = c(0.3,0.8),
        legend.text = element_text(size = 15, face = 'bold'),
        legend.title = element_text(size = 17, face = 'bold'),
        legend.background = element_rect(color = 'black', linewidth = 1.25))

plot

svg('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Figures/Kaplan_Meier_number_major_defects_20230530.svg',
    width = 8, height = 6)

plot

dev.off()

# K-M curves for CHD only and CHD+ groups ---------------------------------

fit <- surv_fit(Surv(person.yrs, cancer) ~ chd.plus, data = goback)

plot <- ggsurvplot(fit, 
                   data = goback, 
                   fun = 'event',
                   xlim = c(0,18),
                   break.x.by = 2,
                   ylim = c(0, 0.03),
                   axes.offset = T,
                   conf.int = F,
                   risk.table = F,
                   pval = T,
                   pval.coord = c(2.7, .015),
                   pval.method = T,
                   pval.method.coord = c(1, .015),
                   surv.scale = 'percent',
                   censor = F,
                   legend.labs = c('No birth defect', 'CHD alone (any type)', 'CHD plus extracardiac defects'),
                   legend.title = '',
                   xlab = 'Time (yrs)',
                   ylab = 'Cumulative incidence',
                   ggtheme = theme_bw())

plot <- plot$plot +
  
  scale_color_viridis(discrete = T,
                      limits = c('CHD plus extracardiac defects', 'CHD alone (any type)', 'No birth defect')) + 
  
  theme(axis.text = element_text(size = 15, face = 'bold'),
        axis.title = element_text(size = 17, face = 'bold'),
        
        panel.grid = element_blank(),
        
        legend.position = c(0.3,0.8),
        legend.text = element_text(size = 15, face = 'bold'),
        legend.background = element_rect(color = 'black', linewidth = 1.25))

plot

svg('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Figures/Kaplan_Meier_CHD_with_without_extracardiac_defects_20230531.svg',
    width = 8, height = 6)

plot

dev.off()

# K-M curves for different CHD phenotypes ---------------------------------

defects <- c('major.heart.circulatory.anomaly', 'asd.or.pfo', 
             'ventricular.septal.defect', 'pulm.valve.atresia.stenosis')

labels <- c('Any CHD', 'ASD or PFO', 'VSD', 'Pulmonary atresia')

plots <- vector('list', length(defects))

for (i in seq_along(defects)){
  
  data <- goback %>% 
    select(cancer, person.yrs, all_of(defects[i]))
  
  names(data) <- c('cancer', 'person.yrs', 'anomaly')
  
  fit <- surv_fit(Surv(person.yrs, cancer) ~ anomaly, data = data)
  
  plot <- ggsurvplot(fit, 
                     data = data, 
                     fun = 'event',
                     linetype = c('solid','dashed'),
                     xlim = c(0,18),
                     break.x.by = 2,
                     ylim = c(0, 0.03),
                     axes.offset = T,
                     conf.int = F,
                     risk.table = F,
                     pval = T,
                     pval.coord = c(2.8, .015),
                     pval.method = T,
                     pval.method.coord = c(1, .015),
                     surv.scale = 'percent',
                     censor = F,
                     legend.labs = c('No birth defect', labels[i]),
                     legend.title = '',
                     xlab = 'Time (yrs)',
                     ylab = 'cumulative incidence',
                     ggtheme = theme_bw())
  
  plot <- plot$plot +
    
    scale_color_aaas() +
    
    theme(axis.text = element_text(size = 15, face = 'bold'),
          axis.title = element_text(size = 17, face = 'bold'),
          
          panel.grid = element_blank(),
          
          legend.position = c(0.3,0.8),
          legend.text = element_text(size = 15, face = 'bold'))
  
  plots[[i]] <- plot
  
}

svg('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Figures/CHD_Kaplan_Meier_curves_20230526.svg',
    height = 8, width = 12)

plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
          nrow = 2)

dev.off()
