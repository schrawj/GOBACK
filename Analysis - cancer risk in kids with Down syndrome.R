require(tidyverse); require(survival); require(survminer); require(ggsci)

setwd('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/')

load('Datasets/goback.v20200123.rdata')

goback  %<>% filter(any.birth.defect == 0 | trisomy.21 == 1) %>% 
  mutate(solid.tumor = ifelse(cancer == 1 & leu.any == 0 & lym.any == 0, 1, 0))

# Cumulative incidence plots ----------------------------------------------

cancers <- c('all','aml','solid.tumor')

titles <- c('Acute Lymphoblastic Leukemia', 'Acute Myeloid Leukemia', 'Solid Tumors')

legend.categories = c('Non-Down syndrome', 'Down syndrome')

fits <- list()

plots <- list()

for (i in seq_along(cancers)) {

  index.cancer <- cancers[i]
  
  fit <- survfit(Surv(person.yrs, goback[,index.cancer]) ~ trisomy.21, data = goback)
  
  plot <- ggsurvplot(fit = fit,
                     fun = 'event',
                     xlim = c(0,18),
                     break.x.by = 1,
                     ylim = c(0, 0.02),
                     conf.int = F,
                     risk.table = F,
                     surv.scale = 'percent',
                     censor = F,
                     legend.labs = legend.categories,
                     title = titles[[i]],
                     legend.title = 'Down syndrome',
                     ggtheme = theme_bw())
  
  plot$plot <- plot$plot +
    
    theme(
      plot.title = element_text(face = 'bold', size = 14, hjust = 0.5),
      
      axis.title = element_text(face = 'bold', size = 14),
      axis.text = element_text(face = 'bold', size = 12),
      
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      
      legend.title = element_text(face = 'bold', size = 14), 
      legend.text = element_text(face = 'bold', size = 12),
      legend.position = c(0.15,0.8),
      legend.background = element_rect(color = 'black')
    ) +
    
    scale_color_jco() + 
    
    labs(x = 'Time (yrs)',
         y = 'Cumulative Incidence')
  
  plots[[i]] <- plot

}

for (i in 1:length(plots)) { print(plots[[i]]) }


# Compute hazard ratios ---------------------------------------------------

out <- data.frame()

for (i in seq_along(cancers)){
  
  index.cancer <- cancers[i]
  
  cox <- coxph(Surv(person.yrs, goback[, index.cancer]) ~ trisomy.21 + sex + m.age + state, data = goback)
  
  summary.cox <- summary(cox)
  
  cox.estimates <- data.frame(cancer = index.cancer,
                              aHR = summary.cox$conf.int[1,1],
                              ci.lower = summary.cox$conf.int[1,3],
                              ci.upper = summary.cox$conf.int[1,4])
  
  out <- rbind(out, cox.estimates)
  
}
