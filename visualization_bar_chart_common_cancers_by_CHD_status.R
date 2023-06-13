require(tidyverse); require(ggsci); require(descr)

#' Load current GOBACK file.
goback <- readRDS("//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/goback_20230523.rds")

#' Remove children with genetic or chromosomal syndromes or unknown sex, right censor at 18, and identify more complex cases.
goback <- goback %>% 
  filter(genetic.anomaly == 0 | is.na(genetic.anomaly), sex != 'Unknown') %>% 
  mutate(cancer = ifelse(cancer == 1 & person.yrs > 18, 0, cancer),
         person.yrs = ifelse(person.yrs > 18, 18, person.yrs))

cases <- goback %>% 
  filter(cancer == 1, (major.heart.circulatory.anomaly == 1 | birth.defect == 0) )

tab <- crosstab(cases$first.primary, cases$major.heart.circulatory.anomaly, prop.c = T)$prop.col %>% 
  as_tibble() %>% 
  pivot_wider(id_cols = `cases$first.primary`,
              names_from = `cases$major.heart.circulatory.anomaly`,
              values_from = n) %>% 
  rename(tumor = `cases$first.primary`,
         no.bd = `0`,
         chd = `1`) %>% 
  filter(tumor %in% c('Precursor cell leukemias', 'Astrocytomas', 'Neuroblastoma and ganglioneuroblastoma',
                      'Nephroblastoma', 'Acute myeloid leukemias', 'Retinoblastoma', 'Hepatoblastoma',
                      'Hodgkin lymphomas')) %>% 
  pivot_longer(!tumor, names_to = 'group', values_to = 'freq') %>% 
  mutate(group = factor(group, labels = c('CHD', 'No birth defect')),
         freq = freq*100,
         tumor = ifelse(tumor == 'Neuroblastoma and ganglioneuroblastoma', 'Neuroblastoma', tumor)) %>% 
  arrange(desc(group), desc(freq)) %>% 
  mutate(tumor = factor(tumor, levels = tumor.order))

tumor.order <- tab$tumor[1:8]

#' Bar chart.
plot <- ggplot(data = tab, aes(x = tumor, y = freq, fill = group)) + 
  
  geom_bar(width = 0.7, stat = 'identity', position = 'dodge') +
  
  theme_bw() + 
  
  scale_fill_npg() +
  
  guides(fill = guide_legend(title = 'BD Status')) + 
  
  labs(y = 'Percent of Diagnoses') +
  
  theme(axis.text = element_text(face = 'bold', size = 15),
        axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = 'bold', size = 15, margin = margin(r = 15)),
        axis.ticks = element_blank(),
        
        legend.title = element_text(face = 'bold', size = 15),
        legend.text = element_text(size = 15, face = 'bold'),
        legend.position = c(0.8, 0.9),
        legend.spacing.x = unit(0.5, 'cm'),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'),

        panel.grid.major.y = element_line(linetype = 'solid', color = 'black'),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) 

print(plot)

svg('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Figures/bar_chart_common_cancers_by_CHD_status_20230531.svg',
    height = 8, width = 6)

plot

dev.off()