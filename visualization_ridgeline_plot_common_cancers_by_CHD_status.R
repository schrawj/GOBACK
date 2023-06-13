require(tidyverse); require(ggridges); require(viridis); require(ggsci)

#' Load current GOBACK file.
goback <- readRDS("//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/goback_20230523.rds")

#' A list of some common childhood cancers to plot.
#' Chose a few for which there were noticeable differences in age at diagnosis and a few for which there were not.
common.cancers <- c('Neuroblastoma and ganglioneuroblastoma', 'Osteosarcomas', 'Acute myeloid leukemias', 
                    'Precursor cell leukemias', 'Hepatoblastoma', 'Hodgkin lymphomas')

#' Remove children with genetic or chromosomal syndromes or unknown sex, right censor at 18, and identify more complex cases.
goback <- goback %>% 
  filter(cancer == 1, (major.heart.circulatory.anomaly == 1 | birth.defect == 0), (genetic.anomaly == 0 | is.na(genetic.anomaly)), sex != 'Unknown') %>% 
  mutate(cancer = ifelse(cancer == 1 & person.yrs > 18, 0, cancer),
         person.yrs = ifelse(person.yrs > 18, 18, 
                      ifelse(person.yrs == 0, 0.01, person.yrs)))

#' Create data for plotting.
plot.data <- data.frame(age = goback$person.yrs, 
                        chd.status = factor(goback$major.heart.circulatory.anomaly, labels = c('no birth defect', 'CHD')),
                        tumor.type = goback$first.primary) %>% 
  filter(tumor.type %in% common.cancers) %>% 
  mutate(tumor.type = ifelse(tumor.type == 'Neuroblastoma and ganglioneuroblastoma',
                             'Neuroblastoma', tumor.type))

#' Mean age at diagnosis will determine order of y-axis.
ages <- aggregate(age ~ tumor.type + chd.status, data = plot.data, mean)
ages <- ages %>% filter(chd.status == 'no birth defect') %>% arrange(age)

plot.data$tumor.type <- factor(plot.data$tumor.type,
                                  levels = ages$tumor.type,
                                  labels = ages$tumor.type)

plot.data$group <- paste0(plot.data$tumor.type,', ', plot.data$chd.status)

group.order <- paste0(rep(levels(plot.data$tumor.type), each = 2),', ', rep(c('no BD', 'CHD')))

#' Create a factor variable to determine the order of the y-axis.
plot.data <- plot.data %>% 
  mutate(group = str_replace(group, 'birth defect', 'BD'),
         group.fct = fct_relevel(factor(group), group.order))

#' Generate Ridgeline plot.
plot <- ggplot(plot.data, aes(x = age, y = group.fct, fill = tumor.type)) + 
  
  geom_density_ridges() +
  
  theme_void() +
  
  labs(x = "Age at cancer diagnosis (years)") +

  scale_fill_viridis(discrete = T) +
  scale_color_viridis(discrete = T) +
  
  scale_x_continuous(limits = c(-5,20),
                     breaks = seq(0, 20, by = 5)) +
  
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(face = 'bold', size = 14,
                                    margin = margin(15,0,0,0)),
        axis.text = element_text(face = 'bold', size = 12, hjust = 1),
        
        legend.title = element_blank(),
        legend.text = element_text(face = 'bold', size = 12),
        legend.position = 'none')

plot 

svg('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Figures/ridgeline_plot_age_at_cancer_DX_by_CHD_status_20230601.svg',
    8, 8)

plot

dev.off()

# Scratch paper -----------------------------------------------------------

tmp <- goback %>% 
  split(.$first.primary) %>% 
  map(~aggregate(person.yrs ~ major.heart.circulatory.anomaly, data = .x, mean)) %>% 
  list_rbind(names_to = 'tumor') %>% 
  pivot_wider(id_cols = tumor,
              names_from = 'major.heart.circulatory.anomaly',
              values_from = 'person.yrs') %>% 
  rename(no.bd = `0`, chd = `1`) %>% 
  filter(!is.na(chd), !is.na(no.bd)) %>% 
  mutate(age.diff = abs(no.bd - chd)) %>% 
  arrange(desc(age.diff))

tab <- table(goback$first.primary, goback$major.heart.circulatory.anomaly) %>% 
  as.data.frame() %>% 
  filter(Var2 == 1, Freq >= 10) %>% 
  pull(Var1)

tmp2 <- filter(tmp, tumor %in% tab)
