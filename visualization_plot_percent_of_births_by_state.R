require(tidyverse); require(descr); require(viridis)

goback <- readRDS("//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/goback_20230418.rds")

goback <- goback %>% 
  select(state)

ar <- rep('AR', times = 629120)

ar <- data.frame(state = ar)

goback <- bind_rows(goback, ar)

CrossTable(goback$state, prop.r = T)

plot.data <- data.frame(state = c('AR', 'FL', 'MA', 'MI', 'NC', 'NJ', 'OK', 'SC', 'TX'),
                        prop = c(2.9, 14.8, 6.1, 15.5, 5.7, 15.1, 3.7, 2.8, 33.5)) %>% 
  arrange(desc(prop)) %>% 
  mutate(state = fct_relevel(state, 'TX', 'MI', 'NJ', 'FL', 'MA', 'MI','NC','OK','AR','SC'))

plot <- ggplot(data = plot.data, aes(x = state, y = prop, color = state, fill = state)) + 
  
  geom_col() + 
  
  labs(x = '', y = "% of Total Cohort") +
  
  scale_fill_viridis(discrete = T) + 
  scale_color_viridis(discrete = T) + 
  
  scale_y_continuous(limits = c(0, 35),
                     breaks = seq(0, 35, by = 5)) +
  
  theme_bw() + 
  
  theme(legend.position = 'none',
        
        axis.text = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size = 22, face = 'bold'),
        axis.title.y = element_text(margin = margin(0,20,0,0)))

svg('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Figures/GOBACK_counts_by_state_20230523.svg',
    width = 10, height = 10)

plot

dev.off()
