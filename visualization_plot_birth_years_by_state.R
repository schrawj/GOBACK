require(tidyverse); require(viridis)

plot.data <- data.frame(state = c('AR','FL','MA','MI',
                                  'NJ', 'NC', 'OK', 
                                  'SC', 'TX'),
                        start = c(1995,1998,2000,1992,
                                  1990,2003,1997,
                                  2008,1999),
                        stop = c(2011,2013,2017,2018,
                                 2018,2012,2012,
                                 2018,2017)) %>% 
  mutate(state = as.factor(state)) %>% 
  mutate(state = fct_relevel(state, 'NJ', 'MI','AR','OK','FL','TX','MA','NC','SC'))
  #pivot_longer(!state, names_to = 'timepoint', values_to = 'value')

plot <- ggplot(data = plot.data, aes(x = value, y = state)) + 
  
  geom_segment(aes(y = state, yend = state, x = start, xend = stop, color = state), linewidth = 20) +
  
  scale_color_viridis(discrete = T) + 
  
  scale_x_continuous(limits = c(1990,2020),
                     breaks = seq(1990,2020, by = 5)) + 
  
  scale_y_discrete(limits = rev(levels(plot.data$state))) +
  
  labs(y = '', x = '') + 
  
  theme_classic() + 
  
  theme(legend.position = 'none',
        axis.text = element_text(size = 20, face = 'bold'))

plot

svg('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Figures/GOBACK_birth_years_20230523.svg',
    width = 14, height = 9)

plot

dev.off()
