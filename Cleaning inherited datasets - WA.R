require(tidyverse); require(magrittr); require(haven)

setwd('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Old Datasets/Washington State/')

# Import and save raw data ------------------------------------------------

wash <- read_dta('baylor2020_data_malf_20220111.dta')

saveRDS(wash, 'washington.state.raw.data.rds')
