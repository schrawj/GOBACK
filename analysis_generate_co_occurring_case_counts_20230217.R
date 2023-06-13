require(tidyverse)

goback <- readRDS("//smb-main.ad.bcm.edu/genepi2/GOBACK/Datasets/goback_20230216.rds") %>% 
  filter(birth.defect == 1)

cancers <- goback %>% 
  select(cancer, all:other.tumor) %>% 
  names()
  
defects <- goback %>% 
  select(birth.defect:major.musculoskeletal.anomaly) %>% 
  names()

co.occurring.counts <- data.frame()

#' All subjects with BDs.
for (i in seq_along(cancers)){
  
  index.cancer <- cancers[i]
  
  for (j in seq_along(defects)){
    
    index.defect <- defects[j]
    
    new.counts <- goback %>% 
      select(all_of(index.defect), all_of(index.cancer)) %>% 
      rename(defect = index.defect, cancer = index.cancer) %>% 
      mutate(across(defect:cancer, ~ factor(.x))) %>% 
      count(defect, cancer, .drop = F) %>% 
      filter(defect == 1, cancer == 1) %>% 
      rename(defect.status = defect, 
             cancer.status = cancer) %>% 
      mutate(defect = index.defect,
             cancer = index.cancer,
             phenotype = 'all.bd.cases') %>% 
      select(phenotype, defect, cancer, defect.status, cancer.status, n)
    
    co.occurring.counts <- rbind(co.occurring.counts, new.counts)

  }
  
  saveRDS(co.occurring.counts,
          '//smb-main.ad.bcm.edu/genepi3/JeremySchraw/GOBACK/Datasets/goback_co_occurring_case_counts_20230224.rds')
  
}

#' Non-syndromic BDs.
for (i in seq_along(cancers)){
  
  index.cancer <- cancers[i]
  
  for (j in seq_along(defects)){
    
    index.defect <- defects[j]
    
    new.counts <- goback %>% 
      filter(genetic.anomaly == 0 | is.na(genetic.anomaly)) %>% 
      select(all_of(index.defect), all_of(index.cancer)) %>% 
      rename(defect = index.defect, cancer = index.cancer) %>% 
      mutate(across(defect:cancer, ~ factor(.x))) %>% 
      count(defect, cancer, .drop = F) %>% 
      filter(defect == 1, cancer == 1) %>% 
      rename(defect.status = defect, 
             cancer.status = cancer) %>% 
      mutate(defect = index.defect,
             cancer = index.cancer,
             phenotype = 'non.syndromic.cases') %>% 
      select(phenotype, defect, cancer, defect.status, cancer.status, n)
    
    co.occurring.counts <- rbind(co.occurring.counts, new.counts)
    
  }
  
  saveRDS(co.occurring.counts,
          '//smb-main.ad.bcm.edu/genepi3/JeremySchraw/GOBACK/Datasets/goback_co_occurring_case_counts_20230224.rds')
  
}