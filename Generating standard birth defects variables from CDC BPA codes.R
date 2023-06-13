
# Notes/instructions ------------------------------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Authored: 2019.12.11.
#' 
#' Last updated: 2020.01.17.
#' 
#' This script defines a function compute.bd.vars which computes birth defects
#' variables for the GOBACK dataset based on CDC BPA codes, and based on 
#' the coding rules defined in the required "defects" data frame.
#' 
#' Notes:
#' - For performance reasons, the function is only designed to be 
#' run in subjects with one or more birth defects. Including children with 
#' no birth defects means certain death.
#' 
#' - The input data frame should have one row per subject. The first colum
#' of the data frame should be a unique study id. All birth defects 
#' diagnosed in that subject should be recorded in columns 2:N, and should
#' be 7-digit character strings: 3 digits, a period, then 3 more digits, e.g.
#' "746.700". Do not supply any other columns.
#' 
#' - Supplying data in some other format means certain death.
#' 
#' - The function is written to ignore possible/probable codes (generally, 
#' ones where the last digit is 8). These will NOT be used to diagnose 
#' the resulting birth defects variables.
#' 
#' - Appendix 3.1 from the NBDPN Birth Defects Surveillance Guidelines, 
#' revised 3/2017, was the primary document used to determine the codes that 
#' represent each birth defect. For phenotypes not included in that document,
#' the Texas Birth Defects Registry dx_map crosswalk was used.
#' 
#' - Each row in the "defects" data frame is a birth defect. 'diagnois' is
#' the name of that defect. 'code.range' has one or more regular expressions
#' for matching via stringr::str_detect that indicate who will be flagged as
#' having that defect. "canonical.codes" is a stricter regular expression
#' representing the expected codes based on birth defects registry documents.
#' Thus, two columns will be produced for every row in the defects data frame:
#' One indiates whether the child has the index birth defect (1 = yes, 
#' NA = no, but has one or more other birth defects). The second indicates 
#' whether the child has one of the "canonical" codes from the registry docs.
#' This may be useful for finding data entry errors or special cases.
#' 
#' - The function makes use of the dplyr and stringr packages. These
#' will need to be installed and loaded prior to use.
#' 
#' - All specific phenotypes tabulated by this function, and ONLY those 
#' phenotypes, are considered major birth defects and used to compute the 
#' number of major birth defects. 
#' 
#' - Some variables were renamed to better reflect the phenotypes they 
#' represent. These are atrial septal defect (now ASD or PFO);
#' holoprosencephaly (now reduction deformities of brain); anotia/microtia
#' (now anotia, microtia, or anomalies causing hearing impairment).
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

# Define the function -----------------------------------------------------

compute.bd.vars <- function (codes.data.frame, id.column.as.quoted.string) {
 
  require(dplyr); require(stringr)
  
#' Minor birth defects. Others will be taken to be major. Derived from Rasmussen et al., Birth Defects Research Part A, 2003 
  minor.birth.defects <- c('^742.28','^742.385', '^742.40', '^743.4[345]', '^743.6[05]','^743.63[05]', '^743.8[01]',
                           '^744.1[12]', '^744.2[0283]', '^744.24[56]', '^744.[49]1', '^744.[59]0', '^744.8[23]',
                           '^745.50', '^746.0[28]', '^746.105', '^746.48', '^746.[468]0', '^746.86', '^746.886',
                           '^746.9[09]', '^747.[05]0', '^747.325', '^748.18', '^748.36', '^748.51', '^748.08', '^749.19',
                           '^750.[05]0', '^750.1[128]', '^750.2[4678]', '^751.01', '^751.58', '^751.62', '^752.4[345678]',
                           '^752.5[0-2]', '^752.6[02]5', '^752.621', '^752.8[128]', '^752.86[05]', '^753.20', '^753.33',
                           '^753.7[01]', '^754.0[02347]', '^754.05[05]', '^754.31', '^754.43', '^754.52', '^754.8[01]',
                           '^754.82[05]', '^755.006', '^755.13', '^755.5[0145]', '^755.60[056]', '^755.6[3568]', '^755.610',
                           '^755.616', '^755.64[056]', '^756.08[05]', '^756.[12]0', '^756.79', '^756.86', '^757.20', 
                           '^757.3[19]', '^757.38[056]', '^757.4[58]', '^757.5[1]', '^757.58[05]', '^757.6[3458]', '^759.02',
                           '^759.24', '^759.90', '^216.[0-9]', '^228.0[01]', '^351.00', '^3[67]8.00', '^378.90', '^379.50',
                           '^520.60', '^524.00', '^550.[0-9]', '^553.10', '^608.20', '^685.10', '^767.6', '^777.[16]', '^778.[06]')
  
  print('Computing number of total, major and minor defects in each child.')
  
  birth.defects.results <<- data.frame(studyid = codes.data.frame[ , id.column.as.quoted.string], stringsAsFactors = F)
  
  birth.defects.results$any.birth.defect <<- 0
  
  birth.defects.results$number.defects <<- 0
  
  birth.defects.results$number.major.defects <<- 0
  
  birth.defects.results$number.minor.defects <<- 0
  
#' Computes number of total, major, and minor birth defects.
  for (j in 1:nrow(codes.data.frame)){
    
    index.child.codes <- subset(paste(codes.data.frame[j, 2:ncol(codes.data.frame)], sep = ','),
                               paste(codes.data.frame[j, 2:ncol(codes.data.frame)], sep = ',') != 'NA')
    
    #' Remove possible/probable codes (those ending in 8).
    index.child.codes <- subset(index.child.codes, str_detect(index.child.codes, '[0-9]{3}.[0-9]{2}8$', negate = T))

    index.child.codes <- subset(index.child.codes, 
                                str_detect(index.child.codes, str_c(c('^7[45]','^216','^228.0[01]','^237.7'), collapse = "|")))
    
    birth.defects.results[j, 'any.birth.defect'] <<- ifelse(length(index.child.codes) > 0, 1, birth.defects.results[j, 'any.birth.defect'])
    
    minor.codes <- sum(outer(index.child.codes, paste(minor.birth.defects, sep = '|'), str_count))
    
    major.codes <- length(index.child.codes) - minor.codes
    
    birth.defects.results[j , 'number.defects'] <<- length(index.child.codes)
    
    birth.defects.results[j , 'number.major.defects'] <<- major.codes
    
    birth.defects.results[j , 'number.minor.defects'] <<- minor.codes
    
  }
  
#' Identifies children with each of the specific birth defects phenotypes in "defects."  
  for (i in 1:nrow(defects)){
    
    defect.name <- as.character(defects[i , 'diagnosis'])
    
    print(paste('Identifying children with', defect.name))

    code.range <- stringr::str_split(defects[i, 'code.range'], ',', simplify  = T)
    
    canonical.codes <- str_split(as.character(defects$canonical.codes[i]), ',', simplify  = T)
    
    birth.defects.results[, defect.name] <<- as.numeric(NA)
    
    birth.defects.results[, paste0(defect.name,'.canonical.code')] <<- as.numeric(NA)
    
    for (j in 1:nrow(codes.data.frame)){
      
#' A list of the birth defects codes recorded for the index child, excluding NA values and anything ending in 8.
#' Do any of them fall within the ranges we're presently searching? Are any "canonical"?
      index.child.codes <- subset(paste(codes.data.frame[j, 2:ncol(codes.data.frame)], sep = ','),
                                  paste(codes.data.frame[j, 2:ncol(codes.data.frame)], sep = ',') != 'NA')
      index.child.codes <- subset(index.child.codes, str_detect(index.child.codes, '[0-9]{3}.[0-9]{2}8$', negate = T))

      has.code.in.range <- any(outer(index.child.codes, code.range, str_detect))
      
      has.canonical.code <- any(outer(index.child.codes, canonical.codes, str_detect))

      if (has.code.in.range == T) {
        
        birth.defects.results[j, defect.name] <<- 1
        
        birth.defects.results[j, paste0(defect.name,'.canonical.code')] <<- ifelse(has.canonical.code == T, 1, 0)
        
      }
      
      else {
        
        next
        
      }
      
    }
    
  }
  
#' No spina bifida DX if the child also has anencephalus. No hydrocephalus DX if the child also has spina bifida.
  print('Recoding spina bifida and hydrocephalus cases.')

  recode <- which(birth.defects.results$spina.bifida.without.anencephalus == 1 & birth.defects.results$anencephalus == 1)
  
  birth.defects.results[recode, 'spina.bifida.without.anencephalus'] <- NA
  
  recode <- which(birth.defects.results$hydrocephalus.without.spina.bifida == 1 & birth.defects.results$spina.bifida.without.anencephalus == 1)
  
  birth.defects.results[recode, 'hydrocephalus.without.spina.bifida'] <- NA
  
#' Compute higher-level CHD variables. See Botto et al., Birth Defects Research, 2007.
  print('Computing level 2 and 3 CHD variables.')

   birth.defects.results$rvot.defects <<- 
    ifelse(rowSums(birth.defects.results[c('pulm.valve.atresia.stenosis', 'ebstein.anomaly')], na.rm = T)  > 0, 1, NA)
  
  birth.defects.results$septal.defects <<- 
    ifelse(rowSums(birth.defects.results[c('ventricular.septal.defect','asd.or.pfo')], na.rm = T) > 0, 1, NA)
  
  birth.defects.results$lvot.defects <<- 
    ifelse(rowSums(birth.defects.results[c('hypoplastic.left.heart.syndrome', 'interrupted.aortic.arch', 'coarctation.of.aorta', 'aortic.valve.stenosis')], na.rm = T) > 0, 1, NA)
  
  birth.defects.results$conotruncal.defects <<- 
    ifelse(rowSums(birth.defects.results[c('common.truncus','tetralogy.of.fallot','transposition.of.great.vessels')], na.rm = T) > 0, 1, NA)
  
}

# Call the function -------------------------------------------------------

require(tictoc)

bd.codes.merge <- readRDS("//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.merge.v20210806.rds")

defects <- readRDS('//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/birth.defects.lookup.table.v20220211.rds')

#' Subsets of different sizes, for testing new code.
#' bd.codes.ok.nc.tx <- slice(bd.codes.ok.nc.tx, 1:100)
#' bd.codes.ok.nc.tx$runif <- runif(nrow(bd.codes.ok.nc.tx),0,1)
#' bd.codes.ok.nc.tx <- slice(arrange(bd.codes.ok.nc.tx, runif), 1:50000)

tic(); compute.bd.vars(bd.codes.merge, 'studyid'); toc()

saveRDS(birth.defects.results, 
        file = '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/computed.birth.defects.variabls.cdc.bpa.states.v20220214.rds')


# Scratch paper -----------------------------------------------------------

bd.codes.merge %<>% mutate(runif = runif(nrow(.)))

tmp <- filter(bd.codes.merge, runif < 0.05)

tmp %<>% select(-runif)

tic(); compute.bd.vars(tmp, 'studyid'); toc()
