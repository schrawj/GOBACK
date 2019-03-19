#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.11.21.
#' 
#' Have been asked to compute HRs and 95% CIs for germ cell tumors 
#' according to specific chromosomal defects.
#' 
#' This isn't actually worth doing, for reasons that become clear when you
#' look at the counts of specific GCTs in children with specific defects.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

load('W:/Old_genepi2/Jeremy/GOBACK/Datasets/goback.chrom.v20180829.rdata')

cancercols <- c(151,122:124)
defectcols <- c(97,98,100:105)

out <- data.frame(defect = as.character(),
                  cancer = as.character(),
                  comorbid.count = as.numeric())

for (i in defectcols){
  
  for (j in cancercols){
    
    tab <- table(goback.chrom[,i], goback.chrom[,j])[2,2]
    
    new.out <- data.frame(cancer = names(goback.chrom[j]),
                          defect = names(goback.chrom[i]),
                          comorbid.count = tab)
    
    out <- rbind(out, new.out)
    
  }
  
}

write.csv(out, file = 'W:/Old_genepi2/Jeremy/GOBACK/R outputs/gct.and.chromosomal.anomalies.comorbid.cases.csv', row.names = FALSE)
