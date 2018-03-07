#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Received Marlena's (former student of Beth Mueller's) code for 
#' generating BD-CC heatmap.
#' 
#' This is a translation of that code for use in our project.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------



# prep environment --------------------------------------------------------
install.packages('pacman')
install.packages('ggthemes')
install.packages('ReporteRs')
install.packages('gridExtra')

require(dplyr)
require(reshape2)
require(ggplot2)
require(ggthemes)
require(ReporteRs)
require(gridExtra)

setwd('//discovery2/bcm-pedi-genepi2/Jeremy/GOBACK/R scripts/R scripts from WA/')





# Old data: load in old data ----------------------------------------------
ors <- read.delim('./Table4_Data_OR.txt', header = TRUE)
ors2 <- melt(ors, id.vars = 'Cancer.Type')

cis <- read.delim('Table4_Data_CI.txt', header = TRUE)
cis2 <- melt(cis, id.vars = 'Cancer.Type')

ors2 <- cbind(ors2, cis2$value)
names(ors2) = c("cancer","malf","OR","ci")

# create labels
labels = paste(format(ors2$OR, digits = 2),"\n", format(ors2$ci,digits = 2), sep = "")
labels = ifelse(labels == "   NA\n           ", NA, labels)




# Old data: generate plot -------------------------------------------------

# create tile plot with ggplot
print(ggplot(ors2, aes(x = cancer, y = malf)) +
  geom_tile(aes(fill = OR),color = "white") +
  geom_text(aes(fill = ors2$OR, label = labels),size = 2.5) +
  scale_fill_gradient2(low = "green", mid = "white",high = "indianred1",midpoint = 0,trans = "log",na.value = "grey90",limits = c(0.3,25)) +
  ggtitle("Major Malformations") +
  labs(y = "Malformation Type", x = "Cancer Type") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1)))



# Discuss how to adapt to current data ------------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Generating the plot is straightforward.
#' 
#' The data frame that feeds into the plot is also pretty simple.
#' 
#' 4 columns: cancer (as factor), defect (as factor), OR, CI (as factor).
#' 
#' The challenge here will be in programatically generating a table with 
#' this information.
#' 
#' Could likely adjust the loop I'm running to generate a data frame with
#' this info.
#' 
#' Is there any way to recover it from the text files I've already 
#' generated with the model summaries?
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

tmp <- gmodels::CrossTable(armitx.nochrom$conganomalies.cns, armitx.nochrom$all)

or <- as.numeric(paste(format(
                    ((tmp$t[1,1]*tmp$t[2,2])/(tmp$t[2,1]*tmp$t[1,2])), digits = 3)))


tmp2 <-  filter(armitx.nochrom, runif < 0.05)

rm(armitx.nochrom)

model <- glm(all ~ down.syndrome, data = tmp2, family = binomial(link='logit'))

tmp <- exp(cbind(OR=coef(x),confint(x,level=0.95)))
tmp2 <- summary(x)$coefficients
str(tmp2)
print(tmp)
print(tmp2)


tmp <- exp(cbind(OR=coef(x),confint(x,level=0.95)))
estimates <- data.frame(OR = tmp[2,1], ci.lower = tmp[2,2], ci.upper = tmp[2,3])
write.table(estimates, file = 'C:/Users/schraw/Desktop/goback models/cns defects/estimates for cns defects.csv', sep=',', append = TRUE, row.names = FALSE)


for (i in 22){
  
  tmp <- table(armitx.nochrom[,i], armitx.nochrom$cancer1)
  tmp <- which(tmp[2, ] > 5)
  tmp <- tmp + 117
  
  if (length(tmp) > 0){
    
    for (j in tmp){
      
      z <- names(armitx.nochrom[i])
      y <- names(armitx.nochrom[j])
      
#      sink(file = paste0('C:/Users/schraw/Desktop/goback models/cns defects/',as.character(y),'_by_',as.character(z),'.txt'), append = TRUE)
      
#      print(paste('THIS IS THE REGRESSION MODEL FOR', y, 'BY', z))
      
      x <- glm(armitx.nochrom[,j] ~ armitx.nochrom[,i], data = armitx.nochrom, family = binomial(link = 'logit'))
#      x.summary <- (summary(x))
#      x.summary <- exp(cbind(OR=coef(x),confint(x,level=0.95)))
      x.summary <- summary(x)$coefficients
      estimates <- data.frame(defect = z, cancer = y, OR = exp(x.summary[2,1]), 
                              ci.lower = exp(x.summary[2,1]-(1.96*x.summary[2,2])), 
                              ci.upper = exp(x.summary[2,1]+(1.96*x.summary[2,2])))
      write.table(estimates, file = 'C:/Users/schraw/Desktop/goback models/cns defects/estimates for cns defects.csv', sep=',', append = TRUE, 
                  row.names = FALSE, col.names = FALSE)
#      sink()
    }
    
  }
  
  else{
    
    sink(file = 'C:/Users/schraw/Desktop/goback models/list of untested defects/list of defects with no models.txt')
    
    print(paste('There were no cancers with 5 or more comorbid cases for',names(armitx.nochrom[i])))
    
    sink()
  }
  
}

