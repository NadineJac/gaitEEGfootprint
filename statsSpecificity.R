# statistical comparisons of the specificity analysis
#
# directory:  */derivates/specificity/group/results
# In: erpZ.txt, topoZ.txt and SNRdB.txt  (output of "compareERP.m")
# Out: "statsSpecificity.xlsx" 
#
# developed in R (v4.0.1) by Nadine Jacobsen, nadine.jacobsen@uol.de
# June 2020, last revision July 1, 2020

#PATH <- file.path("E:", "nadine", "test_footprint_BIDS") # add your path to the BIDS dataset here
PATH = "D:/DATA_D/test_footprint_scripts_PC"

# library
library(effsize) # needed for cohens d, if not installed use: install.packages("effsize")
library(xlsx) #save results as excel file, s.a.

setwd(file.path(PATH, "derivates", "specificity", "group", "results"))

stats_tmp <-  data.frame(Feature   = double(), #setup df to store results
                         Mean      = double(),
                         SD        = double(),
                         shapiro.p = double(),
                         T.value   = double(),
                         df        = double(),
                         p         = double(),
                         p.adj     = double(),
                         d         = double())

COND <- c("erpZ", "topoZ","SNRdB")

for ( c in 1:length(COND)){
  
  # load data
  dat <-read.table(paste(COND[c], ".txt", sep=""),T, sep = ",")
  stats <- stats_tmp
  
  if (c == 3) {#SNR requires dependent samples ttest
    
    dat1 <- dat[,c(1,3,5)]
    dat2 <- dat[,c(2,4,6)]
    datDiff <- dat2-dat1
    
    for (var in c(1:3)){
      
      # check whether data is normal distributed --> OK
      testNorm <- shapiro.test(datDiff[,var])
      if (testNorm$p.value < .05){
        print(paste(colnames(dat[var]), "values of", COND[var], "are not normal distributed"))
      }
      
      stats[var,"shapiro.p"] <- round(testNorm$p.value,3)
      
      # descriptives
      stats[var,"Feature"] <- colnames(dat[var*2])
      stats[var,"Mean"] <- round(mean(datDiff[,var]),2)
      stats[var,"SD"]<- round(sd(datDiff[,var]),2)
      
      # t-test
      # one-sample, two-sided student t-test, mean = 0
      StudentModel <- t.test(dat2[,var],dat1[,var], paired = T)
      stats[var,"T.value"] <- round(StudentModel$statistic,2)
      stats[var,"df"] <- StudentModel$parameter
      stats[var,"p"] <- round(StudentModel$p.value,3)
      
      
      # cohen d
      d <- cohen.d(dat2[,var],dat1[,var], paired = T)
      stats[var,"d"]<- round(d$estimate,2)
    }
    
  } else {
    for (var in c(1:3)){
      
      # check whether data is normal distributed --> OK
      testNorm <- shapiro.test(dat[,var])
      if (testNorm$p.value < .05){
        print(c(colnames(dat[var]), "values of", COND[var], "are not normal distributed"))
      }
      
      stats[var,"shapiro.p"] <- round(testNorm$p.value,3)
      
      # descriptives
      stats[var,"Feature"] <- colnames(dat[var])
      stats[var,"Mean"] <- round(mean(dat[,var]),2)
      stats[var,"SD"]<- round(sd(dat[,var]),2)
      
      # t-test
      # one-sample, two-sided student t-test, mean = 0
      StudentModel <- t.test(dat[,var])
      stats[var,"T.value"] <- round(StudentModel$statistic,2)
      stats[var,"df"] <- StudentModel$parameter
      stats[var,"p"] <- round(StudentModel$p.value,3)
      
      # cohen d
      d <- cohen.d(dat[,var],NA)
      stats[var,"d"]<- round(d$estimate,2)
      
    }
  }
  
  # adjust for multiple comparisons
  stats$p.adj <- round(p.adjust(stats$p, method = "holm", n = 9), 3)
  
  write.xlsx(stats, "statsSpecificity.xlsx", sheetName = COND[c], append = T)
  
}

rm(list = ls())
