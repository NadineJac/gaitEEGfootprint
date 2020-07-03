# statistical comparisons of the footprint sensitivity analysis
#
# directory:  */derivates/footprint/group/results
#
# developed in R (v4.0.1) by Nadine Jacobsen, nadine.jacobsen@uol.de
# June 2020, last revision July 1, 2020

#PATH = file.path("E:", "nadine", "test_footprint_BIDS") # add your path to the BIDS dataset here
PATH = "D:/DATA_D/test_footprint_scripts_PC"

# libraries
library(effsize) # needed for cohens d, if not installed use: install.packages("effsize")
library(xlsx) #save results as excel file, s.a.

setwd(file.path(PATH, "derivates", "footprint", "group", "results"))


# footprint distances -----------------------------------------------------
# In: footprintDistances.txt (output of "calcFootprintDistances.m")
# Out: "statsSensitvity.xlsx" sheet "footprintDistances"


# load data
distFoot <- read.csv("footprintDistances.txt", T)

# data frame for storing results
statsDistFoot <-  data.frame(
  Comp      = double(),
  p.shapiro = double(),
  M         = double(),
  SD        = double(),
  T.stat    = double(),
  df        = integer(),
  p         = double(),
  p.adj     = double(),
  d         = double()
)

for (var in 3:4) {
  # col 2 (raw2ASr not reported in manuscript)
  
  # Descriptives
  statsDistFoot[var - 1, "Comp"] <- colnames(distFoot[var])
  statsDistFoot[var - 1, "M"] <- round(mean(distFoot[, var]), 2)
  statsDistFoot[var - 1, "SD"] <- round(sd(distFoot[, var]), 2)
  
  # assess normality of distances
  # tmp <- shapiro.test(distFoot[,var])
  # statsDistFoot[var-1,"Comp"]<-round(tmp$p.value,3)
  
  
  # perform one-sample t-test
  StudentModel <-
    t.test(distFoot[, var], mu = 0, alternative = "greater")
  
  
  statsDistFoot[var - 1, "T.stat"] <-
    round(StudentModel$statistic, 2)
  statsDistFoot[var - 1, "df"] <- StudentModel$parameter
  statsDistFoot[var - 1, "p"] <- round(StudentModel$p.value, 3)
  
  # correct for 2 comparisons (only raw2ICA and ASR2ICA reported)
  statsDistFoot[var - 1, "p.adj"] <-
    round(p.adjust(StudentModel$p.value, method = "holm", n = 2), 3)
  
  # calculate effect size
  CohenD <- cohen.d(distFoot[, var], NA)
  statsDistFoot[var - 1, "d"] <- round(CohenD$estimate, 2)
  
}
# add note
statsDistFoot[var, "Comp"] <-
  "Note. p-value adjusted for 2 comparisons (Bonferroni-Holm)"

# save results
write.xlsx(statsDistFoot,
           "statsSensitivity.xlsx",
           sheetName = "footprintDistances",
           append = T)




# single feature distances ------------------------------------------------
# In: gait_footprint_before/ _after/ _afterASR.txt(output of "calculateFootprint.m")
# Out: several sheets in "statsSensitivity.xlsx" 

## descriptives ________________________
# Out: sheets "features befor", "fetures after", "feautures afterASR"
COND <- c("before", "afterASR", "after")
# set up dfs
desFeatures <- data.frame(
  Feature = double(),
  Mdn     = double(),
  M       = double(),
  SD      = double()
)

for (c in 1:3) {
  # did not report comarison raw2ASR
  
  # load data
  dat1 <-
    read.csv(paste("gait_footprint_", COND[c], ".txt", sep = ""), T)
  
  for (var in 2:8) {
    # descriptives
    # save descriptors
    desFeatures[var-1, "Feature"] <- LETTERS[var-1]
    desFeatures[var-1, "Mdn"]     <- round(median(dat1[, var]), 2)
    desFeatures[var-1, "M"]       <- round(mean(dat1[, var]), 2)
    desFeatures[var-1, "SD"]      <- round(sd(dat1[, var]), 2)
  }
  write.xlsx(
    desFeatures,
    "statsSensitivity.xlsx",
    sheetName = paste("features", COND[c]),
    append = T
  )
}

## stats_______________________________________________________
# Out: sheets "features raw2ICA", "features ASR2ICA"

COND1 <- c("before", "afterASR", "before")
COND2 <- c("after", "after", "afterASR")
CONDname <- c("raw2ICA", "ASR2ICA", "raw2ASR")

statsFeatures <-  data.frame(
  Feature   = double(),
  M.diff    = double(),
  SD.diff   = double(),
  p.shapiro = double(),
  test.stat = double(),
  p         = double(),
  p.adj     = double(),
  eff.size  = double()
)


## since differences are not normal distributed, use wicox signed rank (dep. samples)
for (c in 1:2) {
  # did not report comarison raw2ASR
  
  # load data
  dat1 <-
    read.csv(paste("gait_footprint_", COND1[c], ".txt", sep = ""), T)
  dat2 <-
    read.csv(paste("gait_footprint_", COND2[c], ".txt", sep = ""), T)
  
  
  for (var in 2:8) {
    
    # comparisons
    statsFeatures[var-1, "Feature"] <- LETTERS[var-1]
    statsFeatures[var-1, "M.diff"] <-
      round(mean(dat2[, var] - dat1[, var]), 2)
    statsFeatures[var-1, "SD.diff"] <-
      round(sd(dat2[, var] - dat1[, var]), 2)
    
    # assess normality of differences w shapiro-wilk
    tmp <- shapiro.test(dat1[, var] - dat2[, var])
    statsFeatures[var-1, "p.shapiro"] <- round(tmp$p.value, 3)
    
    if (tmp$p.value < 0.05) {
      # wilcoxon signed rank test, dependent samples
      wilcoxModel <- wilcox.test(dat1[, var], dat2[, var], paired = TRUE)
      Z <- qnorm(wilcoxModel$p.value / 2)
      r <- Z / sqrt(nrow(dat1))
      
      statsFeatures[var-1, "test.stat"] <-
        round(wilcoxModel$statistic, 2)
      statsFeatures[var-1, "p"] <- round(wilcoxModel$p.value, 3)
      statsFeatures[var-1, "eff.size"] <- round(r, 2)
      
    } else {
      # dependent two-sided dependent samples student t-test
      StudentModel <- t.test(dat2[, var], dat1[, var], paired = TRUE)
      
      statsFeatures[var-1, "test.stat"] <-
        round(StudentModel$statistic, 2)
      statsFeatures[var-1, "p"] <- round(StudentModel$p.value, 3)
      
      # cohen d
      d <- cohen.d(dat2[, var], dat1[, var], paired = TRUE)
      statsFeatures[var-1, "eff.size"] <- round(d$estimate, 2)
    }
  }
  # correct for multiple comparisons (n014 since we will report  two comparisons of 7 features each)
  statsFeatures$p.adj <-
    p.adjust(statsFeatures$p, method = "holm", n = 14)
  
  # add note
  statsFeatures[var, "Feature"] <-
    "Note. two-sided, dependent samples t-test, effect size: Cohen's d, if p-shapiro<.05: dependent samples, Wilcoxon signed-rank test, Effect size: R. p-value adjusted for 14 comparisons (Bonferroni-Holm)"
  
  # save results
  write.xlsx(
    statsFeatures,
    "statsSensitivity.xlsx",
    sheetName = paste("features", CONDname[c]),
    append = T
  )
  
}

rm(list = ls())