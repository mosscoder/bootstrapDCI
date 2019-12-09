library(tidyverse)
source('https://raw.githubusercontent.com/mosscoder/bootstrapDCI/master/bootstrapDCI.R')

set.seed(123)
control <- runif(30, min = 3, max = 5)
treatment1 <- runif(30, min = 3, max = 6) #fake data that won't be different than control
treatment2 <- runif(30, min = 4, max = 6) #fake data that will be different than control

#First calculate Hedges' G (unbiased Cohen's D)
effSize1 <- Unbiased.d(treatment1, control)
effSize2 <- Unbiased.d(treatment2, control)

#Calculate a 95% CI for the effect size, 1000 samples
effSizeCI1 <- bootstrapDCI(treatment1, control, B = 1000, alpha = 0.05)
effSizeCI2 <- bootstrapDCI(treatment2, control, B = 1000, alpha = 0.05)

#The output you want is the bias corrected accelerated row (BCA, row 2)
#LCL = lower limit, UCL = upper limit
dat <- data.frame(treatment = c('Treat1','Treat2'),
                  effSize = c(effSize1, effSize2),
                  lcl = c(effSizeCI1$lcl[2],effSizeCI2$lcl[2]),
                  ucl = c(effSizeCI1$ucl[2],effSizeCI2$ucl[2]))

ggplot(data = dat, aes(x = treatment, y = effSize, ymin = lcl, ymax = ucl)) +
  geom_point() +
  geom_linerange() +
  geom_hline(aes(yintercept = 0), linetype = 'dashed') + 
  xlab('Treatment Type') +
  ylab('Effect Size')