# The supporting simulation code for the Statistical Science discussion paper
# on Sparse regression: 
# Scalable algorithms and empirical performance 
# by Bertsimas, Pauphilet and van Parys (hereafter BPvP) 
# and Best Subset, Forward Stepwise, or Lasso? Analysis and recommendations based on extensive comparisons 
# by Hastie, Tibshirani and Tibshirani (here- after HTT)

# The code is written by Yuansi Chen, Armeen Taeb and Peter Bühlmann

# Figure_3 plots the feature selection stability of LASSO and BS

setwd('~/ETH/Writing/STSDiscussion_code')

library(bestsubset)
library(ggplot2)
library(grid)
library(ncvreg)

source("utils.R")

set.seed(12345)

# 1. Plot distributional robustness as a function of hyper parameters

# choose settings 1, 2
setting_switch = 1

if (setting_switch == 1) {
  # Set some overall simulation parameters
  # setting 1 low dimension high SNR
  n = 100; p = 30 # Size of training set, and number of predictors
  nval = n # Size of validation set
  nrep = 10 # Number of repetitions
  seed = 0 # Random number generator seed
  s = 5 # Number of nonzero coefficients
  beta.type = 2 # Coefficient type
  snr = 2.0 # signal to noise ratio
  rho = 0.00 # correlation
} else if (setting_switch == 2) {
  # setting 2 low dimension low SNR
  n = 100; p = 30 # Size of training set, and number of predictors
  nval = n # Size of validation set
  nrep = 10 # Number of repetitions
  seed = 0 # Random number generator seed
  s = 5 # Number of nonzero coefficients
  beta.type = 2 # Coefficient type
  snr = 0.1 # signal to noise ratio
  rho = 0.00 # correlation
}

# repeats of the same i.i.d. experiments
repeats = 100

selectfreqMinValMSE.las = rep(0.0, p)
selectfreqKMinVal.las = rep(0.0, p)
selectfreqMinValMSE.bs = rep(0.0, p)
selectfreqKMinVal.bs = rep(0.0, p)
for (rep in 1:repeats){
  print(rep)
  # generate data
  xy.obj = sim.xy(n,p,nval,rho=rho,s=s,beta.type=beta.type,snr=snr)
  
  x = xy.obj$x
  y = xy.obj$y
  
  xval = xy.obj$xval
  yval = xy.obj$yval
  
  
  nrel = 1 # number of relaxation, 1 means only lasso is considered
  nlam = 100 # number of lambdas
  res.las = lasso(x,y,intercept=FALSE,nlam=nlam,nrel=nrel)
  
  regMinValMSE.las = beta_select(res.las$beta, xval, yval, 5, choice=1)
  betaMinValMSE.las = regMinValMSE.las$beta
  indexMinValMSE.las = regMinValMSE.las$index
  regKMinVal.las = beta_select(res.las$beta, xval, yval, 5, choice=2)
  betaKMinVal.las = regKMinVal.las$beta
  indexKMinVal.las = regKMinVal.las$index
  
  selectfreqMinValMSE.las = selectfreqMinValMSE.las + (betaMinValMSE.las != 0)
  selectfreqKMinVal.las = selectfreqKMinVal.las + (betaKMinVal.las != 0)
  
  res.bs = bs(x,y,intercept=FALSE)
  
  regMinValMSE.bs = beta_select(res.bs$beta, xval, yval, 5, choice=1)
  betaMinValMSE.bs = regMinValMSE.bs$beta
  indexMinValMSE.bs = regMinValMSE.bs$index
  regKMinVal.bs = beta_select(res.bs$beta, xval, yval, 5, choice=2)
  betaKMinVal.bs = regKMinVal.bs$beta
  indexKMinVal.bs = regKMinVal.bs$index
  selectfreqMinValMSE.bs = selectfreqMinValMSE.bs + (betaMinValMSE.bs != 0)
  selectfreqKMinVal.bs = selectfreqKMinVal.bs + (betaKMinVal.bs != 0)
  
}

ppdf.stabMinValMSE =  data.frame(feat_num=1:p,
                                 LASSO=selectfreqMinValMSE.las,
                                 BS=selectfreqMinValMSE.bs
)
ppdf.stabKMinVal =  data.frame(feat_num=1:p,
                               LASSO=selectfreqKMinVal.las,
                               BS=selectfreqKMinVal.bs
)

library(tidyverse)
ppdf.stabMinValMSE %>%
  gather(key, value, -feat_num) %>% 
  ggplot(aes(x=feat_num, y=value, fill = key)) + 
  geom_col(position = "dodge") + 
  theme(legend.position=c(.80,.72), legend.key.size = unit(1,"line"), legend.key.width=unit(1, "cm")) + 
  scale_fill_manual(values = alpha(c(
    'LASSO' = 'blue',
    'BS' = 'red'), 0.6)) + 
  ylim(0, 100) + 
  labs(x="feature index",
       y="feature selection frequency (%)",
       fill="Methods")

if (setting_switch == 1) {
  ggsave("fig/Figure_3_lowhigh.pdf", width = 10, height = 6.7, unit="cm")
} else if (setting_switch == 2) {
  ggsave("fig/Figure_3_lowlow.pdf", width = 10, height = 6.7, unit="cm")
}

ppdf.stabKMinVal %>%
  gather(key, value, -feat_num) %>% 
  ggplot(aes(x=feat_num, y=value, fill = key)) + 
  geom_col(position = "dodge") + 
  theme(legend.position=c(.80,.72), legend.key.size = unit(1,"line"), legend.key.width=unit(1, "cm")) + 
  scale_fill_manual(values = alpha(c(
    'LASSO' = 'blue',
    'BS' = 'red'), 0.6)) + 
  ylim(0, 100) + 
  labs(x="feature index",
       y="feature selection frequency (%)",
       fill="Methods")

if (setting_switch == 1) {
  ggsave("fig/Figure_5_lowhigh.pdf", width = 10, height = 6.7, unit="cm")
} else if (setting_switch == 2) {
  ggsave("fig/Figure_5_lowlow.pdf", width = 10, height = 6.7, unit="cm")
}

