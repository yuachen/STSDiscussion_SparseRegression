# The supporting simulation code for the Statistical Science discussion paper
# on Sparse regression: 
# Scalable algorithms and empirical performance 
# by Bertsimas, Pauphilet and van Parys (hereafter BPvP) 
# and Best Subset, Forward Stepwise, or Lasso? Analysis and recommendations based on extensive comparisons 
# by Hastie, Tibshirani and Tibshirani (here- after HTT)

# The code is written by Yuansi Chen, Armeen Taeb and Peter Bühlmann

# Figure_4 plots the feature selection stability of LASSO and BS
# as a function of regularization parameter choices

setwd('~/ETH/Writing/STSDiscussion_code')

library(bestsubset)
library(ggplot2)
library(grid)
library(ncvreg)

source("utils.R")

set.seed(12345)

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

# repeats of the same i.i.d. experiments
repeats = 100

selectfreqAll.las = matrix(0, nrow = p, ncol=100)
selectfreqAll.bs = matrix(0, nrow = p, ncol=31)

lambda_grid=10^seq(2,-3,length=100) ##get lambda sequence
loglamMinValMSE_avg = 0
loglamKMinVal_avg = 0
kMinValMSE_avg = 0
kKMinVal_avg = 0

for (rep in 1:repeats){
  print(rep)
  # generate data
  xy.obj = sim.xy(n,p,nval,rho=rho,s=s,beta.type=beta.type,snr=snr)
  
  x = xy.obj$x
  y = xy.obj$y
  
  xval = xy.obj$xval
  yval = xy.obj$yval
  
  # LASSO
  nrel = 1 # number of relaxation, 1 means only lasso is considered
  nlam = 100 # number of lambdas
  res.las = lasso(x,y,intercept=FALSE,nlam=nlam,nrel=nrel,lambda = lambda_grid)
  
  regMinValMSE.las = beta_select(res.las$beta, xval, yval, 5, choice=1)
  betaMinValMSE.las = regMinValMSE.las$beta
  indexMinValMSE.las = regMinValMSE.las$index
  regKMinVal.las = beta_select(res.las$beta, xval, yval, 5, choice=2)
  betaKMinVal.las = regKMinVal.las$beta
  indexKMinVal.las = regKMinVal.las$index
  
  loglamMinValMSE_avg = loglamMinValMSE_avg + log(res.las$lambda[indexMinValMSE.las])/repeats
  loglamKMinVal_avg = loglamKMinVal_avg + log(res.las$lambda[indexKMinVal.las])/repeats
  # select freq lasso 
  selectfreqAll.las = selectfreqAll.las + as.matrix(res.las$beta != 0)
  
  # BS
  res.bs = bs(x,y,intercept=FALSE)
  
  regMinValMSE.bs = beta_select(res.bs$beta, xval, yval, 5, choice=1)
  betaMinValMSE.bs = regMinValMSE.bs$beta
  indexMinValMSE.bs = regMinValMSE.bs$index
  regKMinVal.bs = beta_select(res.bs$beta, xval, yval, 5, choice=2)
  betaKMinVal.bs = regKMinVal.bs$beta
  indexKMinVal.bs = regKMinVal.bs$index
  
  kMinValMSE_avg = kMinValMSE_avg + res.bs$k[indexMinValMSE.bs]/repeats
  kKMinVal_avg = kKMinVal_avg + res.bs$k[indexKMinVal.bs]/repeats
  # selectfreq bs
  selectfreqAll.bs = selectfreqAll.bs + as.matrix(res.bs$beta != 0)
  
}

ppdf.stabreg.las =  data.frame(lambda=log(lambda_grid),
                               selectfreq1=selectfreqAll.las[1, ],
                               selectfreq6=selectfreqAll.las[6, ]
)

ggplot(ppdf.stabreg.las, aes(x=lambda)) +
  geom_line(aes(y=selectfreq1, color="x1", linetype="x1"), size=1) + 
  geom_line(aes(y=selectfreq6, color="x6", linetype="x6"), size=1) + 
  geom_vline(xintercept = loglamMinValMSE_avg, 
             linetype="dotted", 
             color = "magenta") + 
  geom_text(aes(x=loglamMinValMSE_avg, y= 75, label = "min val error"), colour = "magenta", vjust = 1, size = 4, angle=90) + 
  geom_vline(xintercept = loglamKMinVal_avg, 
             linetype="dotted", 
             color = "orange") + 
  geom_text(aes(x=loglamKMinVal_avg, y= 50, label = "known s=5"), colour = "orange", vjust = 1, size = 4, angle=90) + 
  theme(legend.position=c(.8,.2), legend.key.size = unit(1,"line"), legend.key.width=unit(1.5, "cm")) + 
  scale_color_manual(values = alpha(c(
    'x1' = 'darkblue',
    'x6' = 'blue'), 0.6)) + 
  scale_linetype_manual(values = c(
    'x1' = 1,
    'x6' = 2)) +
  labs(x="log(lambda)", y="feature selection frequency (%)", color="feature index", 
       linetype="feature index") + 
  scale_x_reverse()


ggsave("fig/Figure_4_las.pdf",width = 10, height = 6.7, unit="cm")

ppdf.stabreg.bs =  data.frame(k=0:30,
                              selectfreq1=selectfreqAll.bs[1, ],
                              selectfreq6=selectfreqAll.bs[6, ]
)

ggplot(ppdf.stabreg.bs, aes(x=k)) +
  geom_line(aes(y=selectfreq1, color="x1", linetype="x1"), size=1) + 
  geom_line(aes(y=selectfreq6, color="x6", linetype="x6"), size=1) + 
  geom_vline(xintercept = kMinValMSE_avg, 
             linetype="dotted", 
             color = "magenta") + 
  geom_text(aes(x=kMinValMSE_avg, y= 75, label = "min val error"), colour = "magenta", vjust = 1, size = 4, angle=90) + 
  geom_vline(xintercept = kKMinVal_avg, 
             linetype="dotted", 
             color = "orange") + 
  geom_text(aes(x=kKMinVal_avg, y= 30, label = "known s=5"), colour = "orange", vjust = 1, size = 4, angle=90) + 
  theme(legend.position=c(.8,.2), legend.key.size = unit(1,"line"), legend.key.width=unit(1.5, "cm")) + 
  scale_color_manual(values = alpha(c(
    'x1' = 'darkred',
    'x6' = 'red'), 0.6)) + 
  scale_linetype_manual(values = c(
    'x1' = 1,
    'x6' = 2)) +
  labs(x="k", y="feature selection frequency (%)", color="feature index", linetype="feature index")

ggsave("fig/Figure_4_bs.pdf",width = 10, height = 6.7, unit="cm")
