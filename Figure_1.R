# The supporting simulation code for the Statistical Science discussion paper
# on Sparse regression: 
# Scalable algorithms and empirical performance 
# by Bertsimas, Pauphilet and van Parys (hereafter BPvP) 
# and Best Subset, Forward Stepwise, or Lasso? Analysis and recommendations based on extensive comparisons 
# by Hastie, Tibshirani and Tibshirani (here- after HTT)

# The code is written by Yuansi Chen, Armeen Taeb and Peter Bühlmann

# Figure_1 plots the distributional robustness difference (DRD)
# as a function of regularization parameters for both the l1 and l0 approach.

setwd('~/ETH/Writing/STSDiscussion_code')

library(bestsubset)
library(ggplot2)
library(grid)
library(ncvreg)
library(latex2exp)

source("utils.R")

set.seed(12345)

# 1. Plot distributional robustness as a function of hyper parameters

# Set some overall simulation parameters
# setting 1
n = 100; p = 30 # Size of training set, and number of predictors
nval = n # Size of validation set
nrep = 10 # Number of repetitions
seed = 0 # Random number generator seed
s = 5 # Number of nonzero coefficients
beta.type = 2 # Coefficient type
snr = 2.0 # signal to noise ratio
rho = 0.00 # correlation

# generate data
xy.obj = sim.xy(n,p,nval,rho=rho,s=s,beta.type=beta.type,snr=snr)
xy.objtest = sim.xy(n,p,nval,rho=rho,s=s,beta.type=beta.type,snr=snr)
x = xy.obj$x
y = xy.obj$y
# get validation data of the same size
xval = xy.obj$xval
yval = xy.obj$yval

nrel = 1 # number of relaxation, 1 means only lasso is considered
nlam = 50 # number of lambdas
res.las = lasso(x,y,intercept=FALSE,nlam=nlam,nrel=nrel)

# two choices of regularization parameters
regMinValMSE.las = beta_select(res.las$beta, xval, yval, 5, choice=1)
betaMinValMSE.las = regMinValMSE.las$beta
indexMinValMSE.las = regMinValMSE.las$index
regKMinVal.las = beta_select(res.las$beta, xval, yval, 5, choice=2)
betaKMinVal.las = regKMinVal.las$beta
indexKMinVal.las = regKMinVal.las$index
# get distributional robust mse with fixed perturbation using exact method or gradient methods 
drd.las = drd_exact(res.las$beta, x, y, perturb=1, nrepperb=20)
# drd.las = drd_gradient(res.las$beta, x, y, perturb=1, nrepperb=20)

ppdf.las = data.frame(lambda=res.las$lambda,
                      loglambda=log(res.las$lambda),
                      l1normbeta=colSums(abs(res.las$beta)),
                      l0normbeta=colSums(res.las$beta != 0),
                      drd=drd.las
                     )

ggplot(ppdf.las,aes_string(x="loglambda")) +
  geom_line(aes_string(y="drd"), color="blue") + 
  geom_vline(xintercept = log(res.las$lambda)[indexMinValMSE.las], 
             linetype="dotted", 
             color = "magenta") + 
  geom_text(aes(x=log(res.las$lambda)[indexMinValMSE.las], y= 2, label = "min val error"), colour = "magenta", vjust = 1, size = 4, angle=90) + 
  geom_vline(xintercept = log(res.las$lambda)[indexKMinVal.las], 
             linetype="dotted", 
             color = "orange") + 
  geom_text(aes(x=log(res.las$lambda)[indexKMinVal.las], y= 8, label = "known s=5"), colour = "orange", vjust = 1, size = 4, angle=90) + 
  geom_point(aes_string(y="drd"), size=1) + theme_minimal() +
  labs(y=TeX("DRD == l1 norm of beta"), x="log(lambda)") + 
  scale_x_reverse()

ggsave("fig/Figure_1_left.pdf", width = 10, height = 6.7, unit="cm")

# now do the same with BS
res.bs = bs(x,y,intercept=FALSE)

regMinValMSE.bs = beta_select(res.bs$beta, xval, yval, 5, choice=1)
betaMinValMSE.bs = regMinValMSE.bs$beta
indexMinValMSE.bs = regMinValMSE.bs$index
regKMinVal.bs = beta_select(res.bs$beta, xval, yval, 5, choice=2)
betaKMinVal.bs = regKMinVal.bs$beta
indexKMinVal.bs = regKMinVal.bs$index
# get distributional robust mse with fixed perturbation using exact method or gradient methods 
drd.bs = drd_exact(res.bs$beta, x, y, perturb=1, nrepperb=20)
# drd.bs = drd_gradient(res.bs$beta, x, y, perturb=1, nrepperb=20)

ppdf.bs = data.frame(k=res.bs$k,
                     l1normbeta=colSums(abs(res.bs$beta)),
                     l0normbeta=colSums(res.bs$beta != 0),
                     drd=drd.bs)

ggplot(ppdf.bs,aes_string(x="k")) +
  geom_line(aes_string(y="drd"), color="red") + 
  geom_vline(xintercept = res.bs$k[indexMinValMSE.bs], 
             linetype="dotted", 
             color = "magenta") + 
  geom_text(aes(x=res.bs$k[indexMinValMSE.bs], y= 2, label = "min val error"), colour = "magenta", vjust = 1, size = 4, angle=90) + 
  geom_vline(xintercept = res.bs$k[indexKMinVal.bs], 
             linetype="dotted", 
             color = "orange") + 
  geom_text(aes(x=res.bs$k[indexKMinVal.bs], y= 8, label = "known s=5"), colour = "orange", vjust = 1, size = 4, angle=90) + 
  geom_point(aes_string(y="drd"), size=1) + theme_minimal() +
  labs(y=TeX("DRD == l1 norm of beta"), x="k") 


ggsave("fig/Figure_1_right.pdf",width = 10, height = 6.7, unit="cm")
