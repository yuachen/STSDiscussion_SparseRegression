# The supporting simulation code for the Statistical Science discussion paper
# on Sparse regression: 
# Scalable algorithms and empirical performance 
# by Bertsimas, Pauphilet and van Parys (hereafter BPvP) 
# and Best Subset, Forward Stepwise, or Lasso? Analysis and recommendations based on extensive comparisons 
# by Hastie, Tibshirani and Tibshirani (here- after HTT)

# The code is written by Yuansi Chen, Armeen Taeb and Peter Bühlmann

# Figure_2 plots the distributional robustness 
# comparison across LASSO, BS, SS models in low/high dimension low/high SNR

setwd('~/ETH/Writing/STSDiscussion_code')

library(bestsubset)
library(ggplot2)
library(grid)
library(ncvreg)

source("utils.R")

library("JuliaCall")
# julia_setup(JULIA_HOME = "/Applications/Julia-1.0.app/Contents/Resources/julia/bin")

set.seed(12345)

# choose settings 1, 2, 3, 4
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
} else if (setting_switch == 3) {
  # setting 3 high dimension high SNR
  # might take more than 4 hours for BS
  n = 50; p = 1000 # Size of training set, and number of predictors
  nval = n # Size of validation set
  nrep = 10 # Number of repetitions
  seed = 0 # Random number generator seed
  s = 5 # Number of nonzero coefficients
  beta.type = 2 # Coefficient type
  snr = 20. # signal to noise ratio
  rho = 0.00 # correlation
} else if (setting_switch == 4) {
  # setting 3 high dimension low SNR
  # might take more than 4 hours for BS
  n = 50; p = 1000 # Size of training set, and number of predictors
  nval = n # Size of validation set
  nrep = 10 # Number of repetitions
  seed = 0 # Random number generator seed
  s = 5 # Number of nonzero coefficients
  beta.type = 2 # Coefficient type
  snr = 1. # signal to noise ratio
  rho = 0.00 # correlation
}


# generate data
xy.obj = sim.xy(n,p,nval,rho=rho,s=s,beta.type=beta.type,snr=snr)
xy.objtest = sim.xy(n,p,nval,rho=rho,s=s,beta.type=beta.type,snr=snr)
x = xy.obj$x
y = xy.obj$y

xval = xy.obj$xval
yval = xy.obj$yval


# solve LASSO
nrel = 1 # number of relaxation, 1 means only lasso is considered
nlam = 100 # number of lambdas
res.las = lasso(x,y,intercept=FALSE,nlam=nlam,nrel=nrel)
# two choices of regularization parameters
regMinValMSE.las = beta_select(res.las$beta, xval, yval, 5, choice=1)
betaMinValMSE.las = regMinValMSE.las$beta
indexMinValMSE.las = regMinValMSE.las$index
regKMinVal.las = beta_select(res.las$beta, xval, yval, 5, choice=2)
betaKMinVal.las = regKMinVal.las$beta
indexKMinVal.las = regKMinVal.las$index


# solve BS
res.bs = bs(x,y,intercept=FALSE, k=0:p, time.limit=1800)
# if in high dimension setting, your program takes too long, try to limit the k choices 
# res.bs = bs(x,y,intercept=FALSE, k=0:10, time.limit=1800)
regMinValMSE.bs = beta_select(res.bs$beta, xval, yval, 5, choice=1)
betaMinValMSE.bs = regMinValMSE.bs$beta
indexMinValMSE.bs = regMinValMSE.bs$index
regKMinVal.bs = beta_select(res.bs$beta, xval, yval, 5, choice=2)
betaKMinVal.bs = regKMinVal.bs$beta
indexKMinVal.bs = regKMinVal.bs$index

# Solve SS from Julia
# load julia package by BPvP
julia_library("SubsetSelection")


betaSSJulia = function(x, y, gamma){
  betaSS = matrix(0,p,p)
  # assign x and y matrix
  julia_assign("x", x)
  julia_assign("y", y)
  
  julia_assign("γ", gamma)
  
  for (k in 1:p) {
    julia_assign("k", k)
    julia_eval("a = SubsetSelection.OLS()")
    julia_eval("result_ss = subsetSelection(a, Constraint(k), y, x, γ=γ, maxIter=200)")
    indexes = julia_eval("result_ss.indices[:]")
    betanozero = julia_eval("result_ss.w")
    betaSS[indexes, k] = betanozero
  }
  return(betaSS)
}

betaSS = betaSSJulia(x, y, 1000.)

regMinValMSE.SS = beta_select(betaSS, xval, yval, 5, choice=1)
betaMinValMSE.SS = regMinValMSE.SS$beta
indexMinValMSE.SS = regMinValMSE.SS$index
regKMinVal.SS = beta_select(betaSS, xval, yval, 5, choice=2)
betaKMinVal.SS = regKMinVal.SS$beta
indexKMinVal.SS = regKMinVal.SS$index

betaSS001 = betaSSJulia(x, y, 0.01)

regMinValMSE.SS001 = beta_select(betaSS001, xval, yval, 5, choice=1)
betaMinValMSE.SS001 = regMinValMSE.SS001$beta
indexMinValMSE.SS001 = regMinValMSE.SS001$index
regKMinVal.SS001 = beta_select(betaSS001, xval, yval, 5, choice=2)
betaKMinVal.SS001 = regKMinVal.SS001$beta
indexKMinVal.SS001 = regKMinVal.SS001$index

perturbs = 1.5^(seq(-3, 2, length.out=30))
drMSE.las = rep(0.0, length(perturbs))
drMSE.lasK = rep(0.0, length(perturbs))
drMSE.bs = rep(0.0, length(perturbs))
drMSE.bsK = rep(0.0, length(perturbs))
drMSE.SS = rep(0.0, length(perturbs))
drMSE.SSK = rep(0.0, length(perturbs))
drMSE.SS001 = rep(0.0, length(perturbs))
drMSE.SS001K = rep(0.0, length(perturbs))

xtest = xy.objtest$x
ytest = xy.objtest$y
for (ind in seq(1, length(perturbs))) {
  drMSE.las[ind] = dr_exact(betaMinValMSE.las, xtest, ytest, perturbs[ind])
  drMSE.lasK[ind] = dr_exact(betaKMinVal.las, xtest, ytest, perturbs[ind])
  drMSE.bs[ind] = dr_exact(betaMinValMSE.bs, xtest, ytest, perturbs[ind])
  drMSE.bsK[ind] = dr_exact(betaKMinVal.bs, xtest, ytest, perturbs[ind])
  drMSE.SS[ind] = dr_exact(betaMinValMSE.SS, xtest, ytest, perturbs[ind])
  drMSE.SSK[ind] = dr_exact(betaKMinVal.SS, xtest, ytest, perturbs[ind])
  drMSE.SS001[ind] = dr_exact(betaMinValMSE.SS001, xtest, ytest, perturbs[ind])
  drMSE.SS001K[ind] = dr_exact(betaKMinVal.SS001, xtest, ytest, perturbs[ind])
}

ppdf.dr =  data.frame(perturbs=perturbs,
                      lasso=drMSE.las,
                      lassoK=drMSE.lasK,
                      bs=drMSE.bs,
                      bsK=drMSE.bsK,
                      SS=drMSE.SS,
                      SSK=drMSE.SSK,
                      SS001=drMSE.SS001,
                      SS001K=drMSE.SS001K
)


ggplot(data=ppdf.dr, aes(x=perturbs)) +
  geom_line(aes(y=lasso, color="LASSO", linetype="LASSO"), size=0.8) +
  geom_line(aes(y=bs, color="BS", linetype="BS"), size=0.8) +
  geom_line(aes(y=SS, color="SS_1", linetype="SS_1"), size=0.8) +
  geom_line(aes(y=SS001, color="SS_2", linetype="SS_2"), size=0.8) +
  theme(legend.position=c(.16,.72), legend.key.size = unit(1,"line"), legend.key.width=unit(1, "cm")) + 
    scale_color_manual(values = alpha(c(
    'LASSO' = 'blue',
    'BS' = 'red',
    'SS_1' = 'darkgreen',
    'SS_2' = 'darkgreen'), 0.6)) + 
  scale_linetype_manual(values = c(
    'LASSO' = 1,
    'BS' = 1,
    'SS_1' = 2,
    'SS_2' = 3)) +
  labs(x = "Perturbation",
       y = "DR",
       color = "Methods", 
       linetype = "Methods") 

if (setting_switch == 1) {
  ggsave("fig/Figure_2_low_high.pdf",width = 10, height = 6.7, unit="cm")
} else if (setting_switch == 2) {
  ggsave("fig/Figure_2_low_low.pdf",width = 10, height = 6.7, unit="cm")
} else if (setting_switch == 3) {
  ggsave("fig/Figure_2_high_high.pdf",width = 10, height = 6.7, unit="cm")
} else if (setting_switch == 4) {
  ggsave("fig/Figure_2_high_low.pdf",width = 10, height = 6.7, unit="cm")
}


