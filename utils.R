# The supporting simulation code for the Statistical Science discussion paper
# on Sparse regression: 
# Scalable algorithms and empirical performance 
# by Bertsimas, Pauphilet and van Parys (hereafter BPvP) 
# and Best Subset, Forward Stepwise, or Lasso? Analysis and recommendations based on extensive comparisons 
# by Hastie, Tibshirani and Tibshirani (here- after HTT)

# The code is written by Yuansi Chen, Armeen Taeb and Peter Bühlmann

# utils provide some utility functions for generating all Figures

beta_select = function(beta, xval, yval, slimit, choice=1) {
  if (choice == 1) {
    # choose the hyperparameter based on the lowest validation error
    ypredval =  xval %*% beta
    valMSE = colMeans((ypredval-yval)^2)
    indexMinValMSE = which.min(valMSE)
    betaMinValMSE = beta[, indexMinValMSE]
    res <- list()
    res$beta = betaMinValMSE
    res$index = indexMinValMSE
    return(res)
  } else if (choice == 2) {
    # choose the hyperparameter with beta size smaller or equal to slimit
    # then with the lowest validation error
    ypredval =  xval %*% beta
    valMSE = colMeans((ypredval-yval)^2)
    indexK = which(colSums(beta != 0) <= slimit)
    indexKMinVal = indexK[which.min(valMSE[indexK])]
    betaK = beta[, indexKMinVal]
    res <- list()
    res$beta = betaK
    res$index = indexKMinVal
    return(res)
  } else {
    print("Invalid choice parameter")
  }
}

# functions to compute the distributional robust difference (DRD)
drd_exact = function(beta, x, y, perturb=1, nrepperb = 20) {
  # the easy implementation using l1 norm of beta
  n = nrow(x)
  drdperb = perturb* colSums(abs(beta))
  return(drdperb)
}
drd_gradient = function(beta, x, y, perturb=1, nrepperb = 20) {
  # the implementation with projected gradient descent
  # we verify that this is the same as outputing eta * |beta|_1
  n = nrow(x)
  num_final_lam = ncol(beta)
  ypred = x %*% beta
  drdperb = colMeans((ypred-y)^2)
  
  for (j in 1:num_final_lam){
    delta = matrix(rnorm(n*p),n,p)
    perturb = 1
    delta = scale(delta, center=F)/sqrt(n-1) * perturb
    for (rep in 1:nrepperb) {
      xperb = x + delta
      ypredperb = xperb %*% as.vector(beta[, j])
      delta = delta - as.vector(y-ypredperb) %o% as.vector(beta[ , j]) 
      delta = scale(delta, center=F)/sqrt(n-1) * perturb
    }
    drdperb[j] = sqrt(sum((ypredperb-y)^2)) - sqrt(sum((ypred[, j]-y)^2))
  }
  
  return(drdperb)
}

# functions to compute the distributional robust (DR)
dr_exact = function(beta, x, y, perturb=1) {
  # the easy implementation using l1 norm of beta
  n = nrow(x)
  drperb = perturb* sum(abs(beta)) +  sqrt(sum((x %*% beta-y)^2))
  return(drperb)
}