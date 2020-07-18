# STSDiscussion_SparseRegression
The supporting simulation code for the Statistical Science discussion paper on 

Sparse regression: Scalable algorithms and empirical performance
 by Bertsimas, Pauphilet and van Parys (hereafter BPvP) and 

Best Subset, Forward Stepwise, or Lasso? Analysis and recommendations based on extensive comparisons
by Hastie, Tibshirani and Tibshirani (here- after HTT).

The code is written by Yuansi Chen, Armeen Taeb and Peter Bühlmann

## User Guide 
* Lasso is solved via the R package [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html) and best subset (BS) is solved via the R package [bestsubset](https://github.com/ryantibs/best-subset), with maximum computing time limit of 30 minutes following HTT. The convex relaxation of the l0 + l2 approach SS (introduced by BPvP) is solved via the Julia package [SubsetSelection](https://github.com/jeanpauphilet/SubsetSelection.jl).

* utils.R contains utility code such as regularization parameter choices for all Figures

* Figure_1.R  plots the distributional robustness difference (DRD) as a function of regularization parameters for both the l1 and l0 approach in Figure 1 of the discussion paper.
* Figure_2.R plots the distributional robustness comparison across Lasso, BS, SS models in low/high dimension low/high SNR  in Figure 2 of the discussion paper (It requires Julia to be installed).
* Figure_3.R plots the feature selection stability of LASSO and BS  in Figure 3 and Figure 5 of the discussion paper.
* Figure_4.R plots the feature selection stability of LASSO and BS as a function of regularization parameter choices in Figure 4 of the discussion paper.

## Citation

If the code helps your research, please cite the discussion paper to be posted.