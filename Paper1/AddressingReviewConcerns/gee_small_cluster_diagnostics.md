# GEE Small-Cluster Diagnostics

The implementation uses `statsmodels` Gaussian GEE with participant as the group, `Autoregressive(grid=True)` working correlation, and the package default robust sandwich covariance. The robust covariance is asymptotic in the number of clusters; N=22 is modest, so p-values and confidence intervals should be interpreted cautiously. No finite-cluster correction was applied because a validated correction is not directly implemented in this workflow. Mixed-effects and 60/150/300-s thinning analyses are reported as robustness checks.
