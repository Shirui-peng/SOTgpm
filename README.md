# SOTgpm.jl
The code in this package processes acoustic and in situ data to infer temperature change in the ocean using Gaussian process modeling. The acoustic data are corrected travel time changes between repeating sound waves, and functions to obatin these pairs are variants of those in [SOT](https://github.com/joernc/SOT). The principal functions are:
1. Maximum likelihood estimator to determine parameters of stochastic prior covariance for in situ data
2. Integrated covariances among pointwise, path-mean, and area-mean temperature anomalies
3. Invert pointwise, path-mean, or area-mean temperature anomalies and trends with uncertainties
