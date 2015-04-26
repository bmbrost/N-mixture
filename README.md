# Nmixture
MCMC algorithms for N-mixture model implementation

- base.N.mixture.mcmc.R:  Fit N-mixture model with varying detection probability (p) and Poisson intensity (lambda) between sites, though neither are modeled as a function of covariates.
- base.N.mixture.mcmc.R: Simulate data according to model specification of base.N.mixture.mcmc.R
- N.mixture.mcmc.R: Fit N-mixture model with varying detection probability (p) and Poisson intensity (lambda) between sites. Detection probability is modeled as a function of covariates (W), but lambda is not modeled as a function of covariates. 
- N.mixture.sim.R: Simulate data according to model specification of N.mixture.mcmc.R
- N.mixture.trend.mcmc.R: Fit N-mixture model with temporal trend. Detection probability (p) is modeled as a function of covariates (W), and N[t] is modeled as a function of N[t-1].
- N.mixture.trend.sim.R:
