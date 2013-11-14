% THE DIAGNOSTIC TOOLS (in the diag-folder):
%
% Covergence diagnostics
%   PSRF     - Potential Scale Reduction Factor
%   CPSRF    - Cumulative Potential Scale Reduction Factor
%   MPSRF    - Multivariate Potential Scale Reduction Factor
%   CMPSRF   - Cumulative Multivariate Potential Scale Reduction Factor
%   IPSRF    - Interval-based Potential Scale Reduction Factor
%   CIPSRF   - Cumulative Interval-based Potential Scale Reduction Factor
%   KSSTAT   - Kolmogorov-Smirnov goodness-of-fit hypothesis test
%   HAIR     - Brooks' hairiness convergence diagnostic
%   CUSUM    - Yu-Mykland convergence diagnostic for MCMC
%   SCORE    - Calculate score-function convergence diagnostic
%   GBINIT   - Initial iterations for Gibbs iteration diagnostic
%   GBITER   - Estimate number of additional Gibbs iterations
%
% Time series analysis
%   ACORR      - Estimate autocorrelation function of time series
%   ACORRTIME  - Estimate autocorrelation evolution of time series (simple)
%   GEYER_ICSE - Compute autocorrelation time tau using Geyer's
%                initial convex sequence estimator
%                (requires Optimization toolbox) 
%   GEYER_IMSE - Compute autocorrelation time tau using Geyer's
%                initial monotone sequence estimator
%
% Kernel density estimation etc.:
%   KERNEL1  - 1D Kernel density estimation of data
%   KERNELS  - Kernel density estimation of independent components of data
%   KERNELP  - 1D Kernel density estimation, with automatic kernel width
%   NDHIST   - Normalized histogram of N-dimensional data
%   HPDI     - Estimates the Bayesian HPD intervals
%
% Manipulation of MCMC chains
%   THIN     - Delete burn-in and thin MCMC-chains
%   JOIN     - Join similar structures of arrays to one structure of arrays
%   BATCH    - Batch MCMC sample chain and evaluate mean/median of batches
%
% Misc:
%   CUSTATS   - Calculate cumulative statistics of data
%   BBPRCTILE - Bayesian bootstrap percentile
%   GRADCHEK  - Checks a user-defined gradient function using finite
%               differences.
%   DERIVATIVECHECK - Compare user-supplied derivatives to
%                     finite-differencing derivatives.
%
