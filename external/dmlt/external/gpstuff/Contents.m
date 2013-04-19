% THE GP TOOLS (in the gp-folder):
% 
%  Gaussian process utilities:
%   GP_SET    Create and modify a Gaussian Process structure. 
%   GP_OPTIM  Optimize paramaters of a Gaussian process 
%   GP_PAK    Combine GP parameters into one vector.
%   GP_UNPAK  Set GP parameters from vector to structure
%   GP_COV    Evaluate covariance matrix between two input vectors. 
%   GP_TRCOV  Evaluate training covariance matrix (gp_cov + noise covariance). 
%   GP_TRVAR  Evaluate training variance vector. 
%   GP_RND    Random draws from the postrior Gaussian process
%
%  Covariance functions:
%   GPCF_CAT          Create a categorigal covariance function
%   GPCF_CONSTANT     Create a constant covariance function 
%   GPCF_EXP          Create a squared exponential covariance function
%   GPCF_LINEAR       Create a linear covariance function
%   GPCF_MASK         Create a mask covariance function
%   GPCF_MATERN32     Create a Matern nu=3/2 covariance function
%   GPCF_MATERN52     Create a Matern nu=5/2 covariance function
%   GPCF_NEURALNETWORK Create a neural network covariance function
%   GPCF_NOISE        Create a independent noise covariance function
%   GPCF_PERIODIC     Create a periodic covariance function
%   GPCF_PPCS0        Create a piece wise polynomial (q=0) covariance function 
%   GPCF_PPCS1        Create a piece wise polynomial (q=1) covariance function 
%   GPCF_PPCS2        Create a piece wise polynomial (q=2) covariance function 
%   GPCF_PPCS3        Create a piece wise polynomial (q=3) covariance function 
%   GPCF_PROD         Create a product form covariance function 
%   GPCF_RQ           Create a rational quadratic covariance function 
%   GPCF_SCALED       Create a scaled covariance function
%   GPCF_SEXP         Create a squared exponential covariance function
%   GPCF_SUM          Create a sum form covariance function
%
%  Mean functions:
%   GPMF_CONSTANT     Create a constant mean function
%   GPMF_LINEAR       Create a linear mean function
%   GPMF_SQUARED      Create a squared mean function
%
%  Likelihood functions:
%   LIK_BINOMIAL    Create a binomial likelihood structure 
%   LIK_GAUSSIAN    Create a Gaussian likelihood structure
%   LIK_GAUSSIANSMT Create a Gaussian scale mixture approximating t
%   LIK_LOGIT       Create a Logit likelihood structure 
%   LIK_NEGBIN      Create a Negbin likelihood structure 
%   LIK_NEGBINZTR   Create a zero-truncated Negbin likelihood structure
%   LIK_POISSON     Create a Poisson likelihood structure 
%   LIK_PROBIT      Create a Probit likelihood structure 
%   LIK_T           Create a Student-t likelihood structure 
%   LIK_WEIBULL     Create a Weibull likelihood structure 
%
% Inference utilities:
%   GP_E          Evaluate energy function (un-normalized negative marginal 
%                 log posterior) 
%   GP_G          Evaluate gradient of energy (GP_E) for Gaussian Process
%   GP_EG         Evaluate both GP_E and GP_G. Useful in optimisation.
%   GP_PRED       Make predictions with Gaussian process 
%   GP_CPRED      Conditional predictions using specific covariates
%   GP_JPRED      Joint predictions with Gaussian process 
%   GP_PREPRCTMU  Percentiles of the distribution of the location parameter
%   GP_PREPRCTY   Percentiles of the predictive distribution at test points
%   GP_MC         Markov chain sampling for Gaussian process models
%   GPMC_PREDS    Conditional predictions with Gaussian Process MCMC
%                 approximation.
%   GP_IA         Integration approximation with grid, Monte Carlo or
%                 CCD integration
%   LGCP          Log Gaussian Cox Process intensity estimate for 1D and 
%                 2D data
%
%  Model assesment and comparison:
%   GP_KFCV       K-fold cross validation for a GP model
%   GP_LOOPRED    Leave-one-out-predictions with Gaussian Process
%   GP_LOOE       Evaluate the leave-one-out predictive density in case of
%                 Gaussian observation model
%   GP_LOOG       Evaluate the gradient of the leave-one-out predictive 
%                 density (GP_LOOE) in case of Gaussian observation model 
%   GP_WAIC       The widely applicable information criterion
%   GP_DIC        The DIC statistics and effective number of parameters
%   GP_PEFF       The efective number of parameters in GP model with focus 
%                 on latent variables.
%   GP_AVPREDCOMP Average predictive comparison for Gaussian process model
%
%  Metrics:
%   METRIC_EUCLIDEAN   An Euclidean distance for Gaussian process models.
%  
%  Misc:
%   LDLROWMODIFY  Function to modify the sparse cholesky factorization 
%                 L*D*L' = C, when a row and column k of C have changed 
%   LDLROWUPDATE  Multiple-rank update or downdate of a sparse LDL' factorization.
%   SPINV         Evaluate the sparsified inverse matrix
%   SCALED_HMC    A scaled hybric Monte Carlo samping for latent values
%   SCALED_MH     A scaled Metropolis Hastings samping for latent values
%   SURROGATE_SLS Markov chain Monte Carlo sampling using Surrogate data Slice Sampling
%   ESLS          Markov chain update for a distribution with a Gaussian "prior" factored out
%   GP_INSTALL    Matlab function to compile all the c-files to mex in the 
%                 GPstuff/gp folder.
%
%  Demonstration programs:
%   DEMO_BINOMIAL1          Demonstration of Gaussian process model with binomial
%                           likelihood
%   DEMO_BINOMIAL2          Demonstration of Gaussian process model with binomial
%                           likelihood
%   DEMO_BINOMIAL_APC       Demonstration for modeling age-period-cohort data
%                           by a binomial model combined with GP prior.
%   DEMO_CLASSIFIC          Classification problem demonstration for 2 classes 
%   DEMO_DERIVATIVEOBS      Regression problem demonstration with derivative 
%                           observations
%   DEMO_HURDLE             Demonstration of Logit Negative-binomial hurdle model
%                           using Gaussian process prior
%   DEMO_LGCP               Demonstration for a log Gaussian Cox process
%                           with inference via EP or Laplace approximation
%   DEMO_MODELASSESMENT1    Demonstration for model assesment with DIC, number 
%                           of effective parameters and ten-fold cross validation
%   DEMO_MODELASSESMENT2    Demonstration for model assesment when the observation 
%                           model is non-Gaussian
%   DEMO_NEURALNETWORKCOV   Demonstration of Gaussian process with a neural
%                           network covariance function
%   DEMO_PERIODIC           Regression problem demonstration for periodic data
%   DEMO_REGRESSION1        Regression problem demonstration for 2-input 
%                           function with Gaussian process
%   DEMO_REGRESSION_ADDITIVE1 Regression demonstration demonstration with additive model
%   DEMO_REGRESSION_ADDITIVE2 Regression demonstration with additive Gaussian
%                           process using linear, squared exponential and
%                           neural network covariance fucntions 
%   DEMO_REGRESSION_HIER    Hierarchical regression demonstration
%   DEMO_REGRESSION_MEANF   Regression problem demonstration for GP model with a
%                           mean function
%   DEMO_REGRESSION_PPCS    Regression problem demonstration for 2-input 
%                           function with Gaussian process using CS covariance
%   DEMO_REGRESSION_ROBUST  A regression demo with Student-t distribution as a 
%                           residual model.
%   DEMO_REGRESSION_SPARSE1 Regression problem demonstration for 2-input 
%                           function with sparse Gaussian processes
%   DEMO_REGRESSION_SPARSE2 Regression demo comparing different sparse
%                           approximations
%   DEMO_SPATIAL1           Demonstration for a disease mapping problem
%                           with Gaussian process prior and Poisson likelihood
%   DEMO_SPATIAL2           Demonstration for a disease mapping problem with 
%                           Gaussian process prior and negative binomial 
%                           observation model
%   DEMO_SURVIVAL_WEIBULL   Survival model using Weibull baseline hazard
%
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
% PROBABILITY DISTRIBUTION FUNCTIONS (in the dist-folder):
%
% Priors 
%  PRIOR_FIXED       Fix parameter to its current value
%  PRIOR_GAMMA       Gamma prior structure     
%  PRIOR_GAUSSIAN    Gaussian prior structure     
%  PRIOR_INVGAMMA    Inverse-gamma prior structure     
%  PRIOR_INVT        Inverse Student-t prior structure
%  PRIOR_INVT        Inverse uniform prior structure
%  PRIOR_LAPLACE     Laplace (double exponential) prior structure
%  PRIOR_LOGGAUSSIAN Log-Gaussian prior structure     
%  PRIOR_LOGLOGUNIF  Uniform prior structure for the log(log(parameter))
%  PRIOR_LOGT        Student-t prior structure for the logarithm of the parameter
%  PRIOR_LOGUNIF     Uniform prior structure for the logarithm of the parameter
%  PRIOR_SINVCHI2    Scaled inverse-chi-square prior structure
%  PRIOR_SQRTT       Student-t prior structure for the square root of the
%                    parameter
%  PRIOR_SQRTUNIF    Uniform prior structure for the square root of the
%                    parameter
%  PRIOR_UNIF        Uniform prior structure     
%  PRIOR_T           Student-t prior structure
%  PRIOR_SQRTUNIF    Uniform prior structure for the square root of the
%                    parameter
%
% Probability density functions
%
%    BETA_LPDF     - Beta log-probability density function (lpdf).
%    BETA_PDF      - Beta probability density function (pdf).
%    DIR_LPDF      - Log probability density function of uniform Dirichlet
%                    distribution
%    DIR_PDF       - Probability density function of uniform Dirichlet
%                    distribution
%    GAM_CDF       - Cumulative of Gamma probability density function (cdf).
%    GAM_LPDF      - Log of Gamma probability density function (lpdf).
%    GAM_PDF       - Gamma probability density function (pdf).
%    GEO_LPDF      - Geometric log probability density function (lpdf).
%    INVGAM_LPDF   - Inverse-Gamma log probability density function.
%    INVGAM_PDF    - Inverse-Gamma probability density function.
%    LAPLACE_LPDF  - Laplace log-probability density function (lpdf).
%    LAPLACE_PDF   - Laplace probability density function (pdf).
%    LOGN_LPDF     - Log normal log-probability density function (lpdf)
%    LOGT_LPDF     - Log probability density function (lpdf) for log Student's T
%    MNORM_LPDF    - Multivariate-Normal log-probability density function (lpdf).
%    MNORM_PDF     - Multivariate-Normal log-probability density function (lpdf).
%    NORM_LPDF     - Normal log-probability density function (lpdf).
%    NORM_PDF      - Normal probability density function (pdf).
%    POISS_LPDF    - Poisson log-probability density function.
%    POISS_PDF     - Poisson probability density function.
%    SINVCHI2_LPDF - Scaled inverse-chi log-probability density function.
%    SINVCHI2_PDF  - Scaled inverse-chi probability density function.
%    T_LPDF        - Student's T log-probability density function (lpdf)
%    T_PDF         - Student's T probability density function (pdf)
%
% Random number generators
%
%    CATRAND       - Random matrices from categorical distribution.
%    DIRRAND       - Uniform dirichlet random vectors
%    EXPRAND       - Random matrices from exponential distribution.
%    GAMRAND       - Random matrices from gamma distribution.
%    INTRAND       - Random matrices from uniform integer distribution.
%    INVGAMRAND    - Random matrices from inverse gamma distribution
%    INVGAMRAND1   - Random matrices from inverse gamma distribution
%    INVWISHRND    - Random matrices from inverse Wishart distribution.
%    NORMLTRAND    - Random draws from a left-truncated normal
%                    distribution, with mean = mu, variance = sigma2
%    NORMRTRAND    - Random draws from a right-truncated normal
%                    distribution, with mean = mu, variance = sigma2
%    NORMTRAND     - Random draws from a normal truncated to interval
%    NORMTZRAND    - Random draws from a normal distribution truncated by zero
%    SINVCHI2RAND  - Random matrices from scaled inverse-chi distribution
%    TRAND         - Random numbers from Student's t-distribution
%    UNIFRAND      - Generate unifrom random numberm from interval [A,B]
%    WISHRND       - Random matrices from Wishart distribution.
%
% Others
%    KERNELP       - Kernel density estimator for one dimensional distribution.
%    HAMMERSLEY    - Hammersley quasi-random sequence
%
% MCMC FUNCTIONS (in the mc-folder):
%
%    BBMEAN        - Bayesian bootstrap mean
%    GIBBS         - Gibbs sampling
%    HMC2          - Hybrid Monte Carlo sampling.
%    HMC2_OPT      - Default options for Hybrid Monte Carlo sampling.
%    HMC_NUTS      - No-U-Turn Sampler (NUTS)
%    HMEAN         - Harmonic mean
%    METROP2       - Markov Chain Monte Carlo sampling with Metropolis algorithm.
%    METROP2_OPT   - Default options for Metropolis sampling.
%    RANDPICK      - Pick element from x randomly
%                    If x is matrix, pick row from x randomly.
%    RESAMPDET     - Deterministic resampling
%    RESAMPRES     - Residual resampling
%    RESAMPSIM     - Simple random resampling
%    RESAMPSTR     - Stratified resampling
%    SLS           - Markov Chain Monte Carlo sampling using Slice Sampling
%    SLS_OPT       - Default options for Slice Sampling
%    SLS1MM        - 1-dimensional fast minmax Slice Sampling
%    SLS1MM_OPT    - Default options for SLS1MM_OPT
%
% MISCELLANEOUS FUNCTIONS (in the misc-folder)
%
%    ADDLOGS      - Add numbers represented by their logarithms.
%    BINSGEQ      - Binary search of sorted vector
%    CVIT         - Create itr and itst indeces for k-fold-cv
%    CVITR        - Create itr and itst indeces for k-fold-cv with ranndom permutation
%    DENORMDATA   - De-normalize normalized data
%    LOGIT        - Logit transformation
%    LOGITINV     - Inverse of the logit transformation
%    M2KML        - Converts GP prediction results to a KML file
%    MAPCOLOR     - Returns a colormap ranging from blue through gray to red
%    MAPCOLOR2    - Create a blue-gray-red colormap.
%    NORMDATA     - Normalize input to zero mean and unit variance
%    SET_PIC      - Set the inducing inputs and blocks for two dimensional input data
%    SOFTMAX2     - Softmax transfer function
%    STR2FUN      - Compatibility wrapper to str2func
%    SUMLOGS      - Sum of vector where numbers are represented by their logarithms.
%    VIOLINPLOT   - Plot a violinplot
%    WMEAN        - Weighted mean
%
% OPTIMIZATION FUNCTIONS (in the optim-folder)
%
%    BSEARCH       - Finds the minimum of a combinatorial function using backward search
%    BSEARCH_OPT   - Default options for backward search
%    FSEARCH       - Finds the minimum of a combinatorial function using forward search
%    FSEARCH_OPT   - Default options for forward search
%    SCGES         - Scaled conjugate gradient optimization with early stopping
%    SCGES_OPT     - Default options for scaled conjugate gradient optimization
%    SCGES         - Scaled conjugate gradient optimization with early stopping (new options structure).
%    SCG2          - Scaled conjugate gradient optimization
%    SCG2_OPT      - Default options for scaled conjugate gradient optimization (scg2) (new options structure).
%    FMINSCG       - Scaled conjugate gradient optimization
%    FMINLBFGS     - (Limited memory) Quasi Newton
%                    Broyden-Fletcher-Goldfarb-Shanno (BFGS)
