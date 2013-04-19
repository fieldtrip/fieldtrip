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
