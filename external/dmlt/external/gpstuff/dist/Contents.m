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
