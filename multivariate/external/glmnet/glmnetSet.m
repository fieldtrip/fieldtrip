function options = glmnetSet(opts)

%--------------------------------------------------------------------------
% glmnetSet creates or alters an options structure for glmnet.m.
%--------------------------------------------------------------------------
%   options = glmnetSet; (with no input arguments)
%   creates a structure with all fields set to their default values.
%   Each field is an option (also called a parameter).
%
%   glmnetSet (with no input or output arguments)
%   displays all options and their default values.
%
%   options = glmnetSet(opts); 
%   creates a structure with all fields set to their default values,
%   except valid fields in the structure "opts" replace the defaults.
%
% options.weights     Observation weights. Can be total counts if responses
%                     are proportion matrices. Default is 1 for each
%                     observation.
% options.alpha       The elasticnet mixing parameter, with 0 < alpha <= 1.
%                     The penalty is defined as
%                           (1-alpha)/2(||beta||_2)^2+alpha||beta||_1.
%                     Default is alpha = 1, which is the lasso penalty;
%                     Currently alpha < 0.01 is not reliable, unless you
%                     supply your own lambda sequence.
% options.nlambda     The number of lambda values - default is 100.
% options.lambda_min  Smallest value for lambda, as a fraction of
%                     lambda_max, the (data derived) entry value (i.e., the
%                     smallest value for which all coefficients are zero).
%                     The default depends on the sample size nobs relative
%                     to the number of variables nvars. If nobs > nvars,
%                     the default is 0.0001, close to zero. If nobs <
%                     nvars, the defaults is 0.05. A very small value of
%                     lambda_min will lead to a saturated fit. This is
%                     undefined for "binomial" and "multinomial" models,
%                     and glmnet will exit gracefully when the percentage
%                     deviance explained is almost 1.
% options.lambda      A user supplied lambda sequence. Typical usage is to
%                     have the program compute its own lambda sequence
%                     based on nlambda and lambda_min. Supplying a value of
%                     lambda override this. Use with care - it is better to
%                     supply a decreasing sequence of lambda values than a
%                     single (small) value.
% options.standardize Logical flag for variable standardization, prior to
%                     fitting the model sequence. The coefficients are
%                     always returned on the original scale. Default is
%                     standardize = true.
% options.thresh      Convergence threshold for coordinate descent. Each 
%                     inner coordinate-descent loop continues until the 
%                     relative change in any coefficient is less than
%                     thresh. Defaults value is 1E-4.
% options.dfmax       Limit the maximum number of variables in the model. 
%                     Useful for very large nvars, if a partial path is
%                     desired. Default is nvars + 1.
% options.pmax        Limit the maximum number of variables ever to be
%                     nonzero. Default is min(dfmax * 1.2, nvars).
% options.exclude     Indices of variables to be excluded from the model. 
%                     Default is none. Equivalent to an infinite penalty
%                     factor (next item).
% options.penalty_factor
%                     Separate penalty factors can be applied to each
%                     coefficient. This is a number that multiplies lambda
%                     to allow differential shrinkage. Can be 0 for some
%                     variables, which implies no shrinkage, and that
%                     variable is always included in the model. Default is
%                     1 for all variables (and implicitly infinity for
%                     variables listed in exclude).
% options.maxit       Maximum number of outer-loop iterations for
%                     'binomial' or 'multinomial' families. Default is 100.
% options.HessianExact 
%                     Only applies to 'binomial' or 'multinomial' families.
%                     If false (the default), an upper-bound approximation
%                     is made to the hessian, which is not recalculated at
%                     each outer loop.
% options.type        Two algorithm types are supported for (only)
%                     family = 'gaussian'. The default type = 'covariance'
%                     saves all inner-products ever computed, and can be
%                     much faster than type = 'naive'. The latter can be
%                     more efficient for p >> N situations.
%
% LICENSE: GPL-2
%
% DATE: 14 Jul 2009
%
% AUTHORS:
%    Algorithm was designed by Jerome Friedman, Trevor Hastie and Rob Tibshirani 
%    Fortran code was written by Jerome Friedman 
%    R wrapper (from which the MATLAB wrapper was adapted) was written by Trevor Hasite
%    MATLAB wrapper was written and maintained by Hui Jiang, jiangh@stanford.edu 
%    Department of Statistics, Stanford University, Stanford, California, USA.
%
% REFERENCES:
%    Friedman, J., Hastie, T. and Tibshirani, R. (2009)
%    Regularization Paths for Generalized Linear Models via Coordinate Descent.
%    Journal of Statistical Software, 33(1), 2010
%
% SEE ALSO:
%    glmnet, glmnetPrint, glmnetPlot, glmnetPredict and glmnetCoef methods.
%
% EXAMPLES:
% 
% DEVELOPMENT: 14 Jul 2009: Original version of glmnet.m written.

% Set default options.
  options.weights        =           [];
  options.alpha          =          1.0;
  options.nlambda        =          100;
  options.lambda_min     =            0;
  options.lambda         =           [];
  options.standardize    =         true;
  options.thresh         =         1E-4;
  options.dfmax          =            0;
  options.pmax           =            0;
  options.exclude        =           [];
  options.penalty_factor =           [];
  options.maxit          =          100;
  options.HessianExact   =        false;
  options.type           = 'covariance';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End default options.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Quick return if no user opts
  if nargin == 0 || isempty(opts)
    if nargout == 0    % Display options.
      disp('pdco default options:')
      disp( options )
    end
    return
  end

% List of valid field names
  vfields = fieldnames( options );

% Grab valid fields from user's opts
  for i = 1:length(vfields)
    field = vfields{i};
    if isfield( opts, field );
      options.(field) = opts.(field);
    end
  end
