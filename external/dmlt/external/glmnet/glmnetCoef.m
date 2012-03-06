function result = glmnetCoef( object, s )

%--------------------------------------------------------------------------
% glmnetCoef.m: print coefficients from a "glmnet" object
%--------------------------------------------------------------------------
%
% USAGE: 
%    glmnetCoef(fit);
%    glmnetCoef(fit, s);
%
% DETAILS:
%    glmnetCoef(fit, s) is equivalent to glmnetPredict(fit, "coefficients", [], s)
%    See glmnetPredict for more details.
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
%    glmnet, glmnetSet, glmnetPrint, glmnetPredict and glmnetPlot methods.
% 
% EXAMPLES:
%    x=randn(100,20);
%    y=randn(100,1);
%    fit1=glmnet(x,y);
%    glmnetCoef(fit1,0.01) % extract coefficients at a single value of lambda
%
% DEVELOPMENT: 14 Jul 2009: Original version of glmnet.m written.

%

if nargin < 2
    s = object.lambda;
end

result = glmnetPredict(object, 'coefficients', [], s);
