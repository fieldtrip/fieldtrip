function glmnetPrint( x )

%--------------------------------------------------------------------------
% glmnet.m: print a glmnet object
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Print a summary of the glmnet path at each step along the path.
%
% USAGE: 
%    glmnetPrint(fit)
%
% INPUT ARGUMENTS:
% fit         fitted glmnet object
%
% DETAILS:
%    Three-column matrix with columns Df, dev and Lambda is printed. The Df
%    column is the number of nonzero coefficients (Df is a reasonable name
%    only for lasso fits). dev is the percent deviance explained (relative
%    to the null deviance).
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
%    glmnet, glmnetSet, glmnetPlot, glmnetPredict and glmnetCoef methods.
% 
% EXAMPLES:
%    x=randn(100,20);
%    y=randn(100,1);
%    fit1=glmnet(x,y);
%    glmnetPrint(fit1);
%
% DEVELOPMENT: 14 Jul 2009: Original version of glmnet.m written.

disp(sprintf('\tDf\t%%Dev\tLambda'));
% disp([x.df, x.dev, x.lambda]);
for i=1:length(x.lambda)
    disp(sprintf('%d\t%d\t%f\t%f', i, x.df(i), x.dev(i), x.lambda(i)));
end
