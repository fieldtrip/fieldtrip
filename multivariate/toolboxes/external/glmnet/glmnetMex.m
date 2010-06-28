%--------------------------------------------------------------------------
% glmnetMex.m: Lasso and elastic-net regularized generalized linear models
%--------------------------------------------------------------------------
%     [a0,ca,ia,nin,rsq,alm,nlp,jerr] = ...
%        glmnetMex(parm,x,y,jd,vp,ne,nx,nlam,flmin,ulam,thr,isd,w,ka)
%     [a0,ca,ia,nin,dev,alm,nlp,jerr] = ...
%        glmnetMex(parm,x,y,jd,vp,ne,nx,nlam,flmin,ulam,thr,isd,nc,maxit,kopt)
%
% Extremely efficient procedures for fitting the entire lasso or
% elastic-net regularization path for linear regression, logistic and
% multinomial regression models. The algorithm uses cyclical coordinate
% descent in a pathwise as described in the paper on the maintainer's
% website.
%
% NOTES: This is a MEX-file wrapper of GLMnet.f for MATLAB. Should be called
% only by glmnet.m. For details about input and output arguments, see
% GLMnet.f.
%
% LICENSE: GPL-2
%
% DATE: 13 Jul 2009
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
%    glmnet, glmnetSet, glmnetPrint, glmnetPlot, glmnetPredict and glmnetCoef methods.
%
% EXAMPLES:
%    parm = 1.0;
%    x = [1 1; 2 2; 3 3];
%    y = [1 3 2]';
%    jd = 0;
%    vp = [1 1];
%    ne = 3;
%    nx = 2;
%    nlam = 100;
%    flmin = 0.0001;
%    ulam = 0;
%    thr = 1.0e-4;
%    isd = 0;
%    w = [1 1 1]';
%    ka = 1;
%    [a0,ca,ia,nin,rsq,alm,nlp,jerr] = glmnetMex(parm,x,y,jd,vp,ne,nx,nlam,flmin,ulam,thr,isd,w,ka);
% 
% DEVELOPMENT: 13 Jul 2009: Original version of glmnetMex.m written.
