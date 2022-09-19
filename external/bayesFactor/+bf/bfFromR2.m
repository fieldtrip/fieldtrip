function [bf10,bfApprox10] = bfFromR2(R2,nrSamples,nrRegressors)
% Estimate the BayesFactor of a linear regression using R-squared
% coefficient of determination (R^2)of the regression, based on the formula 
% derived by Liang et al.
% 
% Liang, F., Paulo, R., Molina, G., Clyde, M. A., & Berger, J. O. (2008). 
% Mixtures of g priors for Bayesian variable selection. 
% Journal of the American Statistical Association, 103(481), 410–423. 
% https://doi.org/10.1198/016214507000001337
%
% The function uses Gaussian quadraure
% INPUT
%  R2 - The R-squared of the linear regression.
%  N -  The number of samples in the regression.
%  p - The number of regressors, excluding the constant term
% OUTPUT
%  bf10 - Bayes factor comparing the model with non-zero slopes to the
%  model in which all slopes are zero.
% bf10Approx - Bayes factor calculated using an approximation that does not
%           require integration.
% BK - Dec 2019.

% Formula from Liang  et al.
liangFun = @(g) (sqrt(nrSamples/2)/gamma(1/2).*(1+g).^((nrSamples-nrRegressors-1)/2).*(1+g.*(1-R2)).^(-(nrSamples-1)/2).*g.^(-3/2).*exp(-nrSamples./(2*g)));
bf10 = integral(liangFun,0,Inf);

% Zellner & Siow also derive an approximation that does not require
% integration. No longer particularly useful now that integration is cheap.
% but for comparison
if nargout >1
    a =sqrt(pi)*gamma((nrRegressors+1)/2);
    v2 = nrSamples-nrRegressors-1;
    bfApprox01 = (a*(v2/2)^(nrRegressors/2)*(1-R2)^((v2-1)/2)); % Z&S (3.12)
    bfApprox10 = 1./bfApprox01;
end