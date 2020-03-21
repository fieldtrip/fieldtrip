function [bf10,lm] = regression(tbl,formula)
% [bf10,lm] = regression(tbl,formula)
% Given a table and a Wilcoxon formula, perform a regression and return 
% the Bayes Factor for the full model versus the intercept only (null)
% model.
% INPUT
%  tbl - Data table 
% formula - Wilcoxon notation for the regression.
% OUTPUT
% bf10 - Bayes Factor for the regression (against the null model).
% lm  - The LinearModel for the regression.
%
% BK - Dec 2019

if iscell(formula)
    nrFormulas = numel(formula);
    bf10 = nan(1,nrFormulas);
    lm =cell(1,nrFormulas);
    for i=1:numel(formula)
        [bf10(i),lm{i}] = bf.regression(tbl,formula{i});
    end
else
% Fit a standard linear model to get R2
lm = fitlm(tbl,formula);
% Then use the Liang et al formula to determine the Bayes Factor. 
bf10=bf.bfFromR2(lm.Rsquared.Ordinary,lm.NumObservations,lm.NumEstimatedCoefficients-1);
end
