function bf10 = bfFromF(F,df1,df2,n)
% Estimate the Bayes Factor from an F value, using the BIC approximation 
% of Falukenberry 2018.
%
% INPUT
% F - The F-value of the ANOVA
% df1 - Treatment degrees of freedom
% df2 - Error degrees of freedom
% n - Sample size.
% 
% OUTPUT
% bf10 - Bayesfactor for the H1 (effect of treatment) over H0 (no effect)

% Equation 6 in Faulkenberry 2018.
bf01 = sqrt((n.^df1)*(1+F.*df1./df2).^(-n));
bf10 = 1./bf01;