function bf10 = bfFromT(T,df2)
% Estimate the Bayes Factor from a T statistic associated with a 
% two independent sample T-test, using the BIC approximation of Falukenberry 2018. 
%
% Note that this is provided for completeness only; the bf.ttest function
% provides a better calculation of the Bayes factor (also based on T, df, and n)
% that does not make the BIC approximation.  The BIC approximation is ok
% for large sample size, but overestimates the Bayes Factor for smaller n
% (<50) - see commonDesigns.m
%
% INPUT
% F - The T-value of the test
% df - Degrees of freedom of the T test (
%
% OUTPUT
% bf10 - Bayesfactor for the H1 (effect of treatment) over H0 (no effect)


n = df2+2; % Total number of samples.
F = T.^2;
df1 =1; % Treatment df=1

% Now call the bfFromF
bf10 = bf.bfFromF(F,df1,df2,n);