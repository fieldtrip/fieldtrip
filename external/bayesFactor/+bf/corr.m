function [bf10,r,p] = corr(arg1,arg2)
% Calculate the Bayes Factor for Pearson correlation between two
% variables.
% INPUT
% (x,y)  - two vectors of equal length.
% OR
% (r,n)  - the correlation and number of samples
%
% OUTPUT
% bf10 = the Bayes Factor for the hypothesis that r is differnt
%           from zero (two-tailed).
% r - the correlation
% p - the tradiational p-value based on Fisher-transformed
%
if isscalar(arg1) && isscalar(arg2)
    r= arg1;
    n= arg2;
else
    x=arg1;y=arg2;
    [r,p] = corr(x,y,'type','pearson');
    n=numel(x);
end

% Code from Sam Schwarzkopf
F = @(g,r,n) exp(((n-2)./2).*log(1+g)+(-(n-1)./2).*log(1+(1-r.^2).*g)+(-3./2).*log(g)+-n./(2.*g));
bf10 = sqrt((n/2)) / gamma(1/2) * integral(@(g) F(g,r,n),0,Inf);

if nargout>1
    % Compute classical stats too
    t = r.*sqrt((n-2)./(1-r.^2));
    p = 2*tcdf(-abs(t),n-2);
end
end
