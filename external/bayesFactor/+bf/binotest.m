function [bf10,p,pHat] = binotest(k,n,p)
% Bayes factor for binomial test with k successes, n trials and base probability p.
% INPUT
%  k - number of successes
%  n - number of draws
%  p - true binomial probabiliy
% OUTPUT
% bf - Bayes Factor representing the evidence that this n/k
% could result from random draws with p (BF>1) or not (BF<1)
% p - p-value of a traditional test
% pHat - esttimae of the binomial probablity

% Code from Sam Schwarzkopf
F = @(q,k,n,p) nchoosek(n,k) .* q.^k .* (1-q).^(n-k);
bf01 = (nchoosek(n,k) .* p.^k .* (1-p).^(n-k)) / integral(@(q) F(q,k,n,p),0,1);
bf10 = 1/bf01;

if nargout>1
    % Traditional tests
    pHat = binofit(k,n,0.05);
    p = 1-binocdf(k,n,p);
end

end
