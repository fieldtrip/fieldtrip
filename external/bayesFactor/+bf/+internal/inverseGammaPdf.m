function y = inverseGammaPdf(x,alpha,beta)
% The inverse Gamma PDF.
% INPUT
% x  (>0)
% alpha - shape parameter
% beta  - scale parameter
% OUTPUT
% y = the pdf [nrX nrBeta]
% BK - 2018

[nrX,nrBetaInX] = size(x);
[nrBeta,nrXInBeta] = size(beta);
if nrXInBeta<nrBetaInX     
    beta =repmat(beta,[1 nrBetaInX]);
end

if nrX==1 && nrBeta>1
    x = repmat(x,[1 nrBeta]);
end

z = x<0;
x(z) =NaN;
y = (beta.^alpha)./gamma(alpha).*(1./x).^(alpha+1).*exp(-beta./x);
y(z) = 0;
end
