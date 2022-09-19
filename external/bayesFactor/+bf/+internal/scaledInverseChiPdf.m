function y = scaledInverseChiPdf(x,df,scale)
% The Scaled Inverse Chi-Squared Distribution.
% INPUT
% x = the parameter value
% df = degrees of freedom
% scale = scale (tau squared).
% OUTPUT
% y = The probaility density [nrX nrScales]
%
% BK - 2018
if nargin <3
    scale =1; % Default to scaled inverse Chi-squared.
end

[nrScalesInX,nrX] = size(x);
[nrScales,mustBeOne] = size(scale);
assert(mustBeOne==1 && (nrScales==1 || nrScales ==nrScalesInX),'The number of scales does not match');
if nrScalesInX>1
    scale = repmat(scale,[1 nrX]);    
end

z = x<0;
x(z) =NaN;

y = bf.internal.inverseGammaPdf(x,df/2,df.*scale/2);
y(z) = 0;
end
