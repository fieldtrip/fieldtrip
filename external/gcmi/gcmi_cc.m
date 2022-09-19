function I = gcmi_cc(x, y)
% GCMI_CC Gaussian-Copula Mutual Information between two continuous variables.
%   I = gcmi_cc(x,y) returns the MI between two (possibly multidimensional)
%   continuous variables, x and y, estimated via a Gaussian copula.
%   If x and/or y are multivariate rows must correspond to samples, columns
%   to dimensions/variables. (Samples first axis) 
%   This provides a lower bound to the true MI value.

% ensure samples first axis for vectors
if isvector(x)
    x = x(:);
end
if isvector(y)
    y = y(:);
end
if ndims(x)~=2 || ndims(y)~=2
    error('gcmi_cc: input arrays should be 2d')
end
Ntrl = size(x,1);
Nvarx = size(x,2);
Nvary = size(y,2);

if size(y,1) ~= Ntrl
    error('gcmi_cc: number of trials do not match')
end

% check for repeated values
for xi=1:Nvarx
    if length(unique(x(:,xi)))./Ntrl < 0.9
        warning('Input x has more than 10% repeated values.')
        break
    end
end
for yi=1:Nvary
    if length(unique(y(:,yi)))./Ntrl < 0.9
        warning('Input y has more than 10% repeated values.')
        break
    end
end

% copula normalisation
cx = copnorm(x);
cy = copnorm(y);
% parametric Gaussian MI
I = mi_gg(cx,cy,true,true);

