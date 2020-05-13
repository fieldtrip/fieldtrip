function I = gcmi_model_cd(x, y, Ym)
% GCMI_MODEL_CD Gaussian-Copula Mutual Information between a continuous and a 
%         discrete variable in bits based on ANOVA style model comparison.
%   I = gcmi_model_gd(x,y,Ym) returns the MI between the (possibly multidimensional)
%   continuous variable x and the discrete variable y.
%   For 1D x this is a lower bound to the mutual information.
%   Rows of x correspond to samples, columns to dimensions/variables. 
%   (Samples first axis)
%   y should contain integer values in the range [0 Ym-1] (inclusive).
%   See also: GCMI_MIXTURE_CD

% ensure samples first axis for vectors
if isvector(x)
    x = x(:);
end
if ndims(x)~=2
    error('gcmi_model_cd: input array should be 2d')
end

if isvector(y)
    y = y(:);
else
    error('gcmi_model_cd: only univariate discrete variable supported');
end

Ntrl = size(x,1);
Nvar = size(x,2);

if size(y,1) ~= Ntrl
    error('gcmi_model_cd: number of trials do not match');
end

% check for repeated values
for xi=1:Nvar
    if length(unique(x(:,xi)))./Ntrl < 0.9
        warning('Input x has more than 10% repeated values.')
        break
    end
end

% check values of discrete variable
if min(y)~=0 || max(y)~=(Ym-1) || any(round(y)~=y)
    error('Values of discrete variable y are not correct')
end

% copula normalisation
cx = copnorm(x);
% parametric Gaussian MI
I = mi_model_gd(cx,y,Ym,true,true);

