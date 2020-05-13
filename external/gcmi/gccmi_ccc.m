function I = gcmi_ccc(x, y, z)
% GCCMI_CCC Gaussian-Copula CMI between three continuous variables.
%   I = gccmi_ccc(x,y,z) returns the CMI between two (possibly multidimensional)
%   continuous variables, x and y, conditioned on a third, z, estimated via a 
%   Gaussian copula.
%   If x and/or y are multivariate rows must correspond to samples, columns
%   to dimensions/variables. (Samples first axis) 

% ensure samples first axis for vectors
if isvector(x)
    x = x(:);
end
if isvector(y)
    y = y(:);
end
if isvector(z)
    z = z(:);
end
if ndims(x)~=2 || ndims(y)~=2 || ndims(z)~=2
    error('gccmi_ccc: input arrays should be 2d')
end
Ntrl = size(x,1);
Nvarx = size(x,2);
Nvary = size(y,2);
Nvarz = size(z,2);

if (size(y,1) ~= Ntrl) || (size(z,1) ~= Ntrl)
    error('gccmi_ggg: number of trials do not match')
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
for zi=1:Nvarz
    if length(unique(z(:,zi)))./Ntrl < 0.9
        warning('Input z has more than 10% repeated values.')
        break
    end
end
% copula normalisation
cx = copnorm(x);
cy = copnorm(y);
cz = copnorm(z);

% parametric Gaussian CMI
I = cmi_ggg(cx,cy,cz,true,true);

