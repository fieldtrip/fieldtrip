function [CMI I] = gcmi_ccd(x, y, z, Zm)
% GCCMI_CCD Gaussian-Copula CMI between 2 continuous variables conditioned on 
%           a discrete variable.
%   [CMI I] = gccmi_ccd(x,y,z,Zm) returns the CMI between two (possibly
%   multidimensional) continuous variables, x and y, conditioned on a third 
%   discrete variable, z, estimated via Gaussian copulae.
%   If x and/or y are multivariate rows must correspond to samples, columns
%   to dimensions/variables. (Samples first axis) 
%   z should contain integer values in the range [0 Zm-1] (inclusive).
%   In addition to the CMI, the variable I is returned which is the pooled MI
%   calculated with the copula transform performed within each class

% ensure samples first axis for vectors
if isvector(x)
    x = x(:);
end
if isvector(y)
    y = y(:);
end
if ndims(x)~=2 || ndims(y)~=2
    error('gccmi_cc: input arrays should be 2d')
end
if isvector(z)
    z = z(:);
else
    error('gccmi_ggd: only univariate discrete variable supported');
end
Ntrl = size(x,1);
Nvarx = size(x,2);
Nvary = size(y,2);

if (size(y,1) ~= Ntrl) || (size(z,1) ~= Ntrl)
    error('gccmi_ggd: number of trials do not match')
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
% check values of discrete variable
if min(z)~=0 || max(z)~=(Zm-1) || any(round(z)~=z)
    error('Values of discrete variable z are not correct')
end

% calculate gcmi for each z value
Icond = zeros(Zm,1);
Pz = zeros(1,Zm);
cx = cell(Zm,1);
cy = cell(Zm,1);
for zi=1:Zm
    idx = z==(zi-1);
    Pz(zi) = sum(idx);
    thsx = copnorm(x(idx,:));
    thsy = copnorm(y(idx,:));
    cx{zi} = thsx;
    cy{zi} = thsy;
    I = mi_gg(thsx,thsy,true,true);
    Icond(zi) = I;
end
Pz = Pz ./ Ntrl;

% conditional mutual information
CMI = Pz*Icond;

% full mutual information with copnorm for each z value
I = mi_gg(cell2mat(cx),cell2mat(cy),true,false);

