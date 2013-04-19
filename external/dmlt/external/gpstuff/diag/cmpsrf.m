function [R,neff] = cmpsrf(varargin)
%CMPSRF Cumulative Multivariate Potential Scale Reduction Factor
%
%   [R,neff] = CMPSRF(X,[n0])
%   [R,neff] = CMPSRF(x1,x2,x3,...[,n0])
%   returns Cumulative Multivariate Potential Scale Reduction
%   Factor for collection of MCMC-simulations. Analysis is first
%   based on MPSRF-analysis of samples 1 to n0. Then for samples
%   1 to n0+1 and so on.
%
%   Default value for parameter n0 is |X|/2.
%
%   See also
%     MPSRF, PSRF, CPSRF

% Copyright (C) 1999 Simo Särkkä
%
% This software is distributed under the GNU General Public 
% Licence (version 3 or later); please refer to the file 
% Licence.txt, included with the software, for details.

% 2004-01-22 Aki.Vehtari@hut.fi Added neff, R^2->R, and cleaning

% Handle the input arguments
if prod(size(varargin{nargin}))==1
  count = nargin-1;
  n0 = varargin{nargin};
else
  count = nargin;
  n0 = [];
end

if count<=1
  X = varargin{1};
else
  X = zeros([size(varargin{1}) count]);
  for i=1:count
    X(:,:,i) = varargin{i};
  end
end
if isempty(n0)
  n0 = floor(size(X,1)/2) + 1;
end
if n0 < 2
  error('n0 should be at least 2');
end

% Calculate the reduction factors
[m,n]=size(X);
R = zeros(m-n0+1,n);
neff = R;
for i=n0:size(X,1)
  [R(i-n0+1,:),neff(i-n0+1,:)]=mpsrf(X(1:i,:,:));
end
