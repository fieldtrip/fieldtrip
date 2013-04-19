function R = cipsrf(varargin)
%CIPSRF Cumulative Interval Potential Scale Reduction Factor
%
%   [R] = CIPSRF(X,[n0]) or
%   [R] = CIPSRF(x1,x2,x3,...[,n0])
%   returns Cumulative Interval Potential Scale Reduction Factor
%   for collection of MCMC-simulations. Analysis is first based on
%   PSRF-analysis of samples 1 to n0. Then for samples 1 to n0+1
%   and so on.
%
%   Default value for parameter n0 is |X|/2.
%
%   See also
%     IPSRF

% Copyright (C) 1999 Simo Särkkä
% Copyright (C) 2004 Aki Vehtari
%
% This software is distributed under the GNU General Public 
% Licence (version 3 or later); please refer to the file 
% Licence.txt, included with the software, for details.


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

[N,D,M]=size(X);
R = zeros(N-n0+1,D);
for i=n0:N
  R(i-n0+1,:)=ipsrf(X(1:i,:,:));
end
