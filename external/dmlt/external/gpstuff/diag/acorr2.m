function c = acorr2(x,maxlag)
%ACORR Estimate autocorrelation function of time series
%
%   C = ACORR(X,MAXLAG) returns normalized autocorrelation
%   sequences for each column of X computing correlation via FFT
%

% Copyright (C) 2000-2008 Aki Vehtari
%
% This software is distributed under the GNU General Public
% Licence (version 3 or later); please refer to the file
% Licence.txt, included with the software, for details.

if nargin < 1
  error('Not enough input arguments.');
end
if nargin < 2
  maxlag=length(x)-1;
end
[m,n]=size(x);
c=zeros(maxlag,n);
for i1=1:n
  xn = x(:,i1)-mean(x(:,i1));
  xf = fft(xn,2^nextpow2(2*m-1));
  cf = ifft(abs(xf).^2);
  c(:,i1) = cf(2:(maxlag+1))./cf(1);
end
