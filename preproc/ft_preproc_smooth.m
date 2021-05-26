function [datsmooth] = ft_preproc_smooth(dat, n, tol)

% FT_PREPROC_SMOOTH performs boxcar smoothing with specified length.
% Edge behavior is improved by implicit padding with the mean over
% half the boxcar length at the edges of the data segment.
%
% Use as
%   datsmooth = ft_preproc_smooth(dat, n)
%
% Where dat is an Nchan x Ntimepoints data matrix, and n the length
% of the boxcar smoothing kernel
%
% If the data contains NaNs, these are ignored for the computation, but
% retained in the output.
%
% See also PREPROC

% Undocumented options:
%  n can also be a vector containing a custom smoothing kernel
%  n can also be 'regsmooth', using a regularised estimate of the first temporal derivative 
%
% Copyright (C) 2010, Stefan Klanke
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% preprocessing fails on channels that contain NaN
if any(isnan(dat(:)))
  ft_warning('FieldTrip:dataContainsNaN', 'data contains NaN values');
end

% create smoothing kernel
regflag = false;
if isequal(n, 'regsmooth')
  regflag = true;
  n = 0;
elseif isscalar(n)
  krn = ones(1,n)/n;
else
  krn = n(:)'./sum(n(:));
  n   = numel(krn);
end

% deal with padding
dat = ft_preproc_padding(dat, 'localmean', ceil(n/2));

% do the smoothing
if regflag
  if nargin<3
    tol = 1e-9;
  end
  datsmooth = smooth_regularised(dat, tol);
elseif n<100
  % heuristic: for large kernel the convolution is faster when done along
  % the columns, weighing against the costs of doing the transposition.
  % the threshold of 100 is a bit ad hoc.
  datsmooth = convn(dat,   krn,   'same');
else
  datsmooth = convn(dat.', krn.', 'same').';
end

% cut the eges
datsmooth = ft_preproc_padding(datsmooth, 'remove', ceil(n/2));

function out = smooth_regularised(dat, tol)

tol = tol.*std(dat(:));
n   = size(dat, 2);

B    = eye(n);

% create Toeplitz matrix F
r = [1 zeros(1,n-1)];
d = [1 zeros(1,n-1)];
m = 2;
for k=1:m
 d = filter([1 -1],1,d);
end
F = toeplitz(d,r);

% create Toeplitz matrix G
c = ones(n,1);
G = toeplitz(c,r);

% compute gamma parameter

% SVD
H       = B*(G/F);
[U,D,V] = svd(H);
diagD2  = diag(D).^2;
diagD   = diag(D);
epsi    = U'*B*dat.';

gammarange = 10.^(-10:0.2:10);
ngamma     = numel(gammarange);
onevec_n   = ones(n,1);
onevec_ngamma = ones(1,ngamma);

nchan = size(dat,1);
ro    = zeros(n,ngamma,nchan);
for k = 1:nchan
  ro(:,:,k) = (gammarange(onevec_n,:).*epsi(:,k.*onevec_ngamma))./(diagD2(:,onevec_ngamma) + gammarange(onevec_n,:));
end
wrss = squeeze(sum(ro.^2));

gamma = zeros(nchan,1);
ni = zeros(n,nchan);
for k = 1:nchan
  gamma(k,1) = gammarange(find(wrss(:,k)<tol,1,'last'));
  ni(:,k) = (diagD.*epsi(:,k))./(diagD2 + gamma(k));
end

ddata_smooth = transpose(F\V*ni);
out          = ddata_smooth*G';
