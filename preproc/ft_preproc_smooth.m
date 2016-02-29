function [datsmooth] = ft_preproc_smooth(dat, n)

% FT_PREPROC_SMOOTH performs boxcar smoothing with specified length.
% Edge behavior is improved by implicit padding with the mean over
% half the boxcar length at the edges of the data segment.
%
% Use as
%   datsmooth = ft_preproc_smooth(dat, n)
%
% Where dat is an Nchan x Ntimepoints data matrix, and n the length
% of the boxcar smoothing kernel

% Undocumented option:
%  n can also be a vector containing a custom smoothing kernel
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

% create smoothing kernel
if isscalar(n)
  krn = ones(1,n)/n;
else
  krn = n(:)'./sum(n(:));
  n   = numel(krn);
end

% deal with padding
dat = ft_preproc_padding(dat, 'localmean', ceil(n/2));


% do the smoothing
if n<100
  % heuristic: for large kernel the convolution is faster when done along
  % the columns, weighing against the costs of doing the transposition.
  % the threshold of 100 is a bit ad hoc.
  datsmooth = convn(dat,   krn,   'same');
else
  datsmooth = convn(dat.', krn.', 'same').';
end

% cut the eges
datsmooth = ft_preproc_padding(datsmooth, 'remove', ceil(n/2));

