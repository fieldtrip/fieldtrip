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

% Copyright (C) 2010, Stefan Klanke
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

% deal with padding
pad       = ceil(n/2);
dat       = ft_preproc_padding(dat, 'localmean', pad);

% create smoothing kernel
krn       = ones(1,n)/n;

% do the smoothing
datsmooth = convn(dat, krn, 'same');

% cut the eges
datsmooth = ft_preproc_padding(datsmooth, 'remove', pad);

