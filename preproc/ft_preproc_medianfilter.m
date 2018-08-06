function dat = ft_preproc_medianfilter(dat, order)

% FT_PREPROC_MEDIANFILTER applies a median filter, which smooths the data with a
% boxcar-like kernel, except that it keeps steps in the data. This function requires
% the MATLAB Signal Processing toolbox.
%
% Use as
%   [dat] = ft_preproc_medianfilter(dat, order)
% where
%   dat        data matrix (Nchans X Ntime)
%   order      number, the length of the median filter kernel (default = 25)
%
% See also PREPROC

% Copyright (C) 2008, Robert Oostenveld
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

% set the default filter order
if nargin<2 || isempty(order)
  ft_error('the order of the median filter is not specified');
end

% preprocessing fails on channels that contain NaN
if any(isnan(dat(:)))
  ft_warning('FieldTrip:dataContainsNaN', 'data contains NaN values');
end

% deal with padding
pad = ceil(order/2);
dat = ft_preproc_padding(dat, 'localmean', pad);

hasfast = exist('fastmedfilt1d', 'file');
if hasfast == 2 || hasfast == 3
  % use fast median filter mex file
  for k = 1:size(dat,1)
    dat(k,:) = fastmedfilt1d(dat(k,:), order);
  end
else
  is_matlab=ft_platform_supports('matlabversion',1,inf);
  if is_matlab
    % use Mathworks slow version
    dat = medfilt1(dat, order, [], 2);
  else
    % use helper function that uses Octave's medfilt1
    dat = medfilt1_rowwise(dat, order);
  end
end

% cut the eges
dat = ft_preproc_padding(dat, 'remove', pad);

%%%%%%%%%%%%%%%%%%%%%%
% Helper function
%%%%%%%%%%%%%%%%%%%%%%
function y = medfilt1_rowwise(x,order)
% this function is compatible with Octave;
% Octave's medfilt1 accepts only two input arguments
    y = zeros(size(x));
    for k = 1:size(x,1)
        y(k,:) = medfilt1(x(k,:),order);
    end
