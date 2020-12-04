function [dat] = ft_preproc_hilbert(dat, option, handlenan, padnan)

% FT_PREPROC_HILBERT computes the Hilbert transpose of the data and optionally
% performs post-processing on the complex representation, e.g. the absolute
% value of the Hilbert transform of a band-pass filtered signal corresponds
% with the amplitude envelope.
%
% Use as
%   [dat] = ft_preproc_hilbert(dat, option)
% where
%   dat        data matrix (Nchans X Ntime)
%   option     string that determines whether and how the Hilbert transform
%              should be post-processed, can be
%                'abs' (default)
%                'complex'
%                'real'
%                'imag'
%                'absreal'
%                'absimag'
%                'angle'
%
% If the data contains NaNs, the output of the affected channel(s) will be
% all(NaN).
%
% See also PREPROC

% Undocumented, and insufficiently tested, options: handlenan and padnan
% intend to make the function nan-aware, FIXME: needs to be tested more.
%
%   handlenan  boolean, can be false (default) or true
%   padnan     scalar, number of samples to pad the edges of the NaN
%              samples, to remove the ringing, (default = 0)

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

% set the defaults if option is not specified
if nargin<2 || isempty(option)
  option = 'abs';
end

if nargin<3 || isempty(handlenan)
  handlenan = false; % FIXME: consider making default true
end

if nargin<4 || isempty(padnan)
  padnan = 0;
end

% preprocessing fails on channels that contain NaN
if any(isnan(dat(:)))
  ft_warning('FieldTrip:dataContainsNaN', 'data contains NaN values');
end

if handlenan
  nonfinite = ~isfinite(dat);
  dat(nonfinite) = 0;
end

% use the non-conjugate transpose to be sure
dat = transpose(hilbert(transpose(dat)));

% do postprocessing of the complex representation
switch option
    case {'yes' 'abs'}
        dat = abs(dat);   % this is the default if 'yes' is specified
    case {'no' 'complex'}
        dat = dat;        % this is the default if 'no' is specified
    case 'real'
        dat = real(dat);
    case 'imag'
        dat = imag(dat);
    case 'absreal'
        dat = abs(real(dat));
    case 'absimag'
        dat = abs(imag(dat));
    case 'angle'
        dat = (angle(dat./abs(dat)));
    case 'unwrap_angle'
        dat = unwrap(angle(dat./abs(dat)),[],2);
    otherwise
        ft_error('incorrect specification of the optional input argument');
end

if handlenan
  if padnan ~= 0
    nonfinite = convn(double(nonfinite), ones(1,padnan), 'same') >0;
  end
  dat(nonfinite) = nan;
end
