function [FM, xf] = ft_preproc_online_filter_apply(FM, x)

% function [FM, xf] = ft_preproc_online_filter_apply(FM, x)
%
% Passes signal x (channels times samples) through the filter,
% returns updated filter model (delay states) and filtered signal.
%
% See also FT_PREPROC_ONLINE_FILTER_INIT

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

[dimX, numX] = size(x);

if dimX > numX
  % explicit algorithm - faster for many channels, 1 sample (=fMRI)
  xf = zeros(size(x));
  
  for k=1:numX;
    z_old = FM.z;
    z0 = x(:,k) - z_old*FM.A2;
    xf(:,k) = z_old * FM.B2 + z0*FM.B1;
    FM.z(:,2:end) = z_old(:,1:end-1);
    FM.z(:,1) = z0;
  end
else
  % use built-in MATLAB stuff - faster for many samples, few channels
  [xf, z] = filter(FM.B, FM.A, x, FM.z',2);
  FM.z = z';
end
