function [state, xf] = ft_preproc_online_filter_apply(state, x)

% FT_PREPROC_ONLINE_FILTER_APPLY passes a signal through the online filter and
% returns the updated filter model (delay states) and the filtered signal.
%
% Use as
%   [state, dat] = ft_preproc_online_filter_apply(state, dat)
% where
%   dat   = Nchan x Ntime
%   state = filter state, see FT_PREPROC_ONLINE_FILTER_INIT
%
% See also PREPROC

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
    z_old = state.z;
    z0 = x(:,k) - z_old*state.A2;
    xf(:,k) = z_old * state.B2 + z0*state.B1;
    state.z(:,2:end) = z_old(:,1:end-1);
    state.z(:,1) = z0;
  end
else
  % use built-in MATLAB stuff - faster for many samples, few channels
  [xf, z] = filter(state.B, state.A, x, state.z',2);
  state.z = z';
end
