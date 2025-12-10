function browse_topoplotVAR(cfg, data)

% BROWSE_TOPOPLOTVAR is a simple helper function for FT_DATABROWSER that
% computes the variance of band-pass filtered data and makes a topographic
% plot. It serves to make a quick-and-dirty power topography.
%
% See also BROWSE_MOVIEPLOTER, BROWSE_TOPOPLOTER, BROWSE_MULTIPLOTER, BROWSE_TOPOPLOTVAR, BROWSE_SIMPLEFFT

% Copyright (C) 2009, Robert Oostenveld
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

% compute the variance, i.e. the broad-band power
timelock        = [];
timelock.label  = data.label;
timelock.time   = mean(data.time{1});
timelock.avg    = sum(ft_preproc_baselinecorrect(data.trial{1}).^2, 2);
timelock.dimord = 'chan_time';

if isfield(data, 'grad')
  timelock.grad = data.grad;
end
if isfield(data, 'elec')
  timelock.elec = data.elec;
end

default             = [];
default.markers     = 'labels';
default.interactive = 'no';
cfg = mergestruct(cfg, default);

figure;
ft_topoplotER(cfg, timelock);
