function browse_topoplotER(cfg, data)

% BROWSE_TOPOPLOTER is a simple helper function for FT_DATABROWSER and shows
% the average topography of the selected data.
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

% convert to an ERP
timelock = ft_timelockanalysis([], data);

default             = [];
default.xlim        = [min(timelock.time) max(timelock.time)];
default.marker      = 'on';
default.interactive = 'no';
cfg = mergestruct(cfg, default);

figure;
ft_topoplotER(cfg, timelock);
