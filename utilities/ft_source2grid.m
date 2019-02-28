function [grid] = ft_source2grid(source)

% FT_SOURCE2GRID removes the fields from a source structure that are
% not necessary to reuse the dipole grid in another source estimation.
%
% Use as
%   [grid] = ft_source2grid(source);
%
% The resulting grid can be used in the configuration of another
% run of FT_SOURCEANALYSIS.
%
% See also FTSOURCE2SPARSE, FT_SOURCE2FULL

% Copyright (C) 2004, Robert Oostenveld
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

ft_defaults

grid = keepfields(source, {'pos', 'tri', 'inside', 'outside', 'xgrid', 'ygrid', 'zgrid', 'dim'});

if ~isfield(grid, 'dim') && isfield(grid, 'xgrid') && isfield(grid, 'ygrid') && isfield(grid, 'zgrid')
  sourcemodel.dim = [length(sourcemodel.xgrid) length(sourcemodel.ygrid) length(sourcemodel.zgrid)];
end

if issubfield(source, 'filter')
  sourcemodel.filter = source.filter;
elseif issubfield(source, 'avg.filter')
  sourcemodel.filter = source.avg.filter;
elseif issubfield(source, 'trial.filter')
  ft_error('single trial filters are not supported here');
end

if issubfield(source, 'leadfield')
  sourcemodel.leadfield = source.leadfield;
elseif issubfield(source, 'avg.leadfield')
  sourcemodel.leadfield = source.avg.leadfield;
end
