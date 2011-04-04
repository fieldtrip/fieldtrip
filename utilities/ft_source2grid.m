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
% See also SOURCE2SPARSE, SOURCE2FULL

% Copyright (C) 2004, Robert Oostenveld
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

ft_defaults

% these are always supposed to be present
grid.pos     = source.pos;
grid.inside  = source.inside;
grid.outside = source.outside;

% these are optional
try, grid.xgrid   = source.xgrid; end
try, grid.ygrid   = source.ygrid; end
try, grid.zgrid   = source.zgrid; end
try, grid.dim     = source.dim;   end
try, grid.tri     = source.tri;   end % only in case of a tesselated/triangulated cortical sheet source model

if ~isfield(grid, 'dim') && isfield(grid, 'xgrid') && isfield(grid, 'ygrid') && isfield(grid, 'zgrid') 
  grid.dim = [length(grid.xgrid) length(grid.ygrid) length(grid.zgrid)];
end

if issubfield(source, 'filter')
  grid.filter = source.filter;
elseif issubfield(source, 'avg.filter')
  grid.filter = source.avg.filter;
elseif issubfield(source, 'trial.filter')
  error('single trial filters are not supported here');
end

if isfield(source, 'leadfield')
  grid.leadfield = source.leadfield;
end
