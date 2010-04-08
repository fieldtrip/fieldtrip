function [grid] = source2grid(source)

% SOURCE2GRID removes the fields from a source structure that are
% not neccessary to reuse the dipole grid in another source estimation.
%
% Use as
%   [grid] = source2grid(source);
%
% The resulting grid can be used in the configuration of another
% run of SOURCANALYSIS.
%
% See also SOURCE2SPARSE, SOURCE2FULL

% Copyright (C) 2004, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fieldtripdefs

% these are always supposed to be present
grid.pos     = source.pos;
grid.inside  = source.inside;
grid.outside = source.outside;

% these are optional
try, grid.xgrid   = source.xgrid; end
try, grid.ygrid   = source.ygrid; end
try, grid.zgrid   = source.zgrid; end
try, grid.dim     = source.dim;   end

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
