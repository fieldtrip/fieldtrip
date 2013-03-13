function [cfg, varargout] = rollback_provenance(cfg, varargin)

% ROLLBACK_PROVENANCE rolls the provenance one step back and should
% be used whenever a FT function calls another FT function without
% the user being (or having to be) aware of this.
%
% Some examples for use
%
%   tmpcfg            = [];
%   tmpcfg.downsample = cfg.downsample;  % copy over
%   tmpcfg.smooth     = 'no';            % override the default
%   mri = ft_volumedownsample(tmpcfg, mri);
%   [cfg, mri] = rollback_provenance(tmpcfg, mri);
%
%   tmpcfg = [];
%   tmpcfg.parameter = cfg.parameter;
%   [varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
%   [cfg, varargin{:}] = rollback_provenance(cfg, varargin{:});
%
% See also FT_PREAMBLE, FT_POSTAMBLE

% Copyright (C) 2013, Robert Oostenveld
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

for i=1:(nargin-1)
  
  fn0 = fieldnames(cfg);
  fn1 = fieldnames(varargin{i}.cfg);
  % only work on the fields that are explicitly present in the cfg
  fn = intersect(fn0, fn1);
  
  % ignore the provenance fields themselves
  fn = setdiff(fn, { ...
    'callinfo'
    'checkconfig'
    'checksize'
    'debug'
    'showcallinfo'
    'trackcallinfo'
    'trackconfig'
    'trackdatainfo'
    'trackparaminfo'
    'version'
    });
  
  for j=1:length(fn)
    fprintf('updating cfg.%s\n', fn{j});
    cfg.(fn{j}) = varargin{i}.cfg.(fn{j});
  end % for all fields that overlap
  
  if isfield(varargin{i}.cfg, 'previous')
    % take it one step back
    varargin{i}.cfg = varargin{i}.cfg.previous;
  else
    % there is no previous information
    varargin{i} = rmfield(varargin{i}, 'cfg');
  end
  
end % for all input data structures

% return the updated data structures
varargout = varargin;
