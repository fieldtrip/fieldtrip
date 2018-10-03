function [cfg, varargout] = rollback_provenance(cfg, varargin)

% ROLLBACK_PROVENANCE rolls the provenance one step back and should
% be used whenever a FT function calls another FT function without
% the user being (or having to be) aware of this.
%
% Some examples for use
%
%   tmpcfg            = [];
%   tmpcfg.downsample = cfg.downsample;  % simply copy this option
%   tmpcfg.smooth     = 'no';            % override the default for this option
%   mri = ft_volumedownsample(tmpcfg, mri);
%   [cfg, mri] = rollback_provenance(cfg, mri);
%
%   tmpcfg           = [];
%   tmpcfg.parameter = cfg.parameter;
%   [varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
%   [cfg, varargin{:}] = rollback_provenance(cfg, varargin{:});
%
% See also FT_PREAMBLE, FT_POSTAMBLE

% Copyright (C) 2013-2014, Robert Oostenveld
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

for i=1:(nargin-1)
  
  if ~isfield(varargin{i}, 'cfg')
    % nothing to do
    continue
  end
  
  if isempty(cfg)
    % allow for [] as the input cfg
    fn0 = {};
  else
    % get the input cfg
    fn0 = fieldnames(cfg);
  end
  
  if ~isfield(varargin{i}, 'cfg')
    % input does not contain cfg, so no rollback to be performed
    continue;
  else
    % get the data cfg
    fn1 = fieldnames(varargin{i}.cfg);
  end
  
  % work on the fields that are present in both input cfg and data cfg
  fn = intersect(fn0, fn1);
  
  % some of the fields should not be rolled back
  fn = setdiff(fn, ignorefields('rollback'));
  
  for j=1:length(fn)
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
