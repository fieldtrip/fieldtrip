function [montage] = ft_prepare_montage(cfg, data)

% FT_PREPARE_MONTAGE creates a referencing scheme based on the input
% configuration options and the channels in the data structure. The
% resulting montage can be given as input to ft_apply_montage, or as
% cfg.montage to ft_preprocessing.
%
% Use as
%   [montage] = ft_prepare_montage(cfg, data)
%
% Referencing options:
%   cfg.reref         = 'no' or 'yes' (default = 'yes')
%   cfg.refchannel    = cell-array with new EEG reference channel(s), this can be 'all' for a common average reference
%   cfg.implicitref   = 'label' or empty, add the implicit EEG reference as zeros (default = [])
%
% See also FT_APPLY_MONTAGE and FT_PREPROCESSING

% Copyright (C) 2017, Arjen Stolk & Robert Oostenveld
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


% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble loadvar    data
ft_preamble provenance data

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% set the defaults
cfg.reref          = ft_getopt(cfg, 'reref', 'yes');
cfg.refchannel     = ft_getopt(cfg, 'refchannel', 'all');
cfg.implicitref    = ft_getopt(cfg, 'implicitref', []);

% check of input arguments
if strcmp(cfg.reref, 'no')
  error('cfg.reref should be set to ''yes'' in order to create a montage')
end

% create the refchannel-dependent montage
montage.labelold = data.label;
montage.labelnew = data.label;
montage.tra      = eye(numel(data.label)); % a simple identity matrix
if ~strmp(cfgrefchannel, 'all')
  montage.tra(match_str(data.label, cfg.refchannel),:) = []; % fill the refchannel column(s) with -1
end

% append implicitref if specified
if ~isempty(cfg.implicitref)  
  match_str(data.label, cfg.implicitref)
  montage.tra(match_str(data.label, cfg.implicitref),:) = 0; % fill the implicitref row with 0
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble provenance
ft_postamble previous data
ft_postamble history montage