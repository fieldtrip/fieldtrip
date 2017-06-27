function [timelock] = ft_appendtimelock(cfg, varargin)

% FT_APPENDTIMELOCK concatenates multiple timelock (ERP/ERF) data structures that
% have been processed seperately. If the input data structures contain different
% channels, it will be concatenated along the channel direction. If the channels are
% identical in the input data structures, the data will be concatenated along the
% repetition dimension.
%
% Use as
%   combined = ft_appendtimelock(cfg, timelock1, timelock2, ...)
%
% The configuration can optionally contain
%   cfg.appenddim  = string, the dimension to concatenate over which to append,
%                    this can be 'chan' and 'rpt' (default is automatic)
%   cfg.tolerance  = scalar, tolerance to determine how different the time axes
%                    are allowed to still be considered compatible (default = 1e-5)
%
% See also FT_TIMELOCKANALYSIS, FT_DATATYPE_TIMELOCK, FT_APPENDDATA, FT_APPENDFREQ,
% FT_APPENDSENS, FT_APPENDSOURCE

% Copyright (C) 2011-2017, Robert Oostenveld
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
ft_preamble debug
ft_preamble loadvar    varargin
ft_preamble provenance varargin
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
for i=1:length(varargin)
  % FIXME: what about timelock+comp?
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', {'timelock', 'timelock+comp'}, 'feedback', 'yes', 'hassampleinfo', 'ifmakessense');
end

% set the defaults
cfg.channel    = ft_getopt(cfg, 'channel', 'all');
cfg.parameter  = ft_getopt(cfg, 'parameter', []);
cfg.appenddim  = ft_getopt(cfg, 'appenddim', []);
cfg.tolerance  = ft_getopt(cfg, 'tolerance',  1e-5);
cfg.appendsens = ft_getopt(cfg, 'appendsens', 'no');

if isempty(cfg.appenddim) || strcmp(cfg.appenddim, 'auto')
  if checkchan(varargin{:}, 'identical') && checktime(varargin{:}, 'identical', cfg.tolerance)
    cfg.appenddim = 'rpt';
  elseif checktime(varargin{:}, 'identical', cfg.tolerance) && ~checkchan(varargin{:}, 'unique')
    cfg.appenddim = 'rpt';
  elseif checkchan(varargin{:}, 'unique')
    cfg.appenddim = 'chan';
  elseif checktime(varargin{:}, 'unique', cfg.tolerance)
    cfg.appenddim = 'time';
  else
    error('cfg.appenddim should be specified');
  end
end
fprintf('concatenating over the "%s" dimension\n', cfg.appenddim);

if isempty(cfg.parameter)
  fn = fieldnames(varargin{1});
  for i=2:numel(varargin)
    fn = intersect(fn, fieldnames(varargin{i}));
  end
  cfg.parameter = setdiff(fn, ignorefields('appendtimelock'));
elseif ischar(cfg.parameter)
  cfg.parameter = {cfg.parameter};
end
assert(~isempty(cfg.parameter), 'cfg.parameter should be specified');

if any(strcmp(cfg.parameter, 'avg')) && any(strcmp(cfg.parameter, 'trial'))
  warning('appending the individual trials, not the averages');
  % also prevent var and dof from being appended
  cfg.parameter = setdiff(cfg.parameter, {'avg', 'var', 'dof'}); 
end

% use a low-level function that is shared with the other ft_appendxxx functions
timelock = append_common(cfg, varargin{:});

if isfield(timelock, 'avg') && ~isfield(timelock, 'trial')
  warning('renaming the appended averages to "trial"');
  timelock.trial = timelock.avg;
  timelock = rmfield(timelock, 'avg');
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   varargin
ft_postamble provenance timelock
ft_postamble history    timelock
ft_postamble savevar    timelock
