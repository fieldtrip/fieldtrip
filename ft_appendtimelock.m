function [timelock] = ft_appendtimelock(cfg, varargin)

% FT_APPENDTIMELOCK concatenates multiple timelock (ERP/ERF) data
% structures that have been processed seperately. If the input data
% structures contain different channels, it will be concatenated along the
% channel direction. If the channels are identical in the input data
% structures, the data will be concatenated along the repetition dimension.
%
% Use as
%   combined = ft_appendtimelock(cfg, timelock1, timelock2, ...)
%
% See also FT_TIMELOCKANALYSIS, FT_APPENDDATA, FT_APPENDFREQ, FT_APPENDSOURCE

% Copyright (C) 2011, Robert Oostenveld
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble defaults
ft_preamble callinfo
ft_preamble trackconfig
ft_preamble loadvar varargin

% check if the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'timelock', 'feedback', 'yes', 'hassampleinfo', 'ifmakessense');
end

% set the defaults
if ~isfield(cfg, 'inputfile'),    cfg.inputfile  = [];          end
if ~isfield(cfg, 'outputfile'),   cfg.outputfile = [];          end

% use a helper function to select the consistent parts of the data and to concatenate it
timelock = ft_selectdata(varargin{:}, 'param', {'avg' 'trial' 'cov' 'var' 'dof'});

% add version information to the configuration
cfg.version.name = mfilename('fullpath');
cfg.version.id = '$Id$';

% add information about the Matlab version used to the configuration
cfg.callinfo.matlab = version();
  
% add information about the function call to the configuration
cfg.callinfo.proctime = toc(ftFuncTimer);
cfg.callinfo.procmem  = memtoc(ftFuncMem);
cfg.callinfo.calltime = ftFuncClock;
cfg.callinfo.user = getusername();
fprintf('the call to "%s" took %d seconds and an estimated %d MB\n', mfilename, round(cfg.callinfo.proctime), round(cfg.callinfo.procmem/(1024*1024)));

% remember the configuration details of the input data
cfg.previous = cell(1,length(varargin));
for i=1:numel(varargin)
  try, cfg.previous{i} = varargin{i}.cfg; end
end

% remember the exact configuration details in the output
timelock.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'timelock', timelock); % use the variable name "data" in the output file
end

