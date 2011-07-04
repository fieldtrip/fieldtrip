function [freq] = ft_appendfreq(cfg, varargin)

% FT_APPENDFREQ concatenates multiple frequency or time-frequency data
% structures that have been processed seperately. If the input data
% structures contain different channels, it will be concatenated along the
% channel direction. If the channels are identical in the input data
% structures, the data will be concatenated along the repetition dimension.
%
% Use as
%   combined = ft_appendfreq(cfg, freq1, freq2, ...)
%
% See also FT_FREQANALYSIS, FT_APPENDDATA, FT_APPENDTIMELOCK, FT_APPENDSOURCE

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

ft_defaults

% record start time and total processing time
ftFuncTimer = tic();
ftFuncClock = clock();

% set the defaults
if ~isfield(cfg, 'inputfile'),    cfg.inputfile  = [];          end
if ~isfield(cfg, 'outputfile'),   cfg.outputfile = [];          end

hasdata = nargin>1;
if ~isempty(cfg.inputfile) % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  elseif ~iscell(cfg.inputfile)
    error('you should specify cfg.inpoutfile as cell-array with multiple file names');
  else
    for i=1:numel(cfg.inputfile)
      varargin{i} = loadvar(cfg.inputfile{i}, 'freq'); % read datasets from array inputfile
    end
  end
end

% check if the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'freq', 'feedback', 'yes', 'hassampleinfo', 'ifmakessense');
end

% use a helper function to select the consistent parts of the data and to concatenate it
freq = ft_selectdata(varargin{:}, 'param', {'powspctrm' 'crsspctrm' 'fourierspctrm'});

% add version information to the configuration
cfg.version.name = mfilename('fullpath');
cfg.version.id = '$Id$';

% add information about the Matlab version used to the configuration
cfg.callinfo.matlab = version();
  
% add information about the function call to the configuration
cfg.callinfo.proctime = toc(ftFuncTimer);
cfg.callinfo.calltime = ftFuncClock;
cfg.callinfo.user = getusername();

% remember the configuration details of the input data
cfg.previous = cell(1,length(varargin));
for i=1:numel(varargin)
  try, cfg.previous{i} = varargin{i}.cfg; end
end

% remember the exact configuration details in the output
freq.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'freq', freq); % use the variable name "data" in the output file
end

