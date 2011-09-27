function [data] = ft_spiketriggeredinterpolation(cfg, data)

% FT_SPIKETRIGGEREDINTERPOLATION interpolates the data in the LFP channels
% around the spikes that are detected in the spike channels, or replaces
% the LFP around the spike with NaNs.
%
% Use as
%   [data] = ft_spiketriggeredinterpolation(cfg, data)
%
% The input data should be organised in a structure as obtained from the
% FT_PREPROCESSING function. The configuration should be according to
%
%   cfg.method       = string, The interpolation method can be 'nan',
%                     'cubic', 'linear', 'nearest', spline', 'pchip'
%                     (default = 'nan'). See INTERP1 for more details.
%   cfg.timwin       = [begin end], time around each spike (default = [-0.001 0.002])
%   cfg.interptoi    = value, time in seconds used for interpolation, which
%                      must be larger than timwin (default = 0.01)
%   cfg.spikechannel = string, name of single spike channel to trigger on
%   cfg.channel      = Nx1 cell-array with selection of channels (default = 'all'),
%                      see FT_CHANNELSELECTION for details
%   cfg.feedback     = 'no', 'text', 'textbar', 'gui' (default = 'no')
%
% The output will contain all channels of the input, only the data in the
% selected channels will be interpolated or replaced with NaNs.
%
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_SPIKETRIGGEREDSPECTRUM, FT_SPIKETRIGGEREDAVERAGE

% Copyright (C) 2008, Thilo Womelsdorf
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

% record start time and total processing time
ftFuncTimer = tic();
ftFuncClock = clock();
ftFuncMem   = memtic();

% enable configuration tracking
cfg = ft_checkconfig(cfg, 'trackconfig', 'on');

% set the defaults
if ~isfield(cfg, 'timwin'),         cfg.timwin = [-0.001 0.002];    end
if ~isfield(cfg, 'channel'),        cfg.channel = 'all';            end
if ~isfield(cfg, 'method'),         cfg.method = 'nan';             end
if ~isfield(cfg, 'spikechannel'),   cfg.spikechannel = [];          end
if ~isfield(cfg, 'outputexamples'), cfg.outputexamples = 'no';      end
if ~isfield(cfg, 'feedback'),       cfg.feedback = 'no';            end
if ~isfield(cfg, 'inputfile'),      cfg.inputfile = [];             end
if ~isfield(cfg, 'outputfile'),     cfg.outputfile = [];            end
if strcmp(cfg.method, 'nan')
  cfg.interptoi = 0;
else
  cfg.interptoi = 0.010;
end

% load optional given inputfile as data
hasdata = (nargin>1);
if ~isempty(cfg.inputfile)
  % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    data = loadvar(cfg.inputfile, 'data');
  end
end

% autodetect the spike channels
ntrial = length(data.trial);
nchans  = length(data.label);
spikechan = zeros(nchans,1);
for i=1:ntrial
  for j=1:nchans
    spikechan(j) = spikechan(j) + all(data.trial{i}(j,:)==0 | data.trial{i}(j,:)==1 | data.trial{i}(j,:)==2);
  end
end
spikechan = (spikechan==ntrial);

% determine the channels for interpolation
cfg.channel = ft_channelselection(cfg.channel, data.label);
chansel     = match_str(data.label, cfg.channel);

% determine the spike channel on which will be triggered
cfg.spikechannel = ft_channelselection(cfg.spikechannel, data.label);
spikesel         = match_str(data.label, cfg.spikechannel);
nspikesel        = length(cfg.spikechannel);    % number of spike channels

if nspikesel==0
  error('no spike channel selected');
end

if nspikesel>1
  error('only supported for a single spike channel');
end

if ~spikechan(spikesel)
  error('the selected spike channel seems to contain continuous data');
end

% this determines the segment that will be replaced around each spike
begpad = round(cfg.timwin(1)*data.fsample);
endpad = round(cfg.timwin(2)*data.fsample);

% this determines the segment that is used for the inperpolation around each spike
interppad = round( cfg.interptoi*data.fsample);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolate the LFP around the spikes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this is for storing the examples
cnte = 0;
spkexample = {};

for i=1:ntrial
  spikesmp = find(data.trial{i}(spikesel,:));

  fprintf('processing trial %d of %d (%d spikes)\n', i, ntrial, length(spikesmp));

  ft_progress('init', cfg.feedback, 'interpolating spikes');
  for j=1:length(spikesmp)
    ft_progress(i/ntrial, 'interpolating spike %d of %d\n', j, length(spikesmp));
    begsmp = spikesmp(j) + begpad;
    endsmp = spikesmp(j) + endpad;

    begsmp_interp = begsmp - interppad;
    endsmp_interp = endsmp + interppad;

    if begsmp_interp<1
      continue,
    end
    if endsmp_interp>size(data.trial{i},2)
      continue,
    end

    if strcmp(cfg.method,'nan')
      % only replace with NaNs
      data.trial{i}(chansel,begsmp:endsmp) = NaN;

    else
      % interpolate the data around the spike
      xall  = [begsmp_interp          : endsmp_interp];
      x     = [begsmp_interp:begsmp-1   endsmp+1:endsmp_interp];
      y     =  data.trial{i}(chansel,x) ;
      yi    = interp1(x,y,xall,cfg.method);

      % store the interpolated segment back in the data
      data.trial{i}(chansel,xall) = yi;

      if strcmp(cfg.outputexamples, 'yes') && (cnte<100)
        yall = data.trial{i}(chansel,xall);
        cnte = cnte+1;
        spkexample{cnte} = [ xall; yall; yi];
        % plot(x,y,'r.',xall,yall,'bo',xall,yi,'g-d')
      end

    end % if strcmp(cfg.method)

  end % for each spike in this trial
  ft_progress('close');

end % for each trial

if strcmp(cfg.outputexamples, 'yes')
  data.spkexample = spkexample;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the data has been modified on the fly, only update the configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
try, cfg.previous = data.cfg; end

% remember the exact configuration details in the output
data.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'data', data); % use the variable name "data" in the output file
end

