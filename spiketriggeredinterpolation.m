function [data] = spiketriggeredinterpolation(cfg, data)

% SPIKETRIGGEREDINTERPOLATION interpolates the data in the LFP channels
% around the spikes that are detected in the spike channels, or replaces
% the LFP around the spike with NaNs.
%
% Use as
%   [data] = spiketriggeredinterpolation(cfg, data)
%
% The input data should be organised in a structure as obtained from the
% PREPROCESSING function. The configuration should be according to
%
%   cfg.method       = string, The interpolation method can be 'nan',
%                     'cubic', 'linear', 'nearest', spline', 'pchip'
%                     (default = 'nan'). See INTERP1 for more details.
%   cfg.timwin       = [begin end], time around each spike (default = [-0.001 0.002])
%   cfg.interptoi    = value, time in seconds used for interpolation, which
%                      must be larger than timwin (default = 0.01)
%   cfg.spikechannel = string, name of single spike channel to trigger on
%   cfg.channel      = Nx1 cell-array with selection of channels (default = 'all'),
%                      see CHANNELSELECTION for details
%   cfg.feedback     = 'no', 'text', 'textbar', 'gui' (default = 'no')
%
% The output will contain all channels of the input, only the data in the
% selected channels will be interpolated or replaced with NaNs.
%
% See also SPIKETRIGGEREDSPECTRUM, SPIKETRIGGEREDAVERAGE

% Copyright (C) 2008, Thilo Womelsdorf
%
% $Log: spiketriggeredinterpolation.m,v $
% Revision 1.2  2008/09/24 08:04:50  roboos
% added example output, fixed some bugs
%
% Revision 1.1  2008/09/18 09:50:36  roboos
% new implementation, thanks to Thilo
%

% set the defaults
if ~isfield(cfg, 'timwin'),         cfg.timwin = [-0.001 0.002];    end
if ~isfield(cfg, 'channel'),        cfg.channel = 'all';            end
if ~isfield(cfg, 'method'),         cfg.method = 'nan';             end
if ~isfield(cfg, 'spikechannel'),   cfg.spikechannel = [];          end
if ~isfield(cfg, 'outputexamples'), cfg.outputexamples = 'no';      end
if ~isfield(cfg, 'feedback'),       cfg.feedback = 'no';            end

if strcmp(cfg.method, 'nan')
  cfg.interptoi = 0;
else
  cfg.interptoi = 0.010;
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
cfg.channel = channelselection(cfg.channel, data.label);
chansel     = match_str(data.label, cfg.channel);

% determine the spike channel on which will be triggered
cfg.spikechannel = channelselection(cfg.spikechannel, data.label);
spikesel         = match_str(data.label, cfg.spikechannel);
nspikesel        = length(cfg.spikechannel);	% number of spike channels

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

  progress('init', cfg.feedback, 'interpolating spikes');
  for j=1:length(spikesmp)
    progress(i/ntrial, 'interpolating spike %d of %d\n', j, length(spikesmp));
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
  progress('close');

end % for each trial

if strcmp(cfg.outputexamples, 'yes')
  data.spkexample = spkexample;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the data has been modified on the fly, only update the configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: spiketriggeredinterpolation.m,v 1.2 2008/09/24 08:04:50 roboos Exp $';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output
data.cfg = cfg;


