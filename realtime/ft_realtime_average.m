function ft_realtime_average(cfg)

% FT_REALTIME_AVERAGE is an example realtime application for online
% averaging of the data. It should work both for EEG and MEG.
%
% Use as
%   ft_realtime_average(cfg)
% with the following configuration options
%   cfg.channel    = cell-array, see FT_CHANNELSELECTION (default = 'all')
%   cfg.trialfun   = string with the trial function
%
% The source of the data is configured as
%   cfg.dataset       = string
% or alternatively to obtain more low-level control as
%   cfg.datafile      = string
%   cfg.headerfile    = string
%   cfg.eventfile     = string
%   cfg.dataformat    = string, default is determined automatic
%   cfg.headerformat  = string, default is determined automatic
%   cfg.eventformat   = string, default is determined automatic
%
% To stop the realtime function, you have to press Ctrl-C

% Copyright (C) 2009, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% set the default configuration options
if ~isfield(cfg, 'dataformat'),     cfg.dataformat = [];      end % default is detected automatically
if ~isfield(cfg, 'headerformat'),   cfg.headerformat = [];    end % default is detected automatically
if ~isfield(cfg, 'eventformat'),    cfg.eventformat = [];     end % default is detected automatically
if ~isfield(cfg, 'channel'),        cfg.channel = 'all';      end
if ~isfield(cfg, 'bufferdata'),     cfg.bufferdata = 'last';  end % first or last

% translate dataset into datafile+headerfile
cfg = checkconfig(cfg, 'dataset2files', 'yes');
cfg = checkconfig(cfg, 'required', {'datafile' 'headerfile'});

% ensure that the persistent variables related to caching are cleared
clear read_header
% start by reading the header from the realtime buffer
hdr = read_header(cfg.headerfile, 'cache', true);

% define a subset of channels for reading
cfg.channel = channelselection(cfg.channel, hdr.label);
chanindx    = match_str(hdr.label, cfg.channel);
nchan       = length(chanindx);
if nchan==0
  error('no channels were selected');
end

prevSample = 0;
count      = 0;

% initialize the average, it will be filled on the first iteration
avgsum = [];
avgnum = [];
% open a figure in which the average will be plotted
figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while true

  % determine latest header and event information
  event     = read_event(cfg.dataset, 'minsample', prevSample+1);  % only consider events that are later than the data processed sofar
  hdr       = read_header(cfg.dataset, 'cache', true);             % the trialfun might want to use this, but it is not required
  cfg.event = event;                                               % store it in the configuration, so that it can be passed on to the trialfun
  cfg.hdr   = hdr;                                                 % store it in the configuration, so that it can be passed on to the trialfun

  % evaluate the trialfun, note that the trialfun should not re-read the events and header
  fprintf('evaluating ''%s'' based on %d events\n', cfg.trialfun, length(event));
  trl = feval(cfg.trialfun, cfg);

  fprintf('processing %d trials\n', size(trl,1));

  for trllop=1:size(trl,1)

    begsample = trl(trllop,1);
    endsample = trl(trllop,2);
    offset    = trl(trllop,3);

    % remember up to where the data was read
    prevSample  = endsample;
    count       = count + 1;
    fprintf('processing segment %d from sample %d to %d\n', count, begsample, endsample);

    % read data segment from buffer
    dat = read_data(cfg.datafile, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % from here onward it is specific to the processing of the data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % apply some preprocessing options
    dat = preproc_baselinecorrect(dat);

    if isempty(average)
      % initialize the accumulating variables on the first call
      avgsum = dat;
      avgnum = 1;
    else
      avgsum = avgsum + dat;
      avgnum = avgnum + 1;
    end

    % compute the average
    avg = avgsum ./ avgnum;

    % create a time-axis and plot the average
    time = offset2time(offset, hdr.Fs, endsample-begsample+1);
    plot(time, avg);

    % force matlab to redraw the figure
    drawnow

  end % looping over new trials
end % while true

