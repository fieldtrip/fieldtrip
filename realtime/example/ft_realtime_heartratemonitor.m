function ft_realtime_heartratemonitor(cfg)

% FT_REALTIME_HEARTRATEMONITOR is an example realtime application for online
% viewing of the data. It should work both for EEG and MEG.
%
% Use as
%   ft_realtime_heartratemonitor(cfg)
% with the following configuration options
%   cfg.blocksize  = number, size of the blocks/chuncks that are processed (default = 1 second)
%   cfg.channel    = cell-array, see FT_CHANNELSELECTION (default = 'all')
%   cfg.bufferdata = whether to start on the 'first or 'last' data that is available (default = 'last')
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
if ~isfield(cfg, 'blocksize'),      cfg.blocksize = 0.1;      end % in seconds
if ~isfield(cfg, 'threshold'),      cfg.threshold = -4;       end % in seconds
if ~isfield(cfg, 'mindist'),        cfg.mindist = 0.1;        end % in seconds
if ~isfield(cfg, 'overlap'),        cfg.overlap = 0;          end % in seconds
if ~isfield(cfg, 'channel'),        cfg.channel = 'all';      end
if ~isfield(cfg, 'bufferdata'),     cfg.bufferdata = 'last';  end % first or last
if ~isfield(cfg, 'jumptoeof'),      cfg.jumptoeof = 'no';     end % jump to end of file at initialization

% translate dataset into datafile+headerfile
cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');
cfg = ft_checkconfig(cfg, 'required', {'datafile' 'headerfile'});

% ensure that the persistent variables related to caching are cleared
clear ft_read_header
% start by reading the header from the realtime buffer
hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat, 'cache', true, 'retry', true);

% define a subset of channels for reading
cfg.channel = ft_channelselection(cfg.channel, hdr.label);
chanindx    = match_str(hdr.label, cfg.channel);
nchan       = length(chanindx);
if nchan==0
  error('no channels were selected');
elseif nchan>1
  error('this function expects that you select a single channel');
end

% determine the size of blocks to process
blocksize = round(cfg.blocksize * hdr.Fs);
overlap   = round(cfg.overlap*hdr.Fs);

if strcmp(cfg.jumptoeof, 'yes')
  prevSample = hdr.nSamples * hdr.nTrials;
else
  prevSample = 0;
end

prevState = [];
pad = [];
count = 0;

ws_noPeaks = warning('off', 'signal:findpeaks:noPeaks');
ws_PeakHeight = warning('off', 'signal:findpeaks:largeMinPeakHeight');

% start the timer
tic
t0 = toc;
n0 = 0;
t1 = t0;
n1 = n0;

% these are for the feedback
beat = [];
beep = sin(1000*2*pi*(1:800)/8192);
close all
h1 = figure
set(h1, 'Position', [010 300 560 420]);
h2 = figure
set(h2, 'Position', [580 300 560 420]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while true

  % determine number of samples available in buffer
  hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat, 'cache', true);

  % see whether new samples are available
  newsamples = (hdr.nSamples*hdr.nTrials-prevSample);

  if newsamples>=blocksize

    % determine the samples to process
    if strcmp(cfg.bufferdata, 'last')
      begsample  = hdr.nSamples*hdr.nTrials - blocksize + 1;
      endsample  = hdr.nSamples*hdr.nTrials;
    elseif strcmp(cfg.bufferdata, 'first')
      begsample  = prevSample+1;
      endsample  = prevSample+blocksize ;
    else
      error('unsupported value for cfg.bufferdata');
    end

    % this allows overlapping data segments
    if overlap && (begsample>overlap)
      begsample = begsample - overlap;
      endsample = endsample - overlap;
    end

    % remember up to where the data was read
    prevSample  = endsample;
    count       = count + 1;
    fprintf('processing segment %d from sample %d to %d\n', count, begsample, endsample);

    % read data segment from buffer
    dat = ft_read_data(cfg.datafile, 'header', hdr, 'dataformat', cfg.dataformat, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false);

    % keep track of the timing
    t1 = toc;
    n1 = n1 + size(dat,2);
    fprintf('read %d samples in %f seconds, realtime ratio = %f\n', n1-n0, t1-t0, ((n1-n0)/(t1-t0))/hdr.Fs);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % from here onward it is specific to the display of the data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    time   = offset2time(begsample, hdr.Fs, endsample-begsample+1);
    sample = begsample:endsample;
    label  = hdr.label(chanindx);

    % apply some preprocessing options
    [dat, prevState] = ft_preproc_standardize(dat, [], [], prevState);

    if isempty(pad)
      % this only applies to the first segment being processed
      pad = dat(:,1);
    else
      % remember the last sample for padding the next segment
      pad = dat(:,end);
    end

    dat    = [pad dat];
    sample = [sample(1)-1 sample];
    time   = [time(1)-1/hdr.Fs time];

    if cfg.threshold<0
      % detect negative peaks
      [peakval, peakind] = findpeaks(-dat, 'minpeakheight', -cfg.threshold);
      peakval = -peakval;
    else
      [peakval, peakind] = findpeaks(dat, 'minpeakheight', cfg.threshold);
    end
    % the last sample should not be detected as peak but should be carried on to the next segment as padding
    sel = (peakind==size(dat,2));
    peakval(sel) = [];
    peakind(sel) = [];

    if true
      figure(h1)
      % plot the data
      plot(time, dat);
      xlim([time(1) time(end)]);
      ylim([-6 6]);
    end

    for i=1:length(peakval)
      % make a sound for each heart beat
      soundsc(beep);
      % FIXME the beat vector growing is a bad idea
      beat(end+1) = time(peakind(i));
      figure(h2)
      plot(beat(2:end), diff(beat)*60, 'b.');
      title('heartrate');
      xlabel('time (s)');
      ylabel('beats per minute');
    end

    % force Matlab to update the figure
    drawnow

  end % if enough new samples
end % while true

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [time] = offset2time(offset, fsample, nsamples)
offset   = double(offset);
nsamples = double(nsamples);
time = (offset + (0:(nsamples-1)))/fsample;