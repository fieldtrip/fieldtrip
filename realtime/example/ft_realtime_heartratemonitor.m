function ft_realtime_heartratemonitor(cfg)

% FT_REALTIME_HEARTRATEMONITOR is an example realtime application for online
% detection of heart beats. It should work both for EEG and MEG.
%
% Use as
%   ft_realtime_heartratemonitor(cfg)
% with the following configuration options
%   cfg.blocksize  = number, size of the blocks/chuncks that are processed (default = 1 second)
%   cfg.channel    = cell-array, see FT_CHANNELSELECTION (default = 'all')
%   cfg.jumptoeof  = whether to skip to the end of the stream/file at startup (default = 'yes')
%   cfg.bufferdata = whether to start on the 'first or 'last' data that is available (default = 'first')
%   cfg.threshold  = value, after normalization (default = 3)
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

% Copyright (C) 2009-2015, Robert Oostenveld
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

% set the default configuration options
if ~isfield(cfg, 'dataformat'),     cfg.dataformat = [];      end % default is detected automatically
if ~isfield(cfg, 'headerformat'),   cfg.headerformat = [];    end % default is detected automatically
if ~isfield(cfg, 'eventformat'),    cfg.eventformat = [];     end % default is detected automatically
if ~isfield(cfg, 'blocksize'),      cfg.blocksize = 0.1;      end % in seconds
if ~isfield(cfg, 'threshold'),      cfg.threshold = 3;        end % in uV
if ~isfield(cfg, 'mindist'),        cfg.mindist = 0.1;        end % in seconds
if ~isfield(cfg, 'overlap'),        cfg.overlap = 0;          end % in seconds
if ~isfield(cfg, 'channel'),        cfg.channel = 'all';      end
if ~isfield(cfg, 'jumptoeof'),      cfg.jumptoeof = 'yes';    end % jump to end of file at initialization
if ~isfield(cfg, 'bufferdata'),     cfg.bufferdata = 'first'; end % first or last
if ~isfield(cfg, 'demean'),         cfg.demean = 'yes';       end % baseline correction
if ~isfield(cfg, 'detrend'),        cfg.detrend = 'no';       end
if ~isfield(cfg, 'olfilter'),       cfg.olfilter = 'no';      end % continuous online filter
if ~isfield(cfg, 'olfiltord'),      cfg.olfiltord = 4;        end
if ~isfield(cfg, 'olfreq'),         cfg.olfreq = [2 45];      end

if ~isfield(cfg, 'dataset') && ~isfield(cfg, 'header') && ~isfield(cfg, 'datafile')
  cfg.dataset = 'buffer://localhost:1972';
end

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
pad   = [];
pads  = [];
count = 0;

tpl = [];

ws_noPeaks    = warning('off', 'signal:findpeaks:noPeaks');
ws_PeakHeight = warning('off', 'signal:findpeaks:largeMinPeakHeight');

% start the timer
tic
t0 = toc;
n0 = 0;
t1 = t0;
n1 = n0;

% this will keep the heartrate timing
heartrate = [];

% these are for the feedback
close all

h1f = figure;
plot(nan);
h1a = get(h1f, 'children');
h1c = get(h1a, 'children');
set(h1f, 'Position', [010 300 560 420]);
xlabel('time (s)');
ylim([-6 6]);

h2f = figure;
plot(nan, '.');
h2a = get(h2f, 'children');
h2c = get(h2a, 'children');
set(h2f, 'Position', [580 300 560 420]);
title('heartrate');
xlabel('time (s)');
ylabel('beats per minute');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while true
  
  % determine the samples to process
  if strcmp(cfg.bufferdata, 'last')
    % determine number of samples available in buffer
    hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat, 'cache', true);
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
  % fprintf('processing segment %d from sample %d to %d\n', count, begsample, endsample);
  
  % read data segment from buffer
  dat = ft_read_data(cfg.datafile, 'header', hdr, 'dataformat', cfg.dataformat, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false, 'blocking', true);
  dat = double(dat);
  
  % fprintf('time between subsequent reads is %f seconds\n', toc-t1);
  
  % keep track of the timing
  t1 = toc;
  n1 = n1 + size(dat,2);
  
  % fprintf('read %d samples in %f seconds, realtime ratio = %f\n', n1-n0, t1-t0, ((n1-n0)/(t1-t0))/hdr.Fs);
  % fprintf('time lag %6.3f seconds\n', (n1-n0)/hdr.Fs - (t1-t0));
  % fprintf('estimated sampling rate %.2f Hz\n', (n1-n0)/(t1-t0));
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % from here onward it is specific to the display of the data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  time   = offset2time(begsample, hdr.Fs, endsample-begsample+1);
  sample = begsample:endsample;
  
  % apply some preprocessing options
  if strcmp(cfg.demean, 'yes')
    dat = ft_preproc_baselinecorrect(dat);
  end
  if strcmp(cfg.detrend, 'yes')
    dat = ft_preproc_detrend(dat);
  end
  
  if strcmp(cfg.olfilter, 'yes')
    if ~exist('FM', 'var')
      % initialize online filter
      if cfg.olfreq(1)==0
        fprintf('using online low-pass filter\n');
        [B, A] = butter(cfg.olfiltord, cfg.olfreq(2)/hdr.Fs);
      elseif cfg.olfreq(2)>=hdr.Fs/2
        fprintf('using online high-pass filter\n');
        [B, A] = butter(cfg.olfiltord, cfg.olfreq(1)/hdr.Fs, 'high');
      else
        fprintf('using online band-pass filter\n');
        [B, A] = butter(cfg.olfiltord, cfg.olfreq/hdr.Fs);
      end
      % use one sample to initialize
      FM = ft_preproc_online_filter_init(B, A, dat(:,1));
    end
    % apply online filter
    [FM, dat] = ft_preproc_online_filter_apply(FM, dat);
  end
  
  if ~exist('CM', 'var') && numel(heartrate)>10
    % initialize online filter
    tpl = nan(10, round(0.6*hdr.Fs)+1);
    for i=2:11 % skip the first
      begsample = round((heartrate(i)-0.25)*hdr.Fs);
      endsample = round((heartrate(i)+0.35)*hdr.Fs);
      tpl(i-1,:) = ft_read_data(cfg.datafile, 'header', hdr, 'dataformat', cfg.dataformat, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false);
    end
    % apply some preprocessing options
    if strcmp(cfg.demean, 'yes')
      tpl = ft_preproc_baselinecorrect(tpl);
    end
    if strcmp(cfg.detrend, 'yes')
      tpl = ft_preproc_detrend(tpl);
    end
    tpl = median(tpl);
    tpl = tpl-mean(tpl);
    tpl = tpl/sqrt(max(filter(fliplr(tpl), 1, tpl)));
    
    B  = fliplr(tpl);
    A  = 1;
    CM = ft_preproc_online_filter_init(B, A, dat(:,1));
    
    % reset standardization
    prevState = [];
  end
  
  if exist('CM', 'var')
    % apply online filter
    [CM, dat] = ft_preproc_online_filter_apply(CM, dat);
  end
  
  [dat, prevState] = ft_preproc_standardize(dat, [], [], prevState);
  
  if isempty(pad)
    % this only applies to the first segment being processed
  else
    dat    = [pad  dat];
    sample = [pads sample];
    time   = sample./hdr.Fs;
  end
  
  % remember the last few samples, to be used for padding the next segment
  pad  = dat(:,[end-1 end-0]);
  pads = sample([end-1 end-0]);
  
  if cfg.threshold<0
    % detect negative peaks
    [peakval, peakind] = findpeaks(-dat, 'minpeakheight', -cfg.threshold);
    peakval = -peakval;
  else
    % detect positive peaks
    [peakval, peakind] = findpeaks(dat, 'minpeakheight', cfg.threshold);
  end
  
  if numel(peakind)/(blocksize/hdr.Fs)>3
    % heartrate cannot be above 180 bpm
    warning('skipping due to noise');
    peakval = [];
    peakind = [];
  end
  
  % FIXME having the heartrate vector growing is not a very good idea
  if isempty(tpl)
    heartrate = [heartrate time(peakind)];
  else
    % correct for shift due to the crosscorrelation
    heartrate = [heartrate time(peakind)-0.3];
  end
  
  if ishandle(h1f)
    set(h1c, 'xdata', time, 'ydata', dat);
    set(h1a, 'xlim', time([1 end]));
  end
  
  if numel(heartrate)>2 && ishandle(h2f)
    set(h2c, 'xdata', heartrate(2:end), 'ydata', 60./diff(heartrate));
    set(h2a, 'xlim', heartrate([2 end]));
    set(h2a, 'ylim', [0 160]);
  end
  
  % force Matlab to redraw the figures
  drawnow
  
  delay = peakind/hdr.Fs;            % adjust for the timing of the event in the data block
  delay = round(1000*delay)/1000;    % only milisecond precision
  for i=1:length(delay)
    % a heartbeat that is detected at the start should be played back immediately
    % a heartbeat that is detected at the end should be delayed a bit
    % start(timer('TimerFcn', @feedback_beep, 'ExecutionMode', 'singleShot', 'StartDelay', delay(i)));
  end
  
end % while true

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that gives a beep as feedback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function feedback_beep(varargin)
persistent beep time
if isempty(beep)
  beep = sin(1000*2*pi*(1:800)/8192);
end
soundsc(beep);
if isempty(time)
  time = toc;
end
% fprintf('beep interval = %f\n', toc-time);
% disp(varargin{2}.Data);
time = toc;

% delete the timer object that called this function
stop(varargin{1});
delete(varargin{1});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [time] = offset2time(offset, fsample, nsamples)
offset   = double(offset);
nsamples = double(nsamples);
time = (offset + (0:(nsamples-1)))/fsample;
