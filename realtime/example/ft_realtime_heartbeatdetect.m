function ft_realtime_heartbeatdetect(cfg)

% FT_REALTIME_HEARTBEATDETECT is an example realtime application for online
% detection of heart beats. It should work both for EEG and MEG.
%
% Use as
%   ft_realtime_heartbeatdetect(cfg)
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

% set the default configuration options
cfg.dataformat   = ft_getopt(cfg, 'dataformat');          % default is detected automatically
cfg.headerformat = ft_getopt(cfg, 'headerformat');        % default is detected automatically
cfg.eventformat  = ft_getopt(cfg, 'eventformat');         % default is detected automatically
cfg.blocksize    = ft_getopt(cfg, 'blocksize', 0.1);      % in seconds
cfg.threshold    = ft_getopt(cfg, 'threshold', 3);        % after normalization
cfg.mindist      = ft_getopt(cfg, 'mindist', 0.1);        % in seconds
cfg.channel      = ft_getopt(cfg, 'channel', 'all');
cfg.jumptoeof    = ft_getopt(cfg, 'jumptoeof', 'yes');    % jump to end of file at initialization
cfg.bufferdata   = ft_getopt(cfg, 'bufferdata', 'first'); % first or last
cfg.demean       = ft_getopt(cfg, 'demean', 'yes');       % baseline correction
cfg.detrend      = ft_getopt(cfg, 'detrend', 'no');
cfg.olfilter     = ft_getopt(cfg, 'olfilter', 'no');      % continuous online filter
cfg.olfiltord    = ft_getopt(cfg, 'olfiltord',  4);
cfg.olfreq       = ft_getopt(cfg, 'olfreq',  [2 45]);
cfg.dftfilter    = ft_getopt(cfg, 'dftfilter', 'yes');    % filter using discrete Fourier transform
cfg.dftfreq      = ft_getopt(cfg, 'dftfreq',  50);        % line noise frequency

if ~isfield(cfg, 'dataset') && ~isfield(cfg, 'datafile') && ~isfield(cfg, 'headerfile')
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

if strcmp(cfg.jumptoeof, 'yes')
  prevSample = hdr.nSamples * hdr.nTrials;
else
  prevSample = 0;
end

prevState = [];
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

% this will keep the time of each heart beat
heartbeat = [];

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
title('heartbeat');
xlabel('time (s)');
ylabel('beats per minute');

c = onCleanup(@cleanup_cb);

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
    endsample  = min(prevSample+blocksize, hdr.nSamples*hdr.nTrials);
    begsample  = endsample - blocksize + 1;
  else
    error('unsupported value for cfg.bufferdata');
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
  
  sample = begsample:endsample;
  time   = sample./hdr.Fs;
  
  % apply some preprocessing options
  if strcmp(cfg.demean, 'yes')
    dat = ft_preproc_baselinecorrect(dat);
  end
  if strcmp(cfg.detrend, 'yes')
    dat = ft_preproc_detrend(dat);
  end
  if strcmp(cfg.dftfilter, 'yes')
    dat = ft_preproc_dftfilter(dat, hdr.Fs, cfg.dftfreq);
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
  
  [dat, prevState] = ft_preproc_standardize(dat, [], [], prevState);
  
  if cfg.threshold<0
    % detect negative peaks
    [peakval, peakind] = findpeaks(-dat, 'minpeakheight', -cfg.threshold);
    peakval = -peakval;
  else
    % detect positive peaks
    [peakval, peakind] = findpeaks(dat, 'minpeakheight', cfg.threshold);
  end
  
  if numel(peakind)/(blocksize/hdr.Fs)>3
    % heartbeat cannot be above 180 bpm
    warning('skipping due to noise');
    peakval = [];
    peakind = [];
  end
  
  % FIXME having the heartbeat vector growing is not a very good idea
  heartbeat = [heartbeat time(peakind)];
  
  if ishandle(h1f)
    set(h1c, 'xdata', time, 'ydata', dat);
    set(h1a, 'xlim', time([1 end]));
  end
  
  if numel(heartbeat)>5 && ishandle(h2f)
    % skip the first heartbeat for the axes
    set(h2c, 'xdata', heartbeat(2:end), 'ydata', 60./diff(heartbeat));
    set(h2a, 'xlim', heartbeat([2 end]) + [0 1]);
    set(h2a, 'ylim', [0 160]);
  end
  
  %   if numel(heartbeat)>3
  %     event.type = 'heartrate';
  %     event.value = heartbeat(end) - heartbeat(end-1);
  %     event.sample = [];
  %     event.offset = 0;
  %     event.duration = 0;
  %     ft_write_event(cfg.dataset, event);
  %   end
  
  % force Matlab to redraw the figures
  drawnow
  
end % while true


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that gives a beep as feedback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function feedback_beep(varargin)
beep = audioplayer(0.05*sin(1000*2*pi*(1:1024)/8192), 8192);
play(beep);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cleanup_cb(varargin)
delete(timerfindall)

