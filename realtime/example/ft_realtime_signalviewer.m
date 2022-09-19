function ft_realtime_signalviewer(cfg)

% FT_REALTIME_SIGNALVIEWER is an example realtime application for online viewing of
% the data. It should work both for EEG and MEG.
%
% Use as
%   ft_realtime_signalviewer(cfg)
% with the following configuration options
%   cfg.blocksize  = number, size of the blocks/chuncks that are processed (default = 1 second)
%   cfg.channel    = cell-array, see FT_CHANNELSELECTION (default = 'all')
%   cfg.jumptoeof  = whether to skip to the end of the stream/file at startup (default = 'yes')
%   cfg.bufferdata = whether to start on the 'first or 'last' data that is available (default = 'first')
%   cfg.readevent  = whether or not to copy events (default = 'no')
%   cfg.demean     = 'no' or 'yes', whether to apply baseline correction (default = 'yes')
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
% Some notes about skipping data and catching up with the data stream:
%
% cfg.jumptoeof='yes' causes the realtime function to jump to the end when the
% function _starts_. It causes all data acquired prior to starting the realtime
% function to be skipped.
%
% cfg.bufferdata='last' causes the realtime function to jump to the last available data
% while _running_. If the realtime loop is not fast enough, it causes some data to be
% dropped.
%
% If you want to skip all data that was acquired before you start the RT function,
% but don't want to miss any data that was acquired while the realtime function is
% started, then you should use jumptoeof=yes and bufferdata=first. If you want to
% analyze data from a file, then you should use jumptoeof=no and bufferdata=first.
%
% To stop this realtime function, you have to press Ctrl-C

% Copyright (C) 2008, Robert Oostenveld
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
cfg.dataformat   = ft_getopt(cfg, 'dataformat',   []);      % default is detected automatically
cfg.headerformat = ft_getopt(cfg, 'headerformat', []);      % default is detected automatically
cfg.eventformat  = ft_getopt(cfg, 'eventformat',  []);      % default is detected automatically
cfg.blocksize    = ft_getopt(cfg, 'blocksize',    1);       % in seconds
cfg.overlap      = ft_getopt(cfg, 'overlap',      0);       % in seconds
cfg.channel      = ft_getopt(cfg, 'channel',      'all');
cfg.readevent    = ft_getopt(cfg, 'readevent',    'no');    % capture events?
cfg.bufferdata   = ft_getopt(cfg, 'bufferdata',   'first'); % first or last
cfg.jumptoeof    = ft_getopt(cfg, 'jumptoeof',    'yes');   % jump to end of file at initialization
cfg.demean       = ft_getopt(cfg, 'demean',       'yes');   % baseline correction
cfg.detrend      = ft_getopt(cfg, 'detrend',      'no');
cfg.olfilter     = ft_getopt(cfg, 'olfilter',     'no');    % continuous online filter
cfg.olfiltord    = ft_getopt(cfg, 'olfiltord',    4);
cfg.olfreq       = ft_getopt(cfg, 'olfreq',       [2 45]);
cfg.offset       = ft_getopt(cfg, 'offset',       []);      % in units of the data, e.g. uV for the OpenBCI board
cfg.dftfilter    = ft_getopt(cfg, 'dftfilter',    'no');
cfg.dftfreq      = ft_getopt(cfg, 'dftfreq',      [50 100 150]);
cfg.ylim         = ft_getopt(cfg, 'ylim',         []);

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
  ft_error('no channels were selected');
end

if numel(cfg.offset)==0
  % it will be determined on the first data segment
elseif numel(cfg.offset)==1
  cfg.offset = repmat(cfg.offset, size(cfg.channel));
end

% determine the size of blocks to process
blocksize = round(cfg.blocksize * hdr.Fs);
overlap   = round(cfg.overlap*hdr.Fs);

if strcmp(cfg.jumptoeof, 'yes')
  prevSample = hdr.nSamples * hdr.nTrials;
else
  prevSample = 0;
end
count = 0;

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
    endsample  = prevSample+blocksize;
  else
    ft_error('unsupported value for cfg.bufferdata');
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
  dat = ft_read_data(cfg.datafile, 'header', hdr, 'dataformat', cfg.dataformat, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false, 'blocking', true);

  % make a matching time axis
  time = ((begsample:endsample)-1)/hdr.Fs;
  
  % it only makes sense to read those events associated with the currently processed data
  if strcmp(cfg.readevent, 'yes')
    evt = ft_read_event(cfg.eventfile, 'header', hdr, 'minsample', begsample, 'maxsample', endsample);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % from here onward it is specific to the display of the data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % convert the data to a FieldTrip-like raw structure
  %   data          = [];
  %   data.trial{1} = double(dat);
  %   data.time{1}  = time;
  %   data.label    = hdr.label(chanindx);
  %   data.hdr      = hdr;
  %   data.fsample  = hdr.Fs;
  
  % apply some preprocessing options
  if strcmp(cfg.demean, 'yes')
    % demean using the first sample
    dat = ft_preproc_baselinecorrect(dat, 1, 1);
  end
  if strcmp(cfg.detrend, 'yes')
    dat = ft_preproc_detrend(dat);
  end
  
  if strcmp(cfg.dftfilter, 'yes')
    dat = ft_preproc_dftfilter(dat, hdr.Fs, cfg.dftfreq);
  end
    
  if strcmp(cfg.olfilter, 'yes')
    if count==1
      if cfg.olfreq(1)==0
        fprintf('using online low-pass filter\n');
        [B, A] = butter(cfg.olfiltord, cfg.olfreq(2)/hdr.Fs);
      elseif cfg.olfreq(2)>=hdr.Fs/2
        fprintf('using online high-pass filter\n');
        [B, A] = butter(cfg.olfiltord, cfg.olfreq(1)/hdr.Fs, 'high');
      else
        fprintf('using online band-pass filter\n');
        [B, A] = butter(cfg.olfiltord, cfg.olfreq/hdr.Fs);
      end      % use one sample to initialize
      FM = ft_preproc_online_filter_init(B, A, dat(:,1));
    end
    [FM, dat] = ft_preproc_online_filter_apply(FM, dat);
  end
  
  if isempty(cfg.offset)
    cfg.offset = ((1:nchan)-1) .* mean(max(abs(dat),[],2));
  end
  
  % shift each of the channels with a given offset
  nchan = size(dat,1);
  for i=1:nchan
    dat(i,:) = dat(i,:) + cfg.offset(i);
  end
  
  % plot the data
  plot(time, dat);
  xlim([time(1) time(end)]);
  if ~isempty(cfg.ylim)
    ylim(cfg.ylim);
  end
  
  if strcmp(cfg.readevent, 'yes')
    for i=1:length(evt)
      % draw a line and some text to indicate the event
      time = offset2time(evt(i).sample, hdr.Fs, 1);
      if ischar(evt(i).type) && isempty(evt(i).type)
        description = sprintf('%s', evt(i).type);
      elseif ischar(evt(i).type) && ischar(evt(i).type)
        description = sprintf('%s %s', evt(i).type, evt(i).value);
      elseif ischar(evt(i).type) && isnumeric(evt(i).type)
        description = sprintf('%s %s', evt(i).type, num2str(evt(i).value));
      else
        description = 'event';
      end
      
      h = line([time time], ylim);
      set(h, 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k');
      y = ylim; y = y(1);
      h = text(time, y, description, 'VerticalAlignment', 'bottom');
    end
  end
  
  % force Matlab to update the figure
  drawnow
  
end % while true

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [time] = offset2time(offset, fsample, nsamples)
offset   = double(offset);
nsamples = double(nsamples);
time = (offset + (0:(nsamples-1)))/fsample;
