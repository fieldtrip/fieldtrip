function ft_realtime_ouunpod(cfg)

% FT_REALTIME_OUUNPOD is an example realtime application for online power
% estimation and visualisation. It is designed for use with the OuUnPod, an
% OpenEEG based low cost EEG system with two channels, but in principle
% should work for any EEG or MEG system.
%
% Use as
%   ft_realtime_ouunpod(cfg)
% with the following configuration options
%   cfg.channel    = cell-array, see FT_CHANNELSELECTION (default = 'all')
%   cfg.foilim     = [Flow Fhigh] (default = [1 45])
%   cfg.blocksize  = number, size of the blocks/chuncks that are processed (default = 1 second)
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
%
% See also http://ouunpod.blogspot.com

% Copyright (C) 2008-2012, Robert Oostenveld
% Copyright (C) 2012-2014, Stephen Whitmarsh
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
if ~isfield(cfg, 'dataformat'),       cfg.dataformat      = [];       end % default is detected automatically
if ~isfield(cfg, 'headerformat'),     cfg.headerformat    = [];       end % default is detected automatically
if ~isfield(cfg, 'eventformat'),      cfg.eventformat     = [];       end % default is detected automatically
if ~isfield(cfg, 'blocksize'),        cfg.blocksize       = 0.05;     end % stepsize, in seconds
if ~isfield(cfg, 'channel'),          cfg.channel         = 'all';    end
if ~isfield(cfg, 'bufferdata'),       cfg.bufferdata      = 'last';   end % first or last
if ~isfield(cfg, 'dataset'),          cfg.dataset         = 'buffer:\\localhost:1972'; end
if ~isfield(cfg, 'foilim'),           cfg.foilim          = [1 45];   end
if ~isfield(cfg, 'windowsize'),       cfg.windowsize      = 2;        end % length of sliding window, in seconds
if ~isfield(cfg, 'scale'),            cfg.scale           = 1;        end % can be used to fix the calibration
if ~isfield(cfg, 'feedback'),         cfg.feedback        = 'no';     end % use neurofeedback with MIDI, yes or no

% translate dataset into datafile+headerfile
cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');
cfg = ft_checkconfig(cfg, 'required', {'datafile' 'headerfile'});

if strcmp(cfg.feedback, 'yes')
  % setup MIDI, see http://en.wikipedia.org/wiki/General_MIDI
  beatdrum = true;
  
  m = midiOut;     % Microsoft GS Wavetable Synth = device number 2
  midiOut('O', 2); % o for output; 2 for device nr 2
  midiOut('.', 1); % all off
  midiOut('.', 2); % all off
  
  % midiOut('P', 2, 20); organ
  midiOut('P', 1, 53);
  midiOut('P', 2, 53);
  midiOut('+', 1, [64 65], [127 127]);        % command, channelnr, key, velocity
  midiOut('+', 2, [64 65 67], [127 127 127]); % command, channelnr, key, velocity
  midiOut(uint8([175+1, 7, 0]));              % change volume
  midiOut(uint8([175+2, 7, 0]));
end % if MIDI feedback

% these are used by the GUI callbacks
clear global vaxis hdr chanindx
global vaxis hdr chanindx

% this specifies the vertical axis for each of the 6 subplots
vaxis = [
  -300 300
  -300 300
  0 1000
  0 1000
  0 1000
  0 1000
  ];

b2clicked = false;

% schemerlamp = Lamp('com9');

% ensure that the persistent variables related to caching are cleared
clear ft_read_header

% start by reading the header from the realtime buffer
hdr = ft_read_header(cfg.headerfile, 'cache', true, 'retry', true);

% define a subset of channels for reading
cfg.channel = ft_channelselection(cfg.channel, hdr.label);
chanindx = match_str(hdr.label, cfg.channel);
nchan = length(chanindx);
if nchan>2
  chanindx = [1 2];
  nchan = 2;
  ft_warning('exactly two channels should be selected');
end
if nchan<2
  ft_error('exactly two channels should be selected');
end

nhistory = 100;

% determine the size of blocks to process
blocksize = round(cfg.blocksize * hdr.Fs);

prevSample = 0;
count = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

left_thresh_ampl = -100;
left_thresh_time = nan; % [cfg.blockmem*blocksize-blocksize*4 cfg.blockmem*blocksize];

% FIXME these are hardcoded, but might be incompatible with the cfg and data settings
right_freq = [40 45];
right_offset = 0.5;
right_mult = 127/0.5;

TFR = zeros(2, (cfg.foilim(2)-cfg.foilim(1)+1), 100);

f1 = nan;

% these are handles used in drawing
u1=[]; u2=[]; u3=[]; u4=[]; u5=[]; u6=[]; p1=[]; p2=[]; p3=[]; p4=[]; p5=[]; p6=[]; c1=[]; c2=[]; b2=[];

while true
  
  if isempty(f1) || ~ishandle(f1)
    close all;
    f1 = figure;
    set(f1, 'resizeFcn', 'u1=[]; u2=[]; u3=[]; u4=[]; u5=[]; u6=[]; p1=[]; p2=[]; p3=[]; p4=[]; p5=[]; p6=[]; c1=[]; c2=[]; b2=[];');
    u1=[]; u2=[]; u3=[]; u4=[]; u5=[]; u6=[]; p1=[]; p2=[]; p3=[]; p4=[]; p5=[]; p6=[]; c1=[]; c2=[]; b2=[];
  end
  
  % determine number of samples available in buffer
  hdr = ft_read_header(cfg.headerfile, 'cache', true);
  
  % see whether new samples are available
  newsamples = (hdr.nSamples*hdr.nTrials-prevSample);
  
  if newsamples>=blocksize && (hdr.nSamples*hdr.nTrials/hdr.Fs)>cfg.windowsize
    
    % determine the samples to process
    if strcmp(cfg.bufferdata, 'last')
      begsample = hdr.nSamples*hdr.nTrials - round(cfg.windowsize*hdr.Fs) + 1;
      endsample = hdr.nSamples*hdr.nTrials;
    elseif strcmp(cfg.bufferdata, 'first')
      begsample = prevSample+1;
      endsample = prevSample+blocksize;
    else
      ft_error('unsupported value for cfg.bufferdata');
    end
    
    % remember up to where the data was read
    prevSample = endsample;
    count = count + 1;
    fprintf('processing segment %d from sample %d to %d\n', count, begsample, endsample);
    
    % read the data segment from buffer
    dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false);
    dat = cfg.scale * dat;
    
    % construct a matching time axis
    time = ((begsample:endsample)-1)/hdr.Fs;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % from here onward it is specific to the power estimation from the data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % apply some preprocessing to the data
    dat = ft_preproc_polyremoval(dat, 1);
    dat = ft_preproc_highpassfilter(dat, hdr.Fs,  3, 1, 'but', 'twopass');
    dat = ft_preproc_lowpassfilter (dat, hdr.Fs, 35, 3, 'but', 'twopass');
    
    if hdr.Fs<11025
      % sampling range is low, assume it is EEG
      if hdr.Fs>110
        % apply line noise filter
        dat = ft_preproc_bandstopfilter(dat, hdr.Fs, [45 55], 4, 'but', 'twopass');
      end
      if hdr.Fs>230
        % apply line noise filter
        dat = ft_preproc_bandstopfilter(dat, hdr.Fs, [95 115], 4, 'but', 'twopass');
      end
      [spec, ntaper, freqoi] = ft_specest_mtmfft(dat, time, 'taper', 'dpss', 'tapsmofrq', 2, 'freqoi', cfg.foilim(1):cfg.foilim(2));

    else
      % sampling range is high, assume it is audio
      [spec, ntaper, freqoi] = ft_specest_mtmfft(dat, time, 'taper', 'hanning', 'freqoi', cfg.foilim(1):cfg.foilim(2));
    end
    
    pow = squeeze(mean(abs(spec.^2), 1)); % compute power, average over tapers
    
    if ~exist('TFR', 'var')
      TFR = nan(length(chanindx), length(freqoi), nhistory);
    end
    
    TFR(:,:,1:nhistory-1) = TFR(:,:,2:nhistory);
    TFR(:,:,end) = pow;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % translate channel 1 into a neurofeedback command
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if strcmp(cfg.feedback, 'yes')
      % compute the average power in the specified frequency range
      fbeg = nearest(freqoi, cfg.feedback1.foilim(1));
      fend = nearest(freqoi, cfg.feedback1.foilim(2));
      value1 = mean(pow(1, fbeg:fend));
      % scale the value between 0 and 1
      historicalmean = mean(nanmean(TFR(1,fbeg:fend,:),3),2);
      historicalmin  = min (nanmin (TFR(1,fbeg:fend,:),3),2);
      historicalmax  = max (nanmax (TFR(1,fbeg:fend,:),3),2);
      % the value can be larger than expected from the history
      value1 = (value1 - historicalmin) ./ (historicalmax - historicalmin);
      
      controlfunction(cfg.feedback1, value1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % translate channel 2 into a neurofeedback command
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if strcmp(cfg.feedback, 'yes')
      %       if value2>threshold
      %         controlfunction(cfg.feedback2);
      %       end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make the GUI elements
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
      
      if isempty(c1) || ~ishandle(c1)
        pos = [0.25 0.95 0.1 0.05];
        c1 = uicontrol('style', 'edit', 'units', 'normalized', 'callback', @update_channel, 'BackgroundColor', 'white');
        set(c1, 'position', pos);
        set(c1, 'string', chanindx(1));
        set(c1, 'tag', 'c1');
      end
      
      if isempty(c2) || ~ishandle(c2)
        pos = [0.70 0.95 0.1 0.05];
        c2 = uicontrol('style', 'edit', 'units', 'normalized', 'callback', @update_channel, 'BackgroundColor', 'white');
        set(c2, 'position', pos);
        set(c2, 'string', chanindx(2));
        set(c2, 'tag', 'c2');
      end
      
      if isempty(u1) || ~ishandle(u1)
        pos = get(p1, 'position'); % link the position to the subplot
        pos(1) = pos(1)-0.1;
        pos(2) = pos(2)-0.05;
        pos(3) = 0.1;
        pos(4) = 0.05;
        u1 = uicontrol('style', 'edit', 'units', 'normalized', 'callback', @update_axis, 'BackgroundColor', 'white');
        set(u1, 'position', pos);
        set(u1, 'string', num2str(vaxis(1,2)));
        set(u1, 'tag', 'u1');
      end
      
      if isempty(u2) || ~ishandle(u2)
        pos = get(p2, 'position'); % link the position to the subplot
        pos(1) = pos(1)-0.1;
        pos(2) = pos(2)-0.05;
        pos(3) = 0.1;
        pos(4) = 0.05;
        u2 = uicontrol('style', 'edit', 'units', 'normalized', 'callback', @update_axis, 'BackgroundColor', 'white');
        set(u2, 'position', pos);
        set(u2, 'string', num2str(vaxis(2,2)));
        set(u2, 'tag', 'u2');
      end
      
      if isempty(u3) || ~ishandle(u3)
        pos = get(p3, 'position'); % link the position to the subplot
        pos(1) = pos(1)-0.1;
        pos(2) = pos(2)-0.05;
        pos(3) = 0.1;
        pos(4) = 0.05;
        u3 = uicontrol('style', 'edit', 'units', 'normalized', 'callback', @update_axis, 'BackgroundColor', 'white');
        set(u3, 'position', pos);
        set(u3, 'position', pos);
        set(u3, 'string', num2str(vaxis(3,2)));
        set(u3, 'tag', 'u3');
      end
      
      if isempty(u4) || ~ishandle(u4)
        pos = get(p4, 'position'); % link the position to the subplot
        pos(1) = pos(1)-0.1;
        pos(2) = pos(2)-0.05;
        pos(3) = 0.1;
        pos(4) = 0.05;
        u4 = uicontrol('style', 'edit', 'units', 'normalized', 'callback', @update_axis, 'BackgroundColor', 'white');
        set(u4, 'position', pos);
        set(u4, 'string', num2str(vaxis(4,2)));
        set(u4, 'tag', 'u4');
      end
      
      if isempty(u5) || ~ishandle(u5)
        pos = get(p5, 'position'); % link the position to the subplot
        pos(1) = pos(1)-0.1;
        pos(2) = pos(2)-0.05;
        pos(3) = 0.1;
        pos(4) = 0.05;
        u5 = uicontrol('style', 'edit', 'units', 'normalized', 'callback', @update_axis, 'BackgroundColor', 'white');
        set(u5, 'position', pos);
        set(u5, 'string', num2str(vaxis(5,2)));
        set(u5, 'tag', 'u5');
      end
      
      if isempty(u6) || ~ishandle(u6)
        pos = get(p6, 'position'); % link the position to the subplot
        pos(1) = pos(1)-0.1;
        pos(2) = pos(2)-0.05;
        pos(3) = 0.1;
        pos(4) = 0.05;
        u6 = uicontrol('style', 'edit', 'units', 'normalized', 'callback', @update_axis, 'BackgroundColor', 'white');
        set(u6, 'position', pos);
        set(u6, 'string', num2str(vaxis(6,2)));
        set(u6, 'tag', 'u6');
      end
      
      if isempty(b2) || ~ishandle(b2)
        pos = [0.88 0.01 0.1 0.05];
        b2 = uicontrol('style', 'pushbutton', 'units', 'normalized', 'callback', 'evalin(''caller'', ''b2clicked = true;'')');
        set(b2, 'position', pos);
        set(b2, 'string', 'quit');
        set(b2, 'tag', 'b2');
      end
      
    end % try
    
    if b2clicked
      close all
      return
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % visualize the data in 2*3 subplots
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
      
      if isempty(p1) || ~ishandle(p1)
        p1 = subplot(3, 2, 1);
      else
        
        subplot(p1);
        h1 = plot(time, dat(1, :));
        axis([min(time) max(time) vaxis(1, 1) vaxis(1, 2)]);
        set(p1, 'XTickLabel', []);
        ylabel('amplitude (uV)');
        xlabel(sprintf('time: %d seconds', cfg.windowsize));
        grid on
        if strcmp(cfg.feedback, 'yes')
          ax = axis;
          line([ax(1) ax(2)], [left_thresh_ampl left_thresh_ampl], 'color', 'red');
          line([ax(1) + left_thresh_time(1)/hdr.Fs ax(1) + left_thresh_time(1)/hdr.Fs], [-300 300], 'color', 'green');
          line([ax(1) + left_thresh_time(2)/hdr.Fs ax(1) + left_thresh_time(2)/hdr.Fs], [-300 300], 'color', 'green');
        end
        
      end
      
      if isempty(p2) || ~ishandle(p2)
        p2 = subplot(3, 2, 2);
      else
        subplot(p2);
        h2 = plot(time, dat(2, :));
        axis([min(time) max(time) vaxis(2, 1) vaxis(2, 2)]);
        set(p2, 'XTickLabel', []);
        ylabel('amplitude (uV)');
        xlabel(sprintf('time: %d seconds', cfg.windowsize));
        grid on
      end
      
      if isempty(p3) || ~ishandle(p3)
        p3 = subplot(3, 2, 3);
      else
        subplot(p3)
        h3 = bar(1:length(freqoi), pow(1, :), 0.5);
        % plot(pow(1).Frequencies, pow(1).Data);
        % bar(pow(1).Frequencies, pow(1).Data);
        axis([cfg.foilim(1) cfg.foilim(2) vaxis(3, 1) vaxis(3, 2)]);
        % str = sprintf('time = %d s\n', round(mean(time)));
        % title(str);
        xlabel('frequency (Hz)');
        ylabel('power');
        
      end
      
      if isempty(p4) || ~ishandle(p4)
        p4 = subplot(3, 2, 4);
      else
        subplot(p4)
        h4 = bar(1:length(freqoi), pow(2, :), 0.5);
        % plot(pow(2).Frequencies, pow(2).Data);
        % bar(pow(2).Frequencies, pow(2).Data);
        ax = axis;
        axis([cfg.foilim(1) cfg.foilim(2) vaxis(4, 1) vaxis(4, 2)]);
        if strcmp(cfg.feedback, 'yes')
          line([right_freq(1) right_freq(1)], [ax(3) ax(4)]);
          line([right_freq(2) right_freq(2)], [ax(3) ax(4)]);
          line([right_freq(1) right_freq(2)], [right_offset right_offset]);
        end
        xlabel('frequency (Hz)');
        ylabel('power');
      end
      
      if isempty(p5) || ~ishandle(p5)
        p5 = subplot(3, 2, 5);
      else
        subplot(p5)
        h5 = surf(squeeze(TFR(1, :, :)));
        axis([1 100 cfg.foilim(1) cfg.foilim(2) vaxis(5, 1) vaxis(5, 2)]);
        view(110, 45);
        xlabel(''); % this is the historical time
        ylabel('frequency (Hz)');
        zlabel('power');
        set(p5, 'XTickLabel', []);
        set(h5, 'EdgeColor', 'none');
        shading interp
        box off
      end
      
      if isempty(p6) || ~ishandle(p6)
        p6 = subplot(3, 2, 6);
      else
        subplot(p6)
        h6 = surf(squeeze(TFR(2, :, :)));
        axis([1 100 cfg.foilim(1) cfg.foilim(2) vaxis(6, 1) vaxis(6, 2)]);
        view(110, 45);
        xlabel(''); % this is the historical time
        ylabel('frequency (Hz)');
        zlabel('power');
        set(p6, 'XTickLabel', []);
        set(h6, 'EdgeColor', 'none');
        shading interp
        box off
      end
      
    end % try
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % present MIDI feedback if the data exceeds the specified limits
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if strcmp(cfg.feedback, 'yes')
      if left_thresh_ampl < 0
        
        if min((dat(1, left_thresh_time(1):left_thresh_time(2)))) < left_thresh_ampl
          if beatdrum == true
            % midiOut('+', 10, 64, 127);
            beatdrum = false;
          else
            % midiOut('+', 10, 31, 127);
            beatdrum = true;
          end
        else
          midiOut('.', 1);
        end
        
      elseif max((dat(left_thresh_time(1):left_thresh_time(2)))) > left_thresh_ampl
        if beatdrum == true
          % midiOut('+', 10, 64, 127);
          beatdrum = false;
        else
          % midiOut('+', 10, 31, 127);
          beatdrum = true;
        end
      else
        % midiOut('.', 1);
      end
      
      % schemerlamp.setLevel(round(TFR(1, 60, end) / mean(TFR(1, 60, :)))*5);
      volume_right = round((mean(TFR(2, right_freq, end)) - right_offset) * right_mult);
      
      % midiOut(uint8([175+1, 7, volume_left]));
      % midiOut(uint8([16*14+1-1, 0, volume_left])); % ptich
      
      midiOut(uint8([175+2, 7, volume_right]));
      midiOut(uint8([16*14+2-1, 0, volume_right])); % ptich
      
      % schemerlamp.setLevel(9);
      
    end % if MIDI feedback
    
    % force an update of the figure
    drawnow
    
  end % if enough new samples
end % while true

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_channel(h, varargin)
global hdr chanindx
val = abs(str2num(get(h, 'string')));
val = max(1, min(val, hdr.nChans));
if ~isempty(val)
  switch get(h, 'tag')
    case 'c1'
      chanindx(1) = val;
      set(h, 'string', num2str(val));
    case 'c2'
      chanindx(2) = val;
      set(h, 'string', num2str(val));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_axis(h, varargin)
global vaxis
val = abs(str2num(get(h, 'string')));
if ~isempty(val)
  switch get(h, 'tag')
    case 'u1'
      vaxis(1,:) = [-val val];
    case 'u2'
      vaxis(2,:) = [-val val];
    case 'u3'
      vaxis(3,:) = [0 val];
    case 'u4'
      vaxis(4,:) = [0 val];
    case 'u5'
      vaxis(5,:) = [0 val];
    case 'u6'
      vaxis(6,:) = [0 val];
  end
end
