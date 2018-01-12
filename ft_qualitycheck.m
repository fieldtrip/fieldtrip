function [varargout] = ft_qualitycheck(cfg)

% FT_QUALITYCHECK performs a quality inspection of a given MEG/EEG dataset,
% stores (.mat), and visualizes the result (.png and .pdf).
%
% This function segments the data into 10-second pieces and performs the
% following analyses:
%  1) reads the properties of the dataset
%  2) computes the headpositions and distance covered from recording onset (CTF only)
%  3) computes the mean, max, min, and range of the signal amplitude
%  4) detects trigger events
%  5) detects jump artifacts
%  6) computes the powerspectrum
%  7) estimates the low-frequency (<2 Hz) and line noise (~50 Hz)
%
% Use as
%   [info, timelock, freq, summary, headpos] = ft_qualitycheck(cfg)
% where info contains the dataset properties, timelock the timelocked data,
% freq the powerspectra, summary the mean descriptives, and headpos the
% headpositions throughout the recording
%
% The configuration should contain:
%   cfg.dataset = a string (e.g. 'dataset.ds')
%
% The following parameters can be used:
%   cfg.analyze   = string, 'yes' or 'no' to analyze the dataset (default = 'yes')
%   cfg.savemat   = string, 'yes' or 'no' to save the analysis (default = 'yes')
%   cfg.matfile   = string, filename (e.g. 'previousoutput.mat'), preferably in combination
%                    with analyze = 'no'
%   cfg.visualize = string, 'yes' or 'no' to visualize the analysis (default = 'yes')
%   cfg.saveplot  = string, 'yes' or 'no' to save the visualization (default = 'yes')
%   cfg.linefreq  = scalar, frequency of power line (default = 50)
%   cfg.plotunit  = scalar, the length of time to be plotted in one panel (default = 3600)
%
% See also FT_PREPROCESSING, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT

% Copyright (C) 2010-2011, Arjen Stolk, Bram Daams, Robert Oostenveld
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
% $Id:%

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble provenance
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% set the defaults
cfg.analyze   = ft_getopt(cfg, 'analyze',   'yes');
cfg.savemat   = ft_getopt(cfg, 'savemat',   'yes');
cfg.matfile   = ft_getopt(cfg, 'matfile',   []);
cfg.visualize = ft_getopt(cfg, 'visualize', 'yes');
cfg.saveplot  = ft_getopt(cfg, 'saveplot',  'yes');
cfg.linefreq  = ft_getopt(cfg, 'linefreq',  50);
cfg.plotunit  = ft_getopt(cfg, 'plotunit',  3600);

%% ANALYSIS
if strcmp(cfg.analyze,'yes')
  tic
  
  % checks
  cfg   = ft_checkconfig(cfg, 'dataset2files', 'yes'); % translate into datafile+headerfile
  
  % these will be replaced by more appropriate values
  info.filename    = cfg.dataset;
  info.datasetname = 'unknown';
  info.starttime   = 'unknown';
  info.startdate   = 'unknown';
  info.stoptime    = 'unknown';
  info.stopdate    = 'unknown';
  
  % the exportname is also used in the cron job
  exportname = qualitycheck_exportname(cfg.dataset);
  [iseeg, ismeg, isctf, fltp] = filetyper(cfg.dataset);
  if isctf
    try
      % update the info fields
      info = read_ctf_hist(cfg.dataset);
    end
  end
  
  % add info
  info.event                  = ft_read_event(cfg.dataset);
  info.hdr                    = ft_read_header(cfg.dataset);
  info.filetype               = fltp;
  
  % trial definition
  cfgdef                      = [];
  cfgdef.dataset              = cfg.dataset;
  cfgdef.trialdef.triallength = 10;
  %cfgdef.trialdef.ntrials     = 3; % for debugging
  cfgdef.continuous           = 'yes';
  cfgdef                      = ft_definetrial(cfgdef);
  ntrials                     = size(cfgdef.trl,1)-1; % remove last trial
  timeunit                    = cfgdef.trialdef.triallength;
  
  % channelselection for jump detection (all) and for FFT (brain)
  if ismeg
    allchans                   = ft_channelselection({'MEG','MEGREF'}, info.hdr.label);
    chans                      = ft_channelselection('MEG', info.hdr.label); % brain
    allchanindx                = match_str(info.hdr.label, allchans);
    chanindx                   = match_str(chans, allchans);
    jumpthreshold              = 1e-10;
  elseif iseeg
    allchans                   = ft_channelselection('EEG', info.hdr.label);
    if isempty(allchans)
      % some EEG systems and data files use non-standard channel names that are not detected automatically
      ft_warning('no EEG channels detected, selecting all channels');
      allchans = info.hdr.label;
    end
    chans                      = allchans;  % brain
    allchanindx                = match_str(info.hdr.label, allchans);
    chanindx                   = match_str(chans, allchans);
    jumpthreshold              = 1e4;
  end
  
  % find headcoil channels
  if isctf % this fails for older CTF data sets
    Nx = strmatch('HLC0011', info.hdr.label); % x nasion coil
    Ny = strmatch('HLC0012', info.hdr.label); % y nasion
    Nz = strmatch('HLC0013', info.hdr.label); % z nasion
    Lx = strmatch('HLC0021', info.hdr.label); % x left coil
    Ly = strmatch('HLC0022', info.hdr.label); % y left
    Lz = strmatch('HLC0023', info.hdr.label); % z left
    Rx = strmatch('HLC0031', info.hdr.label); % x right coil
    Ry = strmatch('HLC0032', info.hdr.label); % y right
    Rz = strmatch('HLC0033', info.hdr.label); % z right
    headpos.dimord = 'chan_time';
    headpos.time   = (timeunit-timeunit/2:timeunit:timeunit*ntrials-timeunit/2);
    headpos.label  = {'Nx';'Ny';'Nz';'Lx';'Ly';'Lz';'Rx';'Ry';'Rz'};
    headpos.avg    = NaN(length(headpos.label), ntrials);
    headpos.grad   = info.hdr.grad;
    
    if numel(cat(1,Nx,Ny,Nz,Lx,Ly,Lz,Rx,Ry,Rz))==9
      hasheadpos = true;
    else
      hasheadpos = false;
    end
    
  end % if
  
  % analysis settings
  cfgredef             = [];
  cfgredef.length      = 1;
  cfgredef.overlap     = 0;
  
  cfgfreq              = [];
  cfgfreq.output       = 'pow';
  cfgfreq.channel      = allchans;
  cfgfreq.method       = 'mtmfft';
  cfgfreq.taper        = 'hanning';
  cfgfreq.keeptrials   = 'no';
  cfgfreq.foilim       = [0 min(info.hdr.Fs/2, 400)];
  
  % output variables
  timelock.dimord = 'chan_time';
  timelock.label  = allchans;
  timelock.time   = (timeunit-timeunit/2:timeunit:timeunit*ntrials-timeunit/2);
  timelock.avg    = NaN(length(allchans), ntrials); % updated in loop
  timelock.median = NaN(length(allchans), ntrials); % updated in loop
  timelock.jumps  = NaN(length(allchans), ntrials); % updated in loop
  timelock.range  = NaN(length(allchans), ntrials); % updated in loop
  timelock.min    = NaN(length(allchans), ntrials); % updated in loop
  timelock.max    = NaN(length(allchans), ntrials); % updated in loop
  
  freq.dimord     = 'chan_freq_time';
  freq.label      = allchans;
  freq.freq       = (cfgfreq.foilim(1):cfgfreq.foilim(2));
  freq.time       = (timeunit-timeunit/2:timeunit:timeunit*ntrials-timeunit/2);
  freq.powspctrm  = NaN(length(allchans), length(freq.freq), ntrials); % updated in loop
  
  summary.dimord  = 'chan_time';
  summary.time    = (timeunit-timeunit/2:timeunit:timeunit*ntrials-timeunit/2);
  summary.label   = {'Mean';'Median';'Min';'Max';'Range';'HmotionN';'HmotionL';'HmotionR';'LowFreqPower';'LineFreqPower';'Jumps'};
  summary.avg     = NaN(length(summary.label), ntrials); % updated in loop
  
  % try add gradiometer info
  if isfield(info.hdr, 'grad')
    timelock.grad = info.hdr.grad;
    freq.grad     = info.hdr.grad;
    summary.grad  = info.hdr.grad;
  end
  
  
  % process trial by trial
  for t = 1:ntrials
    fprintf('analyzing trial %s of %s \n', num2str(t), num2str(ntrials));
    
    % preprocess
    cfgpreproc     = cfgdef;
    cfgpreproc.trl = cfgdef.trl(t,:);
    data           = ft_preprocessing(cfgpreproc); clear cfgpreproc;
    
    % determine headposition
    if isctf && hasheadpos
      headpos.avg(1,t) = mean(data.trial{1,1}(Nx,:) * 100);  % meter to cm
      headpos.avg(2,t) = mean(data.trial{1,1}(Ny,:) * 100);
      headpos.avg(3,t) = mean(data.trial{1,1}(Nz,:) * 100);
      headpos.avg(4,t) = mean(data.trial{1,1}(Lx,:) * 100);
      headpos.avg(5,t) = mean(data.trial{1,1}(Ly,:) * 100);
      headpos.avg(6,t) = mean(data.trial{1,1}(Lz,:) * 100);
      headpos.avg(7,t) = mean(data.trial{1,1}(Rx,:) * 100);
      headpos.avg(8,t) = mean(data.trial{1,1}(Ry,:) * 100);
      headpos.avg(9,t) = mean(data.trial{1,1}(Rz,:) * 100);
    end
    
    % update values
    timelock.avg(:,t)    = mean(data.trial{1}(allchanindx,:),2);
    timelock.median(:,t) = median(data.trial{1}(allchanindx,:),2);
    timelock.range(:,t)  = max(data.trial{1}(allchanindx,:),[],2) - min(data.trial{1}(allchanindx,:),[],2);
    timelock.min(:,t)    = min(data.trial{1}(allchanindx,:),[],2);
    timelock.max(:,t)    = max(data.trial{1}(allchanindx,:),[],2);
    
    % detect jumps
    for c = 1:size(data.trial{1}(allchanindx,:),1)
      timelock.jumps(c,t) = length(find(diff(data.trial{1,1}(allchanindx(c),:)) > jumpthreshold));
    end
    
    % FFT and noise estimation
    redef                 = ft_redefinetrial(cfgredef, data); clear data;
    FFT                   = ft_freqanalysis(cfgfreq, redef); clear redef;
    freq.powspctrm(:,:,t) = FFT.powspctrm;
    summary.avg(9,t)      = mean(mean(findpower(0,  2,  FFT, chanindx))); % Low Freq Power
    summary.avg(10,t)     = mean(mean(findpower(cfg.linefreq-1, cfg.linefreq+1, FFT, chanindx))); clear FFT; % Line Freq Power
    
    toc
  end % end of trial loop
  
  % determine headmotion: distance from initial trial (in cm)
  if isctf && hasheadpos
    summary.avg(6,:) = sqrt(sum((headpos.avg(1:3,:)-repmat(headpos.avg(1:3,1),1,size(headpos.avg,2))).^2,1)); % N
    summary.avg(7,:) = sqrt(sum((headpos.avg(4:6,:)-repmat(headpos.avg(4:6,1),1,size(headpos.avg,2))).^2,1)); % L
    summary.avg(8,:) = sqrt(sum((headpos.avg(7:9,:)-repmat(headpos.avg(7:9,1),1,size(headpos.avg,2))).^2,1)); % R
  end
  
  % summarize/mean and store variables of brain info only
  summary.avg(1,:)   = mean(timelock.avg(chanindx,:),1);
  summary.avg(2,:)   = mean(timelock.median(chanindx,:),1);
  summary.avg(3,:)   = mean(timelock.min(chanindx,:),1);
  summary.avg(4,:)   = mean(timelock.max(chanindx,:),1);
  summary.avg(5,:)   = mean(timelock.range(chanindx,:),1);
  summary.avg(11,:)  = mean(timelock.jumps(chanindx,:),1);
  
  % save to .mat
  if strcmp(cfg.savemat, 'yes')
    if isctf && hasheadpos
      headpos.cfg        = cfg;
      save(exportname, 'info','timelock','freq','summary','headpos');
    else
      save(exportname, 'info','timelock','freq','summary');
    end
  end
  
end % end of analysis

%% VISUALIZATION
if strcmp(cfg.visualize, 'yes')
  
  % load data
  if strcmp(cfg.analyze, 'no')
    if ~isempty(cfg.matfile)
      exportname = cfg.matfile;
    else
      exportname = qualitycheck_exportname(cfg.dataset);
    end
    fprintf('loading %s \n', exportname);
    load(exportname);
  end
  
  % determine number of 1-hour plots to be made
  nplots = ceil(length(freq.time)/(cfg.plotunit/10));
  
  % create GUI-like figure(s)
  for p = 1:nplots
    fprintf('visualizing %s of %s \n', num2str(p), num2str(nplots));
    toi = [p*cfg.plotunit-(cfg.plotunit-5) p*cfg.plotunit-5]; % select 1-hour chunks
    
    tmpcfg.latency = toi;
    temp_timelock  = ft_selectdata(tmpcfg, timelock);
    temp_freq      = ft_selectdata(tmpcfg, freq);
    temp_summary   = ft_selectdata(tmpcfg, summary);
    if exist('headpos','var')
      temp_headpos  = ft_selectdata(tmpcfg, headpos);
      draw_figure(info, temp_timelock, temp_freq, temp_summary, temp_headpos, toi);
      clear temp_timelock; clear temp_freq; clear temp_summary; clear temp_headpos; clear toi;
    else
      draw_figure(info, temp_timelock, temp_freq, temp_summary, toi);
      clear temp_timelock; clear temp_freq; clear temp_summary; clear toi;
    end
    
    % export to .PNG and .PDF
    if strcmp(cfg.saveplot, 'yes')
      [pathstr,name,extr] = fileparts(exportname);
      if p == 1
        exportfilename = name;
      else
        exportfilename = strcat(name,'_pt',num2str(p));
      end
      fprintf('exporting %s of %s \n', num2str(p), num2str(nplots));
      set(gcf, 'PaperType', 'a4');
      print(gcf, '-dpng', strcat(exportfilename,'.png'));
      orient landscape;
      print(gcf, '-dpdf', strcat(exportfilename,'.pdf'));
      close
    end
  end % end of nplots
end % end of visualization

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
ft_postamble history timelock   % add the input cfg to multiple outputs
ft_postamble history freq       % add the input cfg to multiple outputs
ft_postamble history summary    % add the input cfg to multiple outputs

%% VARARGOUT
if nargout>0
  mOutputArgs{1} = info;
  mOutputArgs{2} = timelock;
  mOutputArgs{3} = freq;
  mOutputArgs{4} = summary;
  try
    mOutputArgs{5} = headpos;
  end
  [varargout{1:nargout}] = mOutputArgs{:};
  clearvars -except varargout
else
  clear
end

%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%
function [x] = clipat(x, v, v2)
v = [v v2]; % clip between value v and v2
if length(v) == 1
  x(x>v) = v;
elseif length(v) == 2
  x(x<v(1)) = v(1);
  x(x>v(2)) = v(2);
end

%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%
function [iseeg, ismeg, isctf, fltp] = filetyper(dataset)
fltp = ft_filetype(dataset);
iseeg = ft_filetype(dataset,'brainvision_eeg') | ...
  ft_filetype(dataset,'ns_eeg') | ...
  ft_filetype(dataset,'bci2000_dat') | ...
  ft_filetype(dataset,'neuroprax_eeg') | ...
  ft_filetype(dataset,'egi_sbin') | ...
  ft_filetype(dataset,'biosemi_bdf');
ismeg = ft_filetype(dataset,'ctf_ds') | ...
  ft_filetype(dataset,'4d') | ...
  ft_filetype(dataset,'neuromag_fif') | ...
  ft_filetype(dataset,'itab_raw');
isctf = ft_filetype(dataset, 'ctf_ds');
if ~ismeg && ~iseeg % if none found, try less strict checks
  [p, f, ext] = fileparts(dataset);
  if strcmp(ext, '.eeg')
    fltp = 'brainvision_eeg';
    iseeg = 1;
  elseif strcmp(ext, '.bdf')
    fltp = 'biosemi_bdf';
    iseeg = 1;
  elseif strcmp(ext, '.ds')
    fltp = 'ctf_ds';
    ismeg = 1;
  else % otherwise use eeg settings for stability reasons
    iseeg = 1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%
function [power, freq] = findpower(low, high, freqinput, chans)
% replace value with the index of the nearest bin
xmin  = nearest(getsubfield(freqinput, 'freq'), low);
xmax  = nearest(getsubfield(freqinput, 'freq'), high);
% select the freq range
power = freqinput.powspctrm(chans, xmin:xmax);
freq  = freqinput.freq(:, xmin:xmax);

%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%
function draw_figure(varargin)
% deal with input
if nargin == 6
  info     = varargin{1};
  timelock = varargin{2};
  freq     = varargin{3};
  summary  = varargin{4};
  headpos  = varargin{5};
  toi      = varargin{6};
elseif nargin == 5
  info     = varargin{1};
  timelock = varargin{2};
  freq     = varargin{3};
  summary  = varargin{4};
  toi      = varargin{5};
end

% determine whether it is EEG or MEG
if isfield(info, 'filename') % supported as of January 2018
  filename = info.filename;
elseif isfield(timelock.cfg, 'dataset')
  filename = timelock.cfg.dataset;
elseif isfield(info.hdr.orig, 'FileName')
  filename = info.hdr.orig.FileName;
elseif exist('headpos','var') && isfield(headpos.cfg.previous, 'dataset')
  filename = headpos.cfg.previous.dataset;
else
  error('could not determine the filename');
end
[iseeg, ismeg, isctf, fltp] = filetyper(filename);

if ismeg
  scaling = 1e15; % assuming data is in T and needs to become fT
  powscaling = scaling^2;
  ylab = 'fT';
elseif iseeg
  scaling = 1e0; % assuming data is in muV already
  powscaling = scaling^2;
  ylab = '\muV';
end

% PARENT FIGURE
h.MainFigure = figure(...
  'MenuBar','none',...
  'Name','ft_qualitycheck',...
  'Units','normalized',...
  'color','white',...
  'Position',[0.01 0.01 .99 .99]); % nearly fullscreen

if strcmp(info.startdate,'unknown')
  tmp = 'unknown';
else
  [d,w] = weekday(info.startdate);
  tmp = [w ' ' info.startdate];
end

h.MainText = uicontrol(...
  'Parent',h.MainFigure,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String',tmp,...
  'Backgroundcolor','white',...
  'Position',[.06 .96 .15 .02]);

h.MainText2 = uicontrol(...
  'Parent',h.MainFigure,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String','Jump artefacts',...
  'Backgroundcolor','white',...
  'Position',[.08 .46 .12 .02]);

h.MainText3 = uicontrol(...
  'Parent',h.MainFigure,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String','Mean powerspectrum',...
  'Backgroundcolor','white',...
  'Position',[.4 .3 .15 .02]);

h.MainText4 = uicontrol(...
  'Parent',h.MainFigure,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String','Timecourses',...
  'Backgroundcolor','white',...
  'Position',[.5 .96 .11 .02]);

h.MainText5 = uicontrol(...
  'Parent',h.MainFigure,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String','Events',...
  'Backgroundcolor','white',...
  'Position',[.81 .3 .06 .02]);

% HEADMOTION PANEL
h.HmotionPanel = uipanel(...
  'Parent',h.MainFigure,...
  'Units','normalized',...
  'Backgroundcolor','white',...
  'Position',[.01 .5 .25 .47]);

h.DataText = uicontrol(...
  'Parent',h.HmotionPanel,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String',info.datasetname,...
  'Backgroundcolor','white',...
  'Position',[.01 .85 .99 .1]);

h.TimeText = uicontrol(...
  'Parent',h.HmotionPanel,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String',[info.starttime ' - ' info.stoptime],...
  'Backgroundcolor','white',...
  'Position',[.01 .78 .99 .1]);

if ismeg
  allchans = ft_senslabel(ft_senstype(timelock));
  misschans = setdiff(ft_channelselection('MEG', info.hdr.label), allchans);
  nchans = num2str(size(ft_channelselection('MEG', info.hdr.label),1));
else
  misschans = '';
  nchans = num2str(size(ft_channelselection('EEG', info.hdr.label),1));
end

h.DataText2 = uicontrol(...
  'Parent',h.HmotionPanel,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String',[fltp ', fs: ' num2str(info.hdr.Fs) ', nchans: ' nchans],...
  'Backgroundcolor','white',...
  'Position',[.01 .71 .99 .1]);

h.DataText3 = uicontrol(...
  'Parent',h.HmotionPanel,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String',['missing chans: ' misschans'],...
  'Backgroundcolor','white',...
  'Position',[.01 .64 .99 .1]);

% boxplot headmotion (*10; cm-> mm) per coil
if exist('headpos','var')
  h.HmotionAxes = axes(...
    'Parent',h.HmotionPanel,...
    'Units','normalized',...
    'color','white',...
    'Position',[.05 .08 .9 .52]);
  
  hmotions = ([summary.avg(8,:)' summary.avg(7,:)' summary.avg(6,:)'])*10;
  boxplot(h.HmotionAxes, hmotions, 'orientation', 'horizontal', 'notch', 'on');
  set(h.HmotionAxes,'YTick',1:3);
  set(h.HmotionAxes,'YTickLabel',{'R','L','N'});
  xlim(h.HmotionAxes, [0 10]);
  xlabel(h.HmotionAxes, 'Headmotion from start [mm]');
end

% TIMECOURSE PANEL
h.SignalPanel = uipanel(...
  'Parent',h.MainFigure,...
  'Units','normalized',...
  'Backgroundcolor','white',...
  'Position',[.28 .34 .71 .63]);

h.SignalAxes = axes(...
  'Parent',h.SignalPanel,...
  'Units','normalized',...
  'color','white',...
  'Position',[.08 .36 .89 .3]);

h.LinenoiseAxes = axes(...
  'Parent',h.SignalPanel,...
  'Units','normalized',...
  'color','white',...
  'Position',[.08 .23 .89 .1]);

h.LowfreqnoiseAxes = axes(...
  'Parent',h.SignalPanel,...
  'Units','normalized',...
  'color','white',...
  'Position',[.08 .1 .89 .1]);

% plot hmotion timecourses per coil (*10; cm-> mm)
if exist('headpos','var')
  h.HmotionTimecourseAxes = axes(...
    'Parent',h.SignalPanel,...
    'Units','normalized',...
    'color','white',...
    'Position',[.08 .73 .89 .22]);
  
  plot(h.HmotionTimecourseAxes, summary.time, clipat(summary.avg(6,:)*10, 0, 10), ...
    summary.time, clipat(summary.avg(7,:)*10, 0, 10), ...
    summary.time, clipat(summary.avg(8,:)*10, 0, 10), 'LineWidth',2);
  ylim(h.HmotionTimecourseAxes,[0 10]);
  ylabel(h.HmotionTimecourseAxes, 'Coil distance [mm]');
  xlim(h.HmotionTimecourseAxes,toi);
  grid(h.HmotionTimecourseAxes,'on');
  legend(h.HmotionTimecourseAxes, 'N','L','R');
end

% plot mean and range of the raw signal
plot(h.SignalAxes, summary.time, summary.avg(5,:)*scaling, summary.time, summary.avg(1,:)*scaling, 'LineWidth', 2);
set(h.SignalAxes,'Nextplot','add');
plot(h.SignalAxes, summary.time, summary.avg(3,:)*scaling, summary.time, summary.avg(4,:)*scaling, 'LineWidth', 1, 'Color', [255/255 127/255 39/255]);
grid(h.SignalAxes,'on');
ylabel(h.SignalAxes, ['Amplitude [' ylab ']']);
xlim(h.SignalAxes,toi);
legend(h.SignalAxes,'Range','Mean','Min','Max');
set(h.SignalAxes,'XTickLabel','');

% plot linenoise
semilogy(h.LinenoiseAxes, summary.time, clipat(summary.avg(10,:)*powscaling, 1e2, 1e4), 'LineWidth',2);
grid(h.LinenoiseAxes,'on');
legend(h.LinenoiseAxes, ['LineFreq [' ylab '^2/Hz]']);
set(h.LinenoiseAxes,'XTickLabel','');
xlim(h.LinenoiseAxes,toi);
ylim(h.LinenoiseAxes,[1e2 1e4]); % before april 28th this was 1e0 - 1e3

% plot lowfreqnoise
semilogy(h.LowfreqnoiseAxes, summary.time, clipat(summary.avg(9,:)*powscaling, 1e10, 1e12), 'LineWidth',2);
grid(h.LowfreqnoiseAxes,'on');
xlim(h.LowfreqnoiseAxes,toi);
ylim(h.LowfreqnoiseAxes,[1e10 1e12]);
legend(h.LowfreqnoiseAxes, ['LowFreq [' ylab '^2/Hz]']);
xlabel(h.LowfreqnoiseAxes, 'Time [seconds]'); % before april 28th this was 1e0 - 1e10

% EVENT PANEL
h.EventPanel = uipanel(...
  'Parent',h.MainFigure,...
  'Units','normalized',...
  'Backgroundcolor','white',...
  'Position',[.7 .01 .29 .3]);

% event details
eventtypes    = {};
eventtriggers = {};
eventvalues   = {};

if ~isempty(info.event)
  [a,b,c] = unique({info.event.type});
  for j=1:length(a)
    eventtypes{j,1}    = a{j};
    eventtriggers{j,1} = sum(c==j);
    eventvalues{j,1}   = length(unique([info.event(c==j).value]));
  end
end
if isempty(eventtypes)
  eventtypes{1,1} = 'no triggers found';
end

h.EventText = uicontrol(...
  'Parent',h.EventPanel,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String',['Types'; ' '; eventtypes],...
  'Backgroundcolor','white',...
  'Position',[.05 .05 .4 .85]);

h.EventText2 = uicontrol(...
  'Parent',h.EventPanel,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String',['Triggers'; ' '; eventtriggers],...
  'Backgroundcolor','white',...
  'Position',[.55 .05 .2 .85]);

h.EventText3 = uicontrol(...
  'Parent',h.EventPanel,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String',['Values'; ' '; eventvalues],...
  'Backgroundcolor','white',...
  'Position',[.8 .05 .15 .85]);

% POWERSPECTRUM PANEL
h.SpectrumPanel = uipanel(...
  'Parent',h.MainFigure,...
  'Units','normalized',...
  'Backgroundcolor','white',...
  'Position',[.28 .01 .4 .3]);

h.SpectrumAxes = axes(...
  'Parent',h.SpectrumPanel,...
  'Units','normalized',...
  'color','white',...
  'Position',[.15 .2 .8 .7]);

% plot powerspectrum
loglog(h.SpectrumAxes, freq.freq, mean(mean(freq.powspctrm,1),3)*powscaling,'r','LineWidth',2);
xlabel(h.SpectrumAxes, 'Frequency [Hz]');
ylabel(h.SpectrumAxes, ['Power [' ylab '^2/Hz]']);

% ARTEFACT PANEL
h.JumpPanel = uipanel(...
  'Parent',h.MainFigure,...
  'Units','normalized',...
  'Backgroundcolor','white',...
  'Position',[.01 .01 .25 .46]);

% jump details
jumpchans  = {};
jumpcounts = {};
[jumps,i] = find(timelock.jumps>0); % find all jumps
[a,b,c] = unique(jumps);
for j=1:length(a)
  jumpchans{j,1}  = timelock.label{a(j)};
  jumpcounts{j,1} = sum(c==j);
end
if isempty(jumpchans)
  jumpchans{1,1} = 'no jumps detected';
end

h.JumpText = uicontrol(...
  'Parent',h.JumpPanel,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String',jumpchans,...
  'Backgroundcolor','white',...
  'Position',[.15 .5 .25 .4]);

h.JumpText2 = uicontrol(...
  'Parent',h.JumpPanel,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String',jumpcounts,...
  'Backgroundcolor','white',...
  'Position',[.65 .5 .2 .4]);

% plot jumps on the dewar sensors
if ismeg
  h.TopoMEG = axes(...
    'Parent',h.JumpPanel,...
    'color','white',...
    'Units','normalized',...
    'Position',[0.4 0.05 0.55 0.4]);
  
  MEGchans                 = ft_channelselection('MEG', timelock.label);
  MEGchanindx              = match_str(timelock.label, MEGchans);
  cfgtopo                  = [];
  cfgtopo.marker           = 'off';
  cfgtopo.colorbar         = 'no';
  cfgtopo.comment          = 'no';
  cfgtopo.style            = 'blank';
  cfgtopo.layout           = ft_prepare_layout(timelock);
  cfgtopo.highlight        = 'on';
  cfgtopo.highlightsymbol  = '.';
  cfgtopo.highlightsize    = 14;
  cfgtopo.highlightchannel = find(sum(timelock.jumps(MEGchanindx,:),2)>0);
  data.label               = MEGchans;
  data.powspctrm           = sum(timelock.jumps(MEGchanindx,:),2);
  data.dimord              = 'chan_freq';
  data.freq                = 1;
  axes(h.TopoMEG);
  ft_topoplotTFR(cfgtopo, data); clear data;
end
