function [varargout] = ft_qualitycheck(cfg)

% FT_QUALITYCHECK facilitates quality inspection of a dataset. 1) The data is
% analyzed, quantified, and stored in a .mat file in a timelock- and
% freq- like fashion. 2) The quantifications are visualized and exported to
% a .PNG and .PDF file.
%
% In case the data is MEG data recorded with a CTF system, the output contains
% the headpos variable.
%
% Use as:
%   [info, timelock, freq, summary, headpos] = ft_qualitycheck(cfg)
%
% The configuration should contain:
%   cfg.dataset = a string (e.g. 'dataset.ds')
%
% The following parameters can be used:
%   cfg.analyze = 'yes' or 'no' to analyze the dataset (default = 'yes')
%   cfg.savemat = 'yes' or 'no' to save the analysis (default = 'yes')
%   cfg.matfile = a string (e.g. 'previousoutput.mat'), in combination with
%                   analyze = 'no' only
%   cfg.visualize = 'yes' or 'no' to visualize the analysis (default = 'yes')
%   cfg.saveplot = 'yes' or 'no' to save the visualization (default = 'yes')
%
% Copyright (C) 2010-2011, Arjen Stolk, Bram Daams, Robert Oostenveld
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


ft_defaults

% defaults
if ~isfield(cfg,'analyze'),        cfg.analyze   = 'yes';                         end
if ~isfield(cfg,'savemat'),        cfg.savemat   = 'yes';                         end
if ~isfield(cfg,'matfile'),        cfg.matfile   = [];                            end
if ~isfield(cfg,'visualize'),      cfg.visualize = 'yes';                         end
if ~isfield(cfg,'saveplot'),       cfg.saveplot  = 'yes';                         end

%% ANALYSIS
if strcmp(cfg.analyze,'yes')
    tic
    
    % checks
    iseeg = ft_filetype(cfg.dataset,'brainvision_eeg','ns_eeg','bci2000_dat','neuroprax_eeg','egi_sbin','biosemi_bdf');
    ismeg = ft_filetype(cfg.dataset,'ctf_ds','4d','neuromag_fif','itab_raw');
    isctf = ft_filetype(cfg.dataset, 'ctf_ds');
    cfg   = ft_checkconfig(cfg, 'dataset2files', 'yes'); % translate into datafile+headerfile
    
    % these will be replaced by more appropriate values
    info.datasetname = 'unknown';
    info.starttime   = 'unknown';
    info.startdate   = 'unknown';
    info.stoptime    = 'unknown';
    info.stopdate    = 'unknown';
    
    % the exportname is also used in the cron job
    exportname = qualitycheck_exportname(cfg.dataset);
    
    if isctf
        try
            % update the info fields
            info = read_ctf_hist(cfg.dataset);
        end
    end
    
    % read events
    info.event                  = ft_read_event(cfg.dataset);
    
    % read header file
    info.hdr                    = ft_read_header(cfg.dataset);
    
    % add filetype
    info.filetype               = ft_filetype(cfg.dataset);
    
    % trial definition
    cfgdef                      = [];
    cfgdef.dataset              = cfg.dataset;
    cfgdef.trialdef.triallength = 10;
    %cfgdef.trialdef.ntrials     = 3;
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
        chans                      = allchans;  % brain
        allchanindx                = match_str(info.hdr.label, allchans);
        chanindx                   = match_str(chans, allchans);
        jumpthreshold              = 1e4;
    end
    
    % find headcoil channels
    if isctf % this fails for older CTF data sets
        try
            Nx = strmatch('HLC0011', info.hdr.label); % x nasion coil
            Ny = strmatch('HLC0012', info.hdr.label); % y nasion
            Nz = strmatch('HLC0013', info.hdr.label); % z nasion
            Lx = strmatch('HLC0021', info.hdr.label); % x left coil
            Ly = strmatch('HLC0022', info.hdr.label); % y left
            Lz = strmatch('HLC0023', info.hdr.label); % z left
            Rx = strmatch('HLC0031', info.hdr.label); % x right coil
            Ry = strmatch('HLC0032', info.hdr.label); % y right
            Rz = strmatch('HLC0033', info.hdr.label); % z right
            hasheadpos     = true;
            headpos.dimord = 'chan_time';
            headpos.time   = [timeunit-timeunit/2:timeunit:timeunit*ntrials-timeunit/2];
            headpos.label  = {'Nx';'Ny';'Nz';'Lx';'Ly';'Lz';'Rx';'Ry';'Rz'};
            headpos.avg    = NaN(length(headpos.label), ntrials);
            headpos.grad   = info.hdr.grad;
        end % try
    end % if
    
    % analysis settings
    cfgredef             = [];
    cfgredef.length      = 1;
    cfgredef.overlap     = 0;
    
    cfgfreq              = [];
    cfgfreq.output       = 'pow';
    cfgfreq.channel      = chans;
    cfgfreq.method       = 'mtmfft';
    cfgfreq.taper        = 'hanning';
    cfgfreq.keeptrials   = 'no';
    cfgfreq.foilim       = [0 min(info.hdr.Fs/2, 400)];
    
    % output variables
    timelock.dimord = 'chan_time';
    timelock.label  = chans;
    timelock.time   = [timeunit-timeunit/2:timeunit:timeunit*ntrials-timeunit/2];
    timelock.avg    = NaN(length(chans), ntrials); % updated in loop
    timelock.median = NaN(length(chans), ntrials); % updated in loop
    timelock.jumps  = NaN(length(chans), ntrials); % updated in loop
    timelock.range  = NaN(length(chans), ntrials); % updated in loop
    timelock.min    = NaN(length(chans), ntrials); % updated in loop
    timelock.max    = NaN(length(chans), ntrials); % updated in loop
    
    freq.dimord     = 'chan_freq_time';
    freq.label      = chans;
    freq.freq       = [cfgfreq.foilim(1):cfgfreq.foilim(2)];
    freq.time       = [timeunit-timeunit/2:timeunit:timeunit*ntrials-timeunit/2];
    freq.powspctrm  = NaN(length(chans), length(freq.freq), ntrials); % updated in loop
    
    summary.dimord  = 'chan_time';
    summary.time    = [timeunit-timeunit/2:timeunit:timeunit*ntrials-timeunit/2];
    summary.label   = {'Mean';'Median';'Min';'Max';'Range';'HmotionN';'HmotionL';'HmotionR';'LowFreqPower';'LineFreqPower';'Jumps'};
    summary.avg     = NaN(length(summary.label), ntrials); % updated in loop
    
    % try add gradiometer info
    try
        timelock.grad   = info.hdr.grad;
        freq.grad       = info.hdr.grad;
        summary.grad    = info.hdr.grad;
    end
    
    % process trial by trial
    for t = 1:ntrials
        fprintf('analyzing trial %s of %s \n', num2str(t), num2str(ntrials));
        
        % preprocess
        cfgpreproc                  = cfgdef;
        cfgpreproc.trl              = cfgdef.trl(t,:);
        data                        = ft_preprocessing(cfgpreproc); clear cfgpreproc;
        
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
        timelock.avg(:,t)           = mean(data.trial{1}(allchanindx,:),2);
        timelock.median(:,t)        = median(data.trial{1}(allchanindx,:),2);
        timelock.range(:,t)         = max(data.trial{1}(allchanindx,:),[],2) - min(data.trial{1}(allchanindx,:),[],2);
        timelock.min(:,t)           = min(data.trial{1}(allchanindx,:),[],2);
        timelock.max(:,t)           = max(data.trial{1}(allchanindx,:),[],2);
        
        % detect jumps
        for c = 1:size(data.trial{1}(allchanindx,:),1)
            timelock.jumps(c,t)     = length(find(diff(data.trial{1,1}(allchanindx(c),:)) > jumpthreshold));
        end
        
        % FFT and noise estimation
        redef                       = ft_redefinetrial(cfgredef, data); clear data;
        FFT                         = ft_freqanalysis(cfgfreq, redef); clear redef;
        freq.powspctrm(:,:,t)       = FFT.powspctrm;
        summary.avg(9,t)            = mean(mean(findpower(0, 2, FFT, chanindx))); % Low Freq Power
        summary.avg(10,t)           = mean(mean(findpower(49, 51, FFT, chanindx))); clear FFT; % Line Freq Power
        
        toc
    end % end of trial loop
    
    % determine headmotion: distance from initial trial (in cm)
    if isctf && hasheadpos
        summary.avg(6,:) = sqrt(sum((headpos.avg(1:3,:)-repmat(headpos.avg(1:3,1),1,size(headpos.avg,2))).^2,1)); % N
        summary.avg(7,:) = sqrt(sum((headpos.avg(4:6,:)-repmat(headpos.avg(4:6,1),1,size(headpos.avg,2))).^2,1)); % L
        summary.avg(8,:) = sqrt(sum((headpos.avg(7:9,:)-repmat(headpos.avg(7:9,1),1,size(headpos.avg,2))).^2,1)); % R
    end
    
    % summarize/mean and store variables
    summary.avg(1,:)   = mean(timelock.avg(chanindx,:),1);
    summary.avg(2,:)   = mean(timelock.median(chanindx,:),1);
    summary.avg(3,:)   = mean(timelock.min(chanindx,:),1);
    summary.avg(4,:)   = mean(timelock.max(chanindx,:),1);
    summary.avg(5,:)   = mean(timelock.range(chanindx,:),1);
    summary.avg(11,:)  = mean(timelock.jumps(chanindx,:),1);
    
    % add the version details of this function call to the configuration
    cfg.version.name   = mfilename('fullpath');
    cfg.version.id     = '$Id$';
    cfg.version.matlab = version(); % Matlab version used
    
    % add the cfg to the output variables
    info.cfg           = cfg;
    timelock.cfg       = cfg;
    freq.cfg           = cfg;
    summary.cfg        = cfg;
    
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
    nplots = ceil(length(freq.time)/(3600/10));
    
    % create GUI-like figure(s)
    for p = 1:nplots
        fprintf('visualizing %s of %s \n', num2str(p), num2str(nplots));
        toi = [p*3600-3595 p*3600-5]; % select 1-hour chunks
        
        if exist('headpos','var')
            temp_timelock = ft_selectdata(timelock, 'toilim', toi);
            temp_freq     = ft_selectdata(freq, 'toilim', toi);
            temp_summary  = ft_selectdata(summary, 'toilim', toi);
            temp_headpos  = ft_selectdata(headpos, 'toilim', toi);
            draw_figure(info, temp_timelock, temp_freq, temp_summary, temp_headpos, toi);
            clear temp_timelock; clear temp_freq; clear temp_summary; clear temp_headpos; clear toi;
        else
            temp_timelock = ft_selectdata(timelock, 'toilim', toi);
            temp_freq     = ft_selectdata(freq, 'toilim', toi);
            temp_summary  = ft_selectdata(summary, 'toilim', toi);
            draw_figure(info, temp_timelock, temp_freq, temp_summary, toi);
            clear temp_timelock; clear temp_freq; clear temp_summary; clear toi;
        end
        
        % export to .PNG and .PDF
        if strcmp(cfg.saveplot, 'yes')
            [pathstr,name,extr,versn] = fileparts(exportname);
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


%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [power, freq] = findpower(low, high, freqinput, chans)
% replace value with the index of the nearest bin
xmin  = nearest(getsubfield(freqinput, 'freq'), low);
xmax  = nearest(getsubfield(freqinput, 'freq'), high);
% select the freq range
power = freqinput.powspctrm(chans, xmin:xmax);
freq  = freqinput.freq(:, xmin:xmax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
iseeg = ft_filetype(timelock.cfg.dataset,'brainvision_eeg','ns_eeg','bci2000_dat','neuroprax_eeg','egi_sbin','biosemi_bdf');
ismeg = ft_filetype(timelock.cfg.dataset,'ctf_ds','4d','neuromag_fif','itab_raw');
if ismeg
    scaling = 1e15;
    powscaling = 1e30;
    ylab = 'fT';
elseif iseeg
    scaling = 1e1;
    powscaling = 1e2;
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

h.HmotionAxes = axes(...
    'Parent',h.HmotionPanel,...
    'Units','normalized',...
    'color','white',...
    'Position',[.05 .08 .9 .52]);

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
    allchans = ft_senslabel('ctf275');
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
    'String',['fs: ' num2str(info.hdr.Fs) ', nchans: ' nchans],...
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
    hmotions = ([summary.avg(8,:)' summary.avg(7,:)' summary.avg(6,:)'])*10;
    boxplot(h.HmotionAxes, hmotions, 'orientation', 'horizontal', 'notch', 'on');
    set(h.HmotionAxes,'YTick',[1:3]);
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

h.HmotionTimecourseAxes = axes(...
    'Parent',h.SignalPanel,...
    'Units','normalized',...
    'color','white',...
    'Position',[.08 .73 .89 .22]);

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
    plot(h.HmotionTimecourseAxes, summary.time, summary.avg(6,:)*10, summary.time, summary.avg(7,:)*10, summary.time, summary.avg(8,:)*10, 'LineWidth',2);
    ylim(h.HmotionTimecourseAxes,[0 10]);
    ylabel(h.HmotionTimecourseAxes, 'Coil distance [mm]');
    xlim(h.HmotionTimecourseAxes,[toi]);
    grid(h.HmotionTimecourseAxes,'on');
    legend(h.HmotionTimecourseAxes, 'N','L','R');
end

% plot mean and range of the raw signal
plot(h.SignalAxes, summary.time, summary.avg(5,:)*scaling, summary.time, summary.avg(1,:)*scaling, 'LineWidth', 2);
set(h.SignalAxes,'Nextplot','add');
plot(h.SignalAxes, summary.time, summary.avg(3,:)*scaling, summary.time, summary.avg(4,:)*scaling, 'LineWidth', 1, 'Color', [255/255 127/255 39/255]);
grid(h.SignalAxes,'on');
ylabel(h.SignalAxes, ['Amplitude [' ylab ']']);
xlim(h.SignalAxes,[toi]);
legend(h.SignalAxes,'Range','Mean','Min','Max');
set(h.SignalAxes,'XTickLabel','');

% plot linenoise
semilogy(h.LinenoiseAxes, summary.time, summary.avg(10,:)*powscaling, 'LineWidth',2);
grid(h.LinenoiseAxes,'on');
legend(h.LinenoiseAxes, ['LineFreq [' ylab '^2/Hz]']);
set(h.LinenoiseAxes,'XTickLabel','');
xlim(h.LinenoiseAxes,[toi]);
ylim(h.LinenoiseAxes,[0 1e3]);

% plot lowfreqnoise
semilogy(h.LowfreqnoiseAxes, summary.time, summary.avg(9,:)*powscaling, 'LineWidth',2);
grid(h.LowfreqnoiseAxes,'on');
xlim(h.LowfreqnoiseAxes,[toi]);
ylim(h.LowfreqnoiseAxes,[0 1e10]);
legend(h.LowfreqnoiseAxes, ['LowFreq [' ylab '^2/Hz]']);
xlabel(h.LowfreqnoiseAxes, 'Time [seconds]');

% EVENT PANEL
h.EventPanel = uipanel(...
    'Parent',h.MainFigure,...
    'Units','normalized',...
    'Backgroundcolor','white',...
    'Position',[.7 .01 .29 .3]);

% event details
[a,b,c] = unique({info.event.type});
eventtypes    = {};
eventtriggers = {};
eventvalues   = {};
for j=1:length(a)
    eventtypes{j,1}    = a{j};
    eventtriggers{j,1} = sum(c==j);
    eventvalues{j,1}   = length(unique([info.event(c==j).value]));
end
if isempty(eventtypes)
    eventtypes{1,1} = 'no events found';
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
loglog(h.SpectrumAxes, freq.freq, squeeze(mean(mean(freq.powspctrm,1),3))*powscaling,'r','LineWidth',2);
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
    jumpchans{1,1} = 'no jumps found';
end

h.JumpText = uicontrol(...
    'Parent',h.JumpPanel,...
    'Style','text',...
    'Units','normalized',...
    'FontSize',10,...
    'String',[jumpchans],...
    'Backgroundcolor','white',...
    'Position',[.15 .5 .25 .4]);

h.JumpText2 = uicontrol(...
    'Parent',h.JumpPanel,...
    'Style','text',...
    'Units','normalized',...
    'FontSize',10,...
    'String',[jumpcounts],...
    'Backgroundcolor','white',...
    'Position',[.65 .5 .2 .4]);

h.TopoMEG = axes(...
    'Parent',h.JumpPanel,...
    'color','white',...
    'Units','normalized',...
    'Position',[0.4 0.05 0.55 0.4]);

% plot jumps on the dewar sensors
if exist('headpos','var')
    MEGchans                 = ft_channelselection('MEG', timelock.label);
    MEGchanindx              = match_str(timelock.label, MEGchans);
    cfgtopo                  = [];
    cfgtopo.marker           = 'off';
    cfgtopo.colorbar         = 'no';
    cfgtopo.comment          = 'no';
    cfgtopo.style            = 'blank';
    cfgtopo.layout           = 'CTF275.lay';
    cfgtopo.highlight        = 'on';
    cfgtopo.highlightsymbol  = '.';
    cfgtopo.highlightsize    = [14];
    cfgtopo.highlightchannel = [find(sum(timelock.jumps(MEGchanindx,:),2)>0)];
    data.label               = MEGchans;
    data.powspctrm           = sum(timelock.jumps(MEGchanindx,:),2);
    data.dimord              = 'chan_freq';
    data.freq                = 1;
    axes(h.TopoMEG);
    ft_topoplotTFR(cfgtopo, data); clear data;
end