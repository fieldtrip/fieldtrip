function test_tutorial_spikefield

% performs all the operations mentioned int the spikefield tutorial
% (http://fieldtrip.fcdonders.nl/tutorial/spikefield), but only plots figures
% that are called by external functions (e.g. ft_singleplotTFR).

% it corresponds to the tutorial on 2012-09-25

% TEST test_tutorial_spikefield
% TEST ft_read_spike ft_spike_select ft_preprocessing ft_appendspike 
% TEST ft_definetrial ft_spike_maketrials ft_spiketriggeredinterpolation
% TEST ft_spiketriggeredspectrum ft_spiketriggeredspectrum_stat 
% TEST ft_spiketriggeredaverage ft_singleplotTFR

addpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/spikefield/')

filename         = 'p029_sort_final_01.nex';
spike            = ft_read_spike(filename); 
cfg              = [];
cfg.spikechannel = {'sig002a_wf', 'sig003a_wf'};
spike  = ft_spike_select(cfg, spike);

% get the cfg.trl
cfg = [];
cfg.dataset  = filename;
cfg.trialfun = 'trialfun_stimon_samples';
cfg          = ft_definetrial(cfg);

% read in the data in trials
cfg.channel   = {'AD01', 'AD02', 'AD03', 'AD04'}; % these channels contain the LFP
cfg.padding   = 10; % length to which we pad for filtering
cfg.dftfreq   = [60-1*(1/10):(1/10):60+1*(1/10) ]; % filter out 60 hz line noise
cfg.dftfilter = 'yes';
[data_lfp] = ft_preprocessing(cfg); % read in the LFP

cfg = [];
cfg.dataset   = filename;
cfg.trialfun  = 'trialfun_stimon_samples';
cfg           = ft_definetrial(cfg);
trl           = cfg.trl;
cfg = []; 
cfg.hdr       = data_lfp.hdr; % contains information for conversion of samples to timestamps
cfg.trlunit   = 'samples';
cfg.trl       = trl; % now in samples
spikeTrials   = ft_spike_maketrials(cfg,spike); 

data_all = ft_appendspike([],data_lfp, spike);

cfg = [];
cfg.method = 'nan'; % replace the removed segment with nans
cfg.timwin = [-0.002 0.002]; % remove 4 ms around every spike
cfg.spikechannel = spike.label{1};
cfg.channel = data_lfp.label(2);
data_nan = ft_spiketriggeredinterpolation(cfg, data_all);
cfg.method = 'linear'; % remove the replaced segment with interpolation
data_i = ft_spiketriggeredinterpolation(cfg, data_all);

% figure, 
% plot(data_i.time{1},data_i.trial{1}(2,:),'g-'), hold on, plot(data_i.time{1}, data_i.trial{1}(5,:),'r')
% hold on
% plot(data_nan.time{1},data_nan.trial{1}(2,:),'go')
% hold on
% plot(data_all.time{1},data_all.trial{1}(2,:),'k-')
% xlim([0.9 1])
% xlabel('time (s)')

%% STA

cfg = [];
cfg.timwin = [-0.25 0.25]; % take 400 ms
cfg.spikechannel = spike.label{1}; % first unit
cfg.channel = data_lfp.label(1:4); % first four chans
cfg.latency = [0.3 10];
staPost = ft_spiketriggeredaverage(cfg, data_all);
 
% % plot the sta
% figure, plot(staPost.time, staPost.avg(:,:)')
% legend(data_lfp.label)
% xlabel('time (s)')
% xlim(cfg.timwin)


%% spike phase

%fft
cfg              = [];
cfg.method       = 'mtmfft';
cfg.foilim       = [20 100]; % cfg.timwin determines spacing
cfg.timwin       = [-0.05 0.05]; % time window of 100 msec
cfg.taper        = 'hanning';
cfg.spikechannel = spike.label{1};
cfg.channel      = data_lfp.label;
stsFFT  = ft_spiketriggeredspectrum(cfg, data_all);

%convol
cfg           = [];
cfg.method    = 'convol';
cfg.foi       = 20:10:100;
cfg.t_ftimwin = 5./cfg.foi; % 5 cycles per frequency
cfg.taper     = 'hanning';
stsConvol = ft_spiketriggeredspectrum(cfg, data_all);

%% spike phase stats

% pairwise phase consistency over complete trial

for k = 1:length(stsConvol.label)
 
  % compute the statistics on the phases
  cfg               = [];
  cfg.method        = 'ppc0'; % compute the Pairwise Phase Consistency
  excludeChan       = str2num(stsConvol.label{k}(6)); % exclude the same channel
  chan              = true(1,4);
  chan(excludeChan) = false;
  cfg.spikechannel  = stsConvol.label{k};
  cfg.channel       = stsConvol.lfplabel(chan); % selected LFP channels
  cfg.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
  cfg.timwin        = 'all'; % compute over all available spikes in the window 
  cfg.latency       = [0.3 nanmax(stsConvol.trialtime(:))]; % sustained visual stimulation period
  statSts = ft_spiketriggeredspectrum_stat(cfg,stsConvol);
 
%   % plot the results
%   figure(1)
%   plot(statSts.freq,statSts.ppc0')
%   hold on
%   xlabel('frequency')
%   ylabel('PPC')
end


% pairwise phase consistency as funciton of time
param              = 'plv'; % set the desired parameter
for k = 1:length(stsConvol.label)
  cfg                = [];
  cfg.method         = param;
  excludeChan        = str2num(stsConvol.label{k}(6)); % this gives us the electrode number of the unit
  chan = true(1,4);
  chan(excludeChan)  = false;
  cfg.spikechannel   = stsConvol.label{k};
  cfg.channel        = stsConvol.lfplabel(chan);
  cfg.avgoverchan    = 'unweighted';
  cfg.winstepsize    = 0.01; % step size of the window that we slide over time
  cfg.timwin         = 0.5; % duration of sliding window
  statSts = ft_spiketriggeredspectrum_stat(cfg,stsConvol);
 
  statSts.(param) = permute(conv2(squeeze(statSts.(param)), ones(1,20)./20, 'same'),[3 1 2]); % apply some   smoothing over 0.2 sec.  
 
  figure, 
  cfg            = [];
  cfg.parameter  = param;
  cfg.refchannel = statSts.labelcmb{1,1};
  cfg.channel    = statSts.labelcmb{1,2};
  cfg.xlim       = [-1 2];
  ft_singleplotTFR(cfg, statSts)    
end
