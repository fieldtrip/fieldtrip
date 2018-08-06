function test_old_ft_multiplotER

% MEM 1gb
% WALLTIME 00:10:00


% this script tests the functionality of ft_multiplotER with respect to the 
% different input datatypes. no other functionality is tested.
% the script has been written in order to test a clean up of the code

filename = fullfile(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/eeg/'), 'preproc_neuroscan16');
load(filename);

%there's an unresolved issue with duplicate labels 'FREE'
%FIXME
data.label{1} = 'FREE1';
data.label{2} = 'FREE2';
data.label{3} = 'FREE3';
data.label{4} = 'FREE4';

cfg = [];
cfg.channel = data.label(5:end);
cfg.preproc.demean = 'yes';
cfg.trials = 1:5;
tlck1 = ft_timelockanalysis(cfg, data);
cfg.trials = 6:10;
tlck2 = ft_timelockanalysis(cfg, data);

cfg = [];
cfg.layout = 'biosemi64.lay';
lay = ft_prepare_layout(cfg);

cfg        = [];
cfg.layout = lay;

%plot timelocked-data
ft_multiplotER(cfg, tlck1, tlck2);

%create frequency-data
cfgf = [];
cfgf.method = 'mtmfft';
cfgf.foilim = [0 100];
cfgf.taper  = 'hanning';
cfgf.output = 'pow';
cfgf.trials = 1:5;
cfgf.channel = data.label(5:end); % the first 4 are empty
freq1 = ft_freqanalysis(cfgf, data);
cfgf.trials = 6:10;
freq2 = ft_freqanalysis(cfgf, data);

%plot frequency-data
ft_multiplotER(cfg, freq1, freq2);

%create connectivity-data
cfgf.output = 'fourier';
cfgf.trials = 1:10;
freqx = ft_freqanalysis(cfgf, data);

cfgc2 = [];
cfgc2.method = 'coh';
coh   = ft_connectivityanalysis(cfgc2, freqx);

%plot connectivity-data
cfg.refchannel = 'gui';
cfg.parameter = 'cohspctrm';
ft_multiplotER(cfg, coh); % FIXME this causes a crash when a new reference is selected and the old one is not unselected
% FIXME it also crashes when more than one cohref is selected

%create connectivity-data with sparse linear indexing
cfgc2.channelcmb = [repmat(freqx.label(5),[numel(freqx.label)-1 1]) freqx.label([1:4 6:end])];
coh2  = ft_connectivityanalysis(cfgc2, freqx);

%plot
cfg.refchannel = coh2.labelcmb{1,1};
ft_multiplotER(cfg, coh2);

%create connectivity-data with very sparse linear indexing
cfgc2.channelcmb = cfgc2.channelcmb(1:25,:);
coh3   = ft_connectivityanalysis(cfgc2, freqx);

%plot
ft_multiplotER(cfg, coh3);

%create connectivity-data with even sparser linear indexing
% cfgc2.channelcmb = [repmat(freq2.label(5),[10 1]) freq2.label(21:30)';repmat(freq2.label(10),[10 1]) freq2.label(21:30)'];
% coh4 = ft_connectivityanalysis(cfgc2, freq2);
% 
% %plot: this breaks
% cfg.cohrefchannel = 'gui';
% ft_topoplotER(cfg, coh4);
% 
% %plot: this works
% cfg.cohrefchannel = coh4.labelcmb(1,1);
% ft_topoplotER(cfg, coh4);
% 
% %create connectivity-data with asymmetry
% %the data are probably not full-rank creating a problem for the sf
% %freq3 = ft_selectdata(freq2, 'channel', freq2.label(1:30));
% %%subselection of channels does not help
% freq4 = freq2transfer([], freq2);
% 
% cfgc2 = [];
% cfgc2.method = 'granger';
% granger = ft_connectivityanalysis(cfgc2, freq4);
% 
% cfg.cohrefchannel = 'gui';
% cfg.zparam = 'grangerspctrm';
% ft_topoplotER(cfg, granger);
