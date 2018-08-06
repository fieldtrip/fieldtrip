function test_ft_topoplotER

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_topoplotER ft_topoplotTFR ft_topoplotIC

% this script tests the functionality of ft_topoplotER with respect to the 
% different input datatypes. no other functionality is tested.
% the script has been written in order to test a clean up of the code

pwdir = pwd;

load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/eeg/preproc_neuroscan16'));

%there's an unresolved issue with duplicate labels 'FREE'
%FIXME
data.label{1} = 'FREE1';
data.label{2} = 'FREE2';
data.label{3} = 'FREE3';
data.label{4} = 'FREE4';

cfg = [];
cfg.preproc.demean = 'yes';
tlck = ft_timelockanalysis(cfg, data);

cfg = [];
cfg.layout = 'biosemi64.lay';
lay = ft_prepare_layout(cfg);

cfg        = [];
cfg.layout = lay;
cfg.interactive = 'yes';

%plot timelocked-data
figure;ft_topoplotER(cfg, tlck);drawnow

%reproduce bug 1239
cfgtmp = cfg;
cfgtmp.highlight = 'on';
cfgtmp.highlightchannel = {};
figure;ft_topoplotER(cfgtmp, tlck);drawnow

%plot subplots
cfg.xlim = (-1:0.25:1);
figure;ft_topoplotER(cfg, tlck);drawnow
cfg = rmfield(cfg, 'xlim');

%create component-data
cfgc = [];
cfgc.method = 'pca';
comp = ft_componentanalysis(cfgc, data);

%plot component-data
cfg.interactive = 'no';
cfg.component   = 1:10;
figure;ft_topoplotIC(cfg, comp);drawnow

cfg.component   = 5;
figure;ft_topoplotIC(cfg, comp);drawnow

%create frequency-data
cfgf = [];
cfgf.method = 'mtmfft';
cfgf.foilim = [0 100];
cfgf.taper  = 'hanning';
cfgf.output = 'pow';
freq = ft_freqanalysis(cfgf, data);

%plot frequency-data
cfg.interactive = 'yes';
%cfg = rmfield(cfg, 'component');
figure;ft_topoplotER(cfg, freq);drawnow

%create connectivity-data
cfgf.output = 'fourier';
freq2 = ft_freqanalysis(cfgf, data);

cfgc2 = [];
cfgc2.method = 'coh';
cfgc2.channel = freq2.label(5:end);
coh   = ft_connectivityanalysis(cfgc2, freq2);

%plot connectivity-data
cfg.refchannel = 'gui';
cfg.parameter = 'cohspctrm';
figure;ft_topoplotER(cfg, coh);drawnow % FIXME this causes a crash when a new reference is selected and the old one is not unselected
% FIXME it also crashes when more than one ref is selected

%create connectivity-data with sparse linear indexing
cfgc2.channelcmb = [repmat(freq2.label(5),[numel(freq2.label)-1 1]) freq2.label([1:4 6:end])];
coh2  = ft_connectivityanalysis(cfgc2, freq2);

%plot
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.refchannel = 'gui';
cfg.parameter = 'cohspctrm';
cfg.refchannel = coh2.labelcmb{1,1};
figure;ft_topoplotER(cfg, coh2);drawnow

%create connectivity-data with very sparse linear indexing
cfgc2.channelcmb = cfgc2.channelcmb(1:25,:);
coh3   = ft_connectivityanalysis(cfgc2, freq2);

%plot
figure;ft_topoplotER(cfg, coh3);drawnow

%create connectivity-data with even sparser linear indexing
cfgc2.channelcmb = [repmat(freq2.label(5),[10 1]) freq2.label(21:30);repmat(freq2.label(10),[10 1]) freq2.label(21:30)];
coh4 = ft_connectivityanalysis(cfgc2, freq2);

%plot: this breaks
cfg.refchannel = 'gui';
figure;ft_topoplotER(cfg, coh4);drawnow

%plot: this works
cfg.refchannel = coh4.labelcmb(1,1);
figure;ft_topoplotER(cfg, coh4);drawnow

%create connectivity-data with asymmetry
%the data are probably not full-rank creating a problem for the sf
%freq3 = ft_selectdata(freq2, 'channel', freq2.label(1:30));
%%subselection of channels does not help

cfgc2 = [];
cfgc2.method = 'granger';
granger = ft_connectivityanalysis(cfgc2, freq2);

cfg.refchannel = 'gui';
cfg.parameter = 'grangerspctrm';
figure;ft_topoplotER(cfg, granger);drawnow

%plot a stat structure, containing only 1 freq bin (so without freq-field)
stat = freq;
stat.stat = freq.powspctrm(:,10);
stat = rmfield(stat, 'freq');
if isfield(stat, 'cumtapcnt'), 
  stat = rmfield(stat, 'cumtapcnt');
end
stat = rmfield(stat, 'powspctrm');
stat.dimord = 'chan';

cfg = rmfield(cfg, 'refchannel');
cfg.parameter = 'stat';
cfg.interactive = 'no';
figure;ft_topoplotER(cfg, stat);drawnow

cd(pwdir);
