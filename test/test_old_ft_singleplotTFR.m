function test_old_ft_singleplotTFR

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY

% this script tests the functionality of ft_singleplotTFR with respect to the 
% different input datatypes. no other functionality is tested.
% the script has been written in order to test a clean up of the code

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/eeg/preproc_neuroscan16');
load(filename)

%there's an unresolved issue with duplicate labels 'FREE'
%FIXME
data.label{1} = 'FREE1';
data.label{2} = 'FREE2';
data.label{3} = 'FREE3';
data.label{4} = 'FREE4';

% create frequency-data
cfgf = [];
cfgf.method = 'mtmconvol';
cfgf.taper  = 'hanning';
cfgf.output = 'pow';
cfgf.channel = data.label(5:end);
cfgf.toi    = [0.5 0.6];
cfgf.foi    = [10 20 30];
cfgf.t_ftimwin = [0.5 0.5 0.5];
freq = ft_freqanalysis(cfgf, data);

% plot frequency-data
cfg = [];
cfg.channel = 'Pz';
ft_singleplotTFR(cfg, freq);

cfg = [];
cfg.channel = 'all';
cfg.layout = 'elec1005';
cfg.highlight = 'labels';
cfg.highlightchannel = {'Cz'};
cfg.markers = 'numbers';
ft_topoplotTFR(cfg, freq);


% create connectivity-data
cfgf.output = 'fourier';
cfgf.trials = 1:10;
freqx = ft_freqanalysis(cfgf, data);

cfgc2 = [];
cfgc2.method = 'coh';
coh   = ft_connectivityanalysis(cfgc2, freqx);

% plot connectivity-data
cfg = [];
cfg.channel = 'Pz';
cfg.parameter = 'cohspctrm';
cfg.refchannel = 'Cz';
ft_singleplotTFR(cfg, coh);

cfg = [];
cfg.channel = 'all';
cfg.parameter = 'cohspctrm';
% cfg.refchannel = 'Cz';
cfg.refchannel = 'gui';
cfg.layout = 'elec1005';
ft_multiplotTFR(cfg, coh);


% create connectivity-data with sparse linear indexing
cfgc2.channelcmb = [repmat({'Cz'},[numel(freqx.label)-1 1]) setdiff(freqx.label,'Cz')];
coh2  = ft_connectivityanalysis(cfgc2, freqx);

% plot
cfg = [];
cfg.channel = 'Pz';
cfg.parameter = 'cohspctrm';
cfg.refchannel = 'Cz';
ft_singleplotTFR(cfg, coh2);

cfg = [];
cfg.channel = 'all';
cfg.parameter = 'cohspctrm';
% cfg.refchannel = 'Cz';
cfg.refchannel = 'gui';
cfg.layout = 'elec1005';
ft_multiplotTFR(cfg, coh2);


%%

cfg.layout = 'elec1005'
