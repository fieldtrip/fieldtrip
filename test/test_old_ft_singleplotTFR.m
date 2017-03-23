function test_old_ft_singleplotTFR

% MEM 1gb
% WALLTIME 00:10:00


% this script tests the functionality of ft_singleplotTFR with respect to the 
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
cfg.channel = 'Pz';

%create frequency-data
cfgf = [];
cfgf.method = 'mtmconvol';
cfgf.taper  = 'hanning';
cfgf.output = 'pow';
cfgf.channel = data.label(5:end);
cfgf.toi    = [0.5 0.6];
cfgf.foi    = [10 20 30];
cfgf.t_ftimwin = [0.5 0.5 0.5];
freq = ft_freqanalysis(cfgf, data);

%plot frequency-data
ft_singleplotTFR(cfg, freq);

%create connectivity-data
cfgf.output = 'fourier';
cfgf.trials = 1:10;
freqx = ft_freqanalysis(cfgf, data);

cfgc2 = [];
cfgc2.method = 'coh';
coh   = ft_connectivityanalysis(cfgc2, freqx);

%plot connectivity-data
cfg.parameter = 'cohspctrm';
cfg.refchannel = 'Cz';
figure;ft_singleplotTFR(cfg, coh);

%create connectivity-data with sparse linear indexing
cfgc2.channelcmb = [repmat({'Cz'},[numel(freqx.label)-1 1]) setdiff(freqx.label,'Cz')];
coh2  = ft_connectivityanalysis(cfgc2, freqx);

%plot
figure;ft_singleplotTFR(cfg, coh2);
