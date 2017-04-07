function test_bug2770

% WALLTIME 00:10:00
% MEM 1500mb

% TEST eeglab2fieldtrip

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2770/164_MIST_prac.mat');

% the *.set file is actually a MATLAB file with the EEG structure in it
% but it require EEGLAB to read it together with the ftd (which has the binary data)

% I imported the data in EEGLAB with the GUI, and then saved it to a *.mat file
load(dccnpath(filename))

% the eeglab2fieldtrip function is maintained
ft_hastoolbox('eeglab', 1);

data = eeglab2fieldtrip(EEG, 'preprocessing', 'none');

cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.foilim = [4 7];
cfg.output = 'pow';
freq = ft_freqanalysis(cfg,data); % this failed

cfg = [];
cfg.trials = 'all';
cfg.channel = 'all';
data1 = ft_selectdata(cfg, data);
assert(isfield(data1, 'trial')); % this failed

% this is actually at the core of the problem
assert(numel(data.label)==size(data.trial{1},1));

