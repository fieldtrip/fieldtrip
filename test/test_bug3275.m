function test_bug3275

% WALLTIME 00:10:00
% MEM 2gb

% TEST ft_sourceplot ft_checkdata

mri = ft_read_mri('single_subj_T1_1mm.nii');
elec = ft_read_sens('standard_1020.elc'); % this is in MNI space

%%

timelock = [];
timelock.label = elec.label;
timelock.avg = randn(numel(elec.label),1000);
timelock.time = (1:1000)/1000;
timelock.elec = elec;

source = ft_checkdata(timelock, 'datatype', 'source');

% now convert on the fly
cfg = [];
cfg.method = 'cloud';
cfg.funparameter = 'avg';
ft_sourceplot(cfg, timelock, mri)

%%

freq = [];
freq.label = elec.label;
freq.powspctrm = randn(numel(elec.label),30);
freq.freq = 1:30;
freq.elec = elec;

source = ft_checkdata(freq, 'datatype', 'source');

% now convert on the fly
cfg = [];
cfg.method = 'cloud';
cfg.funparameter = 'powspctrm';
ft_sourceplot(cfg, freq)

%%

freq = [];
freq.label = elec.label;
freq.powspctrm = randn(numel(elec.label),30,100);
freq.freq = 1:30;
freq.time = linspace(0, 1, 100);
freq.elec = elec;

source = ft_checkdata(freq, 'datatype', 'source');

% now convert on the fly
cfg = [];
cfg.method = 'cloud';
cfg.funparameter = 'powspctrm';
ft_sourceplot(cfg, freq)




