function test_bug3275

% WALLTIME 00:10:00
% MEM 2gb

% TEST ft_sourceplot ft_checkdata

mri = ft_read_mri('single_subj_T1_1mm.nii');
elec = ft_read_sens('standard_1020.elc'); % this is in MNI space


%%
% See http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=3275#c2

functional = mri;
functional.pow = rand([181 217 181 10]);
functional.powdimord = 'dim1_dim2_dim3_time';
functional.time = (1:10)/10;
functional.inside = functional.anatomy>10000; % very coarse segmentation


cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'pow';
cfg.funcolorlim = [0 1];
% cfg.latency = 0.5;
ft_sourceplot(cfg, functional)


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
cfg.latency = 0.5;
ft_sourceplot(cfg, timelock)

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
cfg.frequency = 15;
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
cfg.frequency = 15;
cfg.latency = 0.5;
ft_sourceplot(cfg, freq)




