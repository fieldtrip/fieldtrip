function test_bug3275

% WALLTIME 00:10:00
% MEM 2gb

% TEST ft_sourceplot ft_checkdata

mri = ft_read_mri('single_subj_T1.nii');
elec = ft_read_sens('standard_1020.elc'); % this is in MNI space


%%
% See http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=3275#c2

source = ft_checkdata(mri, 'datatype', 'source');

tmp = sqrt(sum(source.pos.^2, 2));
pow = repmat(max(tmp)-tmp, 1, 10);
pow = pow + randn(size(pow))*10;

functional = source;
functional.anatomydimord = 'pos';
functional.pow = pow;
functional.powdimord = 'pos_time';
functional.time = (1:10)/10;
functional.inside = functional.anatomy>10000; % very coarse segmentation

%%

cfg = [];
cfg.funparameter = 'pow';
cfg.funcolorlim = 'zeromax';

cfg.method = 'ortho';
ft_sourceplot(cfg, functional)
% none of the others supports time/freq
% cfg.method = 'slice';
% ft_sourceplot(cfg, functional)
% cfg.method = 'surface';
% ft_sourceplot(cfg, functional)
% cfg.method = 'glassbrain';
% ft_sourceplot(cfg, functional)
% cfg.method = 'vertex';
% ft_sourceplot(cfg, functional)
% cfg.method = 'cloud';
% ft_sourceplot(cfg, functional)

%%

cfg.latency = 0.5;
cfg.method = 'ortho';
ft_sourceplot(cfg, functional)
cfg.method = 'slice';
ft_sourceplot(cfg, functional)
cfg.method = 'surface';
ft_sourceplot(cfg, functional)
cfg.method = 'glassbrain';
ft_sourceplot(cfg, functional)
% these are too inefficient for volumetric data
% cfg.method = 'vertex';
% ft_sourceplot(cfg, functional)
% cfg.method = 'cloud';
% ft_sourceplot(cfg, functional)


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


