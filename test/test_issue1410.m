function test_issue1410

% MEM 6gb
% WALLTIME 00:20:00
% DEPENDENCY

load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf151.mat'), 'data');
grad = data.grad;
clear data

%%

vol = [];
vol.o = [0 0 4];
vol.r = 12;
vol.unit = 'cm';
vol.type = 'singlesphere';

%%

cfg = [];
cfg.dip.pos = [0 0 9];
cfg.dip.mom = [1 0 0];
cfg.dip.unit = 'cm';
cfg.headmodel = vol;
cfg.grad = grad;
data = ft_dipolesimulation(cfg);

cfg = [];
timelock = ft_timelockanalysis(cfg, data);

cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.output = 'fourier';
freq = ft_freqanalysis(cfg, data);

%%

tfdata = timelock;
method = 'lcmv';

cfg = [];
cfg.xgrid = -8:2:8;
cfg.ygrid = -8:2:8;
cfg.zgrid = (-8:2:8) + 4;
cfg.unit = 'cm';
cfg.headmodel = vol;
cfg.method = method;
cfg.keepleadfield = 'yes';
cfg.(method).keepfilter = 'yes';
cfg.(method).lambda = '10%';
cfg.channel = 'MEG';
source1 = ft_sourceanalysis(cfg, tfdata);

%%

% in the next code the labels are missing for the precomputed filters
% this should assume that they were computed with the same channel selection

cfg = [];
cfg.method = method;
cfg.sourcemodel.pos = source1.pos;
cfg.sourcemodel.filter = source1.avg.filter;
cfg.channel = 'MEG';
source2 = ft_sourceanalysis(cfg, tfdata);

clear source1 source2

%%
% this is mostly from http://www.fieldtriptoolbox.org/tutorial/beamformingextended/

clear all
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/sensor_analysis/subjectK.mat'));

data_combined = ft_appenddata([], data_left, data_right);

cfg           = [];
cfg.toilim    = [-0.8 1.1];
cfg.minlength = 'maxperlen'; % this ensures all resulting trials are equal length
data          = ft_redefinetrial(cfg, data_combined);

[ftver, ftdir] = ft_version;

%%
% the first analysis is for visual gamma

cfg        = [];
cfg.toilim = [-0.8 0];
data_bsl   = ft_redefinetrial(cfg, data);

cfg.toilim = [0.3 1.1];
data_exp   = ft_redefinetrial(cfg, data);

cfg      = [];
data_cmb = ft_appenddata(cfg, data_bsl, data_exp);

% give a number to each trial: 0 = baseline, 1 = experimental condition
data_cmb.trialinfo = [zeros(length(data_bsl.trial), 1); ones(length(data_exp.trial), 1)];

cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq  = 15;
cfg.foi        = 55;
freq_cmb       = ft_freqanalysis(cfg, data_cmb);

cfg                = [];
cfg.trials         = freq_cmb.trialinfo == 0;
freq_bsl           = ft_selectdata(cfg, freq_cmb);

cfg.trials         = freq_cmb.trialinfo == 1;
freq_exp           = ft_selectdata(cfg, freq_cmb);

%%
% this is used both for visual gamma and beta coherence

load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer_extended/hdm.mat'))
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer_extended/sourcemodel.mat'))

cfg             = [];
cfg.sourcemodel = sourcemodel;
cfg.headmodel   = hdm;
cfg.channel     = {'MEG'};
cfg.grad        = data.grad;
cfg.singleshell.batchsize = 10;
sourcemodel_lf  = ft_prepare_leadfield(cfg);

templatedir = fullfile(ftdir, 'template', 'sourcemodel');
template = load(fullfile(templatedir, 'standard_sourcemodel3d8mm')); % 8mm spacing grid

templatedir = fullfile(ftdir, 'external', 'spm8', 'templates');
template_mri = ft_read_mri(fullfile(templatedir, 'T1.nii'));
template_mri.coordsys = 'mni'; % we know it's in MNI space

%%

cfg                   = [];
cfg.frequency         = freq_cmb.freq;
cfg.grad              = freq_cmb.grad;
cfg.method            = 'dics';
cfg.keeptrials        = 'yes';
cfg.channel           = 'MEG';
cfg.sourcemodel       = sourcemodel_lf;
cfg.keeptrials        = 'yes';
cfg.dics.lambda       = '5%';
cfg.dics.keepfilter   = 'yes';
cfg.dics.fixedori     = 'yes';
cfg.dics.realfilter   = 'yes';
source                = ft_sourceanalysis(cfg, freq_cmb);

% beam pre- and poststim by using the common filter
cfg.sourcemodel.filter  = source.avg.filter;
cfg.sourcemodel.label   = source.avg.label;
source_bsl       = ft_sourceanalysis(cfg, freq_bsl);
source_exp       = ft_sourceanalysis(cfg, freq_exp);

cfg = [];
cfg.parameter = 'avg.pow';
cfg.operation = '(x1 ./ x2) - 1';
source_diff = ft_math(cfg, source_exp, source_bsl);

source_diff.pos = template.sourcemodel.pos;
source_diff.dim = template.sourcemodel.dim;

cfg              = [];
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
source_diff_int  = ft_sourceinterpolate(cfg, source_diff, template_mri);

cfg               = [];
cfg.method        = 'slice';
cfg.funparameter  = 'pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [0.0 1.2];
cfg.opacitylim    = [0.0 1.2];
cfg.opacitymap    = 'rampup';
ft_sourceplot(cfg, source_diff_int);

clear source_bsl source_exp source_diff source_diff_int

%%
% now for the EMG coherence

cfg                 = [];
cfg.toilim          = [-1 -0.0025];
cfg.minlength       = 'maxperlen'; % this ensures all resulting trials are equal length
data_stim           = ft_redefinetrial(cfg, data);

cfg                 = [];
cfg.output          = 'powandcsd';
cfg.method          = 'mtmfft';
cfg.taper           = 'dpss';
cfg.tapsmofrq       = 5;
cfg.foi             = 20;
cfg.keeptrials      = 'yes';
cfg.channel         = {'MEG' 'EMGlft' 'EMGrgt'};
cfg.channelcmb      = {'MEG' 'MEG'; 'MEG' 'EMGlft'; 'MEG' 'EMGrgt'};
freq_csd            = ft_freqanalysis(cfg, data_stim);

%%

cfg                 = [];
cfg.method          = 'dics';
cfg.refchan         = 'EMGlft';
cfg.channel         = 'MEG';
cfg.frequency       = 20;
% with precomputed leadfields
cfg.sourcemodel     = sourcemodel_lf;
source_coh_lft      = ft_sourceanalysis(cfg, freq_csd);

% without precomputed leadfields
cfg.sourcemodel     = removefields(sourcemodel_lf, {'leadfield', 'leadfielddimord', 'label'});
cfg.headmodel       = hdm;
source_coh          = ft_sourceanalysis(cfg, freq_csd);

figure; hold on
plot(source_coh.avg.coh(:), '.')
plot(source_coh_lft.avg.coh(:), 'ro')

assert(isalmostequal(source_coh_lft.avg.coh, source_coh.avg.coh, 'reltol', 1e-6));

%%

cfg                  = [];
cfg.method           = 'dics';
cfg.refchan          = 'EMGlft';
cfg.channel          = 'MEG';
cfg.frequency        = 20;
cfg.headmodel        = hdm;
cfg.sourcemodel.pos  = [3.6987   -3.1192   11.2426];
cfg.sourcemodel.unit = 'cm';
cfg.reducerank       = 2;
cfg.keepleadfield    = 'yes';
cfg.dics.keepfilter  = 'yes';
cfg.backproject      = 'yes';
source3              = ft_sourceanalysis(cfg, freq_csd);
cfg.backproject      = 'no';
source2              = ft_sourceanalysis(cfg, freq_csd);

assert(isalmostequal(source3.avg.coh, source2.avg.coh, 'reltol', 1e-6));
assert(~isequal(source3.leadfield, source2.leadfield)); % should have different number of columns
assert(~isequal(source3.avg.filter, source2.avg.filter)); % should have different number of rows

clear source3 source2

%%

source_coh_lft.pos = template.sourcemodel.pos;
source_coh_lft.dim = template.sourcemodel.dim;

source_coh.pos = template.sourcemodel.pos;
source_coh.dim = template.sourcemodel.dim;

cfg              = [];
cfg.parameter    = 'coh';
cfg.interpmethod = 'nearest';
cfg.coordsys     = 'mni';
source_coh_int   = ft_sourceinterpolate(cfg, source_coh_lft, template_mri);

cfg               = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'coh';
cfg.maskparameter = 'coh';
cfg.headmodel     = hdm;
ft_sourceplot(cfg, source_coh_int);

cfg               = [];
cfg.method        = 'surface';
cfg.funparameter  = 'coh'; % use it to represent color
cfg.maskparameter = 'coh'; % and use it for transparency
ft_sourceplot(cfg, source_coh_int);

%%
% let's also try dics_refdip

% these are approximately in the motor cortex
diplft = [4  3 10];
diprgt = [4 -3 10];

% take the nearest grid point
d = sqrt(sum((sourcemodel_lf.pos - repmat(diplft, size(sourcemodel_lf.pos,1), 1)).^2, 2));
[mind, indx] = min(d);
diplft = sourcemodel_lf.pos(indx,:);

cfg                 = [];
cfg.method          = 'dics';
cfg.refdip          = diplft;
cfg.channel         = 'MEG';
cfg.frequency       = 20;
% with precomputed leadfields
cfg.sourcemodel     = sourcemodel_lf; % the leadfields will be discarded
cfg.headmodel       = hdm;
source_coh_diplft   = ft_sourceanalysis(cfg, freq_csd);

cfg = [];
cfg.funparameter = 'coh'; % use it to represent color
ft_sourceplot(cfg, source_coh_diplft); % it should be 1 on the location of diplft

