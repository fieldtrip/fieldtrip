function test_tutorial_beamformingextended

% MEM 5gb
% WALLTIME 00:30:00

% TEST test_beamforming_extended
% TEST ft_read_mri ft_redefinetrial ft_freqanalysis ft_volumesegment ft_appenddata ft_selectdata ft_prepare_singleshell ft_sourceanalysis ft_prepare_leadfield ft_prepare_headmodel ft_prepare_sourcemodel ft_plot_vol ft_plot_sens ft_plot_mesh ft_sourceinterpolate ft_sourceplot

datadir = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/sensor_analysis');
mridir = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer_extended');
templatedir  = dccnpath('/home/common/matlab/fieldtrip/template/sourcemodel');

load(fullfile(datadir, 'subjectK.mat'));

% Time windows of interest
data = ft_appenddata([], data_left, data_right);
cfg = [];
cfg.toilim = [-0.8 1.1];
cfg.minlength = 'maxperlen';
data = ft_redefinetrial(cfg, data);


cfg = [];
cfg.toilim = [-0.8 0];
data_bsl = ft_redefinetrial(cfg, data);

cfg.toilim = [0.3 1.1];
data_exp = ft_redefinetrial(cfg, data);

cfg = [];
data_cmb = ft_appenddata(cfg, data_bsl, data_exp);
% code the trial: 0 = baseline, 1 = experimental condition
data_cmb.trialinfo = [zeros(length(data_bsl.trial), 1); ones(length(data_exp.trial), 1)];
%% calculating the cross spectral density matrix
cfg = [];
cfg.method        = 'mtmfft';
cfg.output        = 'fourier'; % add hint: why fourier?
cfg.tapsmofrq     = 15;
cfg.foi           = 55;
cfg.keeptrials    = 'yes';
freq_cmb = ft_freqanalysis(cfg, data_cmb);

cfg = [];
cfg.trials = freq_cmb.trialinfo == 0;
freq_bsl = ft_selectdata(cfg, freq_cmb);
freq_bsl.cumtapcnt = freq_cmb.cumtapcnt(cfg.trials);
freq_bsl.cumsumcnt = freq_cmb.cumsumcnt(cfg.trials);
cfg.trials = freq_cmb.trialinfo == 1;
freq_exp = ft_selectdata(cfg, freq_cmb);
freq_exp.cumtapcnt = freq_cmb.cumtapcnt(cfg.trials);
freq_exp.cumsumcnt = freq_cmb.cumsumcnt(cfg.trials);

%%foward model and lead field

mri = ft_read_mri(fullfile(mridir, 'subjectK.mri'));
% cfg = [];
% [segmentedmri] = ft_volumesegment(cfg, mri);
% 
oldsegmented = load(fullfile(mridir, 'segmentedmri.mat'));
segmentedmri = oldsegmented.segmentedmri;
% 
% % check whether the segmentation gives results of more than 99% consistency
% assert(max(abs(oldsegmented.segmentedmri.gray(:)-segmentedmri.gray(:))) < .01, 'Gray matter segmentation differs from stored data')
% assert(max(abs(oldsegmented.segmentedmri.csf(:)-segmentedmri.csf(:))) < .01, 'CSF segmentation differs from stored data')
% assert(max(abs(oldsegmented.segmentedmri.white(:)-segmentedmri.white(:))) < .01, 'White matter segmentation differs from stored data')
% % transformation should be absolutely identical
% assert(isequal(oldsegmented.segmentedmri.transform, segmentedmri.transform), 'Transform differs from stored data')
% 
% %save segmentedmri segmentedmri
% 
% % add anatomical information to the segmentation
% segmentedmri.transform = mri.transform;
% segmentedmri.anatomy   = mri.anatomy;
% % call ft_sourceplot
% cfg = [];
% cfg.funparameter = 'gray';
% ft_sourceplot(cfg,segmentedmri);

% create ht head model from the segmented brain surface
cfg = [];
cfg.method = 'singleshell';
hdm = ft_prepare_headmodel(cfg, segmentedmri);


template = load(fullfile(templatedir, 'standard_sourcemodel3d8mm'));
% inverse-warp the subject specific grid to the template grid cfg = [];
cfg.grid.warpmni   = 'yes';
cfg.grid.template  = template.sourcemodel;
cfg.grid.nonlinear = 'yes'; % use non-linear normalization
cfg.mri            = mri;
sourcemodel        = ft_prepare_sourcemodel(cfg);


figure;
hold on;
% note that when calling different plotting routines, all objects that we plot
% need to be in the same unit and coordinate space, here, we need to transform
% the head model to 'cm'
ft_plot_vol(ft_convert_units(hdm, freq_cmb.grad.unit), 'edgecolor', 'none');
alpha 0.4;
ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:));
ft_plot_sens(freq_cmb.grad);

cfg         = [];
cfg.grid    = sourcemodel;
cfg.headmodel = hdm;
cfg.channel = {'MEG'};
cfg.grad    = freq_cmb.grad;
sourcemodel_lf     = ft_prepare_leadfield(cfg, freq_cmb);


%% contrasting source activity
cfg                   = [];
cfg.frequency         = freq_cmb.freq;
cfg.grad              = freq_cmb.grad;
cfg.method            = 'dics';
cfg.keeptrials        = 'yes';
cfg.grid              = sourcemodel_lf;
cfg.headmodel         = hdm;
cfg.keeptrials        = 'yes';
cfg.dics.lambda       = '5%';
cfg.dics.keepfilter   = 'yes';
cfg.dics.fixedori     = 'yes';
cfg.dics.realfilter   = 'yes';
source  = ft_sourceanalysis(cfg, freq_cmb);

% beam pre- and poststim by using the common filter
cfg.grid.filter   = source.avg.filter;
source_bsl  = ft_sourceanalysis(cfg, freq_bsl);
source_exp  = ft_sourceanalysis(cfg, freq_exp);

source_diff = source_exp;
source_diff.avg.pow = (source_exp.avg.pow ./ source_bsl.avg.pow) - 1;
source_diff.pos = template.sourcemodel.pos;
source_diff.dim = template.sourcemodel.dim;

% note that the exact directory is user-specific
templatefile = dccnpath('/home/common/matlab/fieldtrip/external/spm8/templates/T1.nii');
template_mri = ft_read_mri(templatefile);
template_mri.coordsys = 'spm';

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
ft_sourceplot(cfg,source_diff_int);

cfg.method = 'ortho';
cfg.atlas           = dccnpath('/home/common/matlab/fieldtrip/template/atlas/aal/ROI_MNI_V4.nii');
ft_sourceplot(cfg,source_diff_int);

cfg.method = 'surface';
cfg.projmethod     = 'nearest';
cfg.surffile       = 'surface_white_both.mat';
cfg.surfdownsample = 10;
ft_sourceplot(cfg,source_diff_int);


%% coherence beaming

data = ft_appenddata([], data_left, data_right);
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

cfg                 = [];
cfg.method          = 'dics';
cfg.refchan         = 'EMGlft';
cfg.frequency       = 20;
cfg.headmodel       = hdm;
cfg.grid            = sourcemodel;
source_coh_lft      = ft_sourceanalysis(cfg, freq_csd);

source_coh_lft.pos = template.sourcemodel.pos;
source_coh_lft.dim = template.sourcemodel.dim;

% note that the exact directory is user-specific
templatefile = dccnpath('/home/common/matlab/fieldtrip/external/spm8/templates/T1.nii');
template_mri = ft_read_mri(templatefile);
template_mri.coordsys = 'spm';

cfg              = [];
cfg.parameter    = 'coh';
cfg.interpmethod = 'nearest';
source_coh_int   = ft_sourceinterpolate(cfg, source_coh_lft, template_mri);

cfg              = [];
cfg.method       = 'ortho';
cfg.funparameter = 'coh';
cfg.funcolormap = 'jet';

cfg.funcolorlim   = [00 .15];
cfg.opacitylim    = [00 .15];

cfg.maskparameter = cfg.funparameter;
cfg.opacitymap    = 'rampup';

cfg.atlas         = dccnpath('/home/common/matlab/fieldtrip/template/atlas/aal/ROI_MNI_V4.nii');

ft_sourceplot(cfg, source_coh_int);

