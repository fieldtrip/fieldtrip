function failed_tutorial_meguk2015

% WALLTIME 01:30:00
% MEM 8gb

% TEST test_tutorial_meguk2015
% TEST ft_preprocessing ft_appenddata ft_timelockanalysis ft_componentanalysis ft_rejectcomponent ft_freqanalysis ft_math ft_sourceanalysis ft_sourceinterpolate ft_volumerealign ft_determine_coordsys ft_sourceplot ft_timelockstatistics ft_multiplotER ft_multiplotTFR ft_singleplotER ft_singleplotTFR

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

% This script is the concatenation of the various parts of the MEG-UK 2015
% workshop demonstration. It has been editted to make it smaller and faster
% and to avoid parts with external dependencies.

datadir = dccnpath('/home/common/matlab/fieldtrip/data/test/original/rikhenson');
subj    = 15;

%% read the data from all separate runs

% this will contain the runs for a single subject
rundata = {};

for run=1 % the original uses 1:6, but that takes too long and too much storage
  trialdef = fullfile(datadir, sprintf('Sub%02d', subj), 'MEEG', 'Trials', sprintf('run_%02d_trldef.txt', run));
  dataset  = fullfile(datadir, sprintf('Sub%02d', subj), 'MEEG', sprintf('run_%02d_sss.fif', run));
  
  [begsample, endsample, offset, trialtype] = textread(trialdef, '%d%d%d%s');
  
  trialcode = nan(size(trialtype));
  trialcode(strcmp(trialtype, 'Famous'))      = 1;
  trialcode(strcmp(trialtype, 'Unfamiliar'))  = 2;
  trialcode(strcmp(trialtype, 'Scrambled'))   = 3;
  
  % construct the trial definition matrix, usually done with FT_DEFINETRIAL
  trl = [begsample(:) endsample(:) offset(:) trialcode(:)];
  
  cfg         = [];
  cfg.dataset = dataset;
  cfg.trl     = trl;
  
  % MEG specific settings
  cfg.channel = 'MEG';
  cfg.demean  = 'yes';
  data_meg = ft_preprocessing(cfg);
  
  % EEG specific settings
  cfg.channel    = 'EEG';
  cfg.demean     = 'yes';
  cfg.reref      = 'yes';
  cfg.refchannel = 'all'; % average reference
  data_eeg = ft_preprocessing(cfg);
  
  % settings for all other channels
  cfg.channel = {'all', '-MEG', '-EEG'};
  cfg.demean  = 'no';
  cfg.reref   = 'no';
  data_other = ft_preprocessing(cfg);
  
  cfg = [];
  cfg.resamplefs = 300;
  data_meg   = ft_resampledata(cfg, data_meg);
  data_eeg   = ft_resampledata(cfg, data_eeg);
  data_other = ft_resampledata(cfg, data_other);
  
  %% append the different channel sets into a single structure
  
  rundata{run} = ft_appenddata(cfg, data_meg, data_eeg, data_other);
  clear data_meg data_eeg data_other
  
end % for each run


%% append the 6 runs into a single structure
% data = ft_appenddata(cfg, rundata{:});
data = rundata{1};

%% compute the overall average and condition-specific averages

cfg = [];
cfg.trials = find(data.trialinfo==1);
avg_Famous = ft_timelockanalysis(cfg, data);
cfg.trials = find(data.trialinfo==2);
avg_Unfamiliar = ft_timelockanalysis(cfg, data);

cfg.trials = find(data.trialinfo==3);
avg_Scrambled = ft_timelockanalysis(cfg, data);

cfg.trials = find(data.trialinfo==1 | data.trialinfo==2);
avg_Faces = ft_timelockanalysis(cfg, data);

cfg = [];
% cfg.layout = 'neuromag306all.lay';
cfg.layout = 'neuromag306mag.lay';
ft_multiplotER(cfg, avg_Faces, avg_Scrambled);


%% compute the difference between faces and

cfg = [];
cfg.parameter = 'avg';
cfg.operation = 'x1-x2';
avg_Faces_vs_Scrambled = ft_math(cfg, avg_Faces, avg_Scrambled);

cfg        = [];
% cfg.layout = 'neuromag306all.lay';
cfg.layout = 'neuromag306mag.lay';
ft_multiplotER(cfg, avg_Faces_vs_Scrambled);

cfg           = [];
cfg.layout    = 'neuromag306mag.lay';
cfg.colorbar  = 'yes';
figure; ft_movieplotER(cfg, avg_Faces_vs_Scrambled);

% for saving to disk
prefix = sprintf('Sub%02d', subj);

%% save the raw data to disk
% save([prefix '_raw'], 'data');

%% save the averages to disk
prefix = sprintf('/tmp/Sub%02d', subj);
% save([prefix '_avg_Famous'],     'avg_Famous');
% save([prefix '_avg_Unfamiliar'], 'avg_Unfamiliar');
% save([prefix '_avg_Scrambled'],  'avg_Scrambled');
% save([prefix '_avg_Faces'],      'avg_Faces');
% save([prefix '_avg_Faces_vs_Scrambled'], 'avg_Faces_vs_Scrambled');

%% look at the analysis history
cfg           = [];
cfg.filename  = [prefix '_avg_Faces_vs_Scrambled.html'];
ft_analysispipeline(cfg, avg_Faces_vs_Scrambled);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% stats_part2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load the raw data from disk
% prefix = sprintf('Sub%02d', subj);
% load([prefix '_raw']);
% load([prefix '_avg_Faces_vs_Scrambled']);

%% reorganize the timelocked data and compute stats

cfg = [];
cfg.channel    = 'MEGMAG';
cfg.keeptrials = 'yes';
timelock = ft_timelockanalysis(cfg, data);

cfg           = [];
cfg.correctm  = 'no';
cfg.method    = 'analytic';
cfg.statistic = 'indepsamplesT';        % this is implemented in ft_statfun_indepsamplesT
cfg.design    = nan(1, size(timelock.trialinfo,1));
cfg.design(timelock.trialinfo==1) = 1;  % Famous faces
cfg.design(timelock.trialinfo==2) = 1;  % Unfamiliar faces
cfg.design(timelock.trialinfo==3) = 2;  % Scrambled
cfg.ivar      = 1;                      % the first (and only) row of the design represents the independent variable
analytic = ft_timelockstatistics(cfg, timelock);


%% do some sanity checks
figure
imagesc(analytic.time, 1:length(analytic.label), -log10(analytic.prob))
colorbar

cfg = [];
cfg.channel = analytic.label;
tmp = ft_selectdata(cfg, avg_Faces_vs_Scrambled);

analytic.avg = tmp.avg;

% analytic.logprob = -log10(analytic.prob);
% analytic.logprob(isnan(analytic.logprob)) = 0;
% analytic.logprob(isinf(analytic.logprob)) = 10;

% save analytic analytic

cfg               = [];
cfg.layout        = 'neuromag306mag.lay';
cfg.parameter     = 'avg';
cfg.maskparameter = 'mask';
ft_multiplotER(cfg, analytic);

%% use montecarlo and correctm=max

cfg                  = [];
cfg.correctm         = 'max';
cfg.method           = 'montecarlo';
cfg.numrandomization = 1000;
cfg.statistic        = 'indepsamplesT'; % this is implemented in ft_statfun_indepsamplesT
cfg.design           = nan(1, size(timelock.trialinfo,1));
cfg.design(timelock.trialinfo==1) = 1; % Famous faces
cfg.design(timelock.trialinfo==2) = 1; % Unfamiliar faces
cfg.design(timelock.trialinfo==3) = 2; % Scrambled
cfg.ivar             = 1; % the first (and only) row of the design represents the independent variable
cfg.latency          = [0.140 0.180];
montecarlo = ft_timelockstatistics(cfg, timelock);

% save montecarlo montecarlo

figure
hist([montecarlo.negdistribution' montecarlo.posdistribution'], 100)

%% compare the observed statistical values to the distributions

negdistribution = sort(montecarlo.negdistribution);
negthreshold    = negdistribution(26)   % why not at 5%, i.e. 51;

posdistribution = sort(montecarlo.posdistribution);
posthreshold    = posdistribution(975)  % why not at 5%, i.e. 950;

figure
subplot(3,1,1)
hist(montecarlo.negdistribution, 50)
ylabel('negdist');
xlim([-10 10]);
subplot(3,1,2)
hist(montecarlo.stat(:), 100)
ylabel('observed stat');
xlim([-10 10]);
subplot(3,1,3)
hist(montecarlo.posdistribution, 50)
ylabel('posdist');
xlim([-10 10]);


%% use your own trialfunction, e.g. spearman rank correlation
%% this part cannot be tested, since the statfun is not available
% cfg                  = [];
% cfg.channel          = 'MEG2021';
% cfg.statistic        = 'statfun_parametric';
% cfg.design           = nan(1, size(timelock.trialinfo,1));
% cfg.design(timelock.trialinfo==1) = 1; % Famous faces
% cfg.design(timelock.trialinfo==2) = 2; % Unfamiliar faces
% cfg.design(timelock.trialinfo==3) = 3; % Scrambled
% cfg.ivar             = 1; % the first (and only) row of the design represents the independent variable
%
% cfg.correctm         = 'no'; % or another method
% cfg.method           = 'analytic';
% analytic2 = ft_timelockstatistics(cfg, timelock);
%
% cfg.correctm         = 'max';
% cfg.method           = 'montecarlo';
% cfg.numrandomization = 1000;
% montecarlo2 = ft_timelockstatistics(cfg, timelock);
%
% figure
% hold on
% plot(analytic2.time,   -log10(analytic2.prob),   'b')
% plot(montecarlo2.time, -log10(montecarlo2.prob), 'r')
% line([montecarlo2.time(1) montecarlo2.time(end)], [1.3 1.3])

% save analytic2 analytic2
% save montecarlo2 montecarlo2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% beamformer_part1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% get data from SPM
%
% D = spm_eeg_load('../spm-source-demo/PapMcbdspmeeg_run_01_sss.mat');
%
% % convert sensors and volume conduction model from SPM
% volsens = spm_eeg_inv_get_vol_sens(D, 1, 'Head', 'inv', 'MEG');
% vol1    = volsens.MEG.vol;
% sens1   = volsens.MEG.sens;
% mri1    = ft_read_mri('../spm-source-demo/mprage.nii');

%% start from scratch data in FieldTrip

% subj = 15;
% prefix = sprintf('Sub%02d', subj);
% load([prefix '_raw']);  % this is called "data" rather than "raw"

sens = data.grad;

% load the original MRI
mri_orig = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/test/original/rikhenson/Sub15/T1/mprage.nii'));

% load the positions of the anatomical fiducials (as provided by Rik)
load(dccnpath('/home/common/matlab/fieldtrip/data/test/original/rikhenson/Sub15/T1/mri_fids.mat'));

headshape = ft_read_headshape(dccnpath('/home/common/matlab/fieldtrip/data/test/original/rikhenson/Sub15/MEEG/run_01_raw.fif'));
headshape = ft_convert_units(headshape, 'mm');

% the MRI is neither expressed in MNI, nor in Neuromag coordinates
ft_determine_coordsys(mri_orig, 'interactive', 'no');
hold on; % add the subsequent objects to the same figure
ft_plot_headshape(headshape);
plot3(mri_fids(1,1), mri_fids(1,2), mri_fids(1,3), 'm*');
plot3(mri_fids(2,1), mri_fids(2,2), mri_fids(2,3), 'm*');
plot3(mri_fids(3,1), mri_fids(3,2), mri_fids(3,3), 'm*');

%% validate the positions of the fiducials that were provided by Rik

cfg = [];
cfg.location = mri_fids(1,:);
ft_sourceplot(cfg, mri_orig);

cfg = [];
cfg.location = mri_fids(2,:);
ft_sourceplot(cfg, mri_orig);

cfg = [];
cfg.location = mri_fids(3,:);
ft_sourceplot(cfg, mri_orig);

%%

% the location of fiducials is expressed in original MRI coordinates
% ft_volumerealign needs them in voxel coordinates
vox_fids = ft_warp_apply(inv(mri_orig.transform), mri_fids);

cfg = [];
cfg.fiducial.nas = vox_fids(1,:);
cfg.fiducial.lpa = vox_fids(2,:);
cfg.fiducial.rpa = vox_fids(3,:);
cfg.coordsys = 'neuromag';
mri_realigned = ft_volumerealign(cfg, mri_orig);

% save mri_realigned mri_realigned

% check that the MRI is consistent after realignment
ft_determine_coordsys(mri_realigned, 'interactive', 'no');
hold on; % add the subsequent objects to the figure
ft_plot_headshape(headshape);

%%

cfg = [];
cfg.output = {'brain' 'scalp' 'skull'};
seg = ft_volumesegment(cfg, mri_realigned);

% save seg seg

%%

cfg = [];
cfg.method = 'projectmesh';
cfg.numvertices = 2000;
cfg.tissue = 'brain';
brain = ft_prepare_mesh(cfg, seg);
cfg.tissue = 'skull';
skull = ft_prepare_mesh(cfg, seg);
cfg.tissue = 'scalp';
scalp = ft_prepare_mesh(cfg, seg);

% save brain brain
% save skull skull
% save scalp scalp

%% make the volume conduction model

cfg = [];
cfg.method = 'singleshell';
vol = ft_prepare_headmodel(cfg, brain);

% save vol vol
% save sens sens

ft_determine_coordsys(mri_realigned, 'interactive', 'no')
hold on; % add the subsequent objects to the same figure
ft_plot_headshape(headshape);
ft_plot_vol(ft_convert_units(vol, 'mm'));

figure
hold on; % add the subsequent objects to the same figure
ft_plot_headshape(headshape);
ft_plot_sens(ft_convert_units(sens, 'mm'), 'coil', 'yes', 'coildiameter', 10);
ft_plot_vol(ft_convert_units(vol, 'mm'));

% figure
% ft_plot_vol(ft_convert_units(vol,  'mm'), 'facecolor', 'r'); % FT
% ft_plot_vol(ft_convert_units(vol1, 'mm'), 'facecolor', 'g'); % SPM
% alpha 0.5


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% beamformer_part2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% get data from SPM
%
% D = spm_eeg_load('../spm-source-demo/PapMcbdspmeeg_run_01_sss.mat');
%
% disp(D.condlist)
%
% % convert data from SPM
% raw      = D.ftraw(D.indchantype('MEGMAG'), D.indsample(-0.1):D.indsample(0.3), D.indtrial(D.condlist{1}, 'GOOD'));
% timelock = D.fttimelock(D.indchantype('MEGMAG'), D.indsample(-0.1):D.indsample(0.3), D.indtrial(D.condlist{1}, 'GOOD'));
%
% % alternative method
% raw      = spm2fieldtrip(D);
% timelock = ft_timelockanalysis([], raw);


%% start from data that was processed by FieldTrip
% subj = 15;
% prefix = sprintf('Sub%02d', subj);
% load([prefix '_raw']);  % this is called "data" rather than "raw"
% load([prefix '_avg_Faces_vs_Scrambled']);
% load([prefix '_avg_Famous']);
% load([prefix '_avg_Unfamiliar']);
% load([prefix '_avg_Scrambled']);
%
% % load the results from part 1
% load vol
% load sens


%% deal with maxfilter

% the data has been maxfiltered and subsequently contatenated
% this results in an ill-conditioned estimate of covariance or CSD

cfg = [];
cfg.method = 'pca';
cfg.updatesens = 'no';
cfg.channel = 'MEGMAG';
comp = ft_componentanalysis(cfg, data);

cfg = [];
cfg.updatesens = 'no';
cfg.component = comp.label(51:end);
data_fix = ft_rejectcomponent(cfg, comp);


%%

cfg = [];
cfg.channel = 'MEGMAG';
cfg.method = 'wavelet';
cfg.output = 'powandcsd';
cfg.foi = 4:2:70;
cfg.toi = -0.200:0.020:1.000;
wavelet = ft_freqanalysis(cfg, data_fix);

% save wavelet wavelet

cfg = [];
cfg.layout = 'neuromag306mag.lay';
cfg.baseline = [-inf 0];
cfg.baselinetype = 'relative';
ft_multiplotTFR(cfg, wavelet)

%%

cfg = [];
cfg.grid.resolution = 7;
% cfg.inwardshift = -7; % allow dipoles 10mm outside the brain, this improves interpolation at the edges
cfg.grid.unit = 'mm';
cfg.vol  = vol;  % from FT
cfg.grad = sens; % from FT
cfg.senstype = 'meg';
cfg.normalize = 'yes';
grid = ft_prepare_leadfield(cfg, wavelet);

% save grid grid


%% perform whole-brain source reconstruction

cfg = [];
cfg.vol       = vol;  % from FT
cfg.grad      = sens; % from FT
cfg.senstype  = 'meg';
cfg.grid      = grid;
cfg.method    = 'dics';

cfg.frequency = [14 18];
cfg.latency   = [0.140 0.160];
sourceA = ft_sourceanalysis(cfg, wavelet);
cfg.latency   = [-0.100 -0.080];
sourceB = ft_sourceanalysis(cfg, wavelet);

% cfg.frequency = [40 65];
% cfg.latency   = [0.090 0.140];
% sourceA = ft_sourceanalysis(cfg, wavelet);
% cfg.latency   = [-0.100 -0.050];
% sourceB = ft_sourceanalysis(cfg, wavelet);

%
% cfg.frequency = [12 20];
% cfg.latency   = [0.090 0.140];
% sourceA = ft_sourceanalysis(cfg, wavelet);
% cfg.latency   = [-0.050 0.000];
% sourceB = ft_sourceanalysis(cfg, wavelet);

% FT_MATH requires the time axis needs to be the same
sourceA.time = 0;
sourceB.time = 0;

cfg = [];
cfg.parameter = 'pow';
cfg.operation = 'log10(x1/x2)'; % sourceA divided by sourceB
sourceR = ft_math(cfg, sourceA, sourceB);

cfg = [];
cfg.funparameter = 'pow';
ft_sourceplot(cfg, sourceR);


%% interpolate and plot on individual anatomical MRI

cfg = [];
cfg.parameter = 'pow';
sourceI = ft_sourceinterpolate(cfg, sourceR, mri_realigned);

cfg = [];
cfg.funparameter = 'pow';
ft_sourceplot(cfg, sourceI);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% beamformer_part3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% start from data that was processed by FieldTrip
%
% subj = 15;
% prefix = sprintf('Sub%02d', subj);
% load([prefix '_raw']);  % this is called "data" rather than "raw"

%% deal with maxfilter

% the data has been maxfiltered and subsequently contatenated
% this results in an ill-conditioned estimate of covariance or CSD

cfg = [];
cfg.method = 'pca';
cfg.updatesens = 'no';
cfg.channel = 'MEGMAG';
comp = ft_componentanalysis(cfg, data);

cfg = [];
cfg.updatesens = 'no';
cfg.component = comp.label(51:end);
data_fix = ft_rejectcomponent(cfg, comp);

%% compute covariance

cfg = [];
cfg.channel = 'MEGMAG';
cfg.covariance = 'yes'; % compute the covariance for single trials, then average
% cfg.preproc.bpfilter = 'yes';
% cfg.preproc.bpfreq = [5 70];
% cfg.preproc.hpfilter = 'yes';
% cfg.preproc.hpfreq = 1;
% cfg.preproc.derivative = 'yes';
cfg.preproc.demean = 'yes';             % the PCA cleanup shifted the baseline
cfg.preproc.baselinewindow = [-inf 0];  % reapply the baseline correction
cfg.keeptrials = 'yes';
timelock1 = ft_timelockanalysis(cfg, data_fix);

cfg = [];
cfg.covariance = 'yes'; % compute the covariance of the averaged ERF
timelock2 = ft_timelockanalysis(cfg, timelock1);

figure
plot(timelock2.time, timelock2.avg)

cfg = [];
cfg.layout = 'neuromag306mag.lay';
ft_multiplotER(cfg, timelock2);


%%

pos = [21 -64 30];

cfg = [];
cfg.grid.pos = pos;
cfg.grid.unit = 'mm';
% cfg.grid = grid;
cfg.vol  = vol;
cfg.grad = sens;
cfg.senstype = 'meg';
cfg.method = 'lcmv';
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.projectmom = 'yes';
source = ft_sourceanalysis(cfg, timelock2);

figure
plot(source.time, source.avg.mom{1})

%% construct single-trial virtual channel data

virtualchannel_raw = [];
virtualchannel_raw.label = {'source'};
virtualchannel_raw.trialinfo = data_fix.trialinfo;
for i=1:length(data_fix.trialinfo)
  % note that this is the non-filtered raw data
  virtualchannel_raw.time{i}       = data_fix.time{i};
  virtualchannel_raw.trial{i}(1,:) = source.avg.filter{1} * data_fix.trial{i}(:,:);
end

%% average the virtual channel ERP

cfg = [];
cfg.keeptrials = 'yes';
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = [-inf 0];
virtualchannel_avg = ft_timelockanalysis(cfg, virtualchannel_raw);
cfg.trials = virtualchannel_raw.trialinfo==1;
virtualchannel_avg1 = ft_timelockanalysis(cfg, virtualchannel_raw);
cfg.trials = virtualchannel_raw.trialinfo==2;
virtualchannel_avg2 = ft_timelockanalysis(cfg, virtualchannel_raw);
cfg.trials = virtualchannel_raw.trialinfo==3;
virtualchannel_avg3 = ft_timelockanalysis(cfg, virtualchannel_raw);

figure
plot(virtualchannel_avg.time, virtualchannel_avg.avg);

figure
plot(virtualchannel_avg.time, [virtualchannel_avg1.avg; virtualchannel_avg2.avg; virtualchannel_avg3.avg]);
legend({'1-Famous', '2-Unfamiliar', '3-Scrambled'})

figure
imagesc(squeeze(virtualchannel_avg.trial))


%% investigate the virtual channel spectrally

cfg = [];
cfg.method = 'wavelet';
cfg.output = 'pow';
cfg.foi = 4:2:70;
cfg.toi = -0.200:0.020:1.000;
virtualchannel_wavelet = ft_freqanalysis(cfg, virtualchannel_raw);

cfg = [];
cfg.baseline = [-inf 0];
cfg.baselinetype = 'relative';
cfg.interactive = 'no';
ft_singleplotTFR(cfg, virtualchannel_wavelet);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% connectivity_part1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% start with data that was preprocessed in FieldTrip
%
% subj = 15;
% prefix = sprintf('Sub%02d', subj);
% load([prefix '_raw']);  % this is called "data" rather than "raw"
%
% % load the results from beamformer_part1
% load vol
% load sens
% load mri_realigned


%% deal with maxfilter

% the data has been maxfiltered and subsequently contatenated
% this results in an ill-conditioned estimate of covariance or CSD

cfg             = [];
cfg.method      = 'pca';
cfg.updatesens  = 'no';
cfg.channel     = 'MEGMAG';
comp = ft_componentanalysis(cfg, data);

cfg             = [];
cfg.updatesens  = 'no';
cfg.component   = comp.label(51:end);
data_fix = ft_rejectcomponent(cfg, comp);


%%

pos1 = [21 -64 30];
pos2 = [0 35 83];

cfg = [];
cfg.location = pos1;
ft_sourceplot(cfg, mri_realigned);
cfg.location = pos2;
ft_sourceplot(cfg, mri_realigned);

%%

cfg             = [];
cfg.vol         = vol;
cfg.grad        = sens;
cfg.senstype    = 'meg';
cfg.method      = 'lcmv';
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.projectmom = 'yes';
cfg.grid.unit   = 'mm';
cfg.grid.pos    = pos1;
source1 = ft_sourceanalysis(cfg, timelock2);

cfg.grid.pos = pos2;
source2 = ft_sourceanalysis(cfg, timelock2);


%% construct single-trial virtual channel representation

virtualchannel_raw = [];
virtualchannel_raw.label = {'pos1'; 'pos2'};
virtualchannel_raw.trialinfo = data_fix.trialinfo;
for i=1:length(data_fix.trialinfo)
  % note that this is the non-filtered raw data
  virtualchannel_raw.time{i}       = data_fix.time{i};
  virtualchannel_raw.trial{i}(1,:) = source1.avg.filter{1} * data_fix.trial{i}(:,:);
  virtualchannel_raw.trial{i}(2,:) = source2.avg.filter{1} * data_fix.trial{i}(:,:);
end

%%

cfg                 = [];
cfg.keeptrials      = 'yes';
cfg.preproc.demean  = 'yes';
cfg.preproc.baselinewindow = [-inf 0];
virtualchannel_avg = ft_timelockanalysis(cfg, virtualchannel_raw);

figure
plot(virtualchannel_avg.time, virtualchannel_avg.avg)
legend(virtualchannel_avg.label);


%%

cfg         = [];
cfg.method  = 'wavelet';
cfg.output  = 'powandcsd';
cfg.foi     = 4:2:70;
cfg.toi     = -0.200:0.020:1.000;
virtualchannel_wavelet = ft_freqanalysis(cfg, virtualchannel_raw);

cfg                 = [];
cfg.baselinetype    = 'relative';
cfg.baseline        = [-inf 0];
cfg.channel         = {'pos1'};
figure; ft_singleplotTFR(cfg, virtualchannel_wavelet);

cfg.channel         = {'pos2'};
figure; ft_singleplotTFR(cfg, virtualchannel_wavelet);


%%

cfg = [];
cfg.method = 'coh';
coherence = ft_connectivityanalysis(cfg, virtualchannel_wavelet);

cfg = [];
cfg.baseline = [-inf 0];
cfg.baselinetype = 'relative';
cfg.channel = {'pos1'};
figure; ft_singleplotTFR(cfg, virtualchannel_wavelet);
cfg.channel = {'pos2'};
figure; ft_singleplotTFR(cfg, virtualchannel_wavelet);

figure
imagesc(coherence.time, coherence.freq, squeeze(coherence.cohspctrm(1,:,:)));
axis xy


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% connectivity_part1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% start from data that was processed by FieldTrip
%
% subj = 15;
% prefix = sprintf('Sub%02d', subj);
% load([prefix '_raw']);  % this is called "data" rather than "raw"
%
% % load the results from beamformer_part1
% load vol
% load sens
% load mri_realigned


%% deal with maxfilter

% the data has been maxfiltered and subsequently contatenated
% this results in an ill-conditioned estimate of covariance or CSD

cfg = [];
cfg.method = 'pca';
cfg.updatesens = 'no';
cfg.channel = 'MEGMAG';
comp = ft_componentanalysis(cfg, data);

cfg = [];
cfg.updatesens = 'no';
cfg.component = comp.label(51:end);
data_fix = ft_rejectcomponent(cfg, comp);

%%

cfg = [];
cfg.channel = 'MEGMAG';
cfg.method = 'wavelet';
cfg.output = 'fourier';
cfg.keeptrials = 'yes';
cfg.foi = 16;
cfg.toi = 0.150;
freq = ft_freqanalysis(cfg, data_fix);


%%

cfg = [];
cfg.grid.resolution = 7;
% cfg.inwardshift = -7; % allow dipoles 10mm outside the brain, this improves interpolation at the edges
cfg.grid.unit = 'mm';
cfg.vol       = vol;  % from FT
cfg.grad      = sens; % from FT
cfg.senstype  = 'meg';
cfg.normalize = 'yes';
grid = ft_prepare_leadfield(cfg, freq);

% save grid grid


%%

cfg           = [];
cfg.vol       = vol;  % from FT
cfg.grad      = sens; % from FT
cfg.senstype  = 'meg';
cfg.grid      = grid;
cfg.method    = 'pcc';
cfg.pcc.fixedori = 'yes';
cfg.latency   = [0.140 0.160];
cfg.frequency = [14 18];

source = ft_sourceanalysis(cfg, freq);

figure
plot(source.avg.mom{source.inside(1)}, '.')
xlabel('real');
ylabel('imag');


%%

pos = [21 -64 30];

% compute the nearest grid location
dif = grid.pos;
dif(:,1) = dif(:,1)-pos(1);
dif(:,2) = dif(:,2)-pos(2);
dif(:,3) = dif(:,3)-pos(3);
dif = sqrt(sum(dif.^2,2));
[distance, refindx] = min(dif);

cfg = [];
cfg.method    = 'coh';
% cfg.complex   = 'abs';
cfg.complex   = 'absimag';
cfg.refindx   = refindx;
conn = ft_connectivityanalysis(cfg, source);

% the output contains both the actual source position, as well as the position of the reference
% this is ugly and will probably change in future FieldTrip versions
orgpos = conn.pos(:,1:3);
refpos = conn.pos(:,4:6);
conn.pos = orgpos;

%% visualise the seed-based connectivity results

cfg               = [];
cfg.funparameter  = 'cohspctrm';
ft_sourceplot(cfg, conn);

cfg             = [];
cfg.parameter   = 'cohspctrm';
sourceI = ft_sourceinterpolate(cfg, conn, mri_realigned);

cfg               = [];
cfg.funparameter  = 'cohspctrm';
ft_sourceplot(cfg, sourceI);

%% look at connectivity difference

cfg         = [];
cfg.trials  = find(freq.trialinfo==1 | freq.trialinfo==2); % 1=Famous, 2=Unfamiliar
source1 = ft_selectdata(cfg, source);

cfg.trials  = find(freq.trialinfo==3); % 3=Scrambled
source2 = ft_selectdata(cfg, source);

cfg         = [];
cfg.method  = 'coh';
cfg.complex = 'absimag';
cfg.refindx = refindx;
conn1 = ft_connectivityanalysis(cfg, source1);
conn2 = ft_connectivityanalysis(cfg, source2);

conn1.pos = conn1.pos(:,1:3);
conn2.pos = conn2.pos(:,1:3);

cfg           = [];
cfg.parameter = 'cohspctrm';
cfg.operation = 'x1-x2';
conn_dif = ft_math(cfg, conn1, conn2);

cfg           = [];
cfg.parameter = 'cohspctrm';
source1int    = ft_sourceinterpolate(cfg, conn1, mri_realigned);
source2int    = ft_sourceinterpolate(cfg, conn2, mri_realigned);
source_difint = ft_sourceinterpolate(cfg, conn_dif, mri_realigned);

cfg               = [];
cfg.funparameter  = 'cohspctrm';
cfg.funcolorlim   = [-0.1 0.1];
cfg.maskparameter = 'cohspctrm';
cfg.opacitylim    = [-0.15 0.15];
ft_sourceplot(cfg, source_difint);

%% look at the analysis history

% for saving to disk
prefix = sprintf('/tmp/Sub%02d', subj);

cfg           = [];
cfg.filename  = [prefix '_source_difint.html'];
ft_analysispipeline(cfg, source_difint);

