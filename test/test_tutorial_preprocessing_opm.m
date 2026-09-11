function test_tutorial_preprocessing_opm(datadir)

% MEM 2gb
% WALLTIME 00:20:00
% DEPENDENCY ft_preprocessing ft_timelockanalysis
% DATA public

% see http://www.fieldtriptoolbox.org/tutorial/preproc/denoising_opm
% this test function corresponds to the version on the wiki at 10 September 2026

if nargin<1
  datadir = dccnpath('/project/3031000.02/external/download/tutorial/preprocessing_opm');
end

pwdir = pwd;
cd(datadir);

%
cfg = [];
cfg.dataset = 'MedianNerve_StimBreakStim2min_Pos1.fif';
cfg.ploteventlabels = 'no';
cfg.preproc.demean = 'yes';
cfg.ylim = [-1 1]*1e-11;
ft_databrowser(cfg)

%
hdr = ft_read_header('MedianNerve_StimBreakStim2min_Pos1.fif');

%
event = ft_read_event('MedianNerve_StimBreakStim2min_Pos1.fif');

%
cfg = [];
cfg.trialdef.eventtype = 'Input-1';
cfg.trialdef.prestim = 0.1;
cfg.trialdef.poststim = 0.3;

cfg.dataset = 'MedianNerve_StimBreakStim2min_Pos1.fif';
cfg = ft_definetrial(cfg);

%
cfg.channel = '00*';
cfg.detrend = 'yes';
cfg.baselinewindow = [-inf 0];
data_pos1 = ft_preprocessing(cfg);

%
cfg = [];
cfg.hpfilter = 'yes';
cfg.hpfreq = 0.5;
cfg.hpfilttype = 'firws';
cfg.usefftfilt ='yes';
cfg.channel = '00*';
cfg.dataset = 'MedianNerve_StimBreakStim2min_Pos1.fif';
data_pos1_continuous = ft_preprocessing(cfg);

cfg = [];
cfg.dataset = 'MedianNerve_StimBreakStim2min_Pos1.fif';
cfg.trialdef.eventtype = 'Input-1';
cfg.trialdef.prestim = 0.1;
cfg.trialdef.poststim = 0.3;
cfg = ft_definetrial(cfg);

% the trl field defines the trials
trl = cfg.trl;

% cut the trials from the continuous data in memory
cfg = [];
cfg.trl = trl;
data_pos1_segmented = ft_redefinetrial(cfg, data_pos1_continuous);

%
cfg = [];
cfg.trialdef.eventtype = 'Input-1';
cfg.trialdef.prestim = 0.1;
cfg.trialdef.poststim = 0.3;

cfg.dataset = 'MedianNerve_StimBreakStim2min_Pos2.fif';
cfg = ft_definetrial(cfg);

cfg.channel = '00*';
cfg.detrend = 'yes';
cfg.baselinewindow = [-inf 0];
data_pos2 = ft_preprocessing(cfg);

cfg = [];
cfg.trialdef.eventtype = 'Input-1';
cfg.trialdef.prestim = 0.1;
cfg.trialdef.poststim = 0.3;

cfg.dataset = 'MedianNerve_StimBreakStim2min_Pos3.fif';
cfg = ft_definetrial(cfg);

cfg.channel = '00*';
cfg.detrend = 'yes';
cfg.baselinewindow = [-inf 0];
data_pos3 = ft_preprocessing(cfg);

%
chassis = {1   2   3   4   5  6   7  8};
sensor  = {338 119 323 111 62 336 22 246};
pos1 = {'FL30'   'FL21' 'FL20' 'FL23' 'FL36' 'FL35' 'FL34' 'FL84'};
pos2 = {'FL30_2' 'FL38' 'FL37' 'FL28' 'FL27' 'FL19' 'FL18' 'FL84_2'};
pos3 = {'FL30_3' 'FL22' 'FL29' 'FL39' 'FL45' 'FL44' 'FL43' 'FL84_3'};

montage_pos1 = [];
montage_pos1.labelold = {
  '00:01-BZ_OL'
  '00:02-BZ_OL'
  '00:03-BZ_OL'
  '00:04-BZ_OL'
  '00:05-BZ_OL'
  '00:06-BZ_OL'
  '00:07-BZ_OL'
  '00:08-BZ_OL'
  }';
montage_pos1.labelnew = pos1; % labels according to the helmet
montage_pos1.tra = eye(8);

%
montage_pos2 = [];
montage_pos2.labelold = {
  '00:01-BZ_OL'
  '00:02-BZ_OL'
  '00:03-BZ_OL'
  '00:04-BZ_OL'
  '00:05-BZ_OL'
  '00:06-BZ_OL'
  '00:07-BZ_OL'
  '00:08-BZ_OL'
  }';
montage_pos2.labelnew = pos2; % labels according to the helmet
montage_pos2.tra = eye(8);

montage_pos3 = [];
montage_pos3.labelold = {
  '00:01-BZ_OL'
  '00:02-BZ_OL'
  '00:03-BZ_OL'
  '00:04-BZ_OL'
  '00:05-BZ_OL'
  '00:06-BZ_OL'
  '00:07-BZ_OL'
  '00:08-BZ_OL'
  }';
montage_pos3.labelnew = pos3; % labels according to the helmet
montage_pos3.tra = eye(8);

%
cfg = [];
cfg.montage = montage_pos1;
data_pos1 = ft_preprocessing(cfg, data_pos1);

cfg.montage = montage_pos2;
data_pos2 = ft_preprocessing(cfg, data_pos2);

cfg.montage = montage_pos3;
data_pos3 = ft_preprocessing(cfg, data_pos3);

% %
% cfg = [];
% cfg.method = 'summary';
% data_pos1_clean = ft_rejectvisual(cfg, data_pos1);
% 
% data_pos2_clean = ft_rejectvisual(cfg, data_pos2);
% data_pos3_clean = ft_rejectvisual(cfg, data_pos3);

%
cfg = [];
cfg.metric = 'std';
cfg.threshold = 1.4e-12;
[cfg, artifact] = ft_badsegment(cfg, data_pos1);

data_pos1_clean = ft_rejectartifact(cfg, data_pos1);

cfg = [];
cfg.metric = 'std';
cfg.threshold = 1.4e-12;
[cfg, artifact] = ft_badsegment(cfg, data_pos2);

data_pos2_clean = ft_rejectartifact(cfg, data_pos2);

cfg = [];
cfg.metric = 'std';
cfg.threshold = 1.4e-12;
[cfg, artifact] = ft_badsegment(cfg, data_pos3);

data_pos3_clean = ft_rejectartifact(cfg, data_pos3);

%
cfg.artifact.badsegment.artifact = artifact;

%
cfg = [];
timelock_pos1 = ft_timelockanalysis(cfg, data_pos1_clean);
timelock_pos2 = ft_timelockanalysis(cfg, data_pos2_clean);
timelock_pos3 = ft_timelockanalysis(cfg, data_pos3_clean);

%
cfg = [];
cfg.appenddim = 'chan';
timelock = ft_appendtimelock(cfg, timelock_pos1, timelock_pos2, timelock_pos3);

% FIXME I am not sure why this is needed
timelock.avg = timelock.trial;
timelock = rmfield(timelock, 'trial');

%
cfg = [];
cfg.layout = 'fieldlinealpha1_helmet.mat';
cfg.showlabels = 'yes';
cfg.ylim = 'maxabs';
ft_multiplotER(cfg, timelock);

% We can make the same topographic figure by calling ft_topoplotER.
cfg = [];
cfg.layout = 'fieldlinealpha1_helmet.mat';
cfg.marker = 'labels';
cfg.xlim = [0.019 0.028];
cfg.zlim = 'maxabs';
ft_topoplotER(cfg, timelock);

%
load fieldlinealpha1_helmet.mat  % this contains the layout
missing = setdiff(layout.label, timelock.label);
missing = missing(startsWith(missing, 'FL'));

% make a copy of the original timelock structure, remove two fields that we don't care about
timelock_full = rmfield(timelock, {'dof', 'var'});

% add the extra slots from the helmet that are missing in the recording
timelock_full.label = cat(1, timelock_full.label(:), missing(:));
timelock_full.avg   = cat(1, timelock_full.avg, nan(numel(missing), numel(timelock_full.time)));

%
cfg = [];
cfg.layout = 'fieldlinealpha1_helmet.mat';
cfg.showlabels = 'yes';
cfg.ylim = 'maxabs';
ft_multiplotER(cfg, timelock_full);

%
cfg = [];
cfg.layout = 'fieldlinealpha1_helmet.mat';
cfg.marker = 'labels';
cfg.xlim = [0.016 0.024];
cfg.zlim = 'maxabs';
ft_topoplotER(cfg, timelock_full);

%
cfg = [];
cfg.layout = 'fieldlinealpha1_helmet.mat';
cfg.marker = 'labels';
cfg.xlim = [0.016 0.024];
cfg.zlim = 'maxabs';
cfg.interpolatenan = 'no';
ft_topoplotER(cfg, timelock_full);

%
cfg = [];
cfg.layout = 'fieldlinealpha1_helmet.mat';
cfg.skipcomnt = 'no';
cfg.skipscale = 'no';
layout = ft_prepare_layout(cfg);

load fieldlinealpha1_helmet.mat  % this contains the layout

%
layout_full = layout;
clear layout
 
 %
cfg = [];
cfg.layout = rmfield(layout_full, 'mask');
cfg.channel = timelock.label;
cfg.mask = 'convex';
layout_trimmed = ft_prepare_layout(cfg);

%
cfg = [];
cfg.layout = layout_trimmed;
cfg.marker = 'labels';
cfg.xlim = [0.016 0.024];
cfg.zlim = 'maxabs';
ft_topoplotER(cfg, timelock);

%
cfg = [];
cfg.layout = rmfield(layout_full, 'mask');
cfg.channel = setdiff(timelock.label, 'FL84');
cfg.mask = 'convex';
layout_trimmed = ft_prepare_layout(cfg);

cfg = [];
cfg.layout = layout_trimmed;
cfg.marker = 'labels';
cfg.xlim = [0.016 0.024];
cfg.zlim = 'maxabs';
ft_topoplotER(cfg, timelock);

%
close all
load fieldlinealpha1 % this contains the fieldlinealpha1 structure, similar to the grad structure in the data

figure
ft_plot_sens(fieldlinealpha1, 'label', 'yes', 'axes', 1, 'orientation', 1, 'chantype', 'megmag', 'fiducial', 1)

%
fieldlinealpha1.chantype( ismember(fieldlinealpha1.label, missing)) = {'missing'};
fieldlinealpha1.chantype(~ismember(fieldlinealpha1.label, missing)) = {'megmag'};

figure
ft_plot_sens(fieldlinealpha1, 'label', 'yes', 'axes', 1, 'orientation', 1, 'chantype', 'megmag', 'fiducial', 0)

cd(pwdir);
