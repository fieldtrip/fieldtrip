function failed_tutorial_natmeg_dipolefitting

% WALLTIME 00:40:00
% MEM 10gb

% this script executes the MATLAB content from
% http://www.fieldtriptoolbox.org/tutorial/natmeg/dipolefitting
%
% it corresponds to the wiki version of 7 October 2014

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

clear all
close all

cd(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/natmeg'));

mrifile = './dicom/00000113.dcm';

mri_orig = ft_read_mri(mrifile);

dataset = 'oddball1_mc_downsampled.fif';

grad    = ft_read_sens(dataset,'senstype','meg');
elec    = ft_read_sens(dataset,'senstype','eeg');
shape   = ft_read_headshape(dataset,'unit','cm');

figure;
ft_plot_headshape(shape);
ft_plot_sens(grad, 'style', '*b');
ft_plot_sens(elec, 'style', '*g');

figure;
cfg = [];
ft_sourceplot(cfg, mri_orig);

if false
  % skip the interactive section
  cfg = [];
  cfg.method = 'interactive';
  cfg.coordsys = 'neuromag';
  [mri_realigned1] = ft_volumerealign(cfg, mri_orig);
  
  cfg = [];
  cfg.method = 'headshape';
  cfg.headshape = shape;
  [mri_realigned2] = ft_volumerealign(cfg, mri_realigned1);
else
  load('dipolefitting/mri_realigned2.mat');
end

cfg = [];
cfg.resolution = 1;
cfg.xrange = [-100 100];
cfg.yrange = [-110 140];
cfg.zrange = [-80 120];
mri_resliced = ft_volumereslice(cfg, mri_realigned2);

figure
ft_sourceplot([], mri_resliced);

% the low level plotting functions do not know how to deal with units,
% so make sure we have the MRI expressed in cm as well
mri_resliced_cm = ft_convert_units(mri_resliced, 'cm');

cfg           = [];
cfg.output    = {'brain', 'skull', 'scalp'};
mri_segmented = ft_volumesegment(cfg, mri_resliced);
% copy the anatomy into the segmented mri
mri_segmented.anatomy = mri_resliced.anatomy;

cfg = [];
cfg.funparameter = 'brain';
ft_sourceplot(cfg, mri_segmented);

cfg.funparameter = 'skull';
ft_sourceplot(cfg, mri_segmented);

cfg.funparameter = 'scalp';
ft_sourceplot(cfg, mri_segmented);

cfg = [];
cfg.method = 'projectmesh';
cfg.tissue = 'brain';
cfg.numvertices = 3000;
mesh_brain = ft_prepare_mesh(cfg, mri_segmented);

cfg = [];
cfg.method = 'projectmesh';
cfg.tissue = 'skull';
cfg.numvertices = 2000;
mesh_skull = ft_prepare_mesh(cfg, mri_segmented);

cfg = [];
cfg.method = 'projectmesh';
cfg.tissue = 'scalp';
cfg.numvertices = 1000;
mesh_scalp = ft_prepare_mesh(cfg, mri_segmented);

cfg = [];
cfg.method = 'isosurface';
cfg.tissue = 'scalp';
cfg.numvertices = 16000;
highres_scalp = ft_prepare_mesh(cfg, mri_segmented);

figure
ft_plot_mesh(mesh_scalp, 'edgecolor', 'none', 'facecolor', 'skin')
material dull
camlight
lighting phong

figure
ft_plot_mesh(highres_scalp, 'edgecolor', 'none', 'facecolor', 'skin')
material dull
camlight
lighting phong

cfg = [];
cfg.method = 'singleshell';
headmodel_meg = ft_prepare_headmodel(cfg, mesh_brain);

headmodel_meg = ft_convert_units(headmodel_meg,'cm');

figure;
hold on
ft_plot_headshape(shape);
ft_plot_sens(grad, 'style', 'ob');
ft_plot_sens(elec, 'style', 'og');
ft_plot_vol(headmodel_meg, 'facealpha', 0.5, 'edgecolor', 'none'); % "lighting phong" does not work with opacity
material dull
camlight

cfg = [];
cfg.dataset = dataset;
cfg.trialdef.prestim        = 0.2;
cfg.trialdef.poststim       = 0.4;
cfg.trialdef.rsp_triggers   = [256 4096];
cfg.trialdef.stim_triggers  = [1 2];
cfg.trialfun                = 'trialfun_oddball_stimlocked';
cfg = ft_definetrial(cfg);

cfg.continuous    = 'yes';
cfg.hpfilter      = 'no';
cfg.detrend       = 'no';
cfg.demean        = 'yes';
cfg.baselinewindow = [-inf 0];
cfg.dftfilter     = 'yes';
cfg.dftfreq       = [50 100];
cfg.lpfilter      = 'yes';
cfg.lpfreq        = 120;
cfg.channel       = 'MEG';
cfg.precision     = 'single';
data_meg = ft_preprocessing(cfg);

if false
  % skip the interactive section
  cfg = [];
  cfg.method = 'summary';
  cfg.channel = 'MEG*1';
  cfg.keepchannel = 'yes';
  data_meg_clean1 = ft_rejectvisual(cfg, data_meg);

  cfg.channel = {'MEG*2', 'MEG*3'};
  data_meg_clean2 = ft_rejectvisual(cfg, data_meg_clean1);
else
  % simply copy the data over
  data_meg_clean2 = data_meg;
end

cfg = [];
timelock_all = ft_timelockanalysis(cfg, data_meg_clean2);
cfg.trials = find(data_meg_clean2.trialinfo==1);
timelock_std = ft_timelockanalysis(cfg, data_meg_clean2);
cfg.trials = find(data_meg_clean2.trialinfo==2);
timelock_dev = ft_timelockanalysis(cfg, data_meg_clean2);

% save some memory
clear data_meg_clean2

cfg = [];
cfg.layout = 'neuromag306all.lay';
cfg.layout = 'neuromag306planar.lay';
cfg.layout = 'neuromag306mag.lay';
% cfg.channel = 'MEG*1';
% cfg.channel = {'MEG*2', 'MEG*3'};
ft_multiplotER(cfg, timelock_std, timelock_dev);

cfg = [];
cfg.parameter = 'avg';
cfg.operation = 'x1 - x2';
timelock_dif = ft_math(cfg, timelock_dev, timelock_std);

cfg = [];
cfg.layout = 'neuromag306all.lay';
cfg.layout = 'neuromag306planar.lay';
cfg.layout = 'neuromag306mag.lay';
% cfg.channel = 'MEG*1';
% cfg.channel = {'MEG*2', 'MEG*3'};
ft_multiplotER(cfg, timelock_dif);

cfg = [];
cfg.latency = [0.080 0.110];
cfg.numdipoles = 2;
cfg.symmetry = 'x';
cfg.grid.resolution = 1;
cfg.grid.unit = 'cm';
cfg.gridsearch = 'yes';
cfg.vol = headmodel_meg;
cfg.senstype = 'meg';
cfg.channel = {'MEG*2', 'MEG*3'};
source_planar = ft_dipolefitting(cfg, timelock_all);

cfg.channel = 'MEG*1';
source_mag = ft_dipolefitting(cfg, timelock_all);

cfg = [];
cfg.location = source_planar.dip.pos(1,:);
ft_sourceplot(cfg, mri_resliced_cm);

figure
hold on
ft_plot_dipole(source_mag.dip.pos(1,:), mean(source_mag.dip.mom(1:3,:),2), 'color', 'r')
ft_plot_dipole(source_mag.dip.pos(2,:), mean(source_mag.dip.mom(4:6,:),2), 'color', 'r')
ft_plot_dipole(source_planar.dip.pos(1,:), mean(source_planar.dip.mom(1:3,:),2), 'color', 'g')
ft_plot_dipole(source_planar.dip.pos(2,:), mean(source_planar.dip.mom(4:6,:),2), 'color', 'g')
pos = mean(source_mag.dip.pos,1);
ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1)
ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1)
ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1)
ft_plot_crosshair(pos, 'color', [1 1 1]/2);
axis tight
axis off
view(12, -10)

cfg = [];
cfg.latency = [0.080 0.110];
cfg.numdipoles = 2;
cfg.symmetry = [];
cfg.gridsearch = 'no';
cfg.dip.pos = source_planar.dip.pos;
cfg.vol = headmodel_meg;
cfg.channel = {'MEG*2', 'MEG*3'};
cfg.senstype = 'meg';
source_planar_nosym = ft_dipolefitting(cfg, timelock_all);

figure;
hold on
ft_plot_dipole(source_planar.dip.pos(1,:), mean(source_planar.dip.mom(1:3,:),2), 'color', 'g')
ft_plot_dipole(source_planar.dip.pos(2,:), mean(source_planar.dip.mom(4:6,:),2), 'color', 'g')
ft_plot_dipole(source_planar_nosym.dip.pos(1,:), mean(source_planar_nosym.dip.mom(1:3,:),2), 'color', 'm')
ft_plot_dipole(source_planar_nosym.dip.pos(2,:), mean(source_planar_nosym.dip.mom(4:6,:),2), 'color', 'm')
pos = mean(source_planar.dip.pos,1);
ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1)
ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1)
ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1)
ft_plot_crosshair(pos, 'color', [1 1 1]/2);
axis tight
axis off
view(12, -10)

cfg = [];
cfg.latency = 'all';
cfg.numdipoles = 2;
cfg.symmetry = [];
cfg.nonlinear = 'no';  % use a fixed position
cfg.gridsearch = 'no';
cfg.dip.pos = source_planar.dip.pos;
cfg.vol = headmodel_meg;
cfg.channel = {'MEG*2', 'MEG*3'};
cfg.senstype = 'meg';
source_all = ft_dipolefitting(cfg, timelock_all); % estimate the amplitude and orientation
source_std = ft_dipolefitting(cfg, timelock_std); % estimate the amplitude and orientation
source_dev = ft_dipolefitting(cfg, timelock_dev); % estimate the amplitude and orientation
source_dif = ft_dipolefitting(cfg, timelock_dif); % estimate the amplitude and orientation

figure
subplot(3,1,1); title('left: std & dev')
hold on
plot(source_std.time, source_std.dip.mom(1:3,:), '-')
legend({'x', 'y', 'z'});
plot(source_dev.time, source_dev.dip.mom(1:3,:), '.-')
axis([-0.1 0.4 -40e-3 40e-3])
grid on
subplot(3,1,2); title('right: std & dev')
hold on
plot(source_std.time, source_std.dip.mom(4:6,:), '-')
legend({'x', 'y', 'z'});
plot(source_dev.time, source_dev.dip.mom(4:6,:), '.-')
axis([-0.1 0.4 -40e-3 40e-3])
grid on
subplot(3,1,3); title('dif = dev - std')
hold on
plot(source_dif.time, source_dif.dip.mom(1:3,:), '-');
legend({'x', 'y', 'z'});
plot(source_dif.time, source_dif.dip.mom(4:6,:), '-');
axis([-0.1 0.4 -40e-3 40e-3])
grid on

cfg = [];
cfg.numdipoles = 2;
cfg.symmetry = 'x';
cfg.gridsearch = 'no';
cfg.dip.pos = source_planar.dip.pos;
cfg.vol = headmodel_meg;
cfg.channel = {'MEG*2', 'MEG*3'};
cfg.senstype = 'meg';
cfg.latency = [0.080 0.100];
source_all = ft_dipolefitting(cfg, timelock_all);
source_std = ft_dipolefitting(cfg, timelock_std);
source_dev = ft_dipolefitting(cfg, timelock_dev);

cfg.latency = [0.150 0.180];
source_dif = ft_dipolefitting(cfg, timelock_dif);

figure
hold on
ft_plot_dipole(source_all.dip.pos(1,:), mean(source_all.dip.mom(1:3,:),2), 'color', 'r')
ft_plot_dipole(source_all.dip.pos(2,:), mean(source_all.dip.mom(4:6,:),2), 'color', 'r')
ft_plot_dipole(source_std.dip.pos(1,:), mean(source_std.dip.mom(1:3,:),2), 'color', 'g')
ft_plot_dipole(source_std.dip.pos(2,:), mean(source_std.dip.mom(4:6,:),2), 'color', 'g')
ft_plot_dipole(source_dev.dip.pos(1,:), mean(source_dev.dip.mom(1:3,:),2), 'color', 'b')
ft_plot_dipole(source_dev.dip.pos(2,:), mean(source_dev.dip.mom(4:6,:),2), 'color', 'b')
ft_plot_dipole(source_dif.dip.pos(1,:), mean(source_dif.dip.mom(1:3,:),2), 'color', 'y')
ft_plot_dipole(source_dif.dip.pos(2,:), mean(source_dif.dip.mom(4:6,:),2), 'color', 'y')
pos = mean(source_std.dip.pos,1);
ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1)
ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1)
ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1)
ft_plot_crosshair(pos, 'color', [1 1 1]/2);
axis tight
axis off

cfg = [];
cfg.model = 'moving'; % default is rotating
cfg.latency = [0.070 0.140];
cfg.numdipoles = 2;
cfg.gridsearch = 'no';
cfg.dip.pos = source_planar.dip.pos;
cfg.vol = headmodel_meg;
cfg.channel = {'MEG*2', 'MEG*3'};
cfg.senstype = 'meg';
source = ft_dipolefitting(cfg, timelock_std);
% copy the time-varying position of the two dipoles into a single matrix for convenience.
for i=1:numel(source.dip)
  pos1(i,:) = source.dip(i).pos(1,:);
  pos2(i,:) = source.dip(i).pos(2,:);
end

if false
  % the following causes MATLAB to crash (hard!) on mac011 with version 2012b and 2014a
  figure
  hold on
  plot3(pos1(:,1), pos1(:,2), pos1(:,3), 'r.')
  plot3(pos2(:,1), pos2(:,2), pos2(:,3), 'g.')
  pos = (mean(pos1, 1) + mean(pos2, 1))/2;
  ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1)
  ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1)
  ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1)
  ft_plot_crosshair(pos, 'color', [1 1 1]/2);
  axis tight
  axis off
end

figure
ft_plot_mesh(mesh_brain, 'edgecolor', 'none', 'facecolor', 'r')
ft_plot_mesh(mesh_skull, 'edgecolor', 'none', 'facecolor', 'g')
ft_plot_mesh(mesh_scalp, 'edgecolor', 'none', 'facecolor', 'b')
alpha 0.3
view(132, 14)

mesh_scalp_infl = mesh_scalp;
mesh_scalp_infl.pnt = 1.10 * mesh_scalp_infl.pnt;
figure
ft_plot_mesh(mesh_brain, 'edgecolor', 'none', 'facecolor', 'r')
ft_plot_mesh(mesh_skull, 'edgecolor', 'none', 'facecolor', 'g')
ft_plot_mesh(mesh_scalp_infl, 'edgecolor', 'none', 'facecolor', 'b')
alpha 0.3
view(132, 14)

binary_brain = mri_segmented.brain;
binary_skull = mri_segmented.skull | binary_brain;
binary_scalp = mri_segmented.scalp | binary_brain | binary_skull;

close all
% using ft_sourceplot I identified the crossection with voxel
% indices [107 100 100] where the problem is visible and I will
% plot that intersection multiple times
figure(1)
tmp = binary_scalp + binary_skull + binary_brain;
imagesc(squeeze(tmp(:,:,100)));

% use IMDILATE to inflate the segmentation
binary_scalp = imdilate(binary_scalp, strel_bol(1));
figure(2)
tmp = binary_scalp + binary_skull + binary_brain;
imagesc(squeeze(tmp(:,:,100)));

% use IMDILATE to inflate the segmentation a bit more
binary_scalp = imdilate(binary_scalp, strel_bol(1));
figure(3)
tmp = binary_scalp + binary_skull + binary_brain;
imagesc(squeeze(tmp(:,:,100)));

% revert to the oriiginal binary_scalp
binary_scalp = mri_segmented.scalp + binary_skull;

% use boolean logic together with IMERODE
binary_skull = binary_skull & imerode(binary_scalp, strel_bol(2)); % fully contained inside eroded scalp
binary_brain = binary_brain & imerode(binary_skull, strel_bol(2)); % fully contained inside eroded skull
figure(4)
tmp = binary_scalp + binary_skull + binary_brain;
imagesc(squeeze(tmp(:,:,100)));

mri_segmented2 = mri_segmented;
% insert the updated binary volumes, taking out the center part for skull and scalp
mri_segmented2.brain    = binary_brain;
mri_segmented2.skull    = binary_skull & ~binary_brain;
mri_segmented2.scalp    = binary_scalp & ~binary_brain & ~binary_skull;
mri_segmented2.combined = binary_scalp + binary_skull + binary_brain; % only for plotting

cfg = [];
cfg.funparameter = 'combined';
cfg.funcolormap = 'jet';
ft_sourceplot(cfg, mri_segmented2);

% this has to be removed, otherwise ft_prepare_mesh gets confused
mri_segmented2 = rmfield(mri_segmented2, 'combined');

cfg = [];
cfg.method = 'projectmesh';
cfg.tissue = 'brain';
cfg.numvertices = 3000;
mesh_eeg(1) = ft_prepare_mesh(cfg, mri_segmented2);

cfg.tissue = 'skull';
cfg.numvertices = 2000;
mesh_eeg(2) = ft_prepare_mesh(cfg, mri_segmented2);

cfg.tissue = 'scalp';
cfg.numvertices = 1000;
mesh_eeg(3) = ft_prepare_mesh(cfg, mri_segmented2);

figure
ft_plot_mesh(mesh_eeg(1), 'edgecolor', 'none', 'facecolor', 'r')
ft_plot_mesh(mesh_eeg(2), 'edgecolor', 'none', 'facecolor', 'g')
ft_plot_mesh(mesh_eeg(3), 'edgecolor', 'none', 'facecolor', 'b')
alpha 0.3

cfg = [];
cfg.method = 'bemcp';
cfg.conductivity = [1 1/20 1]; % brain, skull, scalp
headmodel_eeg = ft_prepare_headmodel(cfg, mesh_eeg);

cfg = [];
cfg.dataset = dataset;
cfg.trialdef.prestim        = 0.2;
cfg.trialdef.poststim       = 0.4;
cfg.trialdef.rsp_triggers   = [256 4096];
cfg.trialdef.stim_triggers  = [1 2];
cfg.trialfun                = 'trialfun_oddball_stimlocked';

cfg = ft_definetrial(cfg);

cfg.continuous    = 'yes';
cfg.hpfilter      = 'no';
cfg.detrend       = 'no';
cfg.demean        = 'yes';
cfg.baselinewindow = [-inf 0];
cfg.dftfilter     = 'yes';
cfg.dftfreq       = [50 100];
cfg.lpfilter      = 'yes';
cfg.lpfreq        = 120;
cfg.channel       = 'EEG';
cfg.precision     = 'single';

data_eeg = ft_preprocessing(cfg);

if false
  % skip the interactive section
  cfg = [];
  cfg.method = 'summary';
  cfg.keepchannel = 'no';
  cfg.preproc.reref = 'yes';
  cfg.preproc.refchannel = 'all';
  data_eeg_clean = ft_rejectvisual(cfg, data_eeg);
else
  % just copy it over
  data_eeg_clean = data_eeg;
end

cfg = [];
cfg.reref = 'yes';
cfg.refchannel = 'all';
data_eeg_reref = ft_preprocessing(cfg, data_eeg_clean);

cfg = [];
timelock_eeg_all = ft_timelockanalysis(cfg, data_eeg_reref);

cfg.trials = find(data_eeg_reref.trialinfo==1);
timelock_eeg_std = ft_timelockanalysis(cfg, data_eeg_reref);

cfg.trials = find(data_eeg_reref.trialinfo==2);
timelock_eeg_dev = ft_timelockanalysis(cfg, data_eeg_reref);

% save some memory
clear data_eeg_reref

cfg = [];
cfg.layout = 'neuromag306eeg1005_natmeg.lay';
ft_multiplotER(cfg, timelock_eeg_std, timelock_eeg_dev);

cfg = [];
cfg.parameter = 'avg';
cfg.operation = 'x1 - x2';
timelock_eeg_dif = ft_math(cfg, timelock_eeg_dev, timelock_eeg_std);

cfg = [];
cfg.layout = 'neuromag306eeg1005_natmeg.lay';
ft_multiplotER(cfg, timelock_eeg_dif);

cfg = [];
cfg.latency = [0.080 0.110];
cfg.numdipoles = 2;
cfg.symmetry = 'x';
cfg.grid.resolution = 1;
cfg.grid.unit = 'cm';
cfg.gridsearch = 'yes';
cfg.vol = headmodel_eeg;
cfg.senstype = 'eeg';
cfg.channel = 'all';
source_eeg = ft_dipolefitting(cfg, timelock_eeg_all);

cfg = [];
cfg.location = source_eeg.dip.pos(1,:);
ft_sourceplot(cfg, mri_resliced_cm);
figure

ft_plot_dipole(source_eeg.dip.pos(1,:), mean(source_eeg.dip.mom(1:3,:),2), 'color', 'b')
ft_plot_dipole(source_eeg.dip.pos(2,:), mean(source_eeg.dip.mom(4:6,:),2), 'color', 'b')

ft_plot_dipole(source_mag.dip.pos(1,:), mean(source_mag.dip.mom(1:3,:),2), 'color', 'r')
ft_plot_dipole(source_mag.dip.pos(2,:), mean(source_mag.dip.mom(4:6,:),2), 'color', 'r')

ft_plot_dipole(source_planar.dip.pos(1,:), mean(source_planar.dip.mom(1:3,:),2), 'color', 'g')
ft_plot_dipole(source_planar.dip.pos(2,:), mean(source_planar.dip.mom(4:6,:),2), 'color', 'g')

pos = mean(source_eeg.dip.pos,1);
% pos = source_eeg.dip.pos(1,:); % use another crossection for the MRI

ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1)
ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1)
ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1)

ft_plot_crosshair(pos, 'color', [1 1 1]/2);

axis tight
axis off

cfg = [];
cfg.latency = [0.080 0.110];
cfg.numdipoles = 2;
cfg.dip.pos = source_eeg.dip.pos;
cfg.gridsearch = 'no';
cfg.nonlinear = 'yes';
cfg.vol = headmodel_eeg;
cfg.senstype = 'eeg';
cfg.channel = 'all';
source_eeg2 = ft_dipolefitting(cfg, timelock_eeg_all);

figure
ft_plot_dipole(source_eeg.dip.pos(1,:), mean(source_eeg.dip.mom(1:3,:),2), 'color', 'b')
ft_plot_dipole(source_eeg.dip.pos(2,:), mean(source_eeg.dip.mom(4:6,:),2), 'color', 'b')

ft_plot_dipole(source_eeg2.dip.pos(1,:), mean(source_eeg2.dip.mom(1:3,:),2), 'color', 'm')
ft_plot_dipole(source_eeg2.dip.pos(2,:), mean(source_eeg2.dip.mom(4:6,:),2), 'color', 'm')

pos = mean(source_eeg.dip.pos,1);
% pos = source_eeg.dip.pos(1,:); % alternative crossection for the MRI

ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1)
ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1)
ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1)

ft_plot_crosshair(pos, 'color', [1 1 1]/2);

axis tight
axis off
