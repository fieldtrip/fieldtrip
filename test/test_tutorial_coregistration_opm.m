function test_tutorial_coregistration_opm(datadir)

% MEM 2gb
% WALLTIME 00:20:00
% DEPENDENCY ft_read_headshape ft_electroderealign ft_defacemesh
% DATA public

% see http://www.fieldtriptoolbox.org/tutorial/preproc/denoising_opm
% this test function corresponds to the version on the wiki at 10 September
% 2026, NOTE that it has been stripped down to the parts that can run
% non-interactively, so that it is a part of the daily test batch. Also the
% function can be used for quick testing (e.g. just prior to a toolkit) to
% check whether the basic functionality works, provided the version of the
% tutorial online has not deviated much since the current m-file has been
% created.

if nargin<1
  datadir = dccnpath('/project/3031000.02/external/download/tutorial/coregistration_opm');
end
pwdir = pwd;
cd(datadir);

headshape = ft_read_headshape('example1_head_markers.pos');
headshape = ft_convert_units(headshape, 'mm');

%% visualization, coordinate axes are initially ALS
figure
ft_plot_headshape(headshape)
ft_plot_axes(headshape)
view([-27 20])

headshape.coordsys = 'ctf';
headshape = ft_convert_coordsys(headshape, 'neuromag');  % this rotates it such that the X-axis points to the right

%% visualization, coordinate axes are now RAS
figure
ft_plot_headshape(headshape)
ft_plot_axes(headshape)
view([114 20])

%% select the reference points on the helmets, with their corresponding label
fid_measured = [];
fid_measured.pos(1,:) = headshape.pos(end-7,:);
fid_measured.pos(2,:) = headshape.pos(end-6,:);
fid_measured.pos(3,:) = headshape.pos(end-5,:);
fid_measured.pos(4,:) = headshape.pos(end-4,:);
fid_measured.pos(5,:) = headshape.pos(end-3,:);
fid_measured.pos(6,:) = headshape.pos(end-2,:);
fid_measured.pos(7,:) = headshape.pos(end-1,:);
fid_measured.pos(8,:) = headshape.pos(end-0,:);
fid_measured.label = {'A5', 'A6', 'A7', 'A8', 'A1', 'A2', 'A3', 'A4'};

% To perform a later comparison, it is convenient to sort them from 1 to 8.
[fid_measured.label, indx] = sort(fid_measured.label);
fid_measured.pos = fid_measured.pos(indx,:);

headshape.fid = fid_measured;

%
fieldlinebeta2 = ft_read_sens('fieldlinebeta2.mat'); % from fieldtrip/template/grad
fieldlinebeta2 = ft_convert_units(fieldlinebeta2, 'mm');
fid_helmet     = fieldlinebeta2.fid;

%% show the misalignment
figure
ft_plot_headshape(headshape)
ft_plot_axes(headshape)
ft_plot_sens(fieldlinebeta2)
view([102 5]);

%
% We will proceed with **[ft_electroderealign](/reference/ft_electroderealign)**, which was originally implemented to align EEG electrode positions to a head surface. As it turns out, it can also be used more general to align two sets of points.
%
%% ## Calculation of the transformation parameters
%
% The alignment parameters can be estimated using the |template| method in **[ft_electroderealign](/reference/ft_electroderealign)**. Since we want to express the OPM sensors' coordinates in the head coordinate system, the Polhemus measured positions will be used as the target.
%
%% estimate the alignment parameters
cfg         = [];
cfg.method  = 'template';
cfg.target  = fid_measured;
cfg.elec    = fid_helmet;
fid_aligned = ft_electroderealign(cfg);

%% ## Apply the transformation to the OPM sensors
%
% The output data structure `fid_aligned` not only contains the aligned fiducials, but also the parameters that were used to align (or transform) them. We can apply the same transformation parameters to the OPM sensors.
%
fieldlinebeta2_head = ft_transform_geometry(fid_aligned.rigidbody, fieldlinebeta2, 'rigidbody');

figure
ft_plot_headshape(headshape)
ft_plot_axes(headshape)
ft_plot_sens(fieldlinebeta2_head)
view([102 5]);

% load in the data
cfg              = [];
cfg.dataset      = 'example2_magneticphantom_HPIplusdipoleset6_raw.fif';
cfg.coilaccuracy = 0;
data_all         = ft_preprocessing(cfg);

% We can visualize the data with **[ft_databrowser](/reference/ft_databrowser)** to see where the sine-wave signals start and end.
%
cfg = [];
cfg.viewmode = 'vertical';
cfg.blocksize = 300; % seconds
ft_databrowser(cfg, data_all);

% this is the time of a single sample
tsample = 1./data_all.fsample;

cfg         = [];
cfg.latency = [0 60-tsample];
cfg.channel = {'all' '-L212_bz' '-R212_bz'};
data        = ft_selectdata(cfg, data_all);

% We cut the data into 10-second segments with 80% overlap and compute the averaged power spectrum over all segments to verify the expected spectral peaks (and their harmonics) at 8, 11 and 14 Hz.
%
cfg            = [];
cfg.length     = 10;
cfg.overlap    = 0.8;
data_segmented = ft_redefinetrial(cfg, data);

cfg           = [];
cfg.method    = 'mtmfft';
cfg.foilim    = [0 40];
cfg.taper     = 'dpss';
cfg.tapsmofrq = 0.2;
cfg.pad       = 10;
freq          = ft_freqanalysis(cfg, data_segmented);

figure
plot(freq.freq, log10(mean(freq.powspctrm)));
xlabel('frequency (Hz)');
ylabel('log_10 power')

%
% To focus on the signals of the specific HPI-coils, we bandpass filter the data in the frequency bands corresponding to each of the coils, and cut off the edges for any potential filter edge artifacts.
%
cfg            = [];
cfg.bpfilter   = 'yes';
cfg.bpfilttype = 'firws';
cfg.usefftfilt = 'yes';
cfg.bpfreq     = [7 9];
data08         = ft_preprocessing(cfg, data); % nas

cfg.bpfreq     = [10 12];
data11         = ft_preprocessing(cfg, data); % rpa

cfg.bpfreq     = [13 15];
data14         = ft_preprocessing(cfg, data); % lpa

cfg            = [];
cfg.latency    = [4 56-1./data.fsample];
data08         = ft_selectdata(cfg, data08);
data11         = ft_selectdata(cfg, data11);
data14         = ft_selectdata(cfg, data14);

%% look at 2 seconds of the data
figure
plot(data08.time{1}, data08.trial{1});
xlim([4 6]);
xlabel('time (s)');
ylabel('magnetic field strength (T)');

%
%% ## Fit dipoles to the sensor topographies
%
% We proceed by performing a principal component analysis (PCA) on the filtered data. The idea is that - given that the signals from the HPI coils are the strongest signals in the measurement, and given that we have bandpass filtered the data - the strongest principal components will represent the 'spatial fingerprints' of each of the HPI coils. Those fingerprints will be used to perform a dipole fit, i.e., find the position of a dipole that optimally explain those principal components.
%
cfg            = [];
cfg.method     = 'pca';
cfg.updatesens = 'no';
comp08 = ft_componentanalysis(cfg, data08);
comp11 = ft_componentanalysis(cfg, data11);
comp14 = ft_componentanalysis(cfg, data14);

%% look at the topographies
cfg              = [];
cfg.component    = 1;
cfg.layout       = 'fieldlinebeta2bz_helmet.mat';
cfg.gridscale    = 150;
cfg.interplimits = 'sensors';
cfg.figure       = subplot('position',[0 0 1/3 1]);
ft_topoplotIC(cfg, comp08);
cfg.figure       = subplot('position',[1/3 0 1/3 1]);
ft_topoplotIC(cfg, comp11);
cfg.figure       = subplot('position',[2/3 0 1/3 1]);
ft_topoplotIC(cfg, comp14);

%
% For the fitting the magnetic dipole positions, we will use a grid search as an initial scan over the whole volume, followed by a iterative non-linear optimization. The grid search is motivated by the fact that a non-linear search of the whole parameter space (i.e., volume of space covered by the helmet) might result in convergence to a local minimum.
%
% The following creates a source model that consists of a regular grid of dipole positions that will be used for the initial grid search.
%
%% create a regular grid of dipole positions bounded by the helmet
fieldlinebeta2 = ft_read_sens('fieldlinebeta2.mat');  % from fieldtrip/template/grad

% make a fake headshape, we use this to make a fake headmodel
fake_headshape      = [];
fake_headshape.pos  = fieldlinebeta2.coilpos;
fake_headshape.unit = 'm';

% create a fake singleshell headmodel, this will act as the boundary for the grid
cfg = [];
cfg.method = 'singleshell';
cfg.headshape = fake_headshape;
cfg.numvertices = 144; % keep the same number of vertices as OPMs
fake_headmodel = ft_prepare_headmodel(cfg);

%% create the grid, grid points outside the fake head will be flagged as such
cfg            = [];
cfg.headmodel  = fake_headmodel;
cfg.resolution = 0.0075;
sourcemodel    = ft_prepare_sourcemodel(cfg);

% this is the real volume conductor model that we want to use
cfg = [];
cfg.method = 'infinite';
headmodel = ft_prepare_headmodel(cfg);

% Now we can perform the dipole fits.
%
cfg             = [];
cfg.headmodel   = headmodel;
cfg.grad        = data.grad;
cfg.component   = 1;
cfg.gridsearch  = 'yes';
cfg.sourcemodel = sourcemodel;
dip08 = ft_dipolefitting(cfg, comp08);
dip11 = ft_dipolefitting(cfg, comp11);
dip14 = ft_dipolefitting(cfg, comp14);

% for verification
disp(norm(dip11.dip.pos - dip14.dip.pos)*100) % in cm
disp(norm(dip08.dip.pos - dip11.dip.pos)*100)
disp(norm(dip08.dip.pos - dip14.dip.pos)*100)

%15.4760
%10.2201
%10.5455

% transform fiducial coordinates to head coordinates (RAS)
fid1 = dip08.dip.pos; % nas
fid2 = dip14.dip.pos; % lpa
fid3 = dip11.dip.pos; % rpa
transform_sens2head = ft_headcoordinates(fid1, fid2, fid3, 'neuromag');
%
fieldlinebeta2_head = ft_transform_geometry(transform_sens2head, fieldlinebeta2);

% We can plot the sensors, which are now in head coordinates
figure
ft_plot_sens(fieldlinebeta2_head)
ft_plot_axes(fieldlinebeta2_head)
view([130 30]);

% and if we transform the dipole positions from helmet to head coordinates, we can also add those to the figure.
ft_plot_dipole(dip08.dip.pos, dip08.dip.mom, 'length', 0.02, 'diameter', 0.01)
ft_plot_dipole(dip11.dip.pos, dip11.dip.mom, 'length', 0.02, 'diameter', 0.01)
ft_plot_dipole(dip14.dip.pos, dip14.dip.mom, 'length', 0.02, 'diameter', 0.01)

%
scan      = ft_read_headshape('example3_face_helmet.obj');
scan.unit = 'm'; % the estimated 'dm' is not correct

figure; hold on;
ft_plot_headshape(scan);
ft_plot_axes(scan);
lighting gouraud
material dull
light

% % % by clicking on 'dummy' nas/lpa/rpa
% % cfg          = [];
% % cfg.method   = 'fiducial';
% % cfg.coordsys = 'neuromag';
% % scan         = ft_meshrealign(cfg, scan); % this does not work when running non-interactively
% % 
% % figure; hold on;
% % ft_plot_headshape(scan);
% % ft_plot_axes(scan);
% % view([125 10]);
% % lighting gouraud
% % material dull
% % light
% % 
% % %
% % load example3_face_helmet_aligned.mat  % this contains the aligned scan
% % % cut off the irrelevant parts
% % cfg         = [];
% % cfg.method  = 'plane';
% % cfg.rotate  = [-40 0 0];
% % cfg.translate = [0 0 -140];
% % scan_head   = ft_defacemesh(cfg, scan); % viewpoint left,  rotate [-30 0 0], translate [0 0 -130];
% % 
% % figure; hold on;
% % ft_plot_headshape(scan_head);
% % ft_plot_axes(scan_head);
% % view([125 10]);
% % lighting gouraud
% % material dull
% % light
% % 
% % %
% % cfg        = [];
% % cfg.method = 'box';
% % cfg.selection = 'inside';
% % scan_face  = ft_defacemesh(cfg, scan_head); % rotate [-30 0 0], scale [0.15 0.20 0.20], translate [3 0 -80]
% % 
% % figure; hold on;
% % ft_plot_headshape(scan_face);
% % ft_plot_axes(scan_face);
% % view([125 10]);
% % lighting gouraud
% % material dull
% % light
% % 
% % %
% % % The surface mesh of the helmet will be extracted by removing the face from the scan. We use the same selection as in the extraction of the face, but now we keep the outside of the box rather than the inside.
% % %
% % cfg           = [];
% % cfg.method    = 'box';
% % cfg.selection = 'outside';
% % scan_helmet   = ft_defacemesh(cfg, scan_head); % rotate [-30 0 0], scale [0.15 0.20 0.20], translate [3 0 -80]
% % 
% % figure; hold on;
% % ft_plot_headshape(scan_helmet);
% % ft_plot_axes(scan_helmet);
% % view([125 10]);
% % lighting gouraud
% % material dull
% % light

%
mri = ft_read_mri('example3_anatomical.nii');

% check the coordinate system
ft_determine_coordsys(mri, 'interactive', 'no');

% % %
% % cfg           = [];
% % cfg.coordsys  = 'neuromag';
% % mri_realigned = ft_volumerealign(cfg, mri);
% % 
% % % check the coordinate system after realignment
% % ft_determine_coordsys(mri_realigned, 'interactive', 'no');

%
load example3_mri_realigned.mat % this contains the realigned mri

% segment the scalp
cfg          = [];
cfg.output   = 'scalp';
seg          = ft_volumesegment(cfg, mri_realigned);

% create a mesh for the scalp
cfg             = [];
cfg.tissue      = 'scalp';
cfg.numvertices = 10000;
mri_face        = ft_prepare_mesh(cfg, seg);
mri_face        = ft_convert_units(mri_face, 'm');

% % % the following values can be specified (without clicking the 'apply' button in between) 
% % % for the rotation: `[13.5 0 -4]`, and for the translation: `[-0.003 0.07 -0.008]`. Note that the units are now expressed in 'm'.
% % %
% % cfg             = [];
% % cfg.method      = 'interactive';
% % cfg.headshape   = mri_face;
% % cfg.meshstyle   = {'edgecolor', 'k', 'facecolor', 'skin'};
% % scan_face_aligned = ft_meshrealign(cfg, scan_face);
% % 
% % figure; hold on;
% % ft_plot_headshape(mri_face, 'facealpha', 0.4);
% % ft_plot_mesh(scan_face_aligned, 'facecolor','skin');
% % view([125 10]);
% % lighting gouraud
% % material dull
% % light

%
helmet_rim = ft_read_headshape('fieldlinebeta2_helmet_rim.mat');
helmet_rim.coordsys = 'ras';

% % % the following values can be specified (without clicking the 'apply' button in between) 
% % % for the rotation: `[20 -2 0]`, and for the translation: `[-0.003 0.065 -0.045]`. Note that the units are expressed in 'm'.
% % %
% % cfg = [];
% % cfg.method = 'interactive';
% % cfg.headshape = helmet_rim;
% % cfg.meshstyle = {'edgecolor', 'none', 'facecolor', [1 0.5 0.5]};
% % scan_helmet_aligned = ft_meshrealign(cfg, scan_helmet);
% % 
% % figure; hold on;
% % ft_plot_mesh(model_helmet_rim, 'edgecolor', 'none', 'facecolor', [0.5 0.5 1], 'facealpha', 0.4);
% % ft_plot_mesh(scan_helmet_aligned, 'edgecolor', 'none', 'facecolor', [1 0.5 0.5], 'facealpha', 0.4);
% % view([145 10]);
% % lighting gouraud
% % material dull
% % light

% % %
% % transform_scan2helmet = scan_helmet_aligned.cfg.transform;
% % transform_scan2face   = scan_face_aligned.cfg.transform;
% % transform_helmet2face = transform_scan2face/transform_scan2helmet;
% % 
% % %
% % fieldlinebeta2 = ft_read_sens('fieldlinebeta2.mat');  % from fieldtrip/template/grad
% % fieldlinebeta2.coordsys = 'ras';
% % 
% % %
% % fieldlinebeta2_head = ft_transform_geometry(transform_helmet2face, fieldlinebeta2);
% % 
% % figure; hold on;
% % ft_plot_sens(fieldlinebeta2_head);
% % ft_plot_headshape(mri_face, 'facecolor', [0.5 0.5 1], 'facealpha', 0.4, 'edgecolor', 'none');
% % view([125 10]);
% % lighting gouraud
% % material dull
% % light

cd(pwdir);
