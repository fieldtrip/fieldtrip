function test_ft_electroderealign

% MEM 2500mb
% WALLTIME 00:10:00

% TEST ft_electroderealign
% TEST ft_read_mri ft_read_sens ft_prepare_mesh ft_warp_apply

%% load mri, segmentation and electrode definition
mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/Subject01.mri'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/headmodel_eeg/segmentedmri'));
elec = ft_read_sens(dccnpath('/home/common/matlab/fieldtrip/template/electrode/standard_1020.elc'));
temp = ft_read_sens(dccnpath('/home/common/matlab/fieldtrip/template/electrode/standard_1005.elc'));

% create a bem and a fem mesh

cfg = [];
cfg.tissue = {'brain', 'skull', 'scalp'};
cfg.numvertices = [3000 2000 1000];
bem = ft_prepare_mesh(cfg, segmentedmri);

cfg = [];
cfg.method = 'hexahedral';
cfg.tissue = {'brain', 'skull', 'scalp'};
fem = ft_prepare_mesh(cfg, segmentedmri);

%% method: fiducial

nas = mri.hdr.fiducial.mri.nas;
lpa = mri.hdr.fiducial.mri.lpa;
rpa = mri.hdr.fiducial.mri.rpa;

% they are in voxels, hence need to be transformed to head coordinates
nas = ft_warp_apply(mri.transform, nas, 'homogenous');
lpa = ft_warp_apply(mri.transform, lpa, 'homogenous');
rpa = ft_warp_apply(mri.transform, rpa, 'homogenous');

fiducials.chanpos  = [nas; lpa; rpa];
fiducials.label    = {'Nz', 'LPA', 'RPA'};
fiducials.unit     = mri.unit;
fiducials.coordsys = mri.coordsys;

cfg = [];
cfg.method = 'fiducial';
cfg.template = fiducials;
cfg.elec = elec;
cfg.channel = 'all';
cfg.fiducial = {'Nz', 'LPA', 'RPA'};
elec_realigned1 = ft_electroderealign(cfg);

figure
ft_plot_sens(elec_realigned1, 'label', 'on');
ft_plot_axes(elec_realigned1, 'fontcolor', 'k');

%% method: template

cfg = [];
cfg.method = 'template';
cfg.template = elec_realigned1;
cfg.elec = elec;
elec_realigned2 = ft_electroderealign(cfg);

figure
ft_plot_sens(elec_realigned2, 'label', 'on');
ft_plot_axes(elec_realigned2, 'fontcolor', 'k');

%% method: interactive

% rotate    [0 0 -90]
% translate [35 0 40]

cfg = [];
cfg.method = 'interactive';
cfg.headshape = bem(3);
cfg.elec = elec;
elec_realigned3 = ft_electroderealign(cfg);

figure
ft_plot_sens(elec_realigned3, 'label', 'on');
ft_plot_axes(elec_realigned3, 'fontcolor', 'k');

%% method: headshape

cfg = [];
cfg.method = 'headshape';
cfg.headshape = bem(1);
cfg.elec = elec_realigned3;
elec_realigned4 = ft_electroderealign(cfg);

figure
ft_plot_sens(elec_realigned4, 'label', 'on');
ft_plot_axes(elec_realigned4, 'fontcolor', 'k');

