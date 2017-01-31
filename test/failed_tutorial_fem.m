% function failed_tutorial_fem

% WALLTIME 01:00:00
% MEM 6gb

% TEST ft_volumesegment ft_prepare_mesh ft_prepare_headmodel ft_headmodel_simbio

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

% this test script is based on the tutorial under development at
% http://www.fieldtriptoolbox.org/development/simbio

cd(dccnpath('/home/common/matlab/fieldtrip/data'))

mri = ft_read_mri('Subject01.mri');

figure
cfg = [];
ft_sourceplot(cfg, mri);

cfg     = [];
cfg.dim = mri.dim;
mri     = ft_volumereslice(cfg, mri);

figure
cfg = [];
ft_sourceplot(cfg,mri);

cfg           = [];
cfg.output    = {'gray','white','csf','skull','scalp'};
segmentedmri  = ft_volumesegment(cfg, mri);


seg_i            = ft_datatype_segmentation(segmentedmri, 'segmentationstyle', 'indexed');

cfg              = [];
cfg.funparameter = 'seg';
cfg.funcolormap  = lines(6);
cfg.location     = 'center';
cfg.atlas        = seg_i;
ft_sourceplot(cfg, seg_i);

cfg        = [];
cfg.shift  = 0.3;
cfg.method = 'hexahedral';
mesh = ft_prepare_mesh(cfg, segmentedmri);

cfg              = [];
cfg.method       = 'simbio';
cfg.conductivity = [0.33 0.14 1.79 0.01 0.43];           % order follows mesh.tissuelabel
vol = ft_prepare_headmodel(cfg, mesh);
% the call to "ft_prepare_headmodel" took 320 seconds and required the additional allocation of an estimated 501 MB

% you may need to specify the full path to the file
elec = ft_read_sens('standard_1020.elc');

figure
hold on
ft_plot_mesh(mesh,'surfaceonly','yes','vertexcolor','none','edgecolor','none','facecolor',[0.5 0.5 0.5],'face alpha',0.5)
camlight
ft_plot_sens(elec,'style', 'sk');

nas = mri.hdr.fiducial.mri.nas;
lpa = mri.hdr.fiducial.mri.lpa;
rpa = mri.hdr.fiducial.mri.rpa;

vox2head = mri.transform;

nas = ft_warp_apply(vox2head, nas, 'homogenous');
lpa = ft_warp_apply(vox2head, lpa, 'homogenous');
rpa = ft_warp_apply(vox2head, rpa, 'homogenous');

% create a structure similar to a template set of electrodes
fid.chanpos       = [nas; lpa; rpa];       % CTF head coordinates of fiducials
fid.label         = {'Nz','LPA','RPA'};    % use the same labels as those in elec
fid.unit          = 'mm';                  % use the same units as those in mri

% alignment
cfg               = [];
cfg.method        = 'fiducial';
cfg.elec          = elec;                  % the electrodes we want to align
cfg.template      = fid;                   % the template we want to align to
cfg.fiducial      = {'Nz', 'LPA', 'RPA'};  % labels of fiducials in fid and in elec
elec_aligned      = ft_electroderealign(cfg);


figure;
hold on;
ft_plot_mesh(mesh,'surfaceonly','yes','vertexcolor','none','edgecolor','none','facecolor',[0.5 0.5 0.5],'face alpha',0.5)
camlight
ft_plot_sens(elec_aligned,'style','sr');

if false
  % this should only run in interactive mode
  cfg = [];
  cfg.method = 'interactive';
  cfg.elec = elec_aligned;
  cfg.headshape = vol;
  elec_aligned2 = ft_electroderealign(cfg);
else
  % this is more or less the same as the interactive alignment
  transform = [
    1 0 0 15
    0 1 0 0
    0 0 1 -10
    0 0 0 1
    ];
  elec_aligned2 = ft_transform_sens(transform, elec_aligned);
end

%%%%%%%%%

cfg = [];
cfg.vol = vol;
cfg.grid.resolution = 10;
cfg.elec = elec_aligned2;
cfg.grid.unit = 'mm';
source = ft_prepare_sourcemodel(cfg);


