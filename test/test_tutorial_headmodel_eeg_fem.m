function test_tutorial_headmodel_eeg_fem

% WALLTIME 08:00:00
% MEM 12gb
% DEPENDENCY ft_prepare_headmodel ft_prepare_mesh ft_datatype_segmentation

mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.mri'));

% this needs to be done before reslicing
nas = mri.hdr.fiducial.mri.nas;
lpa = mri.hdr.fiducial.mri.lpa;
rpa = mri.hdr.fiducial.mri.rpa;

vox2head = mri.transform;
nas = ft_warp_apply(vox2head, nas, 'homogenous');
lpa = ft_warp_apply(vox2head, lpa, 'homogenous');
rpa = ft_warp_apply(vox2head, rpa, 'homogenous');

cfg     = [];
cfg.dim = mri.dim;
mri     = ft_volumereslice(cfg,mri);

cfg           = [];
cfg.output    = {'gray','white','csf','skull','scalp'};
segmentedmri  = ft_volumesegment(cfg, mri);

seg_i = ft_datatype_segmentation(segmentedmri,'segmentationstyle','indexed');

cfg        = [];
cfg.shift  = 0.3;
cfg.method = 'hexahedral';
mesh = ft_prepare_mesh(cfg,seg_i);

cfg        = [];
cfg.method ='simbio';
cfg.conductivity = [0.33 0.14 1.79 0.01 0.43];   % order follows mesh.tissyelabel
headmodel  = ft_prepare_headmodel(cfg, mesh);

ft_plot_mesh(mesh, 'surfaceonly', 'yes');

elec = ft_read_sens('standard_1020.elc');

figure
hold on
ft_plot_mesh(mesh,'surfaceonly','yes','vertexcolor','none','edgecolor','none','facecolor',[0.5 0.5 0.5],'face alpha',0.7)
camlight

ft_plot_sens(elec);


% Then, we determine the translation and rotation that is needed to get the position of the fiducials in the electrode structure (defined with labels 'Nz', 'LPA', 'RPA') to their counterparts in the CTF head coordinate system that we acquired from the anatomical mri (nas, lpa, rpa).
%
% create a structure similar to a template set of electrodes
fid.pos           = [nas; lpa; rpa];       % CTF head coordinates of fiducials
fid.label         = {'Nz','LPA','RPA'};    % use the same labels as those in elec
fid.unit          = 'mm';                  % use the same units as those in mri

% alignment
cfg               = [];
cfg.method        = 'fiducial';
cfg.elec          = elec;                  % the electrodes we want to align
cfg.target        = fid;                   % the template we want to align to
cfg.fiducial      = {'Nz', 'LPA', 'RPA'};  % labels of fiducials in fid and in elec
elec_aligned      = ft_electroderealign(cfg);

% We can check the alignment by plotting together the scalp surface with the electrodes.
%
figure;
hold on;
ft_plot_mesh(mesh,'surfaceonly','yes','vertexcolor','none','edgecolor','none','facecolor',[0.5 0.5 0.5],'face alpha',0.5)
camlight
ft_plot_sens(elec_aligned);

%
%
% The alignment is much better, but not perfect. Some of the electrodes are below the head surface in the front, while the electrodes in the back do not fit tightly to the head. The remaining misalignment is due to the use of different conventions to [define the fiducials](/faq/how_are_the_lpa_and_rpa_points_defined). We can improve the alignment of the electrodes interactively.
%
%% ## Interactive alignment
%
if false
  % this should not run in the non-interactive test environment
  cfg          = [];
  cfg.method   = 'interactive';
  cfg.elec     = elec_aligned;
  cfg.headshape = headmodel;
  elec_aligned = ft_electroderealign(cfg);
end
