function test_pull1745

% WALLTIME 00:10:00
% MEMORY 8gb
% DEPENDENCY ft_datatype_source ft_datatype_volume ft_datatype_parcellation ft_datatype_segmentation ft_sourceparcellate

filename3d = dccnpath('/home/common/matlab/fieldtrip/data/test/original/mri/nifti/sub-01_ses-mri_acq-mprage_T1w.nii');                    % 3D anatomical
filename4d = dccnpath('/home/common/matlab/fieldtrip/data/test/original/mri/nifti/sub-01_ses-mri_task-facerecognition_run-01_bold.nii');  % 4D functional

%%

volume3d = ft_read_mri(filename3d);
volume4d = ft_read_mri(filename4d, 'outputfield', 'functional');

volume4d.avg = mean(volume4d.functional, 4);
volume4d.std = std(volume4d.functional, [], 4);

tr = volume4d.hdr.niftihdr.pixdim(4);

volume3d = ft_checkdata(volume3d, 'datatype', 'volume', 'insidestyle', 'logical');
volume4d = ft_checkdata(volume4d, 'datatype', 'volume', 'insidestyle', 'logical');

%%

figure
ft_plot_ortho(volume3d.anatomy, 'transform', volume3d.transform, 'style', 'intersect');
ft_plot_ortho(volume4d.std, 'transform', volume4d.transform, 'style', 'intersect', 'colormap', 'jet');

%%

source3d = ft_checkdata(volume3d, 'datatype', 'source');
source4d = ft_checkdata(volume4d, 'datatype', 'source'); % note the warning: could not determine dimord of "functional"

% add the time axis that describes the 4th dimension of the functional data
source4d.time = (0:size(volume4d.functional,4)-1) * tr;
% reshape the functional data
source4d = ft_checkdata(source4d, 'datatype', 'source');

%%

cfg = [];
cfg.anaparameter = [];
cfg.funparameter = 'std';
cfg.funcolormap = 'jet';
ft_sourceplot(cfg, source4d);

%%

% the next sequence of steps does
% - define ROIs on basis of the original 3D anatomical data
% - interpolate the ROIs from the 3D onto the 4D representation
% - average the non-interpolated functional data over some parcels

% make a structure with subject-specific ROIs from the original 3D anatomical representation
parcellation3d = rmfield(source3d, 'anatomy');

% make three ROIs or parcels: each is a slab that spans 1/3rd of the whole volume
% this first representation is "probabilistic"
parcellation3d.roi1 = false(parcellation3d.dim); parcellation3d.roi1(1:end,1:end,1:64) = true;
parcellation3d.roi2 = false(parcellation3d.dim); parcellation3d.roi2(1:end,1:end,65:128) = true;
parcellation3d.roi3 = false(parcellation3d.dim); parcellation3d.roi3(1:end,1:end,129:192) = true;

% this converts it to an "indexed" representation
parcellation3d = ft_checkdata(parcellation3d, 'type', 'parcellation', 'parcellationstyle', 'indexed')

%

cfg = [];
cfg.anaparameter = [];
cfg.funparameter = 'tissue';
cfg.funcolormap = 'jet';
ft_sourceplot(cfg, parcellation3d);

%

% interpolate the 3D tissue parcellation onto the voxel positions of the lower resolution 4D representation
cfg = [];
cfg.parameter = 'all';
cfg.interpmethod = 'nearest';
parcellation4d = ft_sourceinterpolate(cfg, parcellation3d, source4d)

%

cfg = [];
cfg.anaparameter = [];
cfg.funparameter = 'tissue';
cfg.funcolormap = 'jet';
ft_sourceplot(cfg, parcellation4d);

%

cfg = [];
cfg.parcellation = 'tissue';
rawdata3d = ft_sourceparcellate(cfg, source3d, parcellation3d);
rawdata4d = ft_sourceparcellate(cfg, source4d, parcellation4d);

%%

[ftver, ftpath] = ft_version;
atlas = ft_read_atlas(fullfile(ftpath, 'template', 'atlas', 'aal', 'ROI_MNI_V4.nii'));

%%

% the next sequence of steps would be to
% - interpolate the 4D functional data onto the 3D anatomical
% - spatially transform both the 3D and the 4D to the MNI template
% - interpolate the atlas onto the 3D and 4D representation
% - average the interpolated and transformed functional data over some parcels

%%

% the next sequence of steps would be to
% - spatially transform the atlas to the 3D anatomical
% - interpolate the atlas onto the 3D and 4D representation
% - average the non-interpolated functional data over some parcels


