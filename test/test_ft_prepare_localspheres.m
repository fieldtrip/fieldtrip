function test_ft_prepare_localspheres

% MEM 2500mb
% WALLTIME 00:10:00

% TEST test_ft_prepare_localspheres
% TEST ft_headmodel_localspheres ft_prepare_headmodel ft_prepare_localspheres

% function to test ft_headmodel_localspheres. this function is called
% by ft_prepare_headmodel

% the function should work either on an input (segmented) mri, or it
% should have a description of the geometry in the input, or it should
% have a hdmfile (string) that specifies which file to read in addition,
% it needs a description of the sensor array

success = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the data which is needed

% read in the gradiometer information
hdr  = ft_read_header('/home/common/matlab/fieldtrip/data/Subject01.ds');
grad = hdr.grad;

% read in the segmented mri
load('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/segmentedmri.mat');

% for consistency with assumed cm input to old ft_prepare_localspheres
segmentedmri = ft_convert_units(segmentedmri,'cm');

% specify the file for the headshape
hdmfile  = '/home/common/matlab/fieldtrip/data/Subject01.shape';

% read in the headshape
shape = ft_read_headshape(hdmfile);

%%%%%%%%%%%%%%%%%%%%%
% do the computations

% with an MRI in the input
cfg         = [];
cfg.method  = 'localspheres';
cfg.grad    = grad;
vol1        = ft_prepare_headmodel(cfg, segmentedmri);
% the following needs to be done to be able to make the comparison
vol1 = rmfield(vol1, 'cfg');

cfg         = [];
cfg.grad    = grad;
vol1b       = ft_prepare_localspheres(cfg, segmentedmri);
vol1b = rmfield(vol1b, 'cfg');
success     = success && isequaln(vol1, vol1b);
if ~success
  error('ft_prepare_localspheres and ft_prepare_headmodel gave different outputs');
end

% with a filename in the input
cfg         = [];
cfg.method  = 'localspheres';
cfg.grad    = grad;
cfg.headshape = hdmfile;
vol2        = ft_prepare_headmodel(cfg);
% the following needs to be done to be able to make the comparison
vol2        = rmfield(vol2, 'cfg');

cfg         = [];
cfg.grad    = grad;
cfg.headshape = hdmfile;
vol2b       = ft_prepare_localspheres(cfg);
vol2b       = rmfield(vol2b, 'cfg');
success     = success && isequaln(vol2, vol2b);
if ~success
  error('ft_prepare_localspheres and ft_prepare_headmodel gave different outputs');
end

% with a point cloud in the input
cfg         = [];
cfg.method  = 'localspheres';
cfg.grad    = grad;
vol3        = ft_prepare_headmodel(cfg,shape);
% the following needs to be done to be able to make the comparison
vol3        = rmfield(vol3, 'cfg');

cfg         = [];
cfg.headshape = shape;
cfg.grad    = grad;
vol3b       = ft_prepare_localspheres(cfg);
vol3b       = rmfield(vol3b, 'cfg');
success     = success && isequaln(vol3, vol3b);
if ~success
  error('ft_prepare_localspheres and ft_prepare_headmodel gave different outputs');
end
