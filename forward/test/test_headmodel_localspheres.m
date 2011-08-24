function test_ft_headmodel_localspheres

% TEST: ft_headmodel_localspheres ft_prepare_headmodel ft_prepare_localspheres

% function to test ft_headmodel_localspheres. this function is called by
% ft_prepare_headmodel
%
% the function should work either on an input (segmented) mri, or it should
% have a description of the geometry in the input, or it should have a
% hdmfile (string) that specifies which file to read
% in addition, it needs a description of the sensor array

success = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the data which is needed

% read in the gradiometer information
cd('/home/common/matlab/fieldtrip/data/');
hdr  = ft_read_header('Subject01.ds');
grad = hdr.grad;

% read in the segmented mri
cd('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer');
load segmentedmri
mri = segmentedmri; clear segmentedmri;

%%%%%%%%%%%%%%%%%%%%%
% do the computations

% with an MRI in the input
cfg         = [];
cfg.method  = 'localspheres';
cfg.grad    = grad;
vol1        = ft_prepare_headmodel(cfg, mri);
cfg         = rmfield(cfg, 'method');
vol1b       = ft_prepare_localspheres(cfg, mri);
success     = success && isequal(vol1, vol1b);

cd('/home/common/matlab/fieldtrip/data/');
hdmfile     = 'Subject01.shape';

% with a filename in the input
cfg         = [];
cfg.method  = 'localspheres';
cfg.grad    = grad;
cfg.hdmfile = hdmfile;
vol2        = ft_prepare_headmodel(cfg);

cfg         = [];
cfg.grad    = grad;
cfg.headshape = hdmfile;
vol2b       = ft_prepare_localspheres(cfg);
success     = success && isequal(vol2, vol2b);

% read in the headshape
shape = ft_read_headshape(hdmfile);

% with a point cloud in the input
cfg         = [];
cfg.method  = 'localspheres';
cfg.grad    = grad;
cfg.geom    = shape;
vol3        = ft_prepare_headmodel(cfg);

cfg         = [];
cfg.headshape = shape;
cfg.grad    = grad;
vol3b       = ft_prepare_localspheres(cfg);
success     = success && isequal(vol3, vol3b);

if ~success
  error('there is a problem in test_headmodel_localspheres');
end
