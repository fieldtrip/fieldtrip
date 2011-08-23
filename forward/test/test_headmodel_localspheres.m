% script to test ft_headmodel_localspheres. this function is called by
% ft_prepare_headmodel

% the function should work either on an input (segmented) mri, or it should
% have a description of the geometry in the input, or it should have a
% hdmfile (string) that specifies which file to read
% in addition, it needs a description of the sensor array

% read in the gradiometer information
cd('/home/common/matlab/fieldtrip/data/');
%cd('/Volumes/home/common/matlab/fieldtrip/data/');
hdr  = ft_read_header('Subject01.ds');
grad = hdr.grad;

% read in the segmented mri
cd('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer');
%cd('/Volumes/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer');
load segmentedmri
mri = segmentedmri; clear segmentedmri;

cfg         = [];
cfg.method  = 'localspheres';
cfg.grad    = grad;
vol1        = ft_prepare_headmodel(cfg, mri);

cd('/home/common/matlab/fieldtrip/data/');
%cd('/Volumes/home/common/matlab/fieldtrip/data/');
hdmfile     = 'Subject01.shape';
cfg         = [];
cfg.method  = 'localspheres';
cfg.grad    = grad;
cfg.hdmfile = hdmfile;
vol2        = ft_prepare_headmodel(cfg);