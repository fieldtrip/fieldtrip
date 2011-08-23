% script to test ft_headmodel_singleshell. this function is called by
% ft_prepare_headmodel

% the function should work either on an input (segmented) mri, or it should
% have a description of the geometry in the input, or it should have a
% hdmfile (string) that specifies which file to read

curdir = pwd;

% read in the segmented mri
cd('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer');
%cd('/Volumes/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer');
load segmentedmri
mri = segmentedmri; clear segmentedmri;

cfg         = [];
cfg.method  = 'singleshell';
vol1        = ft_prepare_headmodel(cfg, mri);

cd('/home/common/matlab/fieldtrip/data/');
%cd('/Volumes/home/common/matlab/fieldtrip/data/');
hdmfile     = 'Subject01.shape';

cfg         = [];
cfg.method  = 'singleshell';
cfg.hdmfile = hdmfile;
vol2        = ft_prepare_headmodel(cfg);

% read in the headshape
shape = ft_read_headshape(cfg.hdmfile);

cfg         = [];
cfg.method  = 'singleshell';
cfg.geom    = shape;
vol3        = ft_prepare_headmodel(cfg);

cd(curdir);