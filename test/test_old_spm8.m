function test_old_spm8

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_old_spm8

% the following m-files contain a line with 'spm_' and have to be dealt with
%
%./external/bemcp/Makefile
%./spmcalls.txt
%./private/tal2mni.m
%./private/prepare_mesh_segmentation.m
%./private/prepare_dipole_grid.m
%./private/mni2tal.m
%./multivariate/toolboxes/gerven/bayesbrain/examples/html/fmri_bn_demo.html
%./multivariate/toolboxes/gerven/bayesbrain/examples/fmri_bn_demo.m
%./testing/test_bug62.m
%
%%done
%./private/avw_img_read.m
%./fileio/private/avw_img_read.m
%./fileio/read_mri.m
%./fileio/ft_read_mri.m
%
%%done and changed
%./ft_volumesegment.m
%./ft_volumenormalise.m
%./ft_volumedownsample.m
%./private/volumewrite_spm.m


mrifile = '/home/common/meg_mri/schoffelen_jm/schoffelen_jm.mri';
mri     = ft_read_mri(mrifile);

%ft_volumenormalise
cfg            = [];
cfg.nonlinear  = 'no';
cfg.coordinates = 'ctf';

restoredefaultpath
clear global;
addpath('/home/coherence/jansch/matlab/fieldtrip/');
fieldtripdefs;

clear fun;
cfg.spmversion = 'spm2';
n2 = ft_volumenormalise(cfg, mri);
rmpath(spm('dir'));

clear fun;
cfg.spmversion = 'spm8';
n8 = ft_volumenormalise(cfg, mri);
rmpath(spm('dir'));

%ft_volumesegment
cfg             = [];
cfg.coordinates = 'ctf';

clear fun;
cfg.spmversion  = 'spm2';
s2 = ft_volumesegment(cfg, mri);
rmpath(spm('dir'));

clear fun;
cfg.spmversion  = 'spm8';
s8 = ft_volumesegment(cfg, mri);
rmpath(spm('dir'));

%ft_volumedownsample
tmp            = randn(256,256,256);
mri.pow        = tmp;
mri.pow(1)     = 0;
cfg            = [];
cfg.downsample = 2;
cfg.spmversion = 'spm2';
cfg.smooth     = 'no';
d2             = ft_volumedownsample(cfg, mri);
cfg.smooth     = 5;
d2s            = ft_volumedownsample(cfg, mri);
rmpath(spm('dir'));

mri.pow        = tmp;
mri.pow(1)     = 0;
cfg            = [];
cfg.downsample = 2;
cfg.spmversion = 'spm8';
cfg.smooth     = 'no';
d8             = ft_volumedownsample(cfg, mri);
cfg.smooth     = 5;
d8s            = ft_volumedownsample(cfg, mri);
rmpath(spm('dir'));

%mni2tal and tal2mni
cd ~/matlab/fieldtrip/private
inpoints  = randn(100,3);
outpoints = mni2tal(inpoints);
rmpath(spm('dir'));

inpoint2s = tal2mni(outpoints);
rmpath(spm('dir'));

%prepare_dipole_grid
%prepare_mesh_segmentation
