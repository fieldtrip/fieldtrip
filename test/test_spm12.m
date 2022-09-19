function test_spm12

% MEM 4gb
% WALLTIME 00:10:00
% DEPENDENCY

% currently (Jan, 2017) SPM12 support for:
% - ft_volumerealign
% - ft_volumedownsample (incl smoothing)
% - private/mni2tal
% - private/tal2mni
% - private/volumesmooth
% - ft_read_mri (nifti_spm, default remains spm8) - based on function deps, not real data
% - private/volumewrite_spm (default remains spm8) - based on function deps, not real data
% - utilities/private/sn2individual (called by ft_warp_apply) - based on function deps, not real data
% - utilities/private/individual2sn (called by ft_warp_apply) - forced use of spm8 since spm_invdef is missing in spm12 (in toolbox/OldSeg)
% to do (function-specific intelligence required):
% - ft_volumenormalise
% - ft_volumesegment


mrifile1 = dccnpath('/home/common/matlab/fieldtrip/data/test/latest/mri/freesurfer/T1.mgz');
mri1     = ft_read_mri(mrifile1);
mri1.coordsys = 'tal'; % I don't think this is technically correct, it should be fsaverage aka mni305

mrifile2 = dccnpath('/home/common/matlab/fieldtrip/data/test/latest/mri/nifti/single_subj_T1.nii');
mri2     = ft_read_mri(mrifile2);
mri2.coordsys = 'mni';

mrifile3 = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.mri');
mri3     = ft_read_mri(mrifile3);
mri3.coordsys = 'ctf';

%----------------------------- SPM12 -----------------------------------

%ft_convert_coordsys: currently uses OldNorm sub-toolbox (i.e. SPM8)
ft_hastoolbox('spm12',1)
c2a = ft_convert_coordsys(mri3, 'tal', 1);
c2b = ft_convert_coordsys(mri3, 'tal', 2);

try
  rmpath(spm('dir'));
end

%ft_volumerealign: coregistration (used in human ecog tutorial)
cfg               = [];
cfg.method        = 'spm';
cfg.spmversion    = 'spm12';
cfg.spm.cost_fun  = 'nmi';
r12 = ft_volumerealign(cfg, mri1, mri2);

try
  rmpath(spm('dir'));
end

%ft_volumedownsample
tmp            = randn(256,256,256);
mri1.pow        = tmp;
mri1.pow(1)     = 0;
cfg            = [];
cfg.downsample = 2;
cfg.spmversion = 'spm12';
cfg.smooth     = 'no';
d12            = ft_volumedownsample(cfg, mri1);
cfg.smooth     = 5;
d12s           = ft_volumedownsample(cfg, mri1);

try
  rmpath(spm('dir'));
end

%mni2tal and tal2mni
[v, p] = ft_version;
cd(fullfile(p, 'private'));

inpoints  = randn(100,3);
outpoints = mni2tal(inpoints);

try
  rmpath(spm('dir'));
end

inpoints = tal2mni(outpoints);

%----------------------------- SPM8 COMPAT -----------------------------

try
  rmpath(spm('dir'));
end

%ft_volumenormalise
cfg            = [];
cfg.nonlinear  = 'no';
cfg.spmversion = 'spm8';
n8 = ft_volumenormalise(cfg, mri1);

%ft_warp_apply with spm12
elecpos = ft_warp_apply(n8.params, [4 4 4; 1 1 1], 'individual2sn');

try
  rmpath(spm('dir'));
end

%ft_volumesegment
cfg             = [];
cfg.spmversion  = 'spm8';
s8 = ft_volumesegment(cfg, mri1);


%-----------------------------------------------------------------------
% Notes from test_old_spm8:
%
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
%------------------------------------------------------------------------
