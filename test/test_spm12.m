function test_spm12

% MEM 2gb
% WALLTIME 00:10:00

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


mrifile = dccnpath('/home/common/matlab/fieldtrip/data/test/latest/mri/freesurfer/T1.mgz');
mri     = ft_read_mri(mrifile);
%mri = ft_read_mri('/Users/arjsto/Projects/Ecog/data/IR29/freesurfer/mri/T1.mgz');
mri.coordsys = 'tal';
mrifile2 = dccnpath('/home/common/matlab/fieldtrip/data/test/latest/mri/nifti/single_subj_T1.nii');
mri2     = ft_read_mri(mrifile2);
%mri2 = ft_read_mri('/Users/arjsto/Projects/Ecog/data/template/single_subj_T1_1mm.nii');
mri2.coordsys = 'mni';


%----------------------------- SPM12 -----------------------------------

%ft_convert_coordsys: currently uses OldNorm sub-toolbox (i.e. SPM8)
ft_hastoolbox('spm12',1)
mri.coordsys = 'ctf';
c2a = ft_convert_coordsys(mri, 'tal', 1);
c2b = ft_convert_coordsys(mri, 'tal', 2);

rmpath(spm('dir'));

%ft_volumerealign: coregistration (used in human ecog tutorial)
mri.coordsys = 'tal';
cfg             = [];
cfg.method      = 'spm';
cfg.spmversion  = 'spm12';
cfg.spm.coostfun = 'nmi';
cfg.coordsys    = 'tal';
r12 = ft_volumerealign(cfg, mri, mri2);

rmpath(spm('dir'));

%ft_volumedownsample
tmp            = randn(256,256,256);
mri.pow        = tmp;
mri.pow(1)     = 0;
cfg            = [];
cfg.downsample = 2;
cfg.spmversion = 'spm12';
cfg.smooth     = 'no';
d12            = ft_volumedownsample(cfg, mri);
cfg.smooth     = 5;
d12s           = ft_volumedownsample(cfg, mri);

rmpath(spm('dir'));

%mni2tal and tal2mni
cd ~/matlab/fieldtrip/private
inpoints  = randn(100,3);
outpoints = mni2tal(inpoints);

rmpath(spm('dir'));

inpoints = tal2mni(outpoints);

%----------------------------- SPM8 COMPAT -----------------------------

rmpath(spm('dir'));

%ft_volumenormalise
cfg            = [];
cfg.nonlinear  = 'no';
cfg.coordinates = 'tal';
cfg.spmversion = 'spm8';
n8 = ft_volumenormalise(cfg, mri);

%ft_warp_apply with spm12
elecpos = ft_warp_apply(n8.params, [4 4 4; 1 1 1], 'individual2sn');

rmpath(spm('dir'));

%ft_volumesegment
cfg             = [];
cfg.coordinates = 'tal';
cfg.spmversion  = 'spm8';
s8 = ft_volumesegment(cfg, mri);


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
