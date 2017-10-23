function test_ft_volumesegment

% MEM 2000mb
% WALLTIME 01:00:00

% TEST ft_volumesegment  ft_read_mri

% initial version by Lilla Magyari 2012

% read in an mri
mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/mri/nifti/single_subj_T1.nii'));

% the following could also be done using ft_determine_coordsys
mri.coordsys = 'spm';

%% output:tpm

% also test the different spm versions
mfile = mfilename('fullpath');
[pathstr, mfile] = fileparts(mfile);
ftpath           = pathstr(1:(strfind(pathstr, 'test')-1));

restoredefaultpath;
addpath(ftpath);
ft_defaults;
cfg        = [];
cfg.output = 'tpm';
cfg.spmversion = 'spm2';
tpm1 = ft_volumesegment(cfg,mri);

restoredefaultpath;
addpath(ftpath);
ft_defaults;
cfg        = [];
cfg.output = 'tpm';
cfg.spmversion = 'spm8';
tpm2 = ft_volumesegment(cfg,mri);

restoredefaultpath;
addpath(ftpath);
ft_defaults;
cfg        = [];
cfg.output = 'tpm';
cfg.spmversion = 'spm12';
cfg.spmmethod  = 'old';
tpm3 = ft_volumesegment(cfg,mri);

restoredefaultpath;
addpath(ftpath);
ft_defaults;
cfg        = [];
cfg.output = 'tpm';
cfg.spmversion = 'spm12';
cfg.spmmethod  = 'new';
tpm4 = ft_volumesegment(cfg,mri);

assert(isfield(tpm1, 'gray')&&isfield(tpm1, 'white')&&isfield(tpm1, 'csf'), 'tissue probability map is missing a field');
assert(isfield(tpm2, 'gray')&&isfield(tpm2, 'white')&&isfield(tpm2, 'csf'), 'tissue probability map is missing a field');
assert(isfield(tpm3, 'gray')&&isfield(tpm3, 'white')&&isfield(tpm3, 'csf'), 'tissue probability map is missing a field');
assert(isfield(tpm4, 'gray')&&isfield(tpm4, 'white')&&isfield(tpm4, 'csf'), 'tissue probability map is missing a field');

%% output:brain
tpm = tpm2; % use one of the segmentations created above for the rest

cfg = [];
cfg.output = 'brain';
cfg.spmversion = 'spm8'; % this is the old default, and needed for the test
% function to run through, using spm12 leads to rng state dependent minor
% differences, causing the isequals below to fail.
brain = ft_volumesegment(cfg,mri);

if ~(isfield(brain,'brain'))
  error('brainmask segmentation is missing the brain');
elseif isfield(brain,'anatomy') || isfield(brain,'tpm')
  error('inconsistent segmentation structure');
end

cfg = [];
cfg.output = 'brain';
cfg.spmversion = 'spm8';
brain2 = ft_volumesegment(cfg,tpm);

if ~(isfield(brain2,'brain'))
  error('brainmask segmentation is missing the brain');
elseif isfield(brain2,'anatomy') || isfield(brain2,'tpm')
  error('inconsistent segmentation structure');
end

if ~(isequal(brain.brain,brain2.brain))
  error('brainmasks are different for mri vs. tpm inputs');
end

%% output: brain skull and scalp
tpm.anatomy = mri.anatomy; % for scalp segmentation

cfg = [];
cfg.output = {'brain' 'skull' 'scalp'};
cfg.spmversion = 'spm8';
seg = ft_volumesegment(cfg,tpm);



if ~(isfield(seg,'brain')) || ~(isfield(seg,'skull')) || ~(isfield(seg,'scalp'))
  error('segmentation is missing a field');
elseif isfield(seg,'anatomy') || isfield(seg,'tpm')
  error('inconsistent segmentation structure');
end

cfg = [];
cfg.output = {'brain' 'skull' 'scalp'};
cfg.spmversion = 'spm8';
seg2 = ft_volumesegment(cfg, mri);

if ~(isfield(seg2,'brain')) || ~(isfield(seg2,'skull')) || ~(isfield(seg2,'scalp'))
  error('segmentation is missing a field');
elseif isfield(seg2,'anatomy') || isfield(seg2,'tpm')
  error('inconsistent segmentation structure');
end

if ~(isequal(seg.brain,seg2.brain, brain.brain, brain2.brain))
  error('brainmasks are different for different inputs (mri vs tpm) or for different outputs (1 layer vs. 3)');
elseif ~(isequal(seg.skull,seg2.skull))
  error('skullmasks are different for mri vs. tpm inputs');
elseif ~(isequal(seg.scalp,seg2.scalp))
  error('skullmasks are different for mri vs. tpm inputs');
end

clear brain;
clear brain2;

%% output: gray, white, csf, skull and scalp

cfg = [];
cfg.output = {'gray' 'white' 'csf' 'skull' 'scalp'};
seg3 = ft_volumesegment(cfg,tpm);

if ~(isfield(seg3,'skull')) || ~(isfield(seg3,'scalp')) || ~(isfield(seg3,'gray')) || ~(isfield(seg3,'white')) || ~(isfield(seg3,'csf'))
  error('segmentation is missing a field');
elseif isfield(seg3,'anatomy') || isfield(seg3,'tpm') || isfield(seg3,'brain')
  error('inconsistent segmentation structure');
end

assert(islogical(seg3.gray) & islogical(seg3.white) & islogical(seg3.csf),'Output is not a binary representation.')
clear seg3;

%% output: scalp
cfg = [];
cfg.output = 'scalp';
scalp = ft_volumesegment(cfg,tpm);

if ~(isfield(scalp,'scalp'))
  error('scalpmask segmentation is missing the scalp');
elseif isfield(scalp,'anatomy') || isfield(scalp,'tpm') || isfield(scalp,'brain') || isfield(scalp, 'skull')
  error('inconsistent segmentation structure');
end

cfg = [];
cfg.output = 'scalp';
scalp2 = ft_volumesegment(cfg,mri);

if ~(isfield(scalp2,'scalp'))
  error('scalpmask segmentation is missing the scalp');
elseif isfield(scalp2,'anatomy') || isfield(scalp2,'tpm') || isfield(scalp2,'brain') || isfield(scalp2, 'skull')
  error('inconsistent segmentation structure');
end

% when scalp is the only output: the representation is "cummulative"
% for comparison, also the seg.scalp is converted to cummulative
seg.scalp(seg.brain == 1) = 1;
seg.scalp(seg.skull == 1) = 1;
seg2.scalp(seg2.brain == 1) = 1;
seg2.scalp(seg2.skull == 1) = 1;

if ~(isequal(scalp.scalp,scalp2.scalp,seg.scalp,seg2.scalp))
  error('skullmasks are different for different inputs (mri vs tpm) or for different outputs (1 layer vs. 3)');
end

clear seg;
clear seg2;
clear scalp2;
%% output: skullstrip

cfg = [];
cfg.output = 'skullstrip';
cfg.spmversion = 'spm8';
skullstr = ft_volumesegment(cfg,tpm);

if ~(isfield(skullstr,'anatomy'))
  error('skullstrip segmentation is missing anatomy');
elseif isfield(skullstr,'tpm') || isfield(skullstr,'brain') || isfield(skullstr, 'skull') || isfield(skullstr, 'scalp')
  error('inconsistent segmentation structure');
end

cfg = [];
cfg.output = 'skullstrip';
cfg.spmversion = 'spm8';
skullstr2 = ft_volumesegment(cfg,mri);

if ~(isfield(skullstr2,'anatomy'))
  error('skullstrip segmentation is missing anatomy');
elseif isfield(skullstr2,'tpm') || isfield(skullstr2,'brain') || isfield(skullstr2, 'skull') || isfield(skullstr2, 'scalp')
  error('inconsistent segmentation structure');
end

if isequal(skullstr.anatomy,tpm.anatomy) || isequal(skullstr2.anatomy,mri.anatomy)
  error('anatomy has not been segmented');
end

if ~(isequal(skullstr.anatomy,skullstr2.anatomy))
  error('skullstripped anatomies are different for different inputs (mri vs tpm)');
end
clear skullstr2;
% with no threshold and downsampling

cfg = [];
cfg.downsample = 2;
cfg.output = 'scalp';
cfg.spmversion = 'spm8';
scalp3 = ft_volumesegment(cfg,tpm);

if isequal(scalp.scalp,scalp3.scalp)
  error('scalpmasks should be different for downsampled tpm');
end
clear scalp3;

% with user-specified threshold
cfg = [];
cfg.scalpthreshold = 0.3;
cfg.scalpsmooth = 6;
cfg.output = 'scalp';
scalp4 = ft_volumesegment(cfg,tpm);

cfg = [];
cfg.scalpthreshold = 'no';
cfg.scalpsmooth = 'no';
cfg.output = 'scalp';
scalp5 = ft_volumesegment(cfg,tpm);


if isequal(scalp.scalp,scalp4.scalp) || isequal(scalp.scalp,scalp5.scalp) || isequal(scalp4.scalp,scalp5.scalp)
  error('scalpmasks should be different for different smooth and threshold values');
end
clear scalp;
clear scalp5;
% old way of specifying threshold and smooth
cfg = [];
cfg.threshold = 0.3;
cfg.smooth = 6;
cfg.output = 'scalp';
cfg.spmversion = 'spm8';
scalp6 = ft_volumesegment(cfg,tpm);

if ~(isequal(scalp4.scalp,scalp6.scalp))
  error('old and new way of specifying smooth and threshold should give the same output');
end

clear scalp4;
clear scalp6;

cfg = [];
cfg.brainthreshold = 0.1;
cfg.brainsmooth  =  6;
cfg.output = 'skullstrip';
cfg.spmversion = 'spm8';
skullstr3 = ft_volumesegment(cfg,tpm);

cfg = [];
cfg.brainthreshold = 'no';
cfg.brainsmooth  =  'no';
cfg.output = 'skullstrip';
skullstr4 = ft_volumesegment(cfg,tpm);

if isequal(skullstr3.anatomy,skullstr.anatomy) || isequal(skullstr.anatomy,skullstr4.anatomy) || isequal(skullstr3.anatomy,skullstr4.anatomy)
  error('anatomies should be different for different smooth and threshold values');
end
clear skullstr;
clear skullstr4;
cfg = [];
cfg.threshold = 0.1;
cfg.smooth  =  6;
cfg.output = 'skullstrip';
cfg.spmversion = 'spm8';
skullstr5 = ft_volumesegment(cfg,tpm);

if ~(isequal(skullstr5.anatomy,skullstr3.anatomy))
  error('old and new way of specifying smooth and threshold should give the same output');
end

clear skullstr3;

cfg = [];
cfg.threshold = 0.1;
cfg.smooth  =  6;
cfg.output =  {'skullstrip' 'brain'};
cfg.spmversion = 'spm8';
skullstr6 = ft_volumesegment(cfg,tpm);

if ~(isequal(skullstr5.anatomy,skullstr6.anatomy))
  error('specifying output only with skullstrip or with other output should give the same output');
end

%% output: skullstrip & brain (for the MNE pipeline)

cfg = [];
cfg.output = {'skullstrip' 'brain'};
cfg.spmversion = 'spm8';
skullstr = ft_volumesegment(cfg,tpm);

assert(isfield(skullstr,'anatomy') & isfield(skullstr,'brain'),'skullstrip segmentation is missing anatomy');

clear all;



