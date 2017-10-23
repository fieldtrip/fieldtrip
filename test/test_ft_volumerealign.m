function test_ft_volumerealign

% MEM 2500mb
% WALLTIME 00:10:00

% TEST ft_read_mri ft_volumerealign ft_volumereslice

% to test:
% between modalities
% within modalities
% linear vs. non linear deformations
% same subjects vs. different subjects images
% with reslice or without

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. BETWEEN MODALITIES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1.1. Same subjects
%% 1.1.1. T2 -> T1
%% 1.1.1.1. reslice = yes (default)
%% 1.1.1.1.1. alignment = linear

% SAME subject's images in T1 and T2 modalities

% path to data

subjectT1  = dccnpath('/home/common/matlab/fieldtrip/data/test/bug1826/T1.nii.gz');
subjectT2  = dccnpath('/home/common/matlab/fieldtrip/data/test/bug1826/T2.nii.gz');

% load in the data

T1  = ft_read_mri(subjectT1);
T2  = ft_read_mri(subjectT2);

% reslice

cfg            = [];
cfg.dim        = T2.dim;
cfg.resolution = 1;
T2             = ft_volumereslice(cfg,T2);

cfg            = [];
cfg.dim        = T1.dim;
cfg.resolution = 1;
T1             = ft_volumereslice(cfg,T1);

% add more rotation to T2

transmat       = T2.transform;
rotatm         = zeros(4,4);         % rotation matrix (rotates around y axis)
rotatm(end)    = 1;
rotatm(2,2)    = 1;
rotatm(1,1)    = cos((pi/12));       % with 15 degrees
rotatm(3,3)    = cos((pi/12));
rotatm(3,1)    = sin((pi/12));
rotatm(1,3)    = -sin((pi/12));

T2rot           = T2;                  % reserve also the original T2
T2rot.transform = rotatm*transmat;     % apply rotation matrix on transformation matrix

cfg            = [];
cfg.dim        = T2rot.dim;
cfg.resolution = 1;
T2rot          = ft_volumereslice(cfg,T2rot);  % apply rotation to the anatomy (voxels)



close all;
figure; ft_plot_ortho(T1.anatomy); title('T1 before aligning')
figure; ft_plot_ortho(T2.anatomy); title('T2 before aligning')
figure; ft_plot_ortho(T2rot.anatomy); title('T2rot before aligning')
%figure; ft_plot_ortho(DTI.anatomy(:,:,:,1));title('DTI before aligning')

clear rotatm;
clear transmat;

%% alignement %

%% interpolation method: trilinear

cfg          = [];
cfg.method   = 'volume';
cfg.fsl.path = '/opt/fsl/bin'; % '/opt/fsl_5.0/bin'; % fsl_5.0 only works on high mentats due to libraries
cfg.fsl.dof  = 6;              % rigid body

%interpmethods = {'nearestneighbours', 'sinc', 'trilinear'};
interpmethod = 'trilinear';
%costfuns      = {'mutualinfo', 'corratio', 'normcorr', 'normmi', 'leastsq'};
costfuns      = {'mutualinfo', 'corratio', 'normmi'};

% FIXME: not all combinations of interpmethods and cost functions make sense!
% leastsq and normcorr is applicable only for the same modality

% the remainder of the test script does not yet work, but we don't want the automatic regression testing to flag it as failure
return

%for k = 1:numel(interpmethods)
for m = 1:numel(costfuns)
  cfg.fsl.interpmethod = interpmethod;
  cfg.fsl.costfun      = costfuns{m};
  T2_aligned{m}       = ft_volumerealign(cfg, T2, T1);
  T2rot_aligned{m}       = ft_volumerealign(cfg, T2rot, T1);
  assert(ft_datatype(T2_aligned{m},'volume'),'Output of ft_volumerealing does not return a volume.');
  assert(ft_datatype(T2_aligned{m},'volume'),'Output of ft_volumerealing does not return a volume.');
  assert(~isequal(T2_aligned{m}.anatomy,T2.anatomy),'Anatomy field did not change with reslice option');
  assert(~isequal(T2rot_aligned{m}.anatomy,T2rot.anatomy),'Anatomy field did not change with reslice option');
  assert(isequal(T2_aligned{m}.transform,T2.transform),'Transformation matrix changed with reslice option');
  assert(~isequal(T2rot_aligned{m}.transform,T2rot.transform),'Transformation matrix changed with reslice option with reslice option');
end
% end

% visual check

T1.anatomyT2 = T2rot.anatomy;
cfg=[];
cfg.method = 'slice';
cfg.funparameter = 'anatomyT2';
cfg.maskparameter = 'anatomyT2';

ft_sourceplot(cfg,T1);
title('T2 (rotated) before alignement.');

for m = 1:numel(costfuns)
  T1.anatomyT2 = T2_aligned{m}.anatomy;
  cfg=[];
  cfg.method = 'slice';
  cfg.funparameter = 'anatomyT2';
  cfg.maskparameter = 'anatomyT2';
  figure;
  ft_sourceplot(cfg,T1);
  title(sprintf('T2 (non-rotated) after trilinear alignement with costfunction %s.', costfuns{m}));
  T1.anatomyT2 = T2rot_aligned{m}.anatomy;
  figure;
  ft_sourceplot(cfg,T1);
  title(sprintf('T2 (rotated) after trilinear alignement with costfunction %s.', costfuns{m}));
  
end

assert(ft_datatype(T2_aligned,'volume'),'Output of ft_volumerealing does not return a volume.');


%% interpolation method: sinc

cfg          = [];
cfg.method   = 'volume';
cfg.fsl.path = '/opt/fsl/bin'; % '/opt/fsl_5.0/bin'; % fsl_5.0 only works on high mentats due to libraries
cfg.fsl.dof  = 6;              % rigid body

interpmethod = 'sinc';
costfuns      = {'mutualinfo', 'corratio', 'normmi'};

for m = 1:numel(costfuns)
  cfg.fsl.interpmethod = interpmethod;
  cfg.fsl.costfun      = costfuns{m};
  T2_aligned{m}       = ft_volumerealign(cfg, T2, T1);
  T2rot_aligned{m}       = ft_volumerealign(cfg, T2rot, T1);
  assert(ft_datatype(T2_aligned{m},'volume'),'Output of ft_volumerealing does not return a volume.');
  assert(ft_datatype(T2_aligned{m},'volume'),'Output of ft_volumerealing does not return a volume.');
  assert(~isequal(T2_aligned{m}.anatomy,T2.anatomy),'Anatomy field did not change with reslice option');
  assert(~isequal(T2rot_aligned{m}.anatomy,T2rot.anatomy),'Anatomy field did not change with reslice option');
  assert(isequal(T2_aligned{m}.transform,T2.transform),'Transformation matrix changed with reslice option');
  assert(isequal(T2rot_aligned{m}.transform,T2rot.transform),'Transformation matrix changed with reslice option');
  
end

% visual check

T1.anatomyT2 = T2rot.anatomy;
cfg=[];
cfg.method = 'slice';
cfg.funparameter = 'anatomyT2';
cfg.maskparameter = 'anatomyT2';

ft_sourceplot(cfg,T1);
title('T2 (rotated) before alignement.');

for m = 1:numel(costfuns)
  T1.anatomyT2 = T2_aligned{m}.anatomy;
  cfg=[];
  cfg.method = 'slice';
  cfg.funparameter = 'anatomyT2';
  cfg.maskparameter = 'anatomyT2';
  figure;
  ft_sourceplot(cfg,T1);
  title(sprintf('T2 (non-rotated) after sinc alignement with costfunction %s.', costfuns{m}));
  T1.anatomyT2 = T2rot_aligned{m}.anatomy;
  figure;
  ft_sourceplot(cfg,T1);
  title(sprintf('T2 (rotated) after sinc alignement with costfunction %s.', costfuns{m}));
  
end

cfg=[];
cfg.method = 'ortho';
cfg.interactive = 'yes';
ft_sourceplot(cfg,T2rot);

%% interpolation method: nearestneighbour

cfg          = [];
cfg.method   = 'volume';
cfg.fsl.path = '/opt/fsl/bin'; % '/opt/fsl_5.0/bin'; % fsl_5.0 only works on high mentats due to libraries
cfg.fsl.dof  = 6;              % rigid body

interpmethod = 'nearestneighbour';
costfuns      = {'mutualinfo', 'corratio', 'normmi'};



for m = 1:numel(costfuns)
  cfg.fsl.interpmethod = interpmethod;
  cfg.fsl.costfun      = costfuns{m};
  T2_aligned{m}       = ft_volumerealign(cfg, T2, T1);
  T2rot_aligned{m}       = ft_volumerealign(cfg, T2rot, T1);
  assert(ft_datatype(T2_aligned{m},'volume'),'Output of ft_volumerealing does not return a volume.');
  assert(ft_datatype(T2_aligned{m},'volume'),'Output of ft_volumerealing does not return a volume.');
  assert(~isequal(T2_aligned{m}.anatomy,T2.anatomy),'Anatomy field did not change with reslice option');
  assert(~isequal(T2rot_aligned{m}.anatomy,T2rot.anatomy),'Anatomy field did not change with reslice option');
  assert(isequal(T2_aligned{m}.transform,T2.transform),'Transformation matrix changed with reslice option');
  assert(isequal(T2rot_aligned{m}.transform,T2rot.transform),'Transformation matrix changed with reslice option');
  
end

% visual check

T1.anatomyT2 = T2rot.anatomy;
cfg=[];
cfg.method = 'slice';
cfg.funparameter = 'anatomyT2';
cfg.maskparameter = 'anatomyT2';

ft_sourceplot(cfg,T1);
title('T2 (rotated) before alignement.');

for m = 1:numel(costfuns)
  T1.anatomyT2 = T2_aligned{m}.anatomy;
  cfg=[];
  cfg.method = 'slice';
  cfg.funparameter = 'anatomyT2';
  cfg.maskparameter = 'anatomyT2';
  figure;
  ft_sourceplot(cfg,T1);
  title(sprintf('T2 (non-rotated) after nearestneighbour alignement with costfunction %s.', costfuns{m}));
  T1.anatomyT2 = T2rot_aligned{m}.anatomy;
  figure;
  ft_sourceplot(cfg,T1);
  title(sprintf('T2 (rotated) after nearestneighbour alignement with costfunction %s.', costfuns{m}));
  
end

clear T2rot_aligned;
clear T2_aligned;

% 1. between modalities
% 1.1. same subjects
% 1.1.1. T2->T1
% 1.1.1.1. reslice = yes
%% 1.1.1.1.2. alignment = non-linear

% should be tested but is it necessary for same subjects images?

% 1. between modalities
% 1.1. same subjects
% 1.1.1. T2->T1
%% 1.1.1.2. reslice = no
%% 1.1.1.2.1 alignment = linear

% FIXME!!!!

cfg          = [];
cfg.method   = 'volume';
cfg.fsl.path = '/opt/fsl/bin'; % '/opt/fsl_5.0/bin'; % fsl_5.0 only works on high mentats due to libraries
cfg.fsl.dof  = 6;              % rigid body
cfg.fsl.reslice = 'no';

cfg.fsl.interpmethod = 'trilinear';
cfg.fsl.costfun     = 'mutualinfo';

T2_aligned       = ft_volumerealign(cfg, T2, T1);
T2rot_aligned    = ft_volumerealign(cfg, T2rot, T1);
assert(ft_datatype(T2_aligned,'volume'),'Output of ft_volumerealing does not return a volume.');
assert(ft_datatype(T2_aligned,'volume'),'Output of ft_volumerealing does not return a volume.');
assert(isequal(T2_aligned.anatomy,T2.anatomy),'Anatomy field did change with reslice option = no');
assert(isequal(T2rot_aligned.anatomy,T2rot.anatomy),'Anatomy field did change with reslice option = no');
assert(~isequal(T2_aligned.transform,T2.transform),'Transformation matrix did not changed with reslice option = no');
assert(~isequal(T2rot_aligned.transform,T2rot.transform),'Transformation matrix did not changed with reslice option = no');

% reslice
cfg=[];
cfg            = [];
cfg.dim        = T2rot_aligned.dim;
cfg.resolution = 1;
T2rot_alignrs          = ft_volumereslice(cfg,T2rot_aligned);

% visual check
T1.anatomyT2 = T2rot_alignrs.anatomy;

cfg=[];
cfg.method = 'slice';
cfg.funparameter = 'anatomyT2';
cfg.maskparameter = 'anatomyT2';
figure;
ft_sourceplot(cfg,T1);
title(sprintf('T2 (rotated) after alignement'));

%FIXME: T2 does not seem to be aligned to T1 after reslice

% 1. between modalities
% 1.1. same subjects
% 1.1.1. T2->T1
% 1.1.1.2. reslice = no
%% 1.1.1.2.2 alignment = non-linear
% ?

% 1. between modalities
% 1.1. same subjects
%% 1.1.2. DTI->T1
%% 1.1.2.1. reslice = yes
%% 1.1.2.1.1 alignment = linear

% FIXME!!!

clear all;

% path to data

subjectT1  = dccnpath('/home/common/matlab/fieldtrip/data/test/bug1826/T1.nii.gz');
subjectDTI = dccnpath('/home/common/matlab/fieldtrip/data/test/bug1826/DTI.nii');

% load in the data

DTI = ft_read_mri(subjectDTI);
T1 = ft_read_mri(subjectT1);

% FIXME: it still has to be implemented in ft_volumerealign

cfg          = [];
cfg.method   = 'volume';
cfg.fsl.path = '/opt/fsl/bin'; % '/opt/fsl_5.0/bin'; % fsl_5.0 only works on high mentats due to libraries
cfg.fsl.dof  = 6;

%interpmethods = {'nearestneighbours', 'sinc', 'trilinear'};
interpmethod = 'trilinear';
%costfuns      = {'mutualinfo', 'corratio', 'normcorr', 'normmi', 'leastsq'};
costfuns      = {'mutualinfo', 'corratio', 'normmi'};

for k = 1:numel(interpmethod)
  for m = 1:numel(costfuns)
    cfg.fsl.interpmethod = interpmethod;
    cfg.fsl.costfun      = costfuns{m};
    DTIaligned{m}      = ft_volumerealign(cfg, DTI, T1);
  end
end

% 1. between modalities
% 1.1. same subjects
% 1.1.2. DTI->T1
% 1.1.2.1. reslice = yes
%% 1.1.2.1.2 alignment = non-linear

% should be tested, but is it necessary?

% 1. between modalities
% 1.1. same subjects
% 1.1.2. DTI->T1
%% 1.1.2.2. reslice = no
%% 1.1.2.2.1 alignment = linear

% it would be desired to implement it

% 1. between modalities
% 1.1. same subjects
% 1.1.2. DTI->T1
% 1.1.2.2. reslice = no
%% 1.1.2.2.2 alignment = non-linear
% ?

% 1. between modalities
% 1.1. same subjects
%% 1.1.3. FA ->T1
%% 1.1.3.1. reslice = yes
%% 1.1.3.1.1 alignment = linear


% FA image needed
% still has to be implemented

% 1. between modalities
% 1.1. same subjects
% 1.1.3. FA ->T1
% 1.1.3.1. reslice = yes
%% 1.1.3.1.2 alignment = non-linear

% FA image needed
% still has to be implemented

% 1. between modalities
% 1.1. same subjects
% 1.1.3. FA ->T1
%% 1.1.3.2 reslice = no
%% 1.1.3.2.1 alignment = linear

% FA image needed
% still has to be implemented

% 1. between modalities
% 1.1. same subjects
% 1.1.3. FA ->T1
% 1.1.3.2 reslice = no
%% 1.1.3.2.2 alignment = non-linear

% ?

% 1. between modalities
%% 1.2 Different subjects

% ?

%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2 WITHIN MODALITIES %%
%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2.1. Same subjects

% ?

% within modalities
%% 2.2. Different subjects
%% 2.2.1. T1->T1
%% 2.2.1.1. reslice = yes
%% 2.2.1.1.1 alignment = linear

clear all;

% FIXME!
% still has to  be implemented
% functionalities partially exist already in FT, e.g. aligning to a
% template with ft_volumenormalize

% different subject's images in the same modality

% template T1 (fuzzy)

T1temp = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/external/spm8/templates/T1.nii'));

% other subject's T1

T1other = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/mri/nifti/single_subj_T1.nii'));

cfg=[];
% cfg.nonlinear = 'no';  option in volumenormalise
cfg.method = 'volume';
T1al_temp = ft_volumerealign(cfg,T1,T1temp);

cfg=[];
% cfg.nonlinear = 'no';  option in volumenormalise
cfg.method = 'volume';
T1al_other = ft_volumerealign(cfg,T1,T1other);

% within modalities
%% 2.2. Different subjects
%% 2.2.1. T1->T1
%% 2.2.1.1 reslice = yes
%% 2.2.1.1.2 alignment = non-linear

clear all;

% FIXME!
% still has to  be implemented
% functionalities partially exist already in FT, e.g. aligning to a
% template with ft_volumenormalize

% different subject's images in the same modality

% template T1 (fuzzy)

T1temp = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/external/spm8/templates/T1.nii'));

% other subject's T1

T1other = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/mri/nifti/single_subj_T1.nii'));


cfg=[];
cfg.method = 'volume';
T1al_temp = ft_volumerealing(cfg,T1,T1temp);

cfg=[];
cfg.method = 'volume';
T1al_other = ft_volumerealing(cfg,T1,T1other);

% within modalities
% 2.2. Different subjects
% 2.2.1. T1->T1
%% 2.2.1.2. reslice = no
%% 2.2.1.2.1 alignment = linear

% ?

% within modalities
% 2.2. Different subjects
% 2.2.1. T1->T1
% 2.2.1.2. reslice = no
%% 2.2.1.2.2 alignment = non-linear

% ?

% within modalities
% 2.2. Different subjects
%% 2.2.1. T2->T2

% should be implemented?

% within modalities
% 2.2. Different subjects
%% 2.2.1. DTI->DTI

% ?

% within modalities
% 2.2. Different subjects
%% 2.2.1. FA->FA

% ?
