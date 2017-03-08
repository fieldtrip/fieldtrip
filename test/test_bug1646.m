function test_bug1646

% MEM 4gb
% WALLTIME 00:15:00

% TEST ft_prepare_mesh ft_datatype_segmentation

% the purpose of this test script is to ensure that the new implementation
% of ft_prepare_mesh, which is a merger between the old ft_prepare_mesh and
% ft_prepare_mesh_new, works as expected on all normal cases. What I have
% done is to search for as many as possible other tets cases and replicated
% them here.

% http://bugzilla.fcdonders.nl/show_bug.cgi?id=105 is obsolete
% http://bugzilla.fcdonders.nl/show_bug.cgi?id=115 is obsolete
% http://bugzilla.fcdonders.nl/show_bug.cgi?id=964 is obsolete

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% it is not supposed to work with a 3-D array as input
try
  cfg = [];
  bnd = ft_prepare_mesh(cfg, randn(30,30,30));
  status = false;
catch
  status = true;
end
if ~status
  error('invalid input was not dealt with')
end

% it is supposed to wotk with a headmodel as input
headmodel = [];
headmodel.r = 12;
headmodel.o = [0 0 4];
headmodel.unit = 'cm';

cfg = [];
cfg.numvertices = 100;
bnd = ft_prepare_mesh(cfg, headmodel);

cfg = [];
bnd = ft_prepare_mesh(cfg, headmodel);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following piece of test code relates loosely to
% http://bugzilla.fcdonders.nl/show_bug.cgi?id=222
% http://bugzilla.fcdonders.nl/show_bug.cgi?id=373
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.headshape = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.shape');
bnd = ft_prepare_mesh(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the same test data as used for bug 1652
% http://bugzilla.fcdonders.nl/show_bug.cgi?id=1652
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datadir = fullfile(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1646'));

load(fullfile(datadir,'seg1'));
load(fullfile(datadir,'seg2'));
load(fullfile(datadir,'seg3'));
load(fullfile(datadir,'seg4'));
load(fullfile(datadir,'seg5'));
load(fullfile(datadir,'seg6'));

atlas = ft_read_atlas(dccnpath('/home/common/matlab/fieldtrip/template/atlas/afni/TTatlas+tlrc.BRIK'));
mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/Subject01.mri'));

assert(ft_datatype(seg1, 'segmentation'));
assert(ft_datatype(seg2, 'segmentation'));
assert(ft_datatype(seg3, 'segmentation'));
assert(ft_datatype(seg4, 'segmentation'));
assert(ft_datatype(seg5, 'segmentation'));
assert(ft_datatype(seg6, 'segmentation'));
assert(ft_datatype(atlas, 'segmentation'));
assert(ft_datatype(mri, 'volume'));

cfg = [];
cfg.numvertices = 642;
cfg.tissue = 'brain';
bnd = ft_prepare_mesh(cfg, seg1);
assert(isequal(bnd(1).unit, seg1.unit));

% this is a weird case, because seg2 is not indexed but tpm
% the output depends on the order of the cfg.tissue
cfg = [];
cfg.numvertices = [1000, 100, 10];
cfg.tissue = [1 2 3];
bndA = ft_prepare_mesh(cfg, seg2);
cfg.numvertices = [10, 100, 1000];
cfg.tissue = [3 2 1];
bndB = ft_prepare_mesh(cfg, seg2);
assert(isequal(bndA(1), bndB(3)));

cfg = [];
cfg.numvertices = [1000, 100, 10];
cfg.tissue = {'brain', 'skull', 'scalp'};
bndA = ft_prepare_mesh(cfg, seg2);
assert(isequal(bndA(1).unit, seg2.unit));
assert(size(bndA(1).pos,1)==1000);
assert(size(bndA(2).pos,1)==100);
assert(size(bndA(3).pos,1)==10);
% check that the output has the desired order (as specified in the cfg)
cfg.numvertices = [10, 100, 1000];
cfg.tissue = {'scalp', 'skull', 'brain'};
bndB = ft_prepare_mesh(cfg, seg2);
assert(isequal(bndA(1), bndB(3)));

cfg = [];
cfg.numvertices = [1000, 100, 10];
cfg.tissue = {'brain', 'skull', 'scalp'};
bnd = ft_prepare_mesh(cfg, seg3);
assert(isequal(bnd(1).unit, seg3.unit));
assert(size(bnd(1).pos,1)==1000);
assert(size(bnd(2).pos,1)==100);
assert(size(bnd(3).pos,1)==10);

cfg = [];
cfg.numvertices = [1000, 100, 10];
cfg.tissue = {'brain', 'skull', 'scalp'};
bnd = ft_prepare_mesh(cfg, seg4);
assert(isequal(bnd(1).unit, seg4.unit));
assert(size(bnd(1).pos,1)==1000);
assert(size(bnd(2).pos,1)==100);
assert(size(bnd(3).pos,1)==10);

cfg = [];
cfg.numvertices = [1000, 100, 10];
cfg.tissue = [3 2 1];
bnd = ft_prepare_mesh(cfg, seg5);
assert(isequal(bnd(1).unit, seg5.unit));
assert(size(bnd(1).pos,1)==1000);
assert(size(bnd(2).pos,1)==100);
assert(size(bnd(3).pos,1)==10);

cfg = [];
cfg.numvertices = 642;
cfg.tissue = 'brain';
bnd = ft_prepare_mesh(cfg, seg6);
assert(isequal(bnd(1).unit, seg6.unit));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addressing further issues raised in bug1646, on cfg.tissue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If the input is an indexed representation and  cfg.tissue is equal to
% setdiff(unique(mri.seg(:)),0) then there is no problem.

% seg5 is indexed segmentation
cfg = [];
cfg.tissue = setdiff(unique(seg5.seg(:)),0);
cfg.numvertices = 1000;
bnd = ft_prepare_mesh(cfg, seg5);

% "If the input is an indexed representation and cfg.tissue is empty, then
% cfg.tissue should be set to unique(mri.seg(:))"

cfg = [];
cfg.numvertices = [1000 1000 1000];
cfg.tissue = [];
bnd = ft_prepare_mesh(cfg, seg5);

cfg = [];
cfg.numvertices = [1000 1000 1000];
bnd = ft_prepare_mesh(cfg, seg5);

cfg = [];
cfg.numvertices = [1000]; % should still be converted to vector Nx3 since seg5 has 3 indices
cfg.tissue = [];
bnd = ft_prepare_mesh(cfg, seg5);

% "If the input is a tpm representation then cfg.tissue should be a string
% pointing to the field."
% Both seg1 and seg2 are tpm, but only seg2 has 'brain'.

cfg = [];
cfg.numvertices = [800];
cfg.tissue = 'brain'; %this is exception in that .brain DNE in seg1
bnd = ft_prepare_mesh(cfg, seg1);

try
  cfg = [];
  cfg.numvertices = [800];
  cfg.tissue = 'randomfieldname';
  bnd = ft_prepare_mesh(cfg, seg1);
  success = true;
catch
  success = false;
end
if success
  error('randomfieldname should not work for tpm type')
end

try
  cfg = [];
  cfg.numvertices = [800];
  cfg.tissue = 'csf'; %this is exception in that .brain DNE in seg1
  bnd = ft_prepare_mesh(cfg, seg2);
  success = true;
catch
  success = false;
end
if success
  error('randomfieldname should not work for tpm type')
end

cfg = [];
cfg.numvertices = [800];
cfg.tissue = 1;
bnd = ft_prepare_mesh(cfg, seg1);

cfg = [];
cfg.numvertices = [800];
cfg.tissue = 1;
bnd = ft_prepare_mesh(cfg, seg2);

try
  cfg = [];
  cfg.numvertices = [800];
  cfg.tissue = 4;
  bnd = ft_prepare_mesh(cfg, seg2);
  success = true;
catch
  success = false;
end
if success
  error('too large of index for cfg.tissue should not work for tpm type')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this shares the test data with bug 1651
% http://bugzilla.fcdonders.nl/show_bug.cgi?id=1651
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1651.mat'));

cfg = [];
cfg.tissue = {'brain', 'skull', 'scalp'};
cfg.numvertices = [1000 2000 3000];
bnd = ft_prepare_mesh(cfg, seg2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% http://bugzilla.fcdonders.nl/show_bug.cgi?id=997
% this one does not state an explicit way of reproducing the problems, it rather refers to 937
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% http://bugzilla.fcdonders.nl/show_bug.cgi?id=937
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug937.mat'));

mri = [];
mri.anatomy = bkgrnd;
mri.sphere1 = bkgrnd==1;
mri.sphere2 = bkgrnd==2;
mri.sphere3 = bkgrnd==3;
mri.transform = eye(4);
mri.dim = size(bkgrnd);

mri.anatomy = mri.anatomy>0;

cfg = [];
cfg.tissue = {'sphere1' 'sphere2' 'sphere3'};
cfg.numvertices = 1000;
bnd = ft_prepare_mesh(cfg, mri);

