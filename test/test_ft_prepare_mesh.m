function test_ft_prepare_mesh

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_prepare_mesh ft_datatype_segmentation ft_plot_mesh

% test ft_prepare_mesh also used for constructing SIMBIO FEM head models
% see also http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=1815

%% segmentations

example.dim = [30 29 31];   % slightly different numbers
example.transform = eye(4);   
example.coordsys = 'ctf';
example.unit = 'mm';
example.seg = zeros(example.dim);

% adjusting transformation matrix: center of the head-coordinates [0 0 0] should
% be the center of volume
% center of volume in voxel coordinates     
x = round(example.dim(1)/2);
y = round(example.dim(2)/2);
z = round(example.dim(3)/2);
    
x = round(x);
y = round(y);
z = round(z);      
             
origin = [x y z];

example.transform(1:4,4)   = [-origin(:); 1];  % head-coordinate [0 0 0] is in the center of 
                                                % the volume (x y z in voxel-coordinates)

% compute position for each voxel in voxelspace and in headspace

[X, Y, Z] = ndgrid(1:example.dim(1), 1:example.dim(2), 1:example.dim(3));
voxelpos = [X(:) Y(:) Z(:)];
headpos = ft_warp_apply(example.transform, voxelpos);

% create 5 spheres
radius1 = 14;
radius2 = 12;
radius3 = 10;
radius4 = 8;
radius5 = 6;

for i=1:size(headpos,1)
    % from small to large
    if norm(headpos(i,:))<radius5
        example.seg(i) = 5;
    elseif norm(headpos(i,:))<radius4
        example.seg(i) = 4;
    elseif norm(headpos(i,:))<radius3
        example.seg(i) = 3;
    elseif norm(headpos(i,:))<radius2
        example.seg(i) = 2;
    elseif norm(headpos(i,:))<radius1
        example.seg(i) = 1;
    end
end
clear X Y Z headpos origin radius1 radius2 radius3  voxelpos x y z 

% indexed segmentation
% 5 tissue-types
seg5 = example;
clear example;
close all;
figure; imagesc(seg5.seg(:,:,15));
% 3 tissue-types
seg3 = seg5;
seg3.seg(seg5.seg(:)==4)=3;
seg3.seg(seg5.seg(:)==5)=3;
figure; imagesc(seg3.seg(:,:,15));
% 1 tissue-types
seg1 = seg3;
seg1.seg(seg3.seg(:)==2)=1;
seg1.seg(seg3.seg(:)==3)=1;
figure; imagesc(seg1.seg(:,:,15));
% probablistic segmentations
seg5p = ft_datatype_segmentation(seg5,'segmentationstyle','probabilistic');
seg3p = ft_datatype_segmentation(seg3,'segmentationstyle','probabilistic');
seg1p = ft_datatype_segmentation(seg1,'segmentationstyle','probabilistic');

%% mesh %%%%%%%%%%%%%%%%%

%% default: triangulation
cfg=[];
cfg.numvertices = 1000;
meshA = ft_prepare_mesh(cfg,seg1p);
meshB = ft_prepare_mesh(cfg,seg1);
assert(isequalwithoutcfg(meshA,meshB),'error: 01');
assert(isfield(meshA,'pos') && isfield(meshA,'tri') && isfield(meshA,'unit'), 'Missing field(s) in mesh structure');
assert((cfg.numvertices == size(meshA.pos,1)) , 'Number of points is not equal to required');
cfg=[];
cfg.numvertices = 1000;
meshA = ft_prepare_mesh(cfg,seg3p);
meshB = ft_prepare_mesh(cfg,seg3);
assert(isequalwithoutcfg(meshA,meshB),'error: 02');
assert(isfield(meshA(1),'pos') && isfield(meshA(1),'tri') && isfield(meshA(1),'unit'), 'Missing field(s) in mesh structure');
assert((cfg.numvertices == size(meshA(1).pos,1)) && (cfg.numvertices == size(meshA(2).pos,1)) && (cfg.numvertices == size(meshA(3).pos,1)), 'Number of points is not equal to required');
cfg=[];
cfg.numvertices = 1000;
meshA = ft_prepare_mesh(cfg,seg5p);
meshB = ft_prepare_mesh(cfg,seg5);
assert(isequalwithoutcfg(meshA,meshB),'error: 03');
assert(isfield(meshA(1),'pos') && isfield(meshA(1),'tri') && isfield(meshA(1),'unit'), 'Missing field(s) in mesh structure');
assert((cfg.numvertices == size(meshA(1).pos,1)) && (cfg.numvertices == size(meshA(2).pos,1)) && (cfg.numvertices == size(meshA(3).pos,1)) && (cfg.numvertices == size(meshA(4).pos,1)) && (cfg.numvertices == size(meshA(5).pos,1)), 'Number of points is not equal to required');
figure; ft_plot_mesh(meshA,'facecolor','none');
%% method: hexahedral
cfg=[];
cfg.method = 'hexahedral';
cfg.numvertices = 1000;
meshA = ft_prepare_mesh(cfg,seg1p);
meshB = ft_prepare_mesh(cfg,seg1);
assert(isequalwithoutcfg(meshA,meshB),'error: 04');
assert(isfield(meshA,'pos') && isfield(meshA,'hex') && isfield(meshA,'unit'), 'Missing field(s) in mesh structure');
cfg=[];
cfg.method = 'hexahedral';
cfg.numvertices = 1000;
meshA = ft_prepare_mesh(cfg,seg3p);
meshB = ft_prepare_mesh(cfg,seg3);
meshA=rmfield(meshA,'tissuelabel');
meshB=rmfield(meshB,'tissuelabel');
assert(isequalwithoutcfg(meshA,meshB),'error: 05');
assert(isfield(meshA,'pos') && isfield(meshA,'hex') && isfield(meshA,'unit'), 'Missing field(s) in mesh structure');
cfg=[];
cfg.method = 'hexahedral';
cfg.numvertices = 1000;
meshA = ft_prepare_mesh(cfg,seg5p);
meshB = ft_prepare_mesh(cfg,seg5);
meshA=rmfield(meshA,'tissuelabel');
meshB=rmfield(meshB,'tissuelabel');
assert(isequalwithoutcfg(meshA,meshB),'error: 06');
assert(isfield(meshA,'pos') && isfield(meshA,'hex') && isfield(meshA,'unit'), 'Missing field(s) in mesh structure');
figure; ft_plot_mesh(meshA,'surfaceonly','yes')
%% tissue specified
cfg=[];
cfg.tissue='tissue_1';
cfg.numvertices=3000;
meshA=ft_prepare_mesh(cfg,seg3);
meshB=ft_prepare_mesh(cfg,seg3p);
assert(isequalwithoutcfg(meshA,meshB),'error: 07');
assert(isfield(meshA,'pos') && isfield(meshA,'tri') && isfield(meshA,'unit'), 'Missing field(s) in mesh structure');
cfg.method='hexahedral';
meshA=ft_prepare_mesh(cfg,seg3);
meshB=ft_prepare_mesh(cfg,seg3p);
assert(isequalwithoutcfg(meshA,meshB),'error: 08');
assert(isfield(meshA,'pos') && isfield(meshA,'hex') && isfield(meshA,'unit'), 'Missing field(s) in mesh structure');
assert(isequalwithoutcfg(meshA.tissuelabel, {'tissue_1'}), 'error:09');
cfg.tissue='tissue_2';
meshB=ft_prepare_mesh(cfg,seg3);
assert(isequalwithoutcfg(meshB.tissuelabel, {'tissue_2'}), 'error:10');
meshA=rmfield(meshA,'tissuelabel');
meshB=rmfield(meshB,'tissuelabel');
assert(~(isequalwithoutcfg(meshA,meshB)),'error: 11');
cfg.tissue={'tissue_2' 'tissue_1'};
meshC=ft_prepare_mesh(cfg,seg3);
assert(isequalwithoutcfg(meshC.tissuelabel, cfg.tissue), 'error:12');
meshC=rmfield(meshC,'tissuelabel');
assert(~(isequalwithoutcfg(meshA,meshC)),'error: 13');
assert(~(isequalwithoutcfg(meshB,meshC)),'error: 14');


% This section evaluates the isosurface/iso2mesh functionality
% create a segmentation

datamatrix = zeros(31,41,25);
datamatrix(16,21,13) = 1;

[ftver, ftdir] = ft_version;

curr_dir = pwd;
cd(fullfile(ftdir,'private'));
datamatrix = volumesmooth(datamatrix, [15 20 12]);
cd(curr_dir);

voxres = 4.5;

mri         = [];
mri.anatomy = datamatrix;
mri.dim     = size(datamatrix);
mri.unit    = 'mm';
mri.transform = inv([eye(3)./voxres (mri.dim(:)+1)./2; 0 0 0 1]);
mri.coordsys = 'acpc';
figure; ft_sourceplot([], mri);

cfg             = [];
cfg.output      = 'scalp';
cfg.scalpsmooth = 'no';
cfg.scalpthreshold = 0.2;
seg             = ft_volumesegment(cfg, mri);

mri.pow = double(seg.scalp);
cfg              = [];
cfg.funparameter ='pow';
ft_sourceplot(cfg, mri);

cfg = [];
cfg.tissue = {'scalp'};
cfg.method = 'isosurface';
bnd1       = ft_prepare_mesh(cfg, seg);
%figure;ft_plot_ortho(double(seg.scalp),'location',[0 0 0],'transform',seg.transform,'intersectmesh',bnd1,'intersectcolor','r');
cfg.method = 'iso2mesh';
cfg.numvertices = 3000;
bnd2       = ft_prepare_mesh(cfg, seg);
%figure;ft_plot_ortho(double(seg.scalp),'location',[0 0 0],'transform',seg.transform,'intersectmesh',{bnd1 bnd2},'intersectcolor','rgm');
cfg.radbound = 0.8;
bnd3       = ft_prepare_mesh(cfg, seg);
figure;ft_plot_ortho(double(seg.scalp),'location',[0 0 0],'transform',seg.transform,'intersectmesh',{bnd1 bnd2 bnd3},'intersectcolor','rgm','interpmethod','linear');
% conclusion: too large a radbound (default = 3) leads to severe
% downsampling

% at this point, I included a half-voxel shift in the iso2mesh method,
% check whether this is correct in general
seg2 = seg;
curr_dir = pwd;
cd(fullfile(ftdir, 'private'));
seg2 = volumeflip(seg2, [1 0 0]);
cd(curr_dir);
cfg = [];
cfg.dim = seg2.dim;
seg2rs = ft_volumereslice(cfg, seg2);

cfg = [];
cfg.tissue = {'scalp'};
cfg.method = 'isosurface';
bnd1b       = ft_prepare_mesh(cfg, seg2);
%figure;ft_plot_ortho(double(seg.scalp),'location',[0 0 0],'transform',seg.transform,'intersectmesh',bnd1,'intersectcolor','r');
cfg.method = 'iso2mesh';
cfg.numvertices = 3000;
bnd2b       = ft_prepare_mesh(cfg, seg2);
%figure;ft_plot_ortho(double(seg.scalp),'location',[0 0 0],'transform',seg.transform,'intersectmesh',{bnd1 bnd2},'intersectcolor','rgm');
cfg.radbound = 0.8;
bnd3b       = ft_prepare_mesh(cfg, seg2);
figure;ft_plot_ortho(double(seg2.scalp),'location',[0 0 0],'transform',seg2.transform,'intersectmesh',{bnd1b bnd1},'intersectcolor','rgm','interpmethod','linear');
figure;ft_plot_ortho(double(seg2.scalp),'location',[0 0 0],'transform',seg2.transform,'intersectmesh',{bnd2b bnd2},'intersectcolor','rgm','interpmethod','linear');
figure;ft_plot_ortho(double(seg2.scalp),'location',[0 0 0],'transform',seg2.transform,'intersectmesh',{bnd3b bnd3},'intersectcolor','rgm','interpmethod','linear');

ft_hastoolbox('fileexchange',1);
[R,T] = icp(bnd1.pos',bnd3.pos');
assert(norm(R-eye(3))<0.01);
assert(norm(T)<0.1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c = isequalwithoutcfg(a, b)
if isfield(a, 'cfg')
  a = rmfield(a, 'cfg');
end
if isfield(b, 'cfg')
  b = rmfield(b, 'cfg');
end
c = isequal(a, b);

