function test_ft_prepare_mesh

% test ft_prepare_mesh also used for constructing SIMBIO FEM head models
% see also http://bugzilla.fcdonders.nl/show_bug.cgi?id=1815

% TEST test_ft_prepare_mesh
% TEST ft_prepare_mesh ft_datatype_segmentation ft_plot_mesh

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
headpos = warp_apply(example.transform, voxelpos);

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

% default: triangulation
cfg=[];
cfg.numvertices = 1000;
meshA = ft_prepare_mesh(cfg,seg1p);
meshB = ft_prepare_mesh(cfg,seg1);
assert(isequal(meshA,meshB),'error: 01');
assert(isfield(meshA,'pnt') && isfield(meshA,'tri') && isfield(meshA,'unit'), 'Missing field(s) in mesh structure');
assert((cfg.numvertices == size(meshA.pnt,1)) , 'Number of points is not equal to required');
cfg=[];
cfg.numvertices = 1000;
meshA = ft_prepare_mesh(cfg,seg3p);
meshB = ft_prepare_mesh(cfg,seg3);
assert(isequal(meshA,meshB),'error: 02');
assert(isfield(meshA(1),'pnt') && isfield(meshA(1),'tri') && isfield(meshA(1),'unit'), 'Missing field(s) in mesh structure');
assert((cfg.numvertices == size(meshA(1).pnt,1)) && (cfg.numvertices == size(meshA(2).pnt,1)) && (cfg.numvertices == size(meshA(3).pnt,1)), 'Number of points is not equal to required');
cfg=[];
cfg.numvertices = 1000;
meshA = ft_prepare_mesh(cfg,seg5p);
meshB = ft_prepare_mesh(cfg,seg5);
assert(isequal(meshA,meshB),'error: 03');
assert(isfield(meshA(1),'pnt') && isfield(meshA(1),'tri') && isfield(meshA(1),'unit'), 'Missing field(s) in mesh structure');
assert((cfg.numvertices == size(meshA(1).pnt,1)) && (cfg.numvertices == size(meshA(2).pnt,1)) && (cfg.numvertices == size(meshA(3).pnt,1)) && (cfg.numvertices == size(meshA(4).pnt,1)) && (cfg.numvertices == size(meshA(5).pnt,1)), 'Number of points is not equal to required');
figure; ft_plot_mesh(meshA,'facecolor','none');
% method: hexahedral
cfg=[];
cfg.method = 'hexahedral';
cfg.numvertices = 1000;
meshA = ft_prepare_mesh(cfg,seg1p);
meshB = ft_prepare_mesh(cfg,seg1);
assert(isequal(meshA,meshB),'error: 04');
assert(isfield(meshA,'pnt') && isfield(meshA,'hex') && isfield(meshA,'unit'), 'Missing field(s) in mesh structure');
cfg=[];
cfg.method = 'hexahedral';
cfg.numvertices = 1000;
meshA = ft_prepare_mesh(cfg,seg3p);
meshB = ft_prepare_mesh(cfg,seg3);
meshA=rmfield(meshA,'tissuelabel');
meshB=rmfield(meshB,'tissuelabel');
assert(isequal(meshA,meshB),'error: 05');
assert(isfield(meshA,'pnt') && isfield(meshA,'hex') && isfield(meshA,'unit'), 'Missing field(s) in mesh structure');
cfg=[];
cfg.method = 'hexahedral';
cfg.numvertices = 1000;
meshA = ft_prepare_mesh(cfg,seg5p);
meshB = ft_prepare_mesh(cfg,seg5);
meshA=rmfield(meshA,'tissuelabel');
meshB=rmfield(meshB,'tissuelabel');
assert(isequal(meshA,meshB),'error: 06');
assert(isfield(meshA,'pnt') && isfield(meshA,'hex') && isfield(meshA,'unit'), 'Missing field(s) in mesh structure');
figure; ft_plot_mesh(meshA,'edgeonly','yes')

