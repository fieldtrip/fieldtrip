function test_ft_prepare_headmodel

return;

% TEST test_ft_prepare_headmodel


% input: mesh - 1, 3 or 5 compartment
%        headshape
%        hdmfile
%        elec
%        nothing

% method: singlesphere       
%   concentricspheres  
%   openmeeg           
%   bemcp              
%   dipoli             
%   asa                
%   simbio             
%   fns                
%   infinite           
%   halfspace          
%
%   singlesphere       
%   localspheres       
%   singleshell        
%   infinite       


%%%%%%%%%%%%%%%%%%%%
%% input: nothing %%
%%%%%%%%%%%%%%%%%%%%

% methods for which it is applicable: infinite, ? halfspace and
% singlesphere

cfg = [];
cfg.method = 'infinite';
vol = ft_prepare_headmodel(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% input: segmentation %%
%%%%%%%%%%%%%%%%%%%%%%%%%

% methods for which it is applicable: singlesphere and singleshell,?
% localspheres (see test_ft_prepare_localspheres)?

% example segmentation 
example.dim = [91 104 111];   % slightly different numbers
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

% create a 3 spheres and a 1 sphere example

example1 = example;
radius1 = 40;
radius2 = 30;
radius3 = 20;

for i=1:size(headpos,1)
    % from small to large
    if norm(headpos(i,:))<radius3
        example.seg(i) = 3;
        example1.seg(i) = 1;
    elseif norm(headpos(i,:))<radius2
        example.seg(i) = 2;
    elseif norm(headpos(i,:))<radius1
        example.seg(i) = 1;
    end
end
%%
cfg=[];
cfg.method = 'singlesphere';
vol = ft_prepare_headmodel(cfg,example1);
cfg.method = 'singleshell';
vol = ft_prepare_headmodel(cfg,example1);

%%%%%%%%%%%%%%%%%
%% input: mesh %%
%%%%%%%%%%%%%%%%%

cfg=[];
cfg.numvertices = 1000;
mesh = ft_prepare_mesh(cfg,example1);
cfg.numvertices = [3000 2000 1000];
bnd = ft_prepare_mesh(cfg,example);
%%

cfg=[];
cfg.method = 'singlesphere';
vol = ft_prepare_headmodel(cfg,mesh);
cfg.method = 'singleshell';
vol = ft_prepare_headmodel(cfg,mesh);
cfg=[];
cfg.method = 'concentricspheres';
vol = ft_prepare_headmodel(cfg,bnd); 
cfg=[];
cfg.method = 'openmeeg';
vol = ft_prepare_headmodel(cfg,bnd);
cfg=[];
cfg.method = 'bemcp';
vol = ft_prepare_headmodel(cfg,bnd);
cfg=[];
cfg.method = 'dipoli';
vol = ft_prepare_headmodel(cfg,bnd);

% read in the gradiometer information
hdr  = ft_read_header('/home/common/matlab/fieldtrip/data/Subject01.ds');
grad = hdr.grad;

% read in the segmented mri
load('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/segmentedmri.mat');
mri = segmentedmri; clear segmentedmri;

cfg=[];
cfg.numvertices = 1000;
mesh2 = ft_prepare_mesh(cfg,mri);

cfg         = [];
cfg.method  = 'localspheres';
cfg.grad    = grad;
vol        = ft_prepare_headmodel(cfg, mesh2);

% ? halfspace
v = mesh.pnt(1,3);
h = mesh.pnt(:,3)';
l = find(h==v);
plain = mesh.pnt(l,:);
k = find(h~=v);
point = mesh.pnt(k(1),:);

%?


%%%%%%%%%
%%

% test_ft_prepare_localspheres

hdr  = ft_read_header('/home/common/matlab/fieldtrip/data/Subject01.ds');
grad = hdr.grad;


% specify the file for the headshape
hdmfile  = '/home/common/matlab/fieldtrip/data/Subject01.shape';

% read in the headshape
shape = ft_read_headshape(hdmfile);



cfg=[];
cfg.method = 'halfspace';
vol3 = ft_prepare_headmodel(cfg,mesh); % doesn't work





