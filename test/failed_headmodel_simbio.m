function failed_headmodel_simbio

% WALLTIME 00:45:00
% MEM 16gb

% TEST test_headmodel_simbio
% TEST ft_prepare_headmodel ft_prepare_mesh ft_prepare_leadfield ft_prepare_sourcemodel 

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

% this function tests that simbio forward model works, comparing the results with a 3 concentric spheres model
% intial version by Lilla Magyari 2013

%% create example segmentation
example.dim = [250 250 250];   % slightly different numbers
example.transform = eye(4);
example.coordsys = 'ctf';
example.unit = 'mm';
example.seg = zeros(example.dim);

x = round(example.dim(1)/2);
y = round(example.dim(2)/2);
z = round(example.dim(3)/2);

x = round(x);
y = round(y);
z = round(z);

origin = [x y z];

example.transform(1:4,4)   = [-origin(:); 1];  
[X, Y, Z] = ndgrid(1:example.dim(1), 1:example.dim(2), 1:example.dim(3));
voxelpos = [X(:) Y(:) Z(:)];
headpos = ft_warp_apply(example.transform, voxelpos);

% create 3 spheres

radius1 = 100;
radius2 = 110;  % "skull" 
radius3 = 112;  % "skin"

for i=1:size(headpos,1)
    if norm(headpos(i,:))<radius1
        example.seg(i) = 1;
    elseif norm(headpos(i,:))<radius2
        example.seg(i) = 2;
    elseif norm(headpos(i,:))<radius3
        example.seg(i) = 3;
    end
end

% Create a set of 42 electrodes on the outer surface
currdir = pwd;
cd ~/matlab/fieldtrip/test/private/;
r = radius3;
[pnt, tri] = icosahedron42; 
sens.pnt   = r * pnt;
sens.label = {};
nsens  = size(sens.pnt,1);
for ii = 1:nsens
  sens.label{ii} = sprintf('vertex%03d', ii);
end

cd(currdir);

%% create mesh
cfg=[];
cfg.method = 'hexahedral';    
hexmesh = ft_prepare_mesh(cfg,example);

clear example;
%% create headmodel
cfg=[];
cfg.method = 'simbio';
cfg.conductivity = [0.33 0.01 0.43];  
vol_hex = ft_prepare_headmodel(cfg,hexmesh);

% check output structure
assert(isequal(hexmesh.pnt,vol_hex.pos), 'Positions of headmodel is not equal to mesh points');
assert(isequal(hexmesh.hex,vol_hex.hex), 'Hexahedrons in headmodel do not correpond to mesh hexahedrons');
assert(isequal(hexmesh.tissue,vol_hex.tissue), 'Tissue field in headmodel do not correspond to tissues in mesh');
assert(isequal(hexmesh.tissuelabel,vol_hex.tissuelabel), 'Labels in headmodel do not correspond to mesh labels');
clear hexmesh;

%% compute leadfield and compare it with concetric spheres

% concetricspheres
volcs.unit ='mm';
volcs.o    = [0 0 0];
volcs.r    = [100 110 112];
volcs.c    = [0.33 0.01 0.43];  
volcs.type = 'concentricspheres';


% sourcemodel
% Define positions of dipoles (along z axis)
zp = linspace(0,100,50)'; % generates 50 points bw 0 and 100 equally spaced
pos = [zeros(size(zp)) zeros(size(zp)) zp];
 
% Define the corresponding spatial grid 
cfg = [];
cfg.inwardshift = 5; % otherwise some sourcepoints will be on boundary or really close to it
cfg.grid.pos    = pos;
cfg.vol = volcs;
cfg.sens = sens;
gridp = ft_prepare_sourcemodel(cfg);

cfg=[];
cfg.vol = vol_hex;
cfg.elec = sens;
cfg.grid = gridp;
lf_hex  = ft_prepare_leadfield(cfg);
clear vol_hex;

cfg=[];
cfg.vol = volcs;
cfg.elec = sens;
cfg.grid = gridp;
lf_cc  = ft_prepare_leadfield(cfg);

grid1 = lf_cc;
clear lf_cc;
grid2 = lf_hex;
clear lf_hex;

num_dipole = size(grid1.leadfield,2)-3; % the last three points fall outside of sourcespace

corrcoeffs = zeros(num_dipole,3);
for i=1:num_dipole
    vec1 = grid1.leadfield{i}(:,1);
    vec2 = grid2.leadfield{i}(:,1);
    coeff = corr(vec1,vec2);
    corrcoeffs(i,1)=coeff;
    vec1 = grid1.leadfield{i}(:,2);
    vec2 = grid2.leadfield{i}(:,2);
    coeff = corr(vec1,vec2);
    corrcoeffs(i,2)=coeff;
    vec1 = grid1.leadfield{i}(:,3);
    vec2 = grid2.leadfield{i}(:,3);
    coeff = corr(vec1,vec2);
    corrcoeffs(i,3)=coeff;
end    
    
assert(all(corrcoeffs(:)>0.9),'Large difference in leadfield values between concentric spheres and simbio headmodel');

