function test_bug1820

% MEM 12gb
% WALLTIME 04:30:00

% TEST ft_prepare_mesh ft_headmodel_simbio ft_prepare_vol_sens ft_compute_leadfield

% See http://bugzilla.fcdonders.nl/show_bug.cgi?id=1820

%% create segmentation 

example.dim = [200 200 200];   % slightly different numbers
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

% create 3 spheres

radius3 = 80;
radius2 = 90;  % "skull" 
radius1 = 92;  % "skin"

for i=1:size(headpos,1)
    % from small to large
%     if norm(headpos(i,:))<radius4
%         example.seg(i) = 4;
    if norm(headpos(i,:))<radius3
        example.seg(i) = 3;
    elseif norm(headpos(i,:))<radius2
        example.seg(i) = 2;
    elseif norm(headpos(i,:))<radius1
        example.seg(i) = 1;
    end
end


%% create electrodes
% Create a set of 42 electrodes on the outer surface
currdir = pwd;
[~,ftpath] = ft_version();
cd([ftpath '/test/private/']);
r = radius1;
[pnt, tri] = icosahedron42; 
sens.pnt   = r * pnt;
sens.label = {};
nsens  = size(sens.pnt,1);
for ii = 1:nsens
  sens.label{ii} = sprintf('vertex%03d', ii);
end

cd(currdir);

%% create mesh
cfg = [];
cfg.method = 'hexahedral';
mesh = ft_prepare_mesh(cfg,example);

pos = [zeros(181,1)'; -90:1:90; zeros(181,1)']';


%% test

% first calculate stiffness matrix, will be stored in test_stiff.stiff
test_stiff = ft_headmodel_simbio(mesh,'conductivity',[0.33,0.0042,0.33]);

% now compute transfer matrix, will be stored in test_stiff.transfer
test_transfer = ft_prepare_vol_sens(test_stiff,sens);

% finally compute the leadfield
lf = ft_compute_leadfield(pos,sens,test_transfer);

% now it would be good to compute leadfields with an analytic concentric
% sphere solution and compare them
