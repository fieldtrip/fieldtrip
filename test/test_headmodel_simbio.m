function test_headmodel_simbio

% TEST test_headmodel_simbio
% TEST ft_prepare_headmodel ft_prepare_mesh ft_prepare_leadfield ft_datatype_headmodel
% TEST ft_datatype_segmentation ft_datatype_sens ft_datatype_parcellation ft_datatype
% this function tests that simbio forward model works, comparing the results with a 3 concentric
% spheres model
% this script is under development, earlier version of test script written by Cristiano is in
% test_headmodel_simbio_old 

return;

ft_hastoolbox('simbio', 1);

% create a 3D volume with 3 spheres
% volume and mesh creation is adapted from test_bug1815.m
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

% create 3 spheres

radius1 = 40;
radius2 = 30;
radius3 = 20;

for i=1:size(headpos,1)
    % from small to large
    if norm(headpos(i,:))<radius3
        example.seg(i) = 3;
    elseif norm(headpos(i,:))<radius2
        example.seg(i) = 2;
    elseif norm(headpos(i,:))<radius1
        example.seg(i) = 1;
    end
end

elec = {};
radius4 = 45;   % electrode positions 5 mm outside of spheres 

sel=zeros(3);
for i=1:size(headpos,1)
    % create also sensor positions, only in positive z
    if headpos(i,3) > 0
       if norm(headpos(i,:)) == radius4
       n = size(sel,1);
       p = n+1;
            sel(p,:) = headpos(i,:);
       end
    end
end


elec.pnt = sel(4:size(sel,1),:); %first 3 rows are only zeros

for i=1:length(elec.pnt)
  elec.label{i} = sprintf('electrode%d', i);
end
elec.unit = 'mm';

elec = ft_datatype_sens(elec);


clear X Y Z headpos origin radius1 radius2 radius3  radius4 voxelpos x y z sel


figure; imagesc(example.seg(:,:,50));



% create meshes of the volume
% hexahedral

% not working yet; see test_bug1815
cfg=[];
cfg.tissue = {'seg'};
cfg.numvertices = 3000;
cfg.method = 'hex';          % option for hexahedral mesh generation
hexmesh = ft_prepare_mesh(cfg,example);

% check mesh
parcellation = ft_datatype_parcellation(hexmesh);
assert(ft_datatype(parcellation,'parcellation'),'the conversion to a parcellation failed');

% tetrahedral
% not working yet; see test_bug1815
cfg=[];
cfg.tissue = {'seg'};
cfg.numvertices = 3000;
cfg.method = 'tet';          % option for hexahedral mesh generation
tetmesh = ft_prepare_mesh(cfg,example);

% check mesh
parcellation = ft_datatype_parcellation(tetmesh);
assert(ft_datatype(parcellation,'parcellation'),'the conversion to a parcellation failed');


cfg=[];
cfg.method = 'simbio';
cfg.conductivity = [0.1 0.4 0.1];       % I just typed in some numbers, but probablt more information
                                        % is necessary on what is the unit of conductivity and if there are
                                        % "default" or suggested values
vol_hex = ft_prepare_headmodel(hexmesh);
vol_tet = ft_prepare_headmodel(tetmesh);

% check output structure

% disp(vol_hex)
%    vol.pos ...
%    vol.hex ...
%    vol.tissue ...
%    vol.tissuelabel  ...
%    vol.cond ...
%    vol.type ... 
%    vol.unit ...
%    vol.cfg ...

% ft_datatype_headmodel could be updated for simbio
assert(ft_datatype_headmodel(vol_hex),'Headmodel does not correspond to datatype headmodel');
assert(ft_datatype_headmodel(vol_tet),'Headmodel does not correspond to datatype headmodel');

assert(isequal(hexmesh.pnt,vol_hex.pos), 'Positions of headmodel is not equal to mesh points');
assert(isequal(hexmesh.hex,vol_hex.hex), 'Hexahedrons in headmodel do not correpond to mesh hexahedrons');
assert(isequal(hexmesh.tissue,vol_hex.tissue), 'Tissue field in headmodel do not correspond to tissues in mesh');
assert(isequal(hexmesh.label,vol_hex.label), 'Labels in headmodel do not correspond to mesh labels');

assert(isequal(tetmesh.pnt,vol_tet.pos), 'Positions of headmodel is not equal to mesh points');
assert(isequal(tetmesh.tet,vol_tet.tet), 'Tetrahedrons in headmodel do not correpond to mesh tetrahedrons');
assert(isequal(tetmesh.tissue,vol_tet.tissue), 'Tissue field in headmodel do not correspond to tissues in mesh');
assert(isequal(tetmesh.label,vol_tet.label), 'Labels in headmodel do not correspond to mesh labels');

%%%%
 
% % compute the lead fields 
cfg=[];
cfg.vol = vol_tet;
cfg.elec = elec;
lf_tet  = ft_prepare_leadfield(cfg);

cfg=[];
cfg.vol = vol_hex;
cfg.elec = elec;
lf_hex  = ft_prepare_leadfield(cfg);


% compare it with concentric spheres

example = ft_datatype_segmentation(example,'segmentationstyle','probabilistic');
cfg=[];
cfg.tissue = {'tissue_1', 'tissue_2', 'tissue_3'};
cfg.numvertices = [3000 2000 1000];
cfg.method = 'tri';                       % default
bnd = ft_prepare_mesh(cfg,example); 

cfg=[];
cfg.method = 'concentricspheres';
cfg.conductivity = [0.1 0.4 0.1];

vol_cc = ft_prepare_headmodel(cfg,bnd);

cfg=[];
cfg.vol = vol_cc;
cfg.elec = elec;
lf_cc  = ft_prepare_leadfield(cfg);

% plot leadfields

