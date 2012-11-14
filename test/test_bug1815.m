function test_bug1815

% test the function that generates meshes also used for constructing SIMBIO FEM head models
% see http://bugzilla.fcdonders.nl/show_bug.cgi?id=1815

% TEST test_bug1815
% TEST ft_plot_mesh ft_sourceplot ft_prepare_mesh

return;

% path fieldtrip private directory
private_dir = '/home/common/matlab/fieldtrip/private/';

%%%% creating a segmentation structure %%%%%%%%%%%%



% An example segmentation 
example.dim = [98 104 101];   % slightly different numbers
example.transform = eye(4);   
example.coordsys = 'ctf';
example.unit = 'mm';

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
% adding segmentation: creating a sphere
current_dir = pwd;
cd(private_dir);
sphere = strel_bol(30);  % creates a 3D binary "sphere" with radius 30

% adjusting sphere to have the same dimensionality as segmentation
sphere = extend(sphere,example); % subfunction

if size(sphere) ~= example.dim
    error('oops, sphere dimension does not correspond to dim of segmentation')
end
example.seg = sphere;
clear sphere;
% 
cfg=[];
%cfg.interactive = 'yes';
cfg.funparameter = 'seg';
figure;
ft_sourceplot(cfg,example);
%close all;

cd(current_dir);


%%%%%%%%%%%% creating a mesh %%%%%%%%%%%%%%%%

% this is how it is implemented now (output: triangulated mesh)
cfg=[];
cfg.tissue = {'seg'};
cfg.numvertices = 3000;
bnd = ft_prepare_mesh(cfg,example);

assert(cfg.numvertices == size(bnd.pnt,1), 'Number of points is not equal to required');
assert(isfield(bnd,'pnt') || isfield(bnd,'tri') || isfield(bnd,'unit'), 'Missing field(s) in mesh structure');

figure;
ft_plot_mesh(bnd);
%close all;


%%%%%%%%%%% new implementation  %%%%%%%%%%%%

% triangle-mesh
cfg=[];
cfg.tissue = {'seg'};
cfg.numvertices = 3000;
cfg.type = 'tri';          % option for triangulated mesh generation (default)
bnd1 = ft_prepare_mesh(cfg,example);

assert(cfg.numvertices == size(bnd1.pnt,1), 'Number of points is not equal to required');
assert(isfield(bnd1,'pnt') || isfield(bnd1,'tri') || isfield(bnd1,'unit'), 'Missing field(s) in mesh structure');

assert(isequal(bnd,bnd1),'Triangulated mesh is different from default.')

figure;
ft_plot_mesh(bnd1);
%close all;

%hexahedral-mesh
cfg=[];
cfg.tissue = {'seg'};
cfg.numvertices = 3000;
cfg.type = 'hex';          % option for hexahedral mesh generation
bnd2 = ft_prepare_mesh(cfg,example);

assert(cfg.numvertices == size(bnd2.pnt,1), 'Number of points is not equal to required');
assert(isfield(bnd2,'pnt') || isfield(bnd2,'hex') || isfield(bnd2,'unit'), 'Missing field(s) in mesh structure');
assert(~(isfield(bnd2,'tri')),'Triangles in hexahedral mesh.')
assert(~(isfield(bnd2,'tet')),'Tetrahedron in hexahedral mesh.')

figure;
ft_plot_mesh(bnd2);
%close all;

% tetrahedral-mesh
cfg=[];
cfg.tissue = {'seg'};
cfg.numvertices = 3000;
cfg.type = 'tet';          % option for tetrahedral mesh generation
bnd3 = ft_prepare_mesh(cfg,example);

assert(cfg.numvertices == size(bnd3.pnt,1), 'Number of points is not equal to required');
assert(isfield(bnd3,'pnt') || isfield(bnd3,'tet') || isfield(bnd3,'unit'), 'Missing field(s) in mesh structure');
assert(~(isfield(bnd3,'tri')),'Triangles in hexahedral mesh.')
assert(~(isfield(bnd3,'hex')),'Hexagonals in tetrahedral mesh.')

figure;
ft_plot_mesh(bnd3);
%close all;
clear bnd bnd1 bnd2 bnd3


%%%%%%%% 3-layered meshes %%%%%%%%%%%%%%%

cd(private_dir);
example = rmfield(example,'seg');
sphere = strel_bol(20);
sphere = extend(sphere,example);
example.seg1 = sphere;
sphere = strel_bol(30); 
sphere = extend(sphere,example); 
example.seg2 = sphere;
sphere = strel_bol(40); 
sphere = extend(sphere,example); 
example.seg3 = sphere;
cd(current_dir);

example = ft_datatype_segmentation(example);  % non-overlapping segmentation

example_i = ft_datatype_segmentation(example,'segmentationstyle','indexed'); %indexed

cfg=[];
%cfg.interactive = 'yes';
cfg.funparameter = 'seg';
figure;
ft_sourceplot(cfg,example_i);

%%%%%%%%%%%% create mesh with 3 surfaces %%%%%%%%%%%%%%%%
cfg=[];
cfg.tissue = {'seg1' 'seg2' 'seg3'};
cfg.numvertices = [3000 2000 1000];
bnd = ft_prepare_mesh(cfg,example);
% 
assert(cfg.numvertices(1) == size(bnd(1).pnt,1) && cfg.numvertices(2) == size(bnd(2).pnt,1) && cfg.numvertices(3) == size(bnd(3).pnt,1) , 'Number of points is not equal to required');
assert(isfield(bnd,'pnt') || isfield(bnd,'tri') || isfield(bnd,'unit'), 'Missing field(s) in mesh structure');
% 
figure;
ft_plot_mesh(bnd, 'facealpha', 0.3);
 %close all;

%%%%%%%%%% new implementation %%%%%%%%%%%%%

% triangle-mesh
cfg=[];
cfg.tissue = {'seg1' 'seg2' 'seg3'};
cfg.numvertices = [3000 2000 1000];
cfg.type = 'tri';                       % default
bnd1 = ft_prepare_mesh(cfg,example); 

 % 
assert(size(bnd1,2) == 3, 'Number of layers is different from 3.')
assert(cfg.numvertices(1) == size(bnd1(1).pnt,1) && cfg.numvertices(2) == size(bnd1(2).pnt,1) && cfg.numvertices(3) == size(bnd1(3).pnt,1) , 'Number of points is not equal to required');
assert(isfield(bnd,'pnt') || isfield(bnd,'tri') || isfield(bnd,'unit'), 'Missing field(s) in mesh structure');
% 
assert(isequal(bnd,bnd1),'Triangulated mesh is different from default.')

% hexahedral-mesh
cfg=[];
cfg.tissue = {'seg1' 'seg2' 'seg3'};
cfg.numvertices = 5000;
cfg.type = 'hex';                      
bnd2 = ft_prepare_mesh(cfg,example); 
 
assert(isfield(bnd2,'pnt') || isfield(bnd2,'hex') || isfield(bnd2,'unit'), 'Missing field(s) in mesh structure');
assert(~(isfield(bnd2,'tri')),'Triangles in hexahedral mesh.')
assert(~(isfield(bnd2,'tet')),'Tetrahedron in hexahedral mesh.')

figure;
ft_plot_mesh(bnd2, 'facealpha', 0.3);

%tetrahedral-mesh
cfg=[];
cfg.tissue = {'seg1' 'seg2' 'seg3'};
cfg.numvertices = 5000;
cfg.type = 'tet';                      
bnd3 = ft_prepare_mesh(cfg,example); 

assert(isfield(bnd3,'pnt') || isfield(bnd3,'tet') || isfield(bnd3,'unit'), 'Missing field(s) in mesh structure');
assert(~(isfield(bnd3,'tri')),'Triangles in hexahedral mesh.')
assert(~(isfield(bnd3,'hex')),'Hexagonals in tetrahedral mesh.')

figure;
ft_plot_mesh(bnd3, 'facealpha', 0.3);
 
%%%%%%%% helper function %%%%%%%%%%%%%%%%%

function sphere = extend(sphere,segmentation)
for i = 1: 3
    a = size(sphere,1);
    b = size(sphere,2);
    c = size(sphere,3);
    
    if i == 1
        a = (segmentation.dim(1) - a)/2;
        if round(a) == a
            sphere = cat(1,sphere,zeros(a,b,c));
            sphere = cat(1,zeros(a,b,c),sphere);
        else
            sphere = cat(1,sphere,zeros(round(a),b,c));
            sphere = cat(1,zeros((round(a)-1),b,c),sphere);
        end
    elseif i == 2
        b = (segmentation.dim(2) - b)/2;
        if round(b) == b
            sphere = cat(2,sphere,zeros(a,b,c));
            sphere = cat(2,zeros(a,b,c),sphere);
        else
            sphere = cat(2,sphere,zeros(a,round(b),c));
            sphere = cat(2,zeros(a,(round(b)-1),c),sphere);
        end
    elseif i == 3
        c = (segmentation.dim(3) - c)/2;
        if round(c) == c
            sphere = cat(3,sphere,zeros(a,b,c));
            sphere = cat(3,zeros(a,b,c),sphere);
        else
            sphere = cat(3,sphere,zeros(a,b,round(c)));
            sphere = cat(3,zeros(a,b,(round(c)-1)),sphere);
        end
    end
end